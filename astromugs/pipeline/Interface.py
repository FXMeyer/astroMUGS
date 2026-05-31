import sys
import os
import inspect
import numpy as np
import xarray as xr

from astromugs.pipeline.Pipeline import Pipeline
from astromugs import nautilus
from astromugs.modeling.Star import Star
from astromugs.modeling.Disk import Disk
from astromugs.modeling.Envelope import Envelope
from astromugs.modeling.InterstellarRadFields import InterstellarRadFields
from astromugs.constants.constants import autocm, M_sun, R_sun, L_sun


class Interface(Pipeline):
    """High-level interface for assembling radiative transfer models.

    Provides methods to add physical components (star, disk, envelope,
    interstellar radiation field) and chemical models to a grid, building
    on the base ``Pipeline`` class.
    """

    #def __init__(self):
        #self.params = StructureParams()
        #self.envelope_params = self.params.envelope
        #self.disk_params = self.params.disk

    def add_chemical_path(self, chemicalpath):
        """Set the path to the chemical data directory.

        Parameters
        ----------
        chemicalpath : str
            Path to the directory containing chemical model files.
        """
        self.chempath = chemicalpath

    def add_thermal_path(self, thermalpath):
        """Set the path to the thermal dust opacity data directory.

        Parameters
        ----------
        thermalpath : str
            Path to the directory containing thermal dust opacity files.
        """
        self.thermalpath = thermalpath

    def add_star(self):
        """Add a central star to the grid using thermal parameters.

        Reads mass, luminosity, temperature, and position from
        ``self.thermalparams.star`` and creates a ``Star`` instance on
        the grid.
        """
        mass = self.thermalparams.star.mass
        luminosity = self.thermalparams.star.luminosity
        temperature = self.thermalparams.star.temperature
        x = self.thermalparams.star.x
        y = self.thermalparams.star.y
        z = self.thermalparams.star.z
        self.grid.add_star(Star(mass=mass, luminosity=luminosity, \
                temperature=temperature, x=x, y=y, z=z))

    def add_isrf(self, cut=2.e-1, d78=True, vdb82=True):
        """Add an interstellar radiation field (ISRF) to the grid.

        Parameters
        ----------
        cut : float, optional
            Short-wavelength cutoff for the ISRF, in microns. Default is
            0.2 microns.
        d78 : bool, optional
            Include the Draine (1978) UV component. Default is True.
        vdb82 : bool, optional
            Include the van Dishoeck & Black (1982) component. Default
            is True.
        """
        self.isrf = InterstellarRadFields(cut, d78, vdb82)
        self.grid.add_isrf(self.isrf.create_isrf(self.grid.lam))

    def add_disk(self, dust=None, **kwargs):
        """Add a disk component and its dust density to the grid.

        Parameters
        ----------
        dust : object or None, optional
            Dust opacity model. If it has a ``path`` attribute set to
            ``'thermal/'``, the path is replaced by ``self.thermalpath``.
        **kwargs
            Keyword arguments forwarded as overrides to the default disk
            parameters (e.g., ``rin``, ``rout``, ``mass``).
        """
        # Create a copy of the defaults
        params = self.params.disk

        # Apply overrides given by user
        for key, val in kwargs.items():
            setattr(params, key, val)

        # Propagate thermalpath only if the user didn't set a custom path
        if dust is not None and hasattr(dust, 'path') and dust.path == 'thermal/':
            dust.path = self.thermalpath

        self.disk = Disk(params=params, dust=dust)

        if (len(self.grid.dust) > 0):
            self.grid.add_dustdensity(self.disk.density_d(self.grid.r, self.grid.theta, self.grid.phi))
        else:
            print('WARNING: no dust model as input. add your dust model in add_diskt and in the grid (grid.add_dust(dust)) \
                  before creating the disk structure.')

    def add_internalheating(self, **kwargs):
        """Add viscous accretion heating to the grid.

        Parameters
        ----------
        **kwargs
            Keyword arguments forwarded as overrides to the default disk
            parameters (e.g., ``mdot``, ``alpha``).
        """
        # Create a copy of the defaults
        params = self.params.disk
        # Apply overrides given by user
        for key, val in kwargs.items():
            setattr(params, key, val)

        self.grid.add_accretionheating(self.disk.viscous_accretion_heating(self.grid.r, self.grid.theta, self.grid.phi))

    def add_chemdisk(self, dust=None, **kwargs):
        """Add a disk component for chemical modeling on the chemistry grid.

        Populates the chemistry grid with dust number density, gas number
        density, gas temperature, scale height, and vertical visual
        extinction.

        Parameters
        ----------
        dust : object or None, optional
            Dust opacity model. If it has a ``path`` attribute set to
            ``'thermal/'``, the path is replaced by ``self.thermalpath``.
        **kwargs
            Keyword arguments forwarded as overrides to the default disk
            parameters.
        """
        # Create a copy of the defaults
        params = self.params.disk

        # Apply overrides given by user
        for key, val in kwargs.items():
            setattr(params, key, val)

        if dust is not None and hasattr(dust, 'path') and dust.path == 'thermal/':
            dust.path = self.thermalpath

        self.chemdisk = Disk(params=params, dust=dust)

        self.grid.add_dustdensity_chem(self.chemdisk.numberdensity_d(self.grid.rchem, self.grid.zchem))
        self.grid.add_dustdensity_single_chem(self.chemdisk.numberdensity_d_single(self.grid.rchem, self.grid.zchem))
        self.grid.add_gasdensity_chem(self.chemdisk.numberdensity(self.grid.rchem, self.grid.zchem))
        self.grid.add_gastemperature_chem(self.chemdisk.temp_altitude(self.grid.rchem, self.grid.zchem))
        self.grid.add_hg_chem(self.chemdisk.scaleheight(self.grid.rchem))
        self.grid.add_avz(self.chemdisk.av_z(self.grid.lam, self.grid.dustdensity_chem[0], self.grid.rchem, self.grid.zchem))


    def add_envelope(self, dust=None, **kwargs):
        """Add an envelope component and its dust density to the grid.

        Parameters
        ----------
        dust : object or None, optional
            Dust opacity model. If it has a ``path`` attribute set to
            ``'thermal/'``, the path is replaced by ``self.thermalpath``.
            The dust density is only added to the grid when ``dust`` is
            not None.
        **kwargs
            Keyword arguments forwarded as overrides to the default
            envelope parameters.
        """
        # Create a copy of the defaults
        params = self.params.envelope

        # Apply overrides given by user
        for key, val in kwargs.items():
            setattr(params, key, val)

        if dust is not None and hasattr(dust, 'path') and dust.path == 'thermal/':
            dust.path = self.thermalpath

        self.envelope = Envelope(params, dust=dust)

        if dust:
            self.grid.add_dustdensity(
                self.envelope.density_d(self.grid.r, self.grid.theta, self.grid.phi)
            )


    def add_chemistry(self, chempath=None, itime=-1):
        """Read chemistry output and store results on the model.

        Auto-detects whether *chempath* is a single chemistry run (0-D or
        single 1-D) or a disk model composed of multiple radial folders
        ending with ``AU/``.

        In **single mode**, ``self.chemistry`` is a dict with an xarray
        ``abundances`` DataArray containing all species and timesteps.

        In **disk mode**, ``self.chemistry`` is a dict keyed by radius
        (int, in AU), where each value has the same structure.
        ``self.grid.chemmodel`` is also populated so that
        ``convert_nautilus2radmc()`` works out of the box.

        Parameters
        ----------
        chempath : str or None, optional
            Path to the chemistry output directory. Defaults to
            ``self.chempath``.
        itime : int, optional
            Timestep index used to build the grid-level chemistry model
            (``self.grid.chemmodel``) for disk mode. Negative indices
            follow NumPy convention (``-1`` = last timestep). Default
            is ``-1``.
        """
        if chempath is None:
            chempath = self.chempath

        # Detect disk mode: look for subfolders ending with "AU"
        au_folders = sorted(
            [d for d in os.listdir(chempath)
             if d.endswith('AU') and os.path.isdir(os.path.join(chempath, d))],
            key=lambda x: int(x.replace('AU', ''))
        )

        if au_folders:
            self._read_chemistry_disk(chempath, au_folders, itime)
        else:
            self._read_chemistry_single(chempath)

    def _read_chemistry_single(self, chempath):
        """Read a single chemistry folder (0-D or 1-D run)."""
        self.chemistry = nautilus.read.abundances_binary(
            os.path.join(chempath, 'abundances.out')
        )
        species_names = nautilus.read.species_names(chempath)
        self.chemistry['species'] = species_names
        self.chemistry['abundances'] = xr.DataArray(
            self.chemistry['abundances'],
            dims=['time', 'species', 'spatial'],
            coords={
                'time':    self.chemistry['time'],
                'species': species_names,
            }
        )

    def _read_chemistry_disk(self, chempath, au_folders, itime):
        """Read a disk chemistry model (multiple *AU/ subfolders)."""
        radii = [int(d.replace('AU', '')) for d in au_folders]
        self.grid.add_existingchemradii(radii)

        # Read species names from the first folder that has completed output
        species_names = None
        for folder in au_folders:
            candidate = os.path.join(chempath, folder)
            if (os.path.isfile(os.path.join(candidate, 'species.out')) and
                    os.path.isfile(os.path.join(candidate, 'abundances.out'))):
                species_names = nautilus.read.species_names(candidate)
                break

        if species_names is None:
            print('[add_chemistry] no completed radius found (species.out / abundances.out missing in all folders).')
            self.chemistry = {}
            return

        self.chemistry = {}

        for radius, folder in zip(radii, au_folders):
            folder_path = os.path.join(chempath, folder)
            ab_path = os.path.join(folder_path, 'abundances.out')

            if not os.path.isfile(ab_path):
                print(f'  [add_chemistry] skipping {folder}: abundances.out not found')
                continue

            # Read binary output
            chem = nautilus.read.abundances_binary(ab_path)
            chem['species'] = species_names
            chem['abundances'] = xr.DataArray(
                chem['abundances'],
                dims=['time', 'species', 'spatial'],
                coords={
                    'time':    chem['time'],
                    'species': species_names,
                }
            )
            self.chemistry[radius] = chem

            # Read vertical grid from 1D_static.dat
            static_data = np.loadtxt(
                os.path.join(folder_path, '1D_static.dat'), comments='!'
            )
            z = static_data[:, 0]

            # Build self.grid.chemmodel for convert_nautilus2radmc
            nH = chem['H_number_density'][itime]            # (spatial,)
            abundances_t = chem['abundances'].isel(time=itime).values  # (nb_species, spatial)

            for idx_sp, sp in enumerate(species_names):
                if sp not in self.grid.chemmodel:
                    self.grid.chemmodel[sp] = {}
                self.grid.chemmodel[sp][radius] = {
                    'z': z,
                    'nH': nH,
                    'numberdens_species': nH * abundances_t[idx_sp],
                }
