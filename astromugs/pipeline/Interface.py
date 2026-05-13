import sys
import os
import inspect
import numpy as np

from astromugs.pipeline.Pipeline import Pipeline
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


    def add_chemmodel(self, chempath="chemistry/", itime=0, species='CO', reader=None):
        """Load an existing chemistry model onto the grid.

        Parameters
        ----------
        chempath : str, optional
            Path to the chemistry model directory. Default is
            ``'chemistry/'``.
        itime : int, optional
            Index of the time snapshot to use from the chemistry output.
            Default is 0.
        species : str, optional
            Chemical species of interest (e.g., ``'CO'``, ``'H2O'``).
            Default is ``'CO'``.
        reader : object or None, optional
            Reader responsible for parsing chemistry data. Defaults to
            ``self.nautilus.read`` when None.

        Raises
        ------
        FileNotFoundError
            If required files are missing in *chempath*.
        KeyError
            If required parameters are missing in the chemistry data.
        RuntimeError
            For any other errors encountered during processing.
        """

        if reader is None:
            reader = self.nautilus.read  # Default to the current reader if none is provided

        try:
            radii = reader.radii(chempath=chempath)
            self.grid.add_existingchemradii(radii)

            parameters = reader.parameters(self.grid.chemradii, chempath=chempath)
            if 'nb_grains_1D' not in parameters:
                raise KeyError("The parameter 'nb_grains_1D' is missing in the chemistry parameters.")
            self.grid.add_existingchemparam(parameters)

            grid = reader.grid(radlist=self.grid.chemradii, nb_sizes=int(parameters['nb_grains_1D']), itime=itime, species=species, chempath=chempath)
            self.grid.add_existingchemmodel(grid, species)
        except FileNotFoundError as e:
            raise FileNotFoundError(f"Required file not found in {chempath}: {e}")
        except KeyError as e:
            raise KeyError(f"Missing required parameter in chemistry data: {e}")
        except Exception as e:
            raise RuntimeError(f"An error occurred while adding the chemistry model: {e}")
