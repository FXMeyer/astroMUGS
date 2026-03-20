"""
_____________________________________________________________________________________________________________
file name: Grid
@author: P.Sheehan. Adapted by S. Gavino for chemistry codes.
last update: Aug 2022
language: PYTHON 3.8
short description:  class Grid for young stellar objects modeling. 
_____________________________________________________________________________________________________________
"""



from __future__ import annotations

import numpy as np

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from astromugs.utils.struct import DiskParams #only need this for annotations
    from astromugs.utils.thermal import WaveParams #only need this for annotations

    

class Grid:
    """Spatial and wavelength grid for young stellar object modeling.

    Holds the coordinate system, physical quantities (density, temperature,
    radiation field), and optional chemistry sub-grids used by radiative
    transfer and chemical codes.

    Attributes
    ----------
    params : DiskParams
        Disk/envelope structural parameters (radii, resolution, etc.).
    wave : WaveParams
        Wavelength grid parameters.
    density : list of ndarray
        Total (gas + dust) density components, in g/cm^3.
    dustdensity : list of ndarray
        Dust density components for radiative transfer, in g/cm^3.
    gasdensity_chem : list of ndarray
        Gas density components for chemistry, in g/cm^3.
    dustdensity_chem : list of ndarray
        Dust density components for chemistry, in g/cm^3.
    dustdensity_single_chem : list of ndarray
        Single-grain-population dust density for chemistry, in g/cm^3.
    hg_chem : list of ndarray
        Gas scale height arrays for the chemistry grid, in AU.
    tgas_chem : list of ndarray
        Gas temperature arrays for the chemistry grid, in K.
    temperature : list of ndarray
        Dust temperature components, in K.
    localfield : list of ndarray
        Local radiation field components.
    avz : list of ndarray
        Vertical visual extinction components, in mag.
    stars : list
        Stellar source objects.
    isrf : list
        Interstellar radiation field components.
    dust : list
        Dust opacity/property objects.
    accretionheating : list of ndarray
        Viscous accretion heating rate components.
    chemradii : list of ndarray
        Radial coordinate arrays imported from existing chemistry models, in AU.
    chemparam : list
        Parameter sets imported from existing chemistry models.
    chemmodel : dict
        Chemical abundance grids keyed by species name.
    """

    def __init__(self,
                 params: DiskParams,
                 wave: WaveParams
    ):
        """Initialize the Grid with disk parameters and wavelength settings.

        Parameters
        ----------
        params : DiskParams
            Disk/envelope structural parameters.
        wave : WaveParams
            Wavelength grid parameters.
        """
        self.params = params
        self.wave = wave
        self.density = []
        self.dustdensity = []
        self.gasdensity_chem = []
        self.dustdensity_chem = []
        self.dustdensity_single_chem = []
        self.hg_chem = []
        self.tgas_chem = []
        self.temperature = []
        self.localfield = []
        self.avz = []
        self.stars = []
        self.isrf = []
        self.dust = []
        self.accretionheating = []
        self.chemradii = []
        self.chemparam = []
        self.chemmodel = {}

    #-----
    # ADD/CREATE GRID STRUCTURES THAT EXIST OR THAT HAVE BEEN CREATED.
    #-----

    def add_star(self, star):
        """Append a stellar source to the grid.

        Parameters
        ----------
        star : object
            Stellar source object containing luminosity, temperature, etc.
        """
        self.stars.append(star)

    def add_isrf(self, isrf):
        """Append an interstellar radiation field component.

        Parameters
        ----------
        isrf : object
            Interstellar radiation field specification.
        """
        self.isrf.append(isrf)

    def add_temperature(self, temperature):
        """Append a dust temperature component.

        Parameters
        ----------
        temperature : ndarray
            Dust temperature array, in K.
        """
        self.temperature.append(temperature)

    def add_localfield(self, localfield):
        """Append a local radiation field component.

        Parameters
        ----------
        localfield : ndarray
            Local radiation field array.
        """
        self.localfield.append(localfield)

    def add_density(self, density):
        """Append a total density component.

        Parameters
        ----------
        density : ndarray
            Total (gas + dust) density array, in g/cm^3.
        """
        self.density.append(density)

    def add_dustdensity(self, density):
        """Append a dust density component for radiative transfer.

        Parameters
        ----------
        density : ndarray
            Dust density array, in g/cm^3.
        """
        self.dustdensity.append(density)

    def add_dustdensity_chem(self, density):
        """Append a dust density component for chemistry.

        Parameters
        ----------
        density : ndarray
            Dust density array for the chemistry grid, in g/cm^3.
        """
        self.dustdensity_chem.append(density)

    def add_dustdensity_single_chem(self, density):
        """Append a single-grain-population dust density for chemistry.

        Parameters
        ----------
        density : ndarray
            Single-grain dust density array for the chemistry grid, in g/cm^3.
        """
        self.dustdensity_single_chem.append(density)

    def add_gasdensity_chem(self, density):
        """Append a gas density component for chemistry.

        Parameters
        ----------
        density : ndarray
            Gas density array for the chemistry grid, in g/cm^3.
        """
        self.gasdensity_chem.append(density)

    def add_gastemperature_chem(self, gas_temperature):
        """Append a gas temperature component for chemistry.

        Parameters
        ----------
        gas_temperature : ndarray
            Gas temperature array for the chemistry grid, in K.
        """
        self.tgas_chem.append(gas_temperature)

    def add_hg_chem(self, hg):
        """Append a gas scale height array for chemistry.

        Parameters
        ----------
        hg : ndarray
            Gas scale height array, in AU.
        """
        self.hg_chem.append(hg)

    def add_avz(self, av_z):
        """Append a vertical visual extinction component.

        Parameters
        ----------
        av_z : ndarray
            Vertical visual extinction array, in mag.
        """
        self.avz.append(av_z)

    def add_dust(self, dust):
        """Append a dust opacity/property object.

        Parameters
        ----------
        dust : object
            Dust opacity or property specification.
        """
        self.dust.append(dust)

    def add_accretionheating(self, q_visc):
        """Append a viscous accretion heating rate component.

        Parameters
        ----------
        q_visc : ndarray
            Viscous heating rate array.
        """
        self.accretionheating.append(q_visc)

    def add_existingchemradii(self,existingchemradii):
        """Append radial coordinates from an existing chemistry model.

        Parameters
        ----------
        existingchemradii : ndarray
            Array of radial and vertical grid points from an existing
            chemistry model, in AU.
        """
        self.chemradii.append(existingchemradii)
    
    def add_existingchemparam(self,existingchemparam):
        """Append parameters from an existing chemistry model.

        Parameters
        ----------
        existingchemparam : object
            Parameter set from an existing chemistry model.
        """
        self.chemparam.append(existingchemparam)

    def add_existingchemmodel(self,existingchemmodel, species):
        """Store abundance data from an existing chemistry model for a species.

        Parameters
        ----------
        existingchemmodel : ndarray
            Abundance grid for the given species.
        species : str
            Chemical species name used as the dictionary key.
        """
        self.chemmodel[species] = existingchemmodel

        

    #-----
    # SET GRID STRUCTURES.
    #-----

    def set_cartesian_grid(self, xmin, xmax, nx):
        """Build a uniform Cartesian grid and return edges and cell centres.

        The same range is used for all three axes (x, y, z).

        Parameters
        ----------
        xmin : float
            Minimum coordinate value, in AU.
        xmax : float
            Maximum coordinate value, in AU.
        nx : int
            Number of grid edges along each axis.

        Returns
        -------
        edges : ndarray, shape (3, nx)
            Cell edge coordinates stacked as (x, y, z).
        centres : ndarray, shape (3, nx-1)
            Cell centre coordinates stacked as (w1, w2, w3).
        """
        self.coordsystem = "cartesian"

        x = np.linspace(xmin, xmax, nx)
        y = np.linspace(xmin, xmax, nx)
        z = np.linspace(xmin, xmax, nx)

        w1 = 0.5*(x[0:x.size-1] + x[1:x.size])
        w2 = 0.5*(y[0:y.size-1] + y[1:y.size])
        w3 = 0.5*(z[0:z.size-1] + z[1:z.size])

        return np.stack((x, y, z)), np.stack((w1, w2, w3))

    def set_spherical_grid(self, r_edge=None, theta_edge=None, phi_edge=None, log=True):
        """Set the spherical grid.

        If edge arrays are provided (e.g. read from an existing amr_grid.inp),
        use them directly. Otherwise, compute edges from the DiskParams
        (rin, rout, nr, ntheta, nphi).

        Parameters
        ----------
        r_edge : array-like, optional
            Radial cell edges in au. If None, computed from params.
        theta_edge : array-like, optional
            Polar cell edges in radians. If None, computed from params.
        phi_edge : array-like, optional
            Azimuthal cell edges in radians. If None, computed from params.
        log : bool
            Use logarithmic spacing for radial edges (only when computing from params).
        """
        self.coordsystem = "spherical"

        if r_edge is not None:
            r_edge = np.asarray(r_edge)
            theta_edge = np.asarray(theta_edge)
            phi_edge = np.asarray(phi_edge)
        else:
            rmin = self.params.rin
            rmax = self.params.rout
            nr = self.params.nr
            ntheta = self.params.ntheta
            nphi = self.params.nphi
            if log:
                r_edge = np.logspace(np.log10(rmin), np.log10(rmax), nr, base=10)
            else:
                r_edge = np.linspace(rmin, rmax, nr)

            theta_edge = np.linspace(0.0, np.pi, ntheta)
            phi_edge = np.linspace(0.0, 2*np.pi, nphi)

        self.r = 0.5*(r_edge[0:r_edge.size-1] + r_edge[1:r_edge.size])
        self.theta = 0.5*(theta_edge[0:theta_edge.size-1] + theta_edge[1:theta_edge.size])
        self.phi = 0.5*(phi_edge[0:phi_edge.size-1] + phi_edge[1:phi_edge.size])

        self.nr = r_edge.size
        self.ntheta = theta_edge.size
        self.nphi = phi_edge.size
        self.r_edge = r_edge
        self.theta_edge = theta_edge
        self.phi_edge = phi_edge

    def set_chemdisk_grid(self, r, max_H=4, nz_chem=64):
        """Build a vertical chemistry grid for a disk model.

        The vertical coordinate runs from ``max_H`` scale heights down
        to the midplane with uniform spacing in normalised units.

        Parameters
        ----------
        r : array_like
            Radial positions, in AU.
        max_H : float, optional
            Upper limit of the grid expressed in gas scale heights.
        nz_chem : int, optional
            Number of vertical grid points (same at all radii).
        """
        #hg = self.disk.scaleheight(np.array(r))
        pts = np.arange(0, nz_chem, 1)
        #zchem = np.ones((len(rchem), nb_points))

        #hh, ptpt = np.meshgrid(hchem, pts)
        z = (1. - (2.*pts/(2.*nz_chem - 1.)))*max_H#*Hg

        self.rchem = np.array(r)
        self.zchem = z
        self.nz_chem = nz_chem


    def set_chem_grid(self, r, z0=0, zmax=None, msize=None, nbcells=70):
        """Build a custom spatial grid for chemistry (disk, envelope, etc.).

        Two modes are available. If *zmax* is given, the vertical extent is
        uniform at all radii. If *msize* is given, the maximum altitude at
        each radius follows a spherical envelope boundary.  Vertical
        coordinates are stored in decreasing order as required by Nautilus.

        Parameters
        ----------
        r : ndarray
            1-D array of radial positions, in AU. Must lie within the
            RADMC-3D model domain.
        z0 : float, optional
            Minimum altitude, in AU.
        zmax : float, optional
            Maximum altitude applied uniformly at all radii, in AU.
        msize : float, optional
            Model size (sphere radius) in AU. When provided, the maximum
            altitude at each radius is computed as sqrt(msize^2 - r^2).
        nbcells : int, optional
            Number of vertical cells (same for all radii).
        """
        self.nz_chem = nbcells
 
        if zmax != None: #if user provides maximum altitude in au, it sets the maximum z at all radii. By default the minimum value z0 is 0. That creates a 2D structure.
            self.rchem = np.array(r)
            z = np.zeros((len(self.rchem), nbcells))
            for idx, rval in enumerate(self.rchem):
                z[idx,:] = np.linspace(z0, zmax, nbcells)
            self.zchem = np.flip(np.array(z), axis=1) #flip because the user gives increasing values and nautilus needs decreasing values.

        if msize != None: #if user wants the altitude max such that the model is a sphere i.e. zmax at each radius follows the spherical structure.
            zmax = np.sqrt(msize**2 - r**2)  # gives max altitude at each radius (polar coordinates).
            zmax = zmax[:-1] # remove the last value because it is zmax = 0 au.

            self.rchem = r[:-1] # because we removed the last zmax.
            z = np.zeros((len(self.rchem), nbcells))

            for idx, zmax_x in enumerate(zmax):
                z[idx,:] = np.linspace(0, zmax_x, nbcells)

            self.zchem = np.flip(z, axis=1) #flip because the user gives increasing values and nautilus needs decreasing values.

            #dz = np.tan(np.pi/nbthetas)*self.rchem #0.0175 is pi/number of thetas.

            # import matplotlib.pyplot as plt
            # fig = plt.figure(figsize=(8, 8.))
            # ax = fig.add_subplot(111)
            # ax.plot(self.rchem, zmax)
            # plt.xlim(0, 5005)
            # plt.ylim(0, 5005)
            # #plt.show()

    def set_wavelength_grid(self, log=True):
        """Build the wavelength grid used by the radiative transfer.

        Parameters
        ----------
        log : bool, optional
            If True, use logarithmic spacing; otherwise linear spacing.
        """
        lmin = self.wave.lmin
        lmax = self.wave.lmax
        nlam = self.wave.nlam
        if log:
            self.lam = np.logspace(np.log10(lmin), np.log10(lmax), \
                    nlam)
        else:
            self.lam = np.linspace(lmin, lmax, nlam)


    def set_mcmonowavelength_grid(self, log=True):
        """Build the monochromatic Monte Carlo wavelength grid.

        Parameters
        ----------
        log : bool, optional
            If True, use logarithmic spacing; otherwise linear spacing.
        """
        lmin_mono = self.wave.lmin_mono
        lmax_mono = self.wave.lmax_mono
        nlam_mono = self.wave.nlam_mono
        if log:
            self.monolam = np.logspace(np.log10(lmin_mono), np.log10(lmax_mono), \
                    nlam_mono)
        else:
            self.monolam = np.linspace(lmin_mono, lmax_mono, nlam_mono)
