from __future__ import annotations

import numpy as np
from dataclasses import fields as dataclass_fields
from typing import TYPE_CHECKING

from .. constants.constants import M_sun, Ggram, autocm

if TYPE_CHECKING:
    from ..utils.thermal import ControlParams  # only needed for annotation


def control(params: ControlParams):
    with open("thermal/radmc3d.inp", "w") as f:
        for fld in dataclass_fields(params):
            value = getattr(params, fld.name)
            if value is not None:
                if isinstance(value, float) and value.is_integer():
                    f.write(f"{fld.name} = {int(value)}\n")
                else:
                    f.write(f"{fld.name} = {value}\n")

def stars(rstar, mstar, lam, xstar, ystar, zstar, tstar=None, fstar=None):
    nstars = len(rstar)
    nlam = len(lam)

    f = open("thermal/stars.inp","w")

    f.write(str(2)+"\n")
    f.write("{0:d}  {1:d}\n".format(nstars, nlam))

    for istar in range (nstars):
        f.write("{0:e}   {1:e}   {2:e}   {3:e}   {4:e}\n".format(rstar[istar], \
                mstar[istar], xstar[istar], ystar[istar], zstar[istar]))

    for ilam in range(nlam):
        f.write("{0:e}\n".format(lam[ilam]))

    for istar in range(nstars):
        if (tstar[istar] != 0):
            f.write("{0:f}\n".format(-tstar[istar]))
        else:
            for i in range(nlam):
                f.write("{0:e}\n".format(fstar[ilam]))

    f.close()

def wavelength_micron(lam):
    nlam = len(lam)
    f = open("thermal/wavelength_micron.inp","w")
    f.write("{0:d}\n".format(nlam))
    for ilam in range(nlam):
        f.write("{0:f}\n".format(lam[ilam]))
    f.close()

def mcmono_wavelength_micron(lam_mono):
    nlam = len(lam_mono)

    f = open("thermal/mcmono_wavelength_micron.inp","w")

    f.write("{0:d}\n".format(nlam))
    for ilam in range(nlam):
        f.write("{0:f}\n".format(lam_mono[ilam]))

    f.close()

def amr_grid(x, y, z, gridstyle="regular", coordsystem="cartesian"):
    nx = x.size-1
    ny = y.size-1
    nz = z.size-1

    incl_x = int(nx > 1)
    incl_y = int(ny > 1)
    incl_z = int(nz > 1)

    f = open("thermal/amr_grid.inp","w")

    f.write(str(1)+"\n")

    if (gridstyle == "regular"):
        f.write("0\n")
    elif (gridstyle == "octtree"):
        f.write("1\n")
    elif (gridstyle == "amr"):
        f.write("10\n")

    if (coordsystem == "cartesian"):
        f.write("0\n")
    elif (coordsystem == "spherical"):
        f.write("100\n")
    elif (coordsystem == "cylindrical"):
        f.write("200\n")

    f.write("0\n")
    f.write("{0:d}  {1:d}  {2:d}\n".format(incl_x, incl_y, incl_z))
    f.write("{0:d}  {1:d}  {2:d}\n".format(nx, ny, nz))

    if (gridstyle == "octtree"):
        print("OctTree grids not yet implemented.")
    elif (gridstyle == "amr"):
        print("Layer-style AMR grids not yet implemented.")

    for i in range(nx+1):
        f.write("{0:12.9e}\n".format(x[i]))
    for i in range(ny+1):
        f.write("{0:12.9e}\n".format(y[i]))
    for i in range(nz+1):
        f.write("{0:12.9e}\n".format(z[i]))

    # Insert extra info for octtree and amr grids here...

    f.close()

def dustopac(opacity):
    '''
    Desc: write dustopac.inp
    Args: opacity
    '''
    nspecies = len(opacity)

    f = open("thermal/dustopac.inp","w")

    f.write("2\n")
    f.write("{0:d}\n".format(nspecies))
    f.write("==============================================================\n")
    for i in range(nspecies):
        filetype = opacity[i].split("/")[1].split("_")[0]
        species = opacity[i].split("_")[1].split(".")[0]

        if (filetype == "dustkappa"):
            f.write("1\n")
        elif (filetype == "dustkapscatmat"):
            f.write("10\n")
        elif (filetype == "dustopac"):
            f.write("-1\n")

        f.write("0\n")
        f.write("{0:s}\n".format(species))
        f.write("----------------------------------------------------------\n")

    f.close()


def dust_density(density, gridstyle="regular"):
    '''
    Desc: write dust_density.inp
    Args: density
    '''
    nstructures = len(density) #disk, envelope...
    print('number of structures (disk, envelope...): ', nstructures, '\n')
    nspecies = 0
    for istruc in range(nstructures):
        nspecies += len(density[istruc]) #nb of species in all structures
        print('number of species in structure {}: '.format(istruc+1), len(density[istruc]))
    print('total number of grain species: ', nspecies)

    if (gridstyle == "regular"):
        nx, ny, nz = density[0][0].shape
        ncells = nx*ny*nz

    f = open("thermal/dust_density.inp","w")
    f.write("1\n")
    f.write("{0:d}\n".format(ncells))
    f.write("{0:d}\n".format(nspecies))

    for istruc in range(nstructures): #loop over structures (disk, envelope...)
        for ispec in range(len(density[istruc])): #loop over the dust species within the given structure.
            if (gridstyle == "regular"):
                for iz in range(nz):
                    for iy in range(ny):
                        for ix in range(nx):
                            f.write("{0:e}\n".format(density[istruc][ispec, ix,iy,iz]))
    f.close()


def external_rad(isrf):
    '''
    Desc: write external_source.inp
    Args: Interstellar radiation field
    '''
    nlam = len(isrf[0])
    factor = 1e0 #artificially multiply by a factor just to see the impact of different ISRF intensities. Will be removed in future updates
    f = open("thermal/external_source.inp","w")

    f.write("2\n")
    f.write("{0:d}\n".format(nlam))
    for ilam in range(nlam):
        f.write("{0:f}\n".format(isrf[0,ilam]*1e-3))
    for ilam in range(nlam):
        f.write("{0:e}\n".format(isrf[1,ilam]*factor))

    f.close()

def accretion_heating(x, y, z, accretionheating, gridstyle="regular"):
    '''
    Desc: Write viscous accretion heating in the file heatsource.inp 
    Args: Viscous accretion heating structure
    '''
    nx = x.size-1
    ny = y.size-1
    nz = z.size-1
    f = open("thermal/heatsource.inp","w")
    f.write("1\n")
    f.write("{0:d}\n".format(nx*ny*nz))

    if (gridstyle == "regular"):
        for iz in range(nz):
            for iy in range(ny):
                for ix in range(nx):
                    f.write("{0:e}\n".format(accretionheating[ix,iy,iz]))

    f.close()


def numberdens_mol(numberdens, species='CO', gridstyle="regular"):
    '''
    Desc: write numberdens_XXX.inp, where XXX is a chemical species
    Args: 
    '''
    if (gridstyle == "regular"):
        nx, ny = numberdens.shape
        ncells = nx*ny

    print('writing numberdens_{}.inp...'.format(species))

    f = open("thermal/numberdens_{}.inp".format(species),"w")
    f.write("1\n")
    f.write("{0:d}\n".format(ncells))
    if (gridstyle == "regular"):
            for iy in range(ny):
                for ix in range(nx):
                    f.write("{0:e}\n".format(numberdens[ix,iy]))
    f.close()

def lines(species='CO', format='leiden'):
    '''
    Desc: write lines.inp
    Args: species, format
    '''
    f = open("thermal/lines.inp","w")
    f.write("2\n")
    f.write("1\n")
    f.write("{} {} 0 0 0".format(species,format))
    f.close()


def gas_velocity(star_mass, r, theta, phi, object="disk"):
    '''
    Desc: write gas_velocity.inp
    Args: species, format
    '''    

    nr, ntheta, nphi = len(r), len(theta), len(phi)
    print(nr, ntheta, nphi)
    # Create 3D grids
    rr, tt, pp = np.meshgrid(r*autocm, theta, phi, indexing='ij')

    # Keplerian azimuthal velocity
    vphi = np.sqrt(Ggram * star_mass*M_sun / rr)

    # Convert to Cartesian
    vx = -vphi * np.sin(pp)
    vy =  np.zeros_like(vx)
    vz = vphi * np.cos(pp)

    # Write to file
    with open('thermal/gas_velocity.inp', 'w') as f:
        f.write("1\n")  # iformat = 1 (ASCII)
        f.write(f"{nr*ntheta*nphi}\n")

        for k in range(nphi):
            for j in range(ntheta):
                for i in range(nr):
                    f.write(f"{vx[i,j,k]:.6e} {vy[i,j,k]:.6e} {vz[i,j,k]:.6e}\n")