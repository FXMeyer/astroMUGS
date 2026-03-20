RADMC3D -- Nautilus coupling
****************************

The ``nautilus.coupling`` module handles the coordinate transformation and interpolation
between the RADMC3D spherical grid and the Nautilus cylindrical chemistry grid.

Overview
========

RADMC3D operates on a spherical grid :math:`(r, \theta, \phi)`, while Nautilus uses
a cylindrical grid :math:`(R, z)`. The coupling functions remap physical quantities
(density, temperature, UV extinction) from one grid to the other.

The main coupling functions are:

- ``dust_density``: remaps dust mass densities and converts to number densities
- ``dust_temperature`` / ``dust_temperature_single``: remaps dust temperatures
- ``av_z``: computes visual extinction from the RADMC3D local radiation field
- ``to_spherical``: remaps Nautilus chemical abundances back to the RADMC3D grid


Coordinate conversion
=====================

For a point at cylindrical coordinates :math:`(R, z)`, the corresponding spherical coordinates are:

.. math::

   d = \sqrt{R^2 + z^2}, \qquad \theta = \arccos\left(\frac{z}{d}\right)

.. important::

   The RADMC3D grid radii (stored in ``grid.r``) are in **AU**, while the chemistry grid
   coordinates (``grid.rchem``, ``grid.zchem``) are passed in **cm**. The coupling functions
   convert internally using ``autocm`` (AU to cm conversion factor).


Density coupling
================

``dust_density`` uses bilinear interpolation (``RegularGridInterpolator``) to remap
RADMC3D dust densities onto the chemistry grid. For each grain size:

1. Extract the 2D density field :math:`\rho(r, \theta)` from the RADMC3D data
2. Interpolate onto the chemistry grid points :math:`(d, \theta)_\mathrm{chem}`
3. Convert mass density to number density: :math:`n_d = \rho_d / m_\mathrm{grain}`

The gas number density is derived either from the dust-to-gas ratio or from a custom
surface density profile.


Temperature coupling
====================

Two variants exist:

- ``dust_temperature`` / ``dust_temperature_single``: use physical :math:`z` coordinates (for ``set_chem_grid``)
- ``dust_temperature_disk`` / ``dust_temperature_single_disk``: use scale-height coordinates (for ``set_chemdisk_grid``)

Both use nearest-neighbor lookup to find the closest RADMC3D cell for each chemistry grid point.


UV extinction coupling
======================

``av_z`` computes the visual extinction by comparing the local UV radiation field
(from RADMC3D's ``local_field``) to the unattenuated blackbody spectrum. The extinction
in magnitudes is:

.. math::

   A_V = -2.5 \log_{10}\left(\frac{F_\mathrm{local}}{F_\mathrm{unattenuated}}\right)


API Reference
=============

.. automodule:: astromugs.nautilus.coupling
   :members:
   :undoc-members:
