RADMC3D -- Nautilus coupling
****************************

The ``nautilus.coupling`` module handles the coordinate transformation and interpolation
between the RADMC3D spherical grid and the Nautilus Cartesian chemistry grid.

Overview
========

RADMC3D operates on a spherical grid :math:`(r, \theta, \phi)`, while Nautilus uses
a Cartesian grid :math:`(R, z)`. The full workflow has two directions:

1. **RADMC3D → Nautilus** (pre-chemistry): physical quantities (density, temperature,
   UV field, visual extinction) are remapped from the spherical RT grid onto the
   Cartesian chemistry columns and written as Nautilus input files by
   ``write_nautilus()``.

2. **Nautilus → RADMC3D** (post-chemistry): chemical abundances from the completed
   Nautilus run are converted back to number densities on the spherical RT grid and
   written as ``numberdens_XXX.inp`` files for line radiative transfer with RADMC3D.

The main coupling functions are:

- ``dust_density``: remaps dust mass densities and converts to number densities
- ``dust_temperature`` / ``dust_temperature_single``: remaps dust temperatures
- ``av_z``: computes visual extinction from the RADMC3D local radiation field
- ``to_spherical``: remaps Nautilus chemical abundances back to the RADMC3D grid


Coordinate conversion
=====================

For a point at Cartesian coordinates :math:`(R, z)`, the corresponding spherical coordinates are:

.. math::

   d = \sqrt{R^2 + z^2}, \qquad \theta = \arccos\left(\frac{z}{d}\right)

.. important::

   The RADMC3D grid radii (stored in ``grid.r``) are in **AU**, while the chemistry grid
   coordinates (``grid.rchem``, ``grid.zchem``) are passed in **cm**. The coupling functions
   convert internally using ``autocm`` (AU to cm conversion factor).


Converting RADMC-3D to Nautilus input files
===========================================

Density coupling
----------------

``dust_density`` uses bilinear interpolation (``RegularGridInterpolator``) to remap
RADMC3D dust densities onto the chemistry grid. For each grain size:

1. Extract the 2D density field :math:`\rho(r, \theta)` from the RADMC3D data
2. Interpolate onto the chemistry grid points :math:`(d, \theta)_\mathrm{chem}`
3. Convert mass density to number density: :math:`n_d = \rho_d / m_\mathrm{grain}`

The gas number density is derived either from the dust-to-gas ratio or from a custom
surface density profile.


Temperature coupling
--------------------

Two variants exist:

- ``dust_temperature`` / ``dust_temperature_single``: use physical :math:`z` coordinates (for ``set_chem_grid``)
- ``dust_temperature_disk`` / ``dust_temperature_single_disk``: use scale-height coordinates (for ``set_chemdisk_grid``)

Both use bilinear interpolation (``RegularGridInterpolator``) similarly to the dust density.


UV extinction coupling
----------------------

``av_z`` computes the visual extinction by comparing the local UV radiation field
(from RADMC3D's ``local_field``) to the unattenuated blackbody spectrum. The extinction
in magnitudes is:

.. math::

   A_V = -2.5 \log_{10}\left(\frac{F_\mathrm{local}}{F_\mathrm{unattenuated}}\right)


Converting Nautilus abundances back to RADMC-3D
================================================

After your Nautilus simulation has completed, use ``add_chemistry()`` followed by
``convert_nautilus2radmc()`` to create the ``numberdens_XXX.inp`` files that RADMC-3D
needs for line radiative transfer.

.. code-block:: python

    import astromugs.pipeline as pipeline

    pipe1 = pipeline.Interface()
    pipe1.add_thermal_path(thermpath)   # path to the RADMC-3D thermal output directory
    pipe1.add_chemical_path(chempath)   # path to the Nautilus chemistry directory

    # Load the chemistry output.
    # itime selects which timestep is used for the conversion (default: -1 = last output).
    pipe1.add_chemistry(itime=-1)

    # Write numberdens_CO.inp into thermpath
    pipe1.convert_nautilus2radmc(species='CO', numberdens=True)

    # Multiple species at once
    pipe1.convert_nautilus2radmc(species=['CO', 'HCO+', 'N2H+'], numberdens=True)

``add_chemistry()`` reads the binary ``abundances.out`` files from all completed
``XXAU/`` sub-folders, converts fractional abundances to number densities
(:math:`n_X = n_\mathrm{H} \times (n_X / n_\mathrm{H})` in cm\ :sup:`-3`), and
stores the result in ``pipe1.grid.chemmodel``, keyed first by species name then by
radius in AU.

``convert_nautilus2radmc()`` then calls ``nautilus.coupling.to_spherical()`` to
remap each chemistry column onto the full RADMC-3D spherical grid by **bilinear
interpolation** (``RegularGridInterpolator``), and writes the result to
``numberdens_XXX.inp`` in ``thermpath``.

The interpolation works in two steps:

1. Each vertical column (one per radius) is resampled onto a common :math:`z`
   grid built from the union of all column :math:`z` arrays. This handles the fact
   that different radii may have different numbers of vertical points after surface
   truncation (controlled by ``min_gas_density`` in ``write_nautilus()``). Points
   above a column's truncation height are set to zero; points below the midplane
   use the midplane value.

2. The resulting regular 2-D field :math:`n_X(R_\mathrm{cyl}, z)` is passed to
   ``RegularGridInterpolator``, which evaluates it at every spherical cell centre
   :math:`(d \sin\theta,\, |d \cos\theta|)` in one vectorised call.

.. note::

   The Nautilus chemistry grid typically has far fewer radial points than the
   RADMC-3D grid. Bilinear interpolation produces smooth abundance maps between
   chemistry columns, avoiding the step discontinuities that nearest-neighbour
   assignment would introduce. Cells inside the inner chemistry radius are set
   to ``1e-20``; cells outside the outermost radius are set to ``0``.

.. note::

   ``add_chemistry()`` silently skips any radius whose ``abundances.out`` is not yet
   present (e.g. still running on a cluster). You can call it again later and then
   repeat ``convert_nautilus2radmc()`` once more radii have finished.


API Reference
=============

.. automodule:: astromugs.nautilus.coupling
   :members:
   :undoc-members:
