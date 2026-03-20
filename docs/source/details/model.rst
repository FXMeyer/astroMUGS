The Model class
***************

``Model`` is the main user-facing class. It orchestrates the full pipeline:
building the disk structure, running RADMC3D, writing Nautilus input files,
and converting results back to RADMC3D format.

Overview
========

A typical workflow:

.. code-block:: python

   from astromugs.modeling.Model import Model

   # 1. Create model with parameters and grid
   model = Model(params, grid)

   # 2. Build disk structure and write RADMC3D inputs
   model.add_disk(dust=dust_model)
   model.write_continuum()

   # 3. Run radiative transfer
   model.run_thermal_radmc3d()

   # 4. Write chemistry inputs (coupling RADMC3D -> Nautilus)
   model.write_nautilus(
       static=True, coupling_dens=True, coupling_temp=True, coupling_av=True,
       multi_grain=True, ...
   )

   # 5. After Nautilus runs, convert results back
   model.convert_nautilus2radmc(species='CO')


Writing Nautilus input files
=============================

The ``write_nautilus`` method is the main entry point for generating Nautilus input files.
It handles:

- **Density coupling** (``coupling_dens=True``): remaps RADMC3D dust densities onto the chemistry grid
  and computes gas number densities, either from the dust-to-gas ratio or from a custom surface density file.

- **Temperature coupling** (``coupling_temp=True``): remaps RADMC3D dust temperatures onto the chemistry grid.

- **Extinction coupling** (``coupling_av=True``): computes :math:`A_V(z)` from the RADMC3D local radiation field.

- **Floor values**: ``min_gas_density`` and ``min_av`` prevent unphysically low values that cause
  solver convergence issues in Nautilus.

- **UV ceiling**: ``max_uv`` caps the UV field to avoid extreme stiffness in low-density regions.

For each radial point on the chemistry grid, the method writes:

- ``1D_static.dat``: gas density, temperatures, :math:`A_V`, UV factor, dust properties
- ``1D_grain_sizes.in``: per-grain-size inverse abundances and temperatures (multi-grain mode)
- ``parameters.in``: Nautilus solver parameters
- Chemical network files, element list, and initial abundances


Converting Nautilus output to RADMC3D
======================================

After Nautilus has run, ``convert_nautilus2radmc`` reads the molecular abundances
and maps them back onto the RADMC3D spherical grid to produce ``numberdens_<species>.inp``
files for line radiative transfer.


API Reference
=============

.. autoclass:: astromugs.modeling.Model.Model
   :members:
   :undoc-members:
   :show-inheritance:
