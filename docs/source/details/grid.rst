The Grid
********

The ``Grid`` object stores all spatial information: coordinates, densities, temperatures, and dust properties.
It is the central data container that connects the disk model to both RADMC3D and Nautilus.

Overview
========

A ``Grid`` is automatically created when you instantiate a ``Model``. You interact with it to:

- Define the spatial domain (spherical grid for RADMC3D, chemistry grid for Nautilus)
- Store dust and gas densities, temperatures, and opacities
- Access coordinates for plotting and analysis

There are two types of grids:

- **Spherical grid** (``set_spherical_grid``): used by RADMC3D for radiative transfer. Defined in :math:`(r, \theta, \phi)`.
- **Chemistry grid** (``set_chem_grid`` or ``set_chemdisk_grid``): used by Nautilus. Defined in :math:`(r, z)` cylindrical coordinates.


Setting up a spherical grid
===========================

The spherical grid is defined by providing cell edges in :math:`r`, :math:`\theta`, and :math:`\phi`:

.. code-block:: python

   model.grid.set_spherical_grid(
       r_edge=r_edges,        # radial cell edges [AU]
       theta_edge=theta_edges, # polar angle cell edges [rad]
       phi_edge=phi_edges      # azimuthal cell edges [rad]
   )

The grid stores both cell edges and cell centers (midpoints).


Setting up a chemistry grid
============================

Two options exist depending on how the vertical coordinate is defined:

**Option 1: Physical altitude** (``set_chem_grid``) -- vertical coordinates are in AU.
Each radius can have a different vertical extent, producing a 2D array for ``zchem``.

.. code-block:: python

   model.grid.set_chem_grid(r=radii, z0=0, zmax=z_max, msize=70, nbcells=70)

**Option 2: Scale-height units** (``set_chemdisk_grid``) -- vertical coordinates are in units of the gas scale height.
The ``zchem`` array is 1D (same normalized grid at every radius).

.. code-block:: python

   model.grid.set_chemdisk_grid(r=radii, max_H=4, nz_chem=70)

.. tip::

   Use ``set_chem_grid`` when coupling with RADMC3D outputs (the coupling functions expect physical coordinates).
   Use ``set_chemdisk_grid`` when using the parametric disk model directly (via ``add_chemdisk``).


API Reference
=============

.. autoclass:: astromugs.modeling.Grid.Grid
   :members:
   :undoc-members:
   :show-inheritance:
