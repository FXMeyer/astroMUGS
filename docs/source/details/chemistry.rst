Gas-grain simulations
************

Write NAUTILUS input files
==========================

Writing Nautilus input files requires three preparatory steps before calling
``write_nautilus()``: adding a dust model to the grid, defining the radial positions
of the chemistry columns, and setting the vertical grid for each column. Each step is
described below.


Step 1 — Add a dust model
--------------------------

.. important::

   A dust model **must** be attached to the pipeline before any chemistry input
   files can be written. Nautilus works with grain number densities and therefore
   needs to know the grain sizes and material density in order to compute
   the dust-to-gas abundance ratios written to ``1D_static.dat`` and
   ``1D_grain_sizes.in``.

Use :class:`~astromugs.dust.CustomDistrib` or :class:`~astromugs.dust.MRNDistrib`to load grain sizes from a file (e.g.
produced by a DustPy simulation) and register the model on the grid:

.. code-block:: python

    from astromugs import dust

    d = dust.CustomDistrib(
        rho_m    = 1.67,                     # grain material density [g cm⁻³]
        filename = "outputs/dust_sizes.inp", # one grain size per line, in microns
        units    = 'microns'
    )
    pipe.grid.add_dust(d)

``rho_m`` is the bulk grain material density in g cm\ :sup:`-3`. Common values are
1.5–3.3 g cm\ :sup:`-3` depending on composition (silicates, carbonaceous grains,
ice mixtures). ``filename`` points to a plain-text file with one grain size per line
in the units specified by ``units`` (``'microns'``, ``'mm'``, ``'cm'``, or ``'m'``).


Step 2 — Define the radial grid
---------------------------------

The radial positions of the chemistry columns are a plain NumPy array of radii in AU.
Choose them to sample the physical gradients of your disk model — logarithmic spacing
works well over a wide radial range:

.. code-block:: python

    import numpy as np

    r_chem = np.array([5, 10, 15, 20, 30, 50, 75, 100, 150, 200, 300])  # AU


Step 3 — Set the vertical chemistry grid
------------------------------------------

Two methods are available to build the vertical grid.

``set_chemdisk_grid`` — scale-height based (recommended for disk models)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The grid runs from ``max_H`` gas scale heights down to the midplane.
The number of vertical points ``nz_chem`` is the same at all radii.

.. code-block:: python

    pipe.grid.set_chemdisk_grid(
        r       = r_chem,  # radial positions [AU]
        max_H   = 4,       # upper boundary in gas scale heights (default: 4)
        nz_chem = 64,      # number of vertical points (default: 64)
        stretch = 1.5      # > 1 concentrates points near the midplane (default: 1.0)
    )

The ``stretch`` parameter redistributes the vertical points without changing their
total number: ``stretch > 1`` packs resolution near the midplane where density
gradients are steepest; ``stretch < 1`` packs resolution near the disk surface.

``set_chem_grid`` — altitude based (envelope or custom geometries)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When the maximum altitude cannot be expressed as a multiple of the scale height
(e.g. envelope models, or when coupling directly to a hydrodynamical grid), use
``set_chem_grid`` instead. Two modes are available:

- **Uniform ceiling** (``zmax`` in AU): the same maximum altitude at every radius.
- **Spherical ceiling** (``msize`` in AU): the maximum altitude at each radius follows
  :math:`z_\mathrm{max}(r) = \sqrt{r_\mathrm{model}^2 - r^2}`, giving a spherical
  outer boundary.

.. code-block:: python

    # Uniform ceiling — same zmax at all radii
    pipe.grid.set_chem_grid(r=r_chem, zmax=200, nbcells=70, stretch=1.0)

    # Spherical ceiling — outer boundary is a sphere of radius msize AU
    pipe.grid.set_chem_grid(r=r_chem, msize=300, nbcells=70, stretch=1.0)

In both cases, ``nbcells`` sets the number of vertical points and ``stretch`` controls
the spacing in the same way as for ``set_chemdisk_grid``.

.. note::

   Vertical coordinates are stored in **decreasing order** (surface → midplane) as
   required by Nautilus. Both grid methods handle this automatically.


Step 4 — Write Nautilus input files
-------------------------------------

With the dust model and grids in place, call ``write_nautilus()`` to generate one
complete set of Nautilus input files per radial column:

.. code-block:: python

    pipe.write_nautilus(
        stop_time  = 1e6,         # chemistry integration time [yr]
        nb_outputs = 64,          # number of output snapshots (log-spaced by default)
        multi_grain = True,       # use NMGC multi-grain mode
        network    = '/path/to/network',  # copy network files from this directory
        abundances = 'atomic',    # initial abundance preset
    )

A sub-directory named ``<radius>AU/`` is created inside the chemistry path for each
radial column, and the following files are written into it:

- ``1D_static.dat`` — vertical physical structure (gas density, temperature, :math:`A_V`, UV, grain properties)
- ``parameters.in`` (or ``parameters_nmgc.in``) — Nautilus run control parameters
- ``1D_grain_sizes.in`` — grain radii, inverse abundances, and temperatures per bin (multi-grain mode)
- ``temperatures.dat`` — grain temperature table (multi-grain mode)
- ``abundances.in`` — initial chemical abundances
- ``element.in`` — elemental composition
- network files — ``gas_species.in``, ``grain_species.in``, ``gas_reactions.in``, etc.

Key arguments
~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 22 15 63

   * - Argument
     - Default
     - Description
   * - ``stop_time``
     - ``3e6``
     - Total chemistry integration time in years.
   * - ``nb_outputs``
     - ``64``
     - Number of output timesteps written to ``abundances.out``. Outputs are log-spaced between ``start_time`` (1 yr) and ``stop_time``.
   * - ``multi_grain``
     - ``True``
     - If ``True``, use NMGC multi-grain mode: one set of grain surface reactions per size bin, with grain parameters written to ``1D_grain_sizes.in``. Set to ``False`` for standard single-grain Nautilus.
   * - ``network``
     - ``True``
     - Network file source. ``True`` copies from the built-in astroMUGS network. A path string (e.g. ``network='/my/network'``) copies from that directory instead. ``False`` skips copying (if the files are already in place).
   * - ``abundances``
     - ``'atomic'``
     - Initial abundance preset. ``'atomic'`` starts from bare atomic gas (solar composition). Pass a file path to use a custom ``abundances.in``.
   * - ``coupling_dens``
     - ``False``
     - If ``True``, gas and dust number densities are derived from the RADMC-3D dust density grid rather than from the analytic disk model. Useful when post-processing hydrodynamical simulations that provide their own density structure.
   * - ``min_gas_density``
     - ``1.0``
     - Minimum gas number density [cm\ :sup:`-3`]. With the default ``cut_cap=True``, vertical grid points below this threshold are removed from ``1D_static.dat`` entirely (see :ref:`safezone`). With ``cut_cap=False``, they are floored to this value instead.
   * - ``exclude_bins``
     - ``None``
     - 1-indexed list of grain size bins to drop from ``1D_grain_sizes.in`` at every radius. Use this to remove bins that are completely depleted by radial drift (e.g. ``exclude_bins=[8, 9, 10]``).

.. note::

   The complete list of arguments and their defaults is documented in the API
   reference under :meth:`~astromugs.pipeline.Pipeline.write_nautilus`. The table
   above covers only the arguments most commonly adjusted in practice.



Best practices
==============

Stay in the Nautilus safe zone
------------------------------

Nautilus can encounter convergence issues when placed in extreme physical conditions.
For example, a very strong UV field with little attenuation can cause the chemical computation
to fail or make the run extremely long.

Nautilus was initially designed for gas-grain simulations in the interstellar medium.
Protoplanetary disk models are particularly prone to pushing Nautilus into stiff chemistry
regimes, especially when post-processing hydrodynamical simulations. Many such simulations
— whether computed from a single core collapse or from a 1-D surface density profile —
can produce extreme conditions in the disk upper atmosphere: very low gas densities, strong
UV fields, and negligible visual extinction. Nautilus is not designed for these environments.

Here are examples where Nautilus will struggle to converge:

- Very low density (:math:`n_\mathrm{H} \lesssim 10^2\,\mathrm{cm}^{-3}`) combined with a strong UV field and low :math:`A_V`: the spread between fast photodissociation rates and slow two-body formation rates (which scale as :math:`n_\mathrm{H}^2`) creates an extremely stiff system.
- Near-balanced UV shielding: when the UV flux and :math:`A_V` partially cancel, small abundance changes shift the shielding, which shifts the photodissociation rates, which shifts the abundances — a feedback loop that prevents convergence.

To mitigate these issues, astroMUGS provides floor and ceiling parameters (``min_gas_density``, ``min_av``, ``max_uv``)
that clamp physical quantities to user-defined bounds. These can be passed as keyword arguments to ``write_nautilus()``. 
Applying these bounds should not impact the results of your simulations, because they typically apply to regions of your model that  
don't correspond to observed physical environments in stellar formation regions and/or that you are not interested in (very high altitudes, very large/small radii, etc.).  

.. warning::

   Floor and ceiling values must be chosen with care. Poorly chosen bounds can
   *create* convergence problems rather than prevent them. Here are two examples:

   - **Av floor too high for a low-UV model.** A high :math:`A_V` floor pushes the chemistry
     into a shielded regime, which is stable when the UV field is strong. But in a low-UV model,
     it can place the chemistry in a transition zone where photodissociation and molecule
     formation rates nearly balance, creating a feedback loop that makes the solver struggle.
     Always consider the strength of your UV field when choosing an :math:`A_V` floor.

   - **Density floor in the wrong regime.** Raising the gas density floor does not always
     help. At intermediate densities (e.g., :math:`n_\mathrm{H} \sim 10^5\,\mathrm{cm}^{-3}`),
     grain-surface adsorption and desorption rates become competitive with gas-phase reactions,
     introducing a new source of chemical stiffness. A density floor can thus shift the
     chemistry from one problematic regime into another.

.. note::

   A future version of astroMUGS will include a scaling law that ties the floor and ceiling
   values of :math:`n_\mathrm{H}`, :math:`A_V`, and the UV factor together, ensuring the
   chemistry always stays in a well-posed regime. This can be enabled with
   ``safezone=True`` in ``write_nautilus()``. However, this option may alter the physical
   conditions of your model, so use it with care — you may prefer to set specific floor or
   ceiling values manually.


Read NAUTILUS output binaries
=============================

Once your Nautilus simulation has completed (or while it is still running — see :ref:`monitoring`),
astroMUGS provides a high-level interface to load and explore the output.

The recommended entry point is :class:`~astromugs.pipeline.Interface`, which auto-detects
whether your chemistry directory contains a single run or a full disk model, and stores
everything in a consistent structure ready for analysis and radiative-transfer post-processing.


Loading chemistry results
--------------------------

.. code-block:: python

    import astromugs.pipeline as pipeline

    pipe = pipeline.Interface()
    pipe.add_chemical_path('path_to_your_chemistry/')  # path to the chemistry directory
    pipe.add_chemistry()                               # reads all completed radii

``add_chemistry()`` scans the chemistry directory for sub-folders whose names end with ``AU``
(e.g. ``5AU/``, ``30AU/``, ``100AU/``). If such folders exist it operates in **disk mode**;
otherwise it falls back to **single-run mode** for 0-D or single-column models.

Any radius whose ``abundances.out`` is not yet present (e.g. still running on a cluster) is
silently skipped and reported in the terminal — you can call ``add_chemistry()`` again later to
pick up newly finished radii.

The optional ``itime`` argument selects the timestep used to populate
``pipe.grid.chemmodel``, which is consumed by ``convert_nautilus2radmc()`` for line
radiative-transfer post-processing (default: ``itime=-1``, the last output).

.. code-block:: python

    pipe.add_chemistry(itime=-1)   # use the last timestep for the RADMC-3D grid


Disk mode: accessing data per radius
--------------------------------------

After ``add_chemistry()``, results are stored in ``pipe.chemistry``, a **dictionary keyed by
radius in AU** (integer):

.. code-block:: python

    # List all radii that were successfully loaded
    print(list(pipe.chemistry.keys()))
    # e.g. [5, 10, 30, 100, 200, 300]

    # Access a single radius
    chem_30 = pipe.chemistry[30]   # dict for the 30 AU column

Each per-radius entry is itself a dictionary with the following keys:

.. list-table::
   :header-rows: 1
   :widths: 25 20 55

   * - Key
     - Shape
     - Description
   * - ``'time'``
     - ``(nb_timesteps,)``
     - Output times **in seconds**. Divide by ``3.156e7`` (or ``365.25*24*3600``) to convert to years.
   * - ``'abundances'``
     - ``(nb_timesteps, nb_species, nz)``
     - Fractional abundances :math:`n_X / n_\mathrm{H}` stored as an **xarray DataArray** with named dimensions and coordinates (see below).
   * - ``'H_number_density'``
     - ``(nb_timesteps, nz)``
     - Hydrogen number density :math:`n_\mathrm{H}` [cm\ :sup:`-3`] as read back from Nautilus.
   * - ``'gas_temperature'``
     - ``(nb_timesteps, nz)``
     - Gas temperature [K].
   * - ``'dust_temperature'``
     - ``(nb_timesteps, nz)``
     - Dust temperature [K] (first grain bin).
   * - ``'visual_extinction'``
     - ``(nb_timesteps, nz)``
     - Visual extinction :math:`A_V` [mag].
   * - ``'species'``
     - list of str
     - Ordered list of species names, matching axis 1 of ``'abundances'``.

.. code-block:: python

    # Time axis in years
    yr = 365.25 * 24 * 3600
    times_yr = pipe.chemistry[30]['time'] / yr
    print(f"Last output time: {times_yr[-1]:.2e} yr")

    # Gas density at the last timestep, all z-points
    nH = pipe.chemistry[30]['H_number_density'][-1]     # array, shape (nz,)


Working with the abundances DataArray (xarray)
-----------------------------------------------

The ``'abundances'`` entry is an :class:`xarray.DataArray` with three named dimensions:

- ``time``    — the output times in seconds (coordinate)
- ``species`` — the species names (coordinate)
- ``spatial`` — the vertical grid index (no coordinate; use the ``'z'`` column of ``1D_static.dat``)

This makes it straightforward to select by species name or timestep index without
having to track integer indices manually.

.. code-block:: python

    ab = pipe.chemistry[30]['abundances']    # DataArray (nb_timesteps, nb_species, nz)

    # Select a species by name — all timesteps and z-points
    co = ab.sel(species='CO')               # DataArray (nb_timesteps, nz)

    # Select the last timestep
    co_last = ab.isel(time=-1).sel(species='CO')   # DataArray (nz,)

    # Select a specific timestep by index
    co_t48  = ab.sel(species='CO').isel(time=48)

    # Convert to a plain NumPy array
    co_arr = co_last.values                 # shape (nz,)

    # Compute number density [cm⁻³]
    nH     = pipe.chemistry[30]['H_number_density'][-1]   # (nz,)
    n_co   = nH * co_arr

    # List all available species
    print(pipe.chemistry[30]['species'])


Plotting a vertical abundance profile
---------------------------------------

.. code-block:: python

    import numpy as np
    import matplotlib.pyplot as plt
    import pandas as pd

    r = 30   # AU
    species = 'CO'

    # Load z-grid from the static file
    static = pd.read_table(
        f'chemistry/{r}AU/1D_static.dat',
        sep=r'\s+', comment='!', header=None, engine='python'
    )
    z = static[0].values   # z [AU], surface → midplane

    ab   = pipe.chemistry[r]['abundances'].isel(time=-1).sel(species=species).values
    nH   = pipe.chemistry[r]['H_number_density'][-1]
    n_X  = nH * ab

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
    ax1.semilogy(ab,  z);  ax1.set_xlabel(f'n({species})/n$_H$')
    ax2.semilogy(n_X, z);  ax2.set_xlabel(f'n({species}) [cm$^{{-3}}$]')
    for ax in (ax1, ax2):
        ax.set_ylabel('z [AU]')
    plt.tight_layout()
    plt.show()


Single-run mode (0-D / single column)
---------------------------------------

If no ``AU`` sub-folders are found, ``add_chemistry()`` reads the ``abundances.out`` file
directly from the chemistry path and stores the result in ``pipe.chemistry`` as a **single
dict** (same keys as above, but ``'abundances'`` has ``nz = 1`` for a 0-D run).

.. code-block:: python

    pipe.add_chemical_path('single_run/')
    pipe.add_chemistry()

    ab = pipe.chemistry['abundances']           # DataArray (nb_timesteps, nb_species, 1)
    co = ab.sel(species='CO').isel(spatial=0)   # time series for CO


.. _monitoring:

Monitoring a running job
--------------------------

The function :func:`~astromugs.nautilus.read.progress` reads only the 8-byte time record
from each timestep in ``abundances.out`` and seeks past the large abundance and physical-state
blocks. It is therefore essentially instantaneous regardless of file size, and safe to call
while Nautilus is still writing to the file.

.. code-block:: python

    from astromugs import nautilus

    # Simple call — just the folder path
    n, t = nautilus.read.progress('chemistry/30AU/')
    # → Timesteps done: 45   last time: 3.162e+05 yr

    # With percentage — reads nb_outputs from parameters.in
    n, t = nautilus.read.progress('chemistry/30AU/', parameters_path='chemistry/30AU/')
    # → Timesteps done: 45 / 64  (70.3 %)   last time: 3.162e+05 yr

    # Loop over all radii to get an overview
    import os, re
    chempath = 'chemistry/'
    radii = sorted([int(d.replace('AU','')) for d in os.listdir(chempath)
                    if re.match(r'^\d+AU$', d)])
    for r in radii:
        nautilus.read.progress(f'{chempath}{r}AU/', parameters_path=f'{chempath}{r}AU/')
