Gas-grain simulations
************

Write NAUTILUS input files
=================


before writing ``your_model.write_chemistry()``, It is important that the user checks and update the parameters (your_model.pararms). Why?
Because the codes needs to know if it deals with a custom physical model (``your_model.params.disk.sigma_compute = 'custom'``) or not.
In case of a custom model, the gas surface density, which is most likely not used/created as RADMC3D input, must be implemented to create NMGC's input files, in. paticular ``1D_static.dat``.
Another required step, is to add the dust model to the grid, if it hasn't been done before (e.g. ``d = dust.CustomDistrib(rho_m=3.3, path="convert/")``, and ``your_model.grid.add_dust(d)``).
This is necessary, because Nautilus considers dust number densities and therefore needs to know which dust grain sizes you consider.  



Best practices
==============

Stay in the Nautilus safe zone
------------------------------

Nautilus can encounter convergence issues when placed in extreme physical conditions.
For example, a very strong UV field with little attenuation can cause the chemical computation
to fail or make the run extremely long.

Protoplanetary disk models are particularly prone to pushing Nautilus into stiff chemistry
regimes, especially when post-processing hydrodynamical simulations. Many such simulations
— whether computed from a single core collapse or from a 1-D surface density profile —
can produce extreme conditions in the disk upper atmosphere: very low gas densities, strong
UV fields, and negligible visual extinction. Nautilus is not designed for these environments,
which are better suited for PDR codes.

Here are examples where Nautilus will struggle to converge:

- Very low density (:math:`n_\mathrm{H} \lesssim 10^2\,\mathrm{cm}^{-3}`) combined with a strong UV field and low :math:`A_V`: the spread between fast photodissociation rates and slow two-body formation rates (which scale as :math:`n_\mathrm{H}^2`) creates an extremely stiff system.
- Near-balanced UV shielding: when the UV flux and :math:`A_V` partially cancel, small abundance changes shift the shielding, which shifts the photodissociation rates, which shifts the abundances — a feedback loop that prevents convergence.

To mitigate these issues, astroMUGS provides floor and ceiling parameters (``min_gas_density``, ``min_av``, ``max_uv``)
that clamp physical quantities to user-defined bounds. These can be passed as keyword arguments to ``write_nautilus()``.

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