Gas-grain simulations
************

Write NAUTILUS input files
=================


before writing ``your_model.write_chemistry()``, It is important that the user checks and update the parameters (your_model.pararms). Why?
Because the codes needs to know if it deals with a custom physical model (``your_model.params.disk.sigma_compute = 'custom'``) or not.
In case of a custom model, the gas surface density, which is most likely not used/created as RADMC3D input, must be implemented to create NMGC's input files, in. paticular ``1D_static.dat``.
Another required step, is to add the dust model to the grid, if it hasn't been done before (e.g. ``d = dust.CustomDistrib(rho_m=3.3, path="convert/")``, and ``your_model.grid.add_dust(d)``).
This is necessary, because Nautilus considers dust number densities and therefore needs to know which dust grain sizes you consider.  
