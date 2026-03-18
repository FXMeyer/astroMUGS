Radiative transfer simulations
************

Write dust continuum radiative transfer RADMC3D input files
=================


you may ask why stellar mass is defined twice (in struct and thermal.py). That is because it happens, sometimes, that the user needs to define the stellar mass differently when building the physical structure and for the thermal computation.