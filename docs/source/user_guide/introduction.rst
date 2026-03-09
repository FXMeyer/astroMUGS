Introduction
************

What is astroMUGS?
=================

AstroMUGS is a software package designed to couple the radiative transfer code RADMC-3D with the astrochemical code NAUTILUS. It provides a flexible framework for building disk and envelope models that can be used consistently in both radiative transfer and chemical simulations.

AstroMUGS can generate disk and envelope structures from scratch for RADMC-3D, or it can import user-defined disk surface density profiles. It allows the user to run the RADMC-3D Monte Carlo routines such as mctherm and mcmono, and to convert the resulting physical structure into a ready-to-use input model for NAUTILUS.

In addition, astroMUGS provides tools to convert NAUTILUS chemical outputs back into RADMC-3D-compatible files, enabling the computation of synthetic observations such as molecular line emission maps through the image module.

As suggested by its name (MUGS – MUlti-Grain Simulations), astroMUGS is specifically designed to handle models that include multiple dust populations. This capability makes it particularly useful for studies involving dust evolution, radiative transfer, and astrochemistry in protoplanetary disks and star-forming environments.

AstroMUGS therefore provides a streamlined and user-friendly workflow for building complete chemistry post-processing pipelines for astrophysical simulations, allowing complex modeling setups to be generated and executed with only a few lines of code.

Note that there is a Graphical User Interface version of astroMUGS that is uses a React flow and fastAPI structure. See astroMUGS-ui for more information. 

Copyright
=========

Copyright (c) 2025, Sacha Gavino

astroMUGS is free software: you can redistribute it and/or modify it under
the terms of the GPL License.

See the LICENSE file distributed with this software for more information.

Citation
========

If you use astroMUGS in your research, please cite:

Gavino et al. (2021), astroMUGS: Multi-Grain Simulations.