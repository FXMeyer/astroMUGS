The Disk model
**************

The ``Disk`` class computes the physical structure of a protoplanetary disk:
surface densities, scale heights, temperatures, densities, and dust settling.

Overview
========

A ``Disk`` is built from a set of parameters (``params.disk``) and a dust model.
It provides methods to compute:

- **Gas structure**: surface density, scale height, midplane/atmosphere temperatures, vertical density profile.
- **Dust structure**: dust surface density per grain size (with settling), dust scale heights, dust number densities.
- **Optical properties**: extinction efficiency, visual extinction :math:`A_V(z)`.

The disk model supports two surface density modes:

- **Parametric** (``sigma_compute='param'``): power-law with exponential taper.
- **Custom** (``sigma_compute='custom'``): interpolated from a user-provided file (e.g., from an MHD simulation).


Key concepts
============

Surface density and settling
-----------------------------

The gas surface density follows:

.. math::

   \Sigma_g(r) = \Sigma_0 \left(\frac{r}{r_0}\right)^{-p} \exp\left[-\left(\frac{r}{r_c}\right)^{2-p}\right]

Dust settling reduces the scale height of larger grains relative to the gas:

.. math::

   H_d(a, r) = H_g(r) \cdot \min\left(1,\; \sqrt{\frac{\alpha}{\mathrm{St}(a, r)}}\right)

where :math:`\mathrm{St}` is the Stokes number and :math:`\alpha` is the turbulence parameter.


Vertical temperature structure
-------------------------------

The temperature transitions from midplane (:math:`T_\mathrm{mid}`) to atmosphere (:math:`T_\mathrm{atm}`)
following a smooth profile controlled by the parameter ``delta``.


API Reference
=============

.. autoclass:: astromugs.modeling.Disk.Disk
   :members:
   :undoc-members:
   :show-inheritance:
