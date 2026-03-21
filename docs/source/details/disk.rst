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

Following the prescription by `Cuzzi et al. (1993) <https://ui.adsabs.harvard.edu/abs/1993Icar..106..102C/abstract>`_,
`Youdin & Lithwick. (2007) <https://ui.adsabs.harvard.edu/abs/2007Icar..192..588Y/abstract>`_,
`Dong et al. (2015) <https://ui.adsabs.harvard.edu/abs/2015ApJ...809...93D>`_, dust settling reduces the scale height of larger grains relative to the gas:

.. math::

   H_d(a, r) = H_g(r) \frac{1}{\sqrt{1 + T_{s,mid} \frac{S_c}{\alpha}}}

where :math:`\mathrm{T_{s,mid}}` is the dimensionless stopping time in the midplane,  :math:`\alpha` is the turbulence parameter, and :math:`S_c` is the Schmidt number.

In params.disk, :math:`\alpha` and :math:`\S_c` correspond to ``alpha`` and ``schmidtnumber``, respectively.

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
