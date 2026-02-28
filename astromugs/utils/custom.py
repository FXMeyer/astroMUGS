from __future__ import annotations

import numpy as np
from dataclasses import fields as dataclass_fields
from typing import TYPE_CHECKING

from .. constants.constants import M_sun, Ggram, autocm

if TYPE_CHECKING:
    from ..utils.thermal import ControlParams  # only needed for annotation


def surfacedensities(filepath='surfacedensities.inp'):
    """
    Read dust and gas surface densities from a custom file.

    File format (space-separated, one header line):
        col 0  : radii [AU]
        col 1  : gas surface density [g.cm-2]
        col 2+ : dust surface density per size bin [g.cm-2], smallest to largest

    Returns
    -------
    sigma_dust : np.ndarray, shape (n_radii, n_dust)
        Dust surface densities for each radius and size bin.
    sigma_gas : np.ndarray, shape (n_radii,)
        Gas surface density at each radius.
    """
    data = np.loadtxt(filepath, skiprows=1)
    radii  = data[:, 0]
    sigma_gas  = data[:, 1]
    sigma_dust = data[:, 2:].T
    return radii, sigma_dust, sigma_gas

