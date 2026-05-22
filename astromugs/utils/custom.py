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

    Two file formats are supported (space-separated, one header line):

    Without temperature:
        col 0  : radii [AU]
        col 1  : gas surface density [g.cm-2]
        col 2+ : dust surface density per size bin [g.cm-2], smallest to largest

    With temperature (header contains 'Tgas'):
        col 0  : radii [AU]
        col 1  : gas surface density [g.cm-2]
        col 2  : gas midplane temperature [K]
        col 3+ : dust surface density per size bin [g.cm-2], smallest to largest

    Returns
    -------
    radii : np.ndarray, shape (n_radii,)
        Radial grid in AU.
    sigma_dust : np.ndarray, shape (n_dust, n_radii)
        Dust surface densities for each size bin and radius.
    sigma_gas : np.ndarray, shape (n_radii,)
        Gas surface density at each radius.
    T_gas : np.ndarray or None, shape (n_radii,)
        Gas midplane temperature in K, or None if not present in the file.
    """
    with open(filepath, 'r') as f:
        header = f.readline().split()
    has_temp = 'Tgas' in header

    data = np.loadtxt(filepath, skiprows=1)
    radii     = data[:, 0]
    sigma_gas = data[:, 1]

    if has_temp:
        T_gas      = data[:, 2]
        sigma_dust = data[:, 3:].T
    else:
        T_gas      = None
        sigma_dust = data[:, 2:].T

    return radii, sigma_dust, sigma_gas, T_gas