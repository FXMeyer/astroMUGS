"""
Tests for astromugs core components.

These tests verify the fundamental building blocks of the pipeline:
data structures, physical models, grid construction, and dust distributions.
They are designed to run without external simulation data (RADMC3D, NAUTILUS).
"""

import numpy as np
import pytest

from astromugs.constants.constants import (
    autocm, M_sun, L_sun, R_sun, kb, h_erg, c, pi, amu, mu, black_body,
)
from astromugs.modeling.Star import Star
from astromugs.modeling.Disk import Disk
from astromugs.dust.MRNDistrib import MRNDistrib
from astromugs.utils.struct import DiskParams, EnvelopeParams, StructureParams
from astromugs.utils.thermal import WaveParams
from astromugs.pipeline.Grid import Grid


# ============================================================
# Constants
# ============================================================

class TestConstants:
    """Verify that physical constants are consistent and reasonable."""

    def test_au_conversion(self):
        """1 AU in cm should be ~1.496e13."""
        assert pytest.approx(autocm, rel=1e-3) == 1.496e13

    def test_solar_mass(self):
        """Solar mass should be ~1.989e33 g."""
        assert pytest.approx(M_sun, rel=1e-2) == 1.989e33

    def test_solar_luminosity(self):
        """Solar luminosity should be ~3.9e33 erg/s."""
        assert pytest.approx(L_sun, rel=1e-1) == 3.828e33

    def test_boltzmann_cgs(self):
        """Boltzmann constant in CGS should be ~1.381e-16 erg/K."""
        assert pytest.approx(kb, rel=1e-2) == 1.381e-16

    def test_blackbody_peak(self):
        """Blackbody at 5780 K should peak near visible wavelengths (Wien's law)."""
        wavelengths = np.linspace(0.1, 5.0, 500)
        intensities = np.array([black_body(5780, lam) for lam in wavelengths])
        peak_wavelength = wavelengths[np.argmax(intensities)]
        # Astropy BlackBody returns B_nu evaluated at given wavelength,
        # so the peak in this representation is shifted relative to Wien's B_lambda law
        assert 0.3 < peak_wavelength < 1.2, (
            f"Blackbody peak at {peak_wavelength} µm, expected in visible/near-IR range"
        )


# ============================================================
# Star model
# ============================================================

class TestStar:
    """Test the Star class initialization and derived quantities."""

    def test_default_star(self):
        """Default star should have reasonable properties."""
        star = Star()
        assert star.mass == 0.5
        assert star.luminosity == 1.0
        assert star.temperature == 4000.0

    def test_solar_star(self):
        """A solar-type star should have radius ~1 R_sun."""
        star = Star(mass=1.0, luminosity=1.0, temperature=5780.0)
        assert pytest.approx(star.radius, rel=0.01) == 1.0

    def test_radius_scales_with_luminosity(self):
        """Higher luminosity at same temperature should give larger radius."""
        star_dim = Star(luminosity=1.0, temperature=5000.0)
        star_bright = Star(luminosity=4.0, temperature=5000.0)
        assert star_bright.radius > star_dim.radius

    def test_radius_scales_with_temperature(self):
        """Higher temperature at same luminosity should give smaller radius."""
        star_cool = Star(luminosity=1.0, temperature=3000.0)
        star_hot = Star(luminosity=1.0, temperature=6000.0)
        assert star_hot.radius < star_cool.radius

    def test_custom_position(self):
        """Star position should be settable."""
        star = Star(x=1.0, y=2.0, z=3.0)
        assert (star.x, star.y, star.z) == (1.0, 2.0, 3.0)


# ============================================================
# DiskParams and dataclass structures
# ============================================================

class TestDiskParams:
    """Test parameter dataclasses for correct defaults and mutability."""

    def test_default_values(self):
        """DiskParams should initialize with physical defaults."""
        dp = DiskParams()
        assert dp.rin == 1.0
        assert dp.rout == 300.0
        assert dp.dtogas == 0.01
        assert dp.nr == 101
        assert dp.ntheta == 181
        assert dp.coordsystem == "spherical"

    def test_custom_values(self):
        """DiskParams should accept custom values."""
        dp = DiskParams(rin=5.0, rout=500.0, nr=201)
        assert dp.rin == 5.0
        assert dp.rout == 500.0
        assert dp.nr == 201

    def test_envelope_defaults(self):
        """EnvelopeParams should have reasonable defaults."""
        ep = EnvelopeParams()
        assert ep.rmin < ep.rmax
        assert ep.dtogas > 0

    def test_structure_params_contains_both(self):
        """StructureParams should bundle disk and envelope."""
        sp = StructureParams()
        assert isinstance(sp.disk, DiskParams)
        assert isinstance(sp.envelope, EnvelopeParams)

    def test_wave_params_defaults(self):
        """WaveParams should have min < max and positive nlam."""
        wp = WaveParams()
        assert wp.lmin < wp.lmax
        assert wp.nlam > 0


# ============================================================
# MRN dust distribution
# ============================================================

class TestMRNDistrib:
    """Test the MRN grain size distribution model."""

    def test_single_grain(self):
        """Single grain mode should return one size."""
        mrn = MRNDistrib(nb_sizes=1, rsingle=0.1)
        sizes = mrn.sizes()
        assert sizes.shape == (1, 1)
        assert sizes[0, 0] == 0.1

    def test_multiple_grain_sizes(self):
        """Multiple grain sizes should span amin to amax."""
        mrn = MRNDistrib(nb_sizes=10, amin=0.005, amax=1000.0)
        sizes = mrn.sizes()
        assert sizes.shape[0] == 3  # (amin, amax, average) per bin
        assert sizes.shape[1] == 10
        # First bin starts at amin, last bin ends at amax
        assert pytest.approx(sizes[0, 0], rel=1e-3) == 0.005
        assert pytest.approx(sizes[1, -1], rel=1e-3) == 1000.0

    def test_average_between_edges(self):
        """Average grain size in each bin should lie between its edges."""
        mrn = MRNDistrib(nb_sizes=5, amin=0.01, amax=100.0)
        sizes = mrn.sizes()
        for i in range(5):
            assert sizes[0, i] < sizes[2, i] < sizes[1, i]

    def test_mass_fractions_sum_to_one(self):
        """Mass fractions across all bins should sum to 1."""
        mrn = MRNDistrib(nb_sizes=10, amin=0.005, amax=1000.0)
        fractions = mrn.massfraction()
        assert pytest.approx(np.sum(fractions), rel=1e-6) == 1.0

    def test_single_grain_fraction_is_one(self):
        """Single grain mode should have mass fraction = 1."""
        mrn = MRNDistrib(nb_sizes=1)
        fractions = mrn.massfraction()
        assert fractions[0] == 1

    def test_grain_mass_positive(self):
        """Grain mass should be positive."""
        mrn = MRNDistrib(nb_sizes=5, amin=0.01, amax=100.0)
        masses = mrn.grainmass()
        assert np.all(masses > 0)

    def test_grain_mass_increases_with_size(self):
        """Larger grains should be more massive."""
        mrn = MRNDistrib(nb_sizes=5, amin=0.01, amax=100.0)
        masses = mrn.grainmass()
        # Masses correspond to average sizes, which increase across bins
        for i in range(len(masses) - 1):
            assert masses[i] < masses[i + 1]


# ============================================================
# Grid construction
# ============================================================

class TestGrid:
    """Test spatial and wavelength grid construction."""

    @pytest.fixture
    def default_grid(self):
        """Create a Grid with default parameters."""
        dp = DiskParams()
        wp = WaveParams()
        return Grid(params=dp, wave=wp)

    def test_spherical_grid_shapes(self, default_grid):
        """Spherical grid should produce correct number of edges and centres."""
        g = default_grid
        g.set_spherical_grid()
        # nr edges = DiskParams.nr = 101, so nr-1 = 100 centres
        assert len(g.r_edge) == 101
        assert len(g.r) == 100
        assert len(g.theta_edge) == 181
        assert len(g.theta) == 180

    def test_spherical_grid_radial_bounds(self, default_grid):
        """Radial edges should span from rin to rout."""
        g = default_grid
        g.set_spherical_grid()
        assert pytest.approx(g.r_edge[0], rel=1e-3) == 1.0   # rin
        assert pytest.approx(g.r_edge[-1], rel=1e-3) == 300.0  # rout

    def test_spherical_grid_theta_bounds(self, default_grid):
        """Theta edges should span 0 to pi."""
        g = default_grid
        g.set_spherical_grid()
        assert pytest.approx(g.theta_edge[0]) == 0.0
        assert pytest.approx(g.theta_edge[-1]) == np.pi

    def test_spherical_grid_log_spacing(self, default_grid):
        """Log spacing should make inner cells smaller than outer cells."""
        g = default_grid
        g.set_spherical_grid(log=True)
        dr_inner = g.r_edge[1] - g.r_edge[0]
        dr_outer = g.r_edge[-1] - g.r_edge[-2]
        assert dr_inner < dr_outer

    def test_cartesian_grid(self, default_grid):
        """Cartesian grid should return edges and centres with correct shapes."""
        g = default_grid
        edges, centres = g.set_cartesian_grid(xmin=-100, xmax=100, nx=50)
        assert edges.shape == (3, 50)
        assert centres.shape == (3, 49)

    def test_wavelength_grid(self, default_grid):
        """Wavelength grid should span lmin to lmax with correct number of points."""
        g = default_grid
        g.set_wavelength_grid()
        assert len(g.lam) == 200  # default nlam
        assert g.lam[0] > 0
        assert g.lam[0] < g.lam[-1]

    def test_wavelength_grid_log_spacing(self, default_grid):
        """Log-spaced wavelength grid should have smaller steps at short wavelengths."""
        g = default_grid
        g.set_wavelength_grid(log=True)
        dl_short = g.lam[1] - g.lam[0]
        dl_long = g.lam[-1] - g.lam[-2]
        assert dl_short < dl_long

    def test_add_star(self, default_grid):
        """Adding a star should append to the stars list."""
        g = default_grid
        star = Star()
        g.add_star(star)
        assert len(g.stars) == 1
        assert g.stars[0].mass == 0.5

    def test_add_density(self, default_grid):
        """Adding density arrays should append correctly."""
        g = default_grid
        rho = np.ones((10, 20, 30))
        g.add_density(rho)
        g.add_dustdensity(rho * 0.01)
        assert len(g.density) == 1
        assert len(g.dustdensity) == 1
        assert g.density[0].shape == (10, 20, 30)

    def test_chemdisk_grid(self, default_grid):
        """Chemistry disk grid should produce correct vertical structure."""
        g = default_grid
        radii = [10.0, 50.0, 100.0, 200.0]
        g.set_chemdisk_grid(r=radii, max_H=4, nz_chem=64)
        assert len(g.rchem) == 4
        assert len(g.zchem) == 64
        # z should go from max_H down toward 0
        assert g.zchem[0] > g.zchem[-1]

    def test_chem_grid_uniform_zmax(self, default_grid):
        """Chemistry grid with uniform zmax should have correct shape."""
        g = default_grid
        radii = np.array([10.0, 50.0, 100.0])
        g.set_chem_grid(r=radii, zmax=50.0, nbcells=30)
        assert g.zchem.shape == (3, 30)
        # Should be flipped (decreasing z, as NAUTILUS expects)
        assert g.zchem[0, 0] > g.zchem[0, -1]

    def test_chem_grid_spherical_envelope(self, default_grid):
        """Chemistry grid with msize should follow spherical boundary."""
        g = default_grid
        radii = np.array([10.0, 50.0, 100.0, 200.0])
        g.set_chem_grid(r=radii, msize=250.0, nbcells=20)
        # Last radius is dropped (zmax=0 at r=msize)
        assert len(g.rchem) == 3
        assert g.zchem.shape == (3, 20)

    def test_chemmodel_storage(self, default_grid):
        """Chemical abundance grids should be stored by species name."""
        g = default_grid
        abundance = np.random.rand(10, 20)
        g.add_existingchemmodel(abundance, "CO")
        g.add_existingchemmodel(abundance * 0.1, "H2O")
        assert "CO" in g.chemmodel
        assert "H2O" in g.chemmodel
        assert g.chemmodel["CO"].shape == (10, 20)