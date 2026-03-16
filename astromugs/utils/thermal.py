from dataclasses import dataclass, fields, is_dataclass, field
import html
from typing import Literal, Any


# =======================================================
# Utility: scientific number formatting
# =======================================================
def fmt_value(v, sci_threshold=1e3):
    if isinstance(v, (float, int)):
        if v != 0 and (abs(v) >= sci_threshold or abs(v) < 1e10):
            return f"{v:.3e}"
    return repr(v)


# =======================================================
# HTML REPR FOR JUPYTER (collapsible panel)
# =======================================================

def _html_repr(obj):
    """Return a collapsible HTML panel with MathJax-enabled table."""
    if not is_dataclass(obj):
        return None

    cls_name = obj.__class__.__name__
    flds = fields(obj)

    rows = ""
    for f in flds:
        desc_raw = f.metadata.get("desc", "")
        # HTML escape text except math (MathJax handles $...$)
        desc = html.escape(desc_raw)
        desc = desc.replace("&dollar;", "$")  # Allow $ in math mode

        val = fmt_value(getattr(obj, f.name))
        typ = f.type.__name__ if hasattr(f.type, "__name__") else str(f.type)

        rows += f"""
        <tr>
            <td><b>{f.name}</b></td>
            <td>{typ}</td>
            <td><code>{val}</code></td>
            <td style="text-align:left; white-space:normal;">{desc}</td>
        </tr>
        """

    # note: MathJax automatically renders inline $...$ or block $$...$$
    html_block = f"""
    <details style="margin:8px 0; padding:6px; border:1px solid #ccc; border-radius:6px;">
      <summary style="font-size:16px; font-weight:bold; cursor:pointer;">{cls_name}</summary>

      <div style="padding:10px;">
        <table style="border-collapse: collapse; margin-top:10px;">
          <thead>
            <tr style="text-align:left; border-bottom:1px solid #888;">
              <th style="padding-right:20px;">Field</th>
              <th style="padding-right:20px;">Type</th>
              <th style="padding-right:20px;">Value</th>
              <th style="text-align:left;">Description</th>
            </tr>
          </thead>
          <tbody>
            {rows}
          </tbody>
        </table>
      </div>

    </details>
    """

    return html_block

# =======================================================
# TEXT FALLBACK (for terminal printing)
# =======================================================
def fancy_repr(obj) -> str:
    """Non-HTML pretty repr for terminal."""
    if not is_dataclass(obj):
        return repr(obj)

    cls_name = obj.__class__.__name__
    flds = fields(obj)

    max_name = max(len(f.name) for f in flds)
    max_type = max(len(f.type.__name__) if hasattr(f.type, "__name__") else len(str(f.type))
                   for f in flds)
    max_desc = max(len(f.metadata.get("desc", "")) for f in flds)

    VALUE_COL_WIDTH = 14

    lines = [f"{cls_name}:"]
    lines.append(
        f"{'Field'.ljust(max_name)}  {'Type'.ljust(max_type)}  "
        f"{'Value'.ljust(VALUE_COL_WIDTH)}  {'Description'.ljust(max_desc)}"
    )
    lines.append(
        f"{'-'*max_name}  {'-'*max_type}  {'-'*VALUE_COL_WIDTH}  {'-'*max_desc}"
    )

    for f in flds:
        name = f.name.ljust(max_name)
        typ = (
            f.type.__name__.ljust(max_type)
            if hasattr(f.type, "__name__")
            else str(f.type).ljust(max_type)
        )
        desc = f.metadata.get("desc", "").ljust(max_desc)
        value = fmt_value(getattr(obj, f.name)).ljust(VALUE_COL_WIDTH)
        lines.append(f"{name}  {typ}  {value}  {desc}")

    return "\n".join(lines)

# =======================================================
# PARAMETER DATACLASSES
# =======================================================

@dataclass
class ControlParams:
    incl_dust: float = field(default=None,
        metadata={'desc': 'Switch to include dust continuum emission (0=off, 1=on)'})
    incl_lines: float = field(default=None,
        metadata={'desc': 'Switch to include line emission (0=off, 1=on)'})
    incl_freefree: float = field(default=None,
        metadata={'desc': 'Switch to include free-free emission (0=off, 1=on)'})
    nphot_therm: float = field(default=None,
        metadata={'desc': 'Number of photon packages for thermal Monte Carlo'})
    nphot_scat: float = field(default=None,
        metadata={'desc': 'Number of photon packages for scattering Monte Carlo'})
    nphot_mono: float = field(default=None,
        metadata={'desc': 'Number of photon packages for monochromatic Monte Carlo'})
    nphot_spec: float = field(default=None,
        metadata={'desc': 'Number of photon packages for spectrum calculation'})
    iseed: int = field(default=None,
        metadata={'desc': 'Random number generator seed'})
    ifast: float = field(default=None,
        metadata={'desc': 'Fast method for thermal emission computation (1=on)'})
    enthres: float = field(default=None,
        metadata={'desc': 'Energy threshold for thermal Monte Carlo convergence'})
    itempdecoup: int = field(default=None,
        metadata={'desc': 'Decouple temperatures of different dust species (1=on)'})
    istar_sphere: int = field(default=None,
        metadata={'desc': 'Treat stars as finite spheres instead of point sources (1=on)'})
    ntemp: int = field(default=None,
        metadata={'desc': 'Number of temperature grid points for the dust emissivity table'})
    temp0: float = field(default=None,
        metadata={'desc': '[K] Minimum temperature in the dust emissivity table'})
    temp1: float = field(default=None,
        metadata={'desc': '[K] Maximum temperature in the dust emissivity table'})
    scattering_mode_max: int = field(default=None,
        metadata={'desc': 'Maximum scattering mode (0=no scattering, 1=isotropic, 2=anisotropic, 5=full polarization)'})
    rto_style: float = field(default=None,
        metadata={'desc': 'Output file format (1=ASCII, 2=F77 unformatted, 3=C binary)'})
    camera_tracemode: float = field(default=None,
        metadata={'desc': 'Ray tracing mode for the camera (1=full ray tracing, -1=emission only)'})
    camera_nrrefine: float = field(default=None,
        metadata={'desc': 'Number of sub-pixel refinements for anti-aliasing in camera images'})
    camera_refine_criterion: float = field(default=None,
        metadata={'desc': 'Intensity contrast criterion triggering sub-pixel refinement'})
    camera_incl_stars: float = field(default=None,
        metadata={'desc': 'Include direct stellar emission in camera images (1=on)'})
    camera_starsphere_nrpix: float = field(default=None,
        metadata={'desc': 'Number of pixels used to render stellar spheres'})
    camera_spher_cavity_relres: float = field(default=None,
        metadata={'desc': 'Relative resolution for the spherical cavity around the star in camera'})
    camera_localobs_projection: float = field(default=None,
        metadata={'desc': 'Projection type for local observer camera mode'})
    camera_min_dangle: float = field(default=None,
        metadata={'desc': 'Minimum angle step for camera ray integration [rad]'})
    camera_max_dangle: float = field(default=None,
        metadata={'desc': 'Maximum angle step for camera ray integration [rad]'})
    camera_min_dr: float = field(default=None,
        metadata={'desc': 'Minimum spatial step for camera ray integration [cm]'})
    camera_diagnostics_subpix: float = field(default=None,
        metadata={'desc': 'Enable sub-pixel diagnostics output for the camera (1=on)'})
    camera_secondorder: float = field(default=None,
        metadata={'desc': 'Use second-order integration scheme for camera ray tracing (1=on)'})
    camera_interpol_jnu: float = field(default=None,
        metadata={'desc': 'Interpolate emissivity instead of opacity along camera rays (1=on)'})
    camera_scatsrc_allfreq: float = field(default=None,
        metadata={'desc': 'Compute scattering source function at all frequencies (1=on)'})
    mc_weighted_photons: float = field(default=None,
        metadata={'desc': 'Use luminosity-weighted photon packages in Monte Carlo (1=on)'})
    optimized_motion: float = field(default=None,
        metadata={'desc': 'Use optimized photon propagation in optically thick regions (1=on)'})
    lines_mode: float = field(default=None,
        metadata={'desc': 'Line transfer mode (1=LTE, 3=LVG/Sobolev, etc.)'})
    lines_maxdoppler: float = field(default=None,
        metadata={'desc': 'Maximum Doppler shift in units of the local line width'})
    lines_partition_ntempint: float = field(default=None,
        metadata={'desc': 'Number of temperature intervals for the partition function table'})
    lines_partition_temp0: float = field(default=None,
        metadata={'desc': '[K] Minimum temperature for the partition function table'})
    lines_partition_temp1: float = field(default=None,
        metadata={'desc': '[K] Maximum temperature for the partition function table'})
    lines_show_pictograms: float = field(default=None,
        metadata={'desc': 'Show ASCII pictograms of line transitions in the terminal (1=on)'})
    tgas_eq_tdust: int = field(default=None,
        metadata={'desc': 'Set gas temperature equal to dust temperature (1=on)'})
    subbox_nx: float = field(default=None,
        metadata={'desc': 'Number of cells in x-direction for subbox data output'})
    subbox_ny: float = field(default=None,
        metadata={'desc': 'Number of cells in y-direction for subbox data output'})
    subbox_nz: float = field(default=None,
        metadata={'desc': 'Number of cells in z-direction for subbox data output'})
    subbox_x0: float = field(default=None,
        metadata={'desc': 'Minimum x coordinate of the subbox [cm]'})
    subbox_x1: float = field(default=None,
        metadata={'desc': 'Maximum x coordinate of the subbox [cm]'})
    subbox_y0: float = field(default=None,
        metadata={'desc': 'Minimum y coordinate of the subbox [cm]'})
    subbox_y1: float = field(default=None,
        metadata={'desc': 'Maximum y coordinate of the subbox [cm]'})
    subbox_z0: float = field(default=None,
        metadata={'desc': 'Minimum z coordinate of the subbox [cm]'})
    subbox_z1: float = field(default=None,
        metadata={'desc': 'Maximum z coordinate of the subbox [cm]'})
    modified_random_walk: float = field(default=None,
        metadata={'desc': 'Enable Modified Random Walk for optically thick cells (1=on)'})
    mrw_gamma: float = field(default=None,
        metadata={'desc': 'Gamma factor controlling the MRW step size in optically thick regions'})
    mrw_tauthres: float = field(default=None,
        metadata={'desc': 'Optical depth threshold above which the Modified Random Walk is triggered'})
    mrw_count_trigger: float = field(default=None,
        metadata={'desc': 'Photon count threshold triggering the Modified Random Walk algorithm'})
    setthreads: float = field(default=None,
        metadata={'desc': 'Number of OpenMP threads for parallel computation'})
    mc_scat_maxtauabs: float = field(default=None,
        metadata={'desc': 'Maximum absorption optical depth allowed per Monte Carlo scattering step'})
    writeimage_unformatted: float = field(default=None,
        metadata={'desc': 'Write image output files in unformatted binary format (1=on)'})

    def __repr__(self):  # terminal
        return fancy_repr(self)

    def _repr_html_(self):  # Jupyter
        return _html_repr(self)

@dataclass
class StarParams:
    mass: float = field(default=1, 
        metadata={'desc': r'[Msun] Stellar mass. Note that this value can be independent from the value used for the physical model.'})
    temperature: float = field(default=4000, 
        metadata={'desc': r'[K] Stellar surface temperature'})
    luminosity: float = field(default=4000, 
        metadata={'desc': r'[Lsun] Stellar luminosity'})
    x: float = field(default=0, 
        metadata={'desc': r'[AU] Stellar position in x'})
    y: float = field(default=0, 
        metadata={'desc': r'[AU] Stellar position in y'})
    z: float = field(default=0, 
        metadata={'desc': r'[AU] Stellar position in z'})


    def __repr__(self):
        # add a blank line before StarParams for readability
        return "\n" + fancy_repr(self)

    def _repr_html_(self):
        return "<br/>" + _html_repr(self)   # blank line in Jupyter
    
@dataclass
class WaveParams:
    lmin: float = field(default=9.12e-2, 
        metadata={'desc': r'[micron] Minimum wavelength for wavelength_micron.inp'})
    lmax: float = field(default=1e4, 
        metadata={'desc': r'[micron] Maximum wavelength for wavelength_micron.inp'})
    nlam: int = field(default=200,
            metadata={'desc': 'Number of wavelengths for wavelength_micron.inp'})
    lmin_mono: float = field(default=9.12e-2, 
        metadata={'desc': r'[micron] Minimum wavelength for mcmono_wavelength_micron.inp'})
    lmax_mono: float = field(default=1e4, 
        metadata={'desc': r'[micron] Maximum wavelength for mcmono_wavelength_micron.inp'})
    nlam_mono: int = field(default=200,
            metadata={'desc': 'Number of wavelengths for mcmono_wavelength_micron.inp'})


    def __repr__(self):
        # add a blank line before StarParams for readability
        return "\n" + fancy_repr(self)

    def _repr_html_(self):
        return "<br/>" + _html_repr(self)   # blank line in Jupyter
    

# ---------------- thermal Params ----------------
@dataclass
class ThermalParams:
    control: ControlParams = field(default_factory=ControlParams)
    star: StarParams = field(default_factory=StarParams)
    wave: WaveParams = field(default_factory=WaveParams)

    def __repr__(self):
        # blank line between dataclasses
        return fancy_repr(self.control) + "\n\n" + fancy_repr(self.star) + "\n\n" + fancy_repr(self.wave)

    def _repr_html_(self):
        # Jupyter collapsible panels with spacing
        return _html_repr(self.control) + "<br/><br/>" + _html_repr(self.star) + "<br/><br/>" + _html_repr(self.wave)