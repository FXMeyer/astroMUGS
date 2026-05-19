# Research: Dust Evolution & Planetesimal Formation Simulations — State of the Art

## Summary

Dust evolution in protoplanetary disks is now modeled with sophisticated 1D codes (DustPy, two-population models) that track growth, fragmentation, and radial drift, providing pebble fluxes as input for planet formation. The streaming instability (SI) is the leading mechanism for planetesimal formation, with recent 3D simulations (2024–2025) refining triggering conditions and the initial mass function (IMF). A key remaining gap is the self-consistent connection between dust evolution outputs and N-body initial conditions for subsequent planet growth.

## Findings

### 1. Dust Growth and Evolution Codes

1. **DustPy (Stammler & Birnstiel 2022)** — A Python package that solves the full Smoluchowski coagulation equation coupled to 1D radial gas/dust transport. It handles viscous gas evolution, dust advection and diffusion, growth, fragmentation, and vertical settling (via effective scale heights). Uses implicit 1st-order Euler integration. Open-source (GPLv3), modular via the Simframe framework. [Stammler & Birnstiel 2022, ApJ 935, 35](https://arxiv.org/abs/2207.00322)

2. **Two-population model (Birnstiel, Klahr & Ercolano 2012)** — A simplified analytical/semi-analytical approach that tracks only the upper end of the size distribution and total dust surface density. Identifies two regimes: fragmentation-limited (inner disk, $a_\mathrm{max} \propto v_\mathrm{frag}^2 / \alpha c_s^2$) and drift-limited (outer disk, $a_\mathrm{max} \propto \Sigma_d / \eta$). Reproduces full coagulation results well and is computationally cheap — widely used in planet formation population synthesis. Code: `two-pop-py`. [Birnstiel et al. 2012, A&A 539, A148](https://www.aanda.org/articles/aa/abs/2012/03/aa18136-11/aa18136-11.html)

3. **TriPoD (Pfeil, Birnstiel & Klahr 2024)** — A tri-population extension for coupling dust coagulation to 2D hydrodynamic simulations (vertically integrated). Adds a third population to better capture the mass distribution peak. [Pfeil et al. 2024, A&A 691, A45](https://ui.adsabs.harvard.edu/abs/2024A%26A...691A..45P)

4. **LA-COMPASS (Li et al.)** — A 2D (r,φ) hydrodynamic code that includes full dust coagulation coupled to disk dynamics. Used to study dust evolution in planet-carved gaps and vortices. Demonstrates that coagulation significantly modifies dust ring morphology and trapping efficiency near planets. [Li et al. 2019, ApJ 878, 39](http://obelix.rice.edu/~ai14/publications/2019_journal/Li_2019_ApJ_878_39.pdf); [Chen & Li 2020](https://www.osti.gov/servlets/purl/1764233)

### 2. Pebble Flux and Pebble Accretion

5. **Pebble accretion mechanism** — Lambrechts & Johansen (2012) showed that pebbles entering a core's gravitational reach lose kinetic energy via gas drag, yielding accretion rates orders of magnitude higher than classical planetesimal accretion. The accretion rate scales as $\dot{M} \propto \epsilon_\mathrm{peb} F_\mathrm{peb}$, where $\epsilon_\mathrm{peb}$ is the accretion efficiency and $F_\mathrm{peb}$ is the radial pebble mass flux. [Lambrechts & Johansen 2012, A&A 544, A32](https://www.aanda.org/articles/aa/full_html/2012/08/aa19127-12/aa19127-12.html)

6. **Pebble flux from dust evolution** — Lambrechts & Johansen (2014) derived an analytical pebble flux model based on a "growth front" moving outward: $F_\mathrm{peb} = 2\pi r_g (dr_g/dt) \Sigma_{d,0}(r_g)$. This gives a time-decaying flux (~$t^{-1/3}$ to $t^{-2/3}$). Drążkowska, Stammler & Birnstiel (2021) showed that coupling detailed DustPy-like evolution to pebble accretion provides more realistic, spatially varying Stokes numbers and fluxes — fragmentation can actually *benefit* growth by keeping pebbles small enough to avoid rapid drift loss. [Lambrechts & Johansen 2014, A&A 572, A107](https://www.aanda.org/articles/aa/full_html/2014/12/aa24343-14/aa24343-14.html); [Drążkowska et al. 2021, A&A 647, A15](https://www.aanda.org/articles/aa/full_html/2021/03/aa39925-20/aa39925-20.html)

7. **Pebble flux determines planetary system architecture** — Lambrechts et al. (2019) showed that high pebble fluxes produce super-Earths (efficient growth → migration → compact systems) while low fluxes yield terrestrial-like systems. The radial pebble flux is the key bifurcation parameter. [Lambrechts et al. 2019, A&A 627, A83](http://ui.adsabs.harvard.edu/abs/2019A&A...627A..83L/abstract)

### 3. Streaming Instability and Planetesimal Formation

8. **SI triggering conditions** — The streaming instability requires enhanced midplane dust-to-gas ratios (ρ_peb ≈ ρ_gas). Key parameters: dimensionless stopping time τ_s and local dust-to-gas ratio Z. Classical threshold: Z ≥ 0.02 for τ_s ~ 0.01–0.1. Lim et al. (2024) showed that external turbulence (α_D ~ 10⁻⁴–10⁻³) significantly raises the required Z threshold, potentially making SI harder to trigger in turbulent disks. [Lim et al. 2024, ApJ 969, 130](https://iopscience.iop.org/article/10.3847/1538-4357/ad47a2)

9. **SI in 3D with small grains** — Lim et al. (2025a, ApJ 981, 160) probed SI with small τ_s ≤ 0.01 and low Z, finding that strong clumping is suppressed for very small particles — the minimum τ_s for effective concentration is ~0.01 in 3D. Lim et al. (2025b, arXiv:2509.18270) extended to full 3D, confirming that axisymmetric 2D results overestimate clumping; 3D simulations show filamentary rather than ring-like structures. [Lim et al. 2025, ApJ 981, 160](https://iopscience.iop.org/article/10.3847/1538-4357/adb311); [Lim et al. 2025, arXiv:2509.18270](https://arxiv.org/abs/2509.18270)

10. **Large-domain SI simulations and the planetesimal IMF** — Schäfer et al. (2024) ran SI simulations in boxes up to 6.4H per side (32× larger than typical), forming ~4000 planetesimals. Key results: (a) filaments have azimuthal extent limited to ~1H (not axisymmetric rings), (b) the high-mass end of the IMF shows steep exponential tapering, (c) intermediate masses follow $dN/dM \propto M^{-p}$ with p ≈ 1.6. [Schäfer et al. 2024, arXiv:2410.08347](https://arxiv.org/abs/2410.08347)

11. **Universality of the planetesimal IMF** — Simon et al. (2016, 2017) found that the IMF is approximately universal (independent of τ_s) with a characteristic mass peaked around ~100 km radius bodies and described by: $N_{\geq}(m)/N_\mathrm{tot} = (m/m_\mathrm{min})^{-p} \exp[(m_\mathrm{min}/m_p)^q - (m/m_p)^q]$ with p ≈ 0.6, q ≈ 0.4. The mass budget is dominated by bodies near the characteristic mass $m_p$. [Simon et al. 2016, ApJ 822, 55](https://iopscience.iop.org/article/10.3847/0004-637X/822/1/55); [Simon et al. 2017, ApJL 847, L12](https://iopscience.iop.org/article/10.3847/2041-8213/aa8c79)

12. **SI in infall-dominated young disks** — Magnin et al. (2025, A&A 696, A162) showed that during the infall phase (Class 0/I), conditions for SI can be met earlier and more easily due to high mass loading, potentially forming the first generation of planetesimals within ~0.1–0.5 Myr. [Magnin et al. 2025, A&A 696, A162](https://www.aanda.org/articles/aa/abs/2025/04/aa52689-24/aa52689-24.html)

### 4. Connecting Dust Evolution to N-body Initial Conditions

13. **Planetesimal formation locations** — Dust evolution + SI models predict planetesimal formation is localized: at ice lines (water, silicate), dead zone edges, pressure bumps, and planet gap edges. Drążkowska & Alibert (2017), Schoonenberg & Ormel (2017), and Lenz et al. (2019) show that pile-up of drifting pebbles in the inner disk (0.3–3 AU) creates steep planetesimal surface density profiles, much steeper than MMSN. [Drążkowska & Alibert 2017; Lenz et al. 2019, ApJ 874, 36](https://www.aanda.org/articles/aa/full_html/2016/10/aa28983-16/aa28983-16.html)

14. **From SI to pebble accretion: growth in planetesimal rings** — Kaufmann et al. (2025, A&A 696, A65) simulated growth of embryos from SI-formed planetesimal rings. Key findings: (a) embryos in inner disk (< 1 AU) grow rapidly to pebble accretion transition mass; (b) at >10 AU, growth from planetesimal accretion alone is negligible; (c) pebble accretion (even in Bondi regime) is essential and cannot be neglected; (d) diffusive widening of the ring is a major growth inhibitor; (e) timing of core formation is strongly dependent on stellar mass and semi-major axis. [Kaufmann et al. 2025, A&A 696, A65](https://www.aanda.org/articles/aa/full_html/2025/04/aa52428-24/aa52428-24.html)

15. **Gap between dust evolution and N-body** — Most planet formation models insert partially assembled cores (0.01 M⊕) as initial conditions without self-consistently modeling the dust-to-planetesimal pathway. The Kaufmann et al. (2025) grid provides transition timings (time for largest SI-formed body to reach efficient pebble accretion mass) as a function of location and stellar mass — a step toward bridging this gap. However, fully self-consistent chains (dust evolution → SI triggering → planetesimal IMF → N-body) remain absent.

### 5. Observational Constraints from ALMA

16. **ALMA substructures** — DSHARP and subsequent surveys reveal ubiquitous rings and gaps in mm-continuum emission. Birnstiel et al. (2018, DSHARP V) showed these are consistent with dust trapping in pressure bumps where particles grow to maximum sizes set by fragmentation. The sharpness of rings constrains α_turb/St ratios. [Birnstiel et al. 2018, ApJL 869, L45](https://ui.adsabs.harvard.edu/abs/2018ApJ...869L..45B)

17. **Birnstiel (2024) ARA&A review** — Comprehensive review of dust growth and evolution. Highlights: (a) ALMA spectral index measurements confirm mm-sized grains in disks; (b) polarization observations constrain maximum grain sizes; (c) dust evolution models reproduce observed disk sizes and fluxes; (d) remaining puzzles include why drift doesn't deplete disks faster than observed (suggesting ubiquitous substructure/trapping). [Birnstiel 2024, ARA&A 62, 157](https://ui.adsabs.harvard.edu/abs/2024ARA%26A..62..157B/abstract)

18. **Substructures as planetesimal factories** — Pressure bumps that create observed rings also trap pebbles, raising local Z above SI thresholds. Stammler et al. (2019), Lau et al. (2022, 2024) showed that planetesimal and planet formation can proceed rapidly in such bumps, with the first planet potentially triggering sequential formation of additional planets at downstream bumps. [Lau et al. 2024, A&A 688, A22; Stammler et al. 2019, ApJ 884, L5]

## Sources

### Kept
- Stammler & Birnstiel 2022 (ApJ 935, 35) — DustPy code paper, primary reference for dust evolution simulations
- Birnstiel, Klahr & Ercolano 2012 (A&A 539, A148) — Two-population model, widely used analytical framework
- Birnstiel 2024 (ARA&A 62, 157) — Authoritative recent review of entire field
- Lim et al. 2024 (ApJ 969, 130) — SI + turbulence conditions, state-of-art
- Schäfer et al. 2024 (arXiv:2410.08347) — Largest SI simulation domain to date, IMF statistics
- Simon et al. 2016/2017 — Foundational planetesimal IMF from SI
- Kaufmann et al. 2025 (A&A 696, A65) — SI-to-pebble-accretion bridge, directly relevant
- Lambrechts & Johansen 2012, 2014 — Foundational pebble accretion papers
- Drążkowska et al. 2021 (A&A 647, A15) — Dust fragmentation + pebble accretion coupling
- Birnstiel et al. 2018 (DSHARP V) — ALMA constraints on dust models

### Dropped
- Generic FARGO dust module references — older, less relevant to current state
- Johansen et al. 2021 (Science Advances, terrestrial planets) — tangential to simulation methodology
- Various early SI papers (pre-2016) — superseded by Simon et al. and Lim et al.

## Gaps

1. **No end-to-end self-consistent pipeline exists** from dust coagulation → SI triggering → planetesimal IMF → N-body evolution in a single simulation framework. Each step uses different codes with hand-off assumptions.

2. **3D SI with realistic size distributions** — Most SI simulations use single-τ_s particles. The effect of a realistic poly-disperse size distribution on clumping and IMF is poorly explored (Krapp et al. 2019 is an exception).

3. **Temporal evolution of planetesimal formation** — When and where planetesimals form over disk lifetime (sequential vs. burst formation) is model-dependent and not observationally constrained.

4. **Fragmentation velocity uncertainties** — The critical parameter v_frag (1–10 m/s for silicates vs. ices) remains poorly constrained experimentally for realistic aggregate compositions, yet determines whether growth is drift- or fragmentation-limited.

5. **Connection to observed exoplanet architectures** — While pebble flux models predict bifurcation between terrestrial and super-Earth systems, quantitative predictions depend sensitively on the dust-to-planetesimal conversion efficiency.

**Suggested next steps:** Look into Drążkowska et al. (2023, PPVII chapter) for the most comprehensive recent synthesis; investigate whether any groups have coupled DustPy outputs directly to SI-triggering criteria in global disk models.
