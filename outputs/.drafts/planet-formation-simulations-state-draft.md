# State of the Art: Planet Formation Simulations
## N-body Dynamics, Dust Evolution, Chemistry — A Postdoc Interview Briefing

**Prepared for:** GREP postdoc interview (Centre for Planetary Habitability, University of Oslo)
**Date:** 2026-05-13

---

## Executive Summary

Planet formation modeling has entered a new era of multi-physics integration. Three historically separate simulation domains — **N-body dynamics**, **dust evolution**, and **disk chemistry** — are now being coupled to predict not only where and when planets form, but *what they are made of*. This briefing synthesizes the current state of each domain and their emerging connections, tailored to the GREP project's goal of predicting exoplanet architectures and compositions testable against PLATO observations.

**Key takeaways:**
- GPU-accelerated N-body codes (GENGA) now simulate from ~4000 planetesimals to planets; hybrid integrators (TRACE in REBOUND) have improved accuracy for close encounters.
- Population synthesis (Bern/NGPPS Generation III) is being statistically validated against RV surveys, identifying four architectural classes of planetary systems.
- Dust evolution is well-modeled by DustPy (full coagulation) and two-population approximations; the streaming instability (SI) is the leading planetesimal formation mechanism, with 3D simulations refining triggering conditions.
- Chemical composition tracking — led by the **chemcomp** framework — can now follow C/O ratios, heavy element content, and volatile-to-refractory ratios through disk evolution and planet accretion.
- The critical finding from Danti, Bitsch & Mah (2023, A&A 679, L7) — the user's reference paper — is that planetesimal formation *dramatically reduces* atmospheric heavy element enrichment compared to pure pebble accretion, and distinguishing formation pathways requires measuring both C/H and O/H, not just C/O.
- **The major gap** is the absence of end-to-end frameworks coupling N-body dynamics with detailed chemistry. This is precisely the frontier the GREP position targets.

---

## 1. N-body Simulations of Planet Formation

### 1.1 Leading Codes

| Code | Key feature | Parallelism | Reference |
|------|-------------|-------------|-----------|
| **REBOUND** | Modular, multiple integrators (WHFast, IAS15, MERCURIUS, TRACE) | CPU/OpenMP | Rein & Liu (2012) |
| **GENGA** | GPU-accelerated hybrid symplectic, handles ~10⁴ bodies | CUDA/GPU | Grimm & Stadel (2014); Grimm et al. (2022) |
| **Mercury6** | Classic hybrid symplectic + Bulirsch-Stoer | CPU | Chambers (1999) |
| **SyMBA/SyMBAp** | Recursive subdivision for encounters; recently parallelized | CPU/OpenMP | Duncan et al. (1998) |

### 1.2 Recent Methodological Advances

**TRACE integrator** (Lu, Hernandez & Rein 2024, MNRAS 533, 3708): A near-time-reversible hybrid integrator for REBOUND that switches between WHFast and high-accuracy methods during close encounters. Outperforms MERCURIUS in accuracy for planet-planet scattering. Extended in 2025 to handle collisional fragmentation — enabling large-N disk simulations with realistic collision outcomes.

**GENGA II** (Grimm et al. 2022, ApJ 932, 124): Extended to non-Newtonian forces (GR, oblateness, tides, gas drag, type I/II migration) and high particle counts (>10⁴). Enables direct simulation from full planetesimal populations to planets on GPU hardware in days–weeks.

**Coupling to disk models:** Modern N-body codes routinely include pebble accretion prescriptions (Lambrechts et al. 2019), gas accretion, and migration with eccentricity/inclination damping. The disk is typically modeled as a 1D evolving α-disk with prescribed surface density, temperature, and accretion rate profiles.

### 1.3 Key Recent Results

**From planetesimals to planets** (Lorek & Lambrechts 2026, A&A 707, A244): GPU N-body simulations (GENGA) starting from streaming-instability-inspired planetesimal distributions between 10–50 AU show that: (a) initial spatial distribution does *not* significantly affect outcomes — narrow rings scatter within ~1 Myr; (b) systems consistently form 1–2 gas giants plus several icy planets >2 M⊕; (c) pebble flux self-regulates final planet numbers regardless of initial planetesimal mass.

**Four architectural classes** (Emsenhuber, Mordasini & Burn 2023, EPJ Plus 138, 181): The Bern population synthesis model identifies four classes: similar-mass tightly-packed (Kepler "peas-in-a-pod"), ordered, anti-ordered, and mixed architectures. Compact low-mass systems emerge naturally from the core accretion paradigm.

**Statistical validation** (Burn et al. 2025, A&A 701, A64): The NGPPS Generation III model is now compared quantitatively against HARPS/Coralie RV survey completeness. The nominal population reproduces observed mass-period distributions well, with systematic differences identified at specific regimes — providing clear targets for model improvement.

**Metallicity effects** (NGPPS VIII, 2025): The same framework predicts how occurrence rates, orbital periods, eccentricities, and the radius valley vary with host star metallicity — directly relevant to PLATO's ability to characterize stars.

### 1.4 The Initial Conditions Problem

Most N-body simulations still begin from pre-placed lunar-mass embryos (0.01 M⊕), skipping the planetesimal→embryo transition. GPU codes (GENGA) now enable starting from ~4000 planetesimals, but sensitivity tests show weak dependence on initial distributions in the outer disk, primarily because pebble accretion dominates early growth and self-regulates planet numbers (Lorek & Lambrechts 2026).

The streaming instability IMF — an exponentially-tapered power law with slope p ≈ 0.6 and characteristic mass corresponding to ~100 km bodies (Simon et al. 2016, 2017) — provides the standard initial condition. However, the IMF is resolution-dependent at the low-mass end and its universality across disk conditions is debated.

### 1.5 Disk Models Feeding N-body

**Standard approach:** Steady-state α-disk with irradiation-dominated outer regions (T ∝ r⁻³/⁷). Accretion rate decreases as a power law with time (Hartmann et al. 2016).

**Frontier: Disk winds.** MHD disk wind models (Tabone et al. 2022; Weder et al. 2023) are replacing simple viscous-α evolution. X-ray photoevaporation rates have been revised downward by ~10× due to enhanced cooling (Sellek et al. 2024). Direct coupling of wind-driven disk evolution into N-body simulations remains at an early stage — an active frontier.

---

## 2. Dust Evolution and Planetesimal Formation

### 2.1 Dust Growth Codes

**DustPy** (Stammler & Birnstiel 2022, ApJ 935, 35): The state-of-the-art open-source Python code that solves the full Smoluchowski coagulation equation coupled to 1D radial gas/dust transport. Handles growth, fragmentation, radial drift, and vertical settling. Modular via the Simframe framework.

**Two-population model** (Birnstiel, Klahr & Ercolano 2012, A&A 539, A148): A widely-used analytical approximation tracking only the peak grain size and total dust density. Identifies two growth regimes: fragmentation-limited (inner disk) and drift-limited (outer disk). Computationally cheap — used inside population synthesis (Bern model) and composition-tracking codes (chemcomp).

**TriPoD** (Pfeil, Birnstiel & Klahr 2024, A&A 691, A45): A tri-population extension enabling coupling of dust coagulation to 2D hydrodynamics.

### 2.2 Pebble Flux and Pebble Accretion

Pebble accretion (Lambrechts & Johansen 2012) yields growth rates orders of magnitude faster than planetesimal accretion for cores above the transition mass (~0.01 M⊕). The pebble mass flux F_peb is the critical parameter — it decays with time as the pebble production front sweeps outward through the disk.

**A key bifurcation result** (Lambrechts et al. 2019): High pebble fluxes produce super-Earths that migrate into compact chains; low fluxes yield terrestrial-like systems. The pebble flux is thus a primary architectural control parameter.

Drążkowska, Stammler & Birnstiel (2021) showed that coupling detailed DustPy-like evolution to pebble accretion provides more realistic, spatially varying Stokes numbers and fluxes — and that fragmentation can actually *benefit* growth by keeping pebbles small enough to avoid rapid drift loss.

### 2.3 Streaming Instability

The streaming instability concentrates pebbles into dense filaments that collapse gravitationally into planetesimals. It remains the leading planetesimal formation mechanism.

**Triggering conditions:** Requires enhanced midplane dust-to-gas ratios (Z ≥ ~0.02 for τ_s ~ 0.01–0.1). External turbulence significantly raises this threshold (Lim et al. 2024, ApJ 969, 130). 3D simulations show filamentary (not ring-like) structures and less efficient clumping than 2D axisymmetric runs overestimate (Lim et al. 2025).

**Planetesimal IMF:** Large-domain simulations (Schäfer et al. 2024, arXiv:2410.08347) forming ~4000 planetesimals confirm an intermediate-mass power law dN/dM ∝ M⁻¹·⁶ with steep exponential tapering at high masses. The mass budget is dominated by ~100 km bodies.

**Early formation in young disks:** Magnin et al. (2025, A&A 696, A162) show that during the infall phase (Class 0/I), conditions for SI can be met within ~0.1–0.5 Myr — potentially forming the first generation of planetesimals before the disk is fully assembled.

### 2.4 The Dust-to-N-body Gap

**This remains the weakest link in the chain.** No end-to-end, self-consistent pipeline exists from dust coagulation → SI triggering → planetesimal IMF → N-body planet formation. Each step uses different codes with manual hand-off. Kaufmann et al. (2025, A&A 696, A65) provide transition timings from SI-formed planetesimal rings to pebble accretion as a function of location and stellar mass — a step toward bridging this gap.

Key uncertainties: the fragmentation velocity v_frag (1–10 m/s for silicates vs ices) remains poorly constrained experimentally, yet determines whether growth is drift- or fragmentation-limited.

---

## 3. Chemical Composition in Planet Formation

### 3.1 The chemcomp Framework

**chemcomp** (Schneider & Bitsch 2021a,b, A&A 654, A71/A72) is the leading open-source code coupling disk chemistry to planet formation. It simulates:
- Viscous disk evolution
- Pebble growth (two-population model)
- Pebble drift and evaporation/condensation at ice lines
- Planet growth via pebble and gas accretion
- Per-species chemical tracking (C, O, N, S, Si, Fe, Mg, etc.)

The code solves advection-diffusion equations per chemical species, following how inward-drifting pebbles evaporate at their respective ice lines (H₂O ~150 K, CO₂ ~70 K, CO ~20 K), enriching the gas that planets subsequently accrete. The C/O ratio in the planetary atmosphere depends on which evaporation fronts the planet crosses during migration.

### 3.2 The Danti et al. (2023) Result — Reference Paper

Danti, Bitsch & Mah (2023, A&A 679, L7) used chemcomp to compare three formation scenarios:

1. **Pure pebble accretion** — highest heavy element enrichment
2. **Pebble accretion + planetesimal formation** (but no planetesimal accretion) — dramatic drop in enrichment because pebbles locked into planetesimals cannot drift inward and enrich disk gas
3. **Combined pebble + planetesimal accretion** — enrichment remains much lower than pure pebble case

**Critical findings:**
- C/O ratio alone is *insufficient* to distinguish formation pathways — it depends on too many parameters (viscosity, starting position, migration history).
- Distinguishing scenarios requires measuring **both** C/H and O/H alongside C/O and volatile/refractory ratios.
- Planets with sub-stellar atmospheric metallicity likely formed in the outer disk and were scattered inward; super-stellar metallicity planets are migration-driven.

This paper exemplifies the kind of chemistry-dynamics coupling the GREP position targets.

### 3.3 Stellar Composition → Planet Composition

A key GREP theme is tracing galactic/stellar chemistry to planet building blocks:

**Bitsch & Battistini (2020, A&A 633, A10):** Stellar [Fe/H] alone is an insufficient proxy for disk composition. Individual element ratios (Mg/Si, C/O, Fe/Si) vary independently among stars (from GALAH survey) and propagate into significantly different planet compositions.

**Thiabaud et al. (2015, A&A 580, A30):** Using Gibbs energy minimization + formation models, they showed planets do not simply inherit stellar ratios — condensation sequences and volatile partitioning modify them, particularly for C/O.

**Adibekyan et al. (2021, Science 374, 330):** Established an empirical correlation between rocky exoplanet iron fractions (from mass+radius) and host star iron fractions — but planets are systematically more iron-rich than simple stellar scaling predicts.

**ECCOplanets** (Timmermann et al. 2023, A&A 676, A52): An open-source equilibrium condensation code linking stellar abundances to rocky planet bulk compositions via Gibbs energy minimization.

**Galactic chemical evolution matters** (Adibekyan et al. 2015): Thin disk, thick disk, and halo stars produce systematically different building blocks — the "galactic recipe" that GREP's name references.

### 3.4 Connecting Chemistry to N-body Dynamics

**This is the frontier the GREP position addresses.** Current approaches:

- **chemcomp** couples chemistry to 1D disk + single-planet growth tracks. It does *not* do N-body — each planet is an independent 1D simulation.
- **Bern model** (NGPPS) includes N-body + population synthesis + disk evolution but uses simplified composition tracking.
- **No existing code** self-consistently couples multi-body N-body dynamics with detailed chemical species tracking across all accretion pathways.

The GREP vision — implementing chemical compositions for forming exoplanets within N-body simulations, using stellar compositions as input — would fill this gap.

---

## 4. The PLATO Mission Context

**PLATO** (ESA M3 mission, launch late 2026) will provide:
- Planet radii to 3% accuracy for Earth-sized planets
- Planet masses to ~10% (with ground-based RV follow-up)
- **Stellar ages via asteroseismology** — unique capability enabling time-resolved formation constraints
- Precise stellar parameters (T_eff, log g, [Fe/H], individual abundances)
- System architectures including multi-planet systems at habitable-zone distances

**Why PLATO matters for formation models:**
- Age information allows testing whether planet composition or architecture evolves with time
- Precise stellar compositions enable the star→disk→planet chemical pipeline
- Detection of smaller, more distant planets will reveal whether "missing" Solar System analogs exist
- Statistical sample size enables population-level model testing (as NGPPS VII already demonstrates with RV data)

The GREP project is explicitly designed to produce predictions testable against PLATO data.

---

## 5. The GREP Project: Where This Position Fits

**GREP** ("Galactic Recipe for Exo-Planets") is led by Prof. Stephanie Werner at the University of Oslo, funded by the Research Council of Norway (12M NOK, announced September 2025). It sits within the Centre for Planetary Habitability (PHAB).

**The position's core task:** Implement chemical composition tracking within N-body planet formation simulations, using PLATO stellar compositions as input. Specifically:
1. Assess element distribution during condensation in the stellar nebula
2. Determine elemental availability for rock-forming and volatile elements across disks
3. Trace chemical composition through dynamical formation models — building blocks and feeding zones for each forming planet
4. Incorporate new disk models for N-body simulations

**How this connects to the state of the art:**
- The position would bridge the gap between chemcomp-style chemical tracking (currently 1D, single-planet) and N-body multi-planet dynamics (currently chemistry-poor).
- It would use the stellar composition → disk composition pipeline (Bitsch & Battistini 2020; ECCOplanets; Thiabaud et al. 2015) as input.
- It would test predictions against PLATO observables.
- This is genuinely novel: no existing framework does full N-body + per-species chemical tracking.

---

## 6. Key Open Questions and Interview Talking Points

### Debates and Unresolved Issues

1. **Pebble vs. planetesimal accretion for composition:** Danti et al. (2023) showed these pathways leave distinct chemical signatures — but measuring them requires atmospheric characterization beyond what current facilities routinely provide. JWST is beginning to constrain C/O and metallicities for hot Jupiters; PLATO + Ariel will extend this.

2. **Where do planetesimals form?** The streaming instability requires enhanced dust-to-gas ratios, achievable at ice lines, dead zone edges, and pressure bumps. The spatial distribution of planetesimal formation directly sets the chemical inventory available to forming planets — and current models produce very different predictions.

3. **Disk wind vs. viscous evolution:** MHD disk winds are replacing α-viscosity as the angular momentum transport mechanism. This fundamentally changes disk structure, pebble drift, and ice line locations — but has not yet been coupled to chemical composition tracking.

4. **The missing Solar System analogs:** Observed exoplanet systems are typically compact; Solar System-like architectures with wide-orbit gas giants + inner terrestrial planets are rare in the data. Is this observational bias, or does it reflect a rare formation pathway?

5. **Initial conditions sensitivity:** Lorek & Lambrechts (2026) suggest weak sensitivity to initial planetesimal distributions (pebble accretion self-regulates), but this has only been tested for outer-disk giant planet formation. Inner disk terrestrial planet formation may be more sensitive.

6. **Fragmentation velocity:** This single poorly-constrained parameter (v_frag ~ 1–10 m/s) determines whether dust growth is fragmentation- or drift-limited, affecting everything downstream.

### What the Candidate Should Know

- The **chemcomp** code and its limitations (1D, single planet)
- The **Bern/NGPPS** model and how it compares with pebble-focused approaches
- How **feeding zones** are computed in N-body simulations and what they mean for composition
- The **condensation sequence** and how it maps stellar abundances to solid/gas compositions
- **PLATO's** unique contribution: ages + precise stellar parameters
- Why **C/O alone is insufficient** — the Danti et al. (2023) result
- The **streaming instability** as the planetesimal formation mechanism and its uncertainties
- **ALMA** disk substructures and what they imply for dust trapping and planetesimal formation

---

## Open Questions

- How will GREP handle the computational cost of coupling per-species chemistry to N-body? Will they extend chemcomp, or build a new module within an existing N-body code?
- Which N-body code will they use? REBOUND (modular, Python-friendly) vs. GENGA (GPU, high particle count)?
- How will condensation sequences be computed — equilibrium (Gibbs minimization as in ECCOplanets) or kinetic chemistry?
- Will the disk model be 1D (computationally tractable) or 2D (needed for non-axisymmetric features)?
- How will PLATO target selection and detection biases be forward-modeled into population predictions?
