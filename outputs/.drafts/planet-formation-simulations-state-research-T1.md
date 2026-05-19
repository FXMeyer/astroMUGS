# Research: State of N-body Planet Formation Simulations

## Summary

N-body simulations of planet formation have advanced substantially in 2022–2026, driven by GPU acceleration (GENGA), new hybrid integrators (TRACE in REBOUND), and increasingly realistic coupling to disk models including pebble accretion and migration. The field is converging on a picture where pebble accretion combined with N-body dynamics naturally forms 1–2 gas giants plus several icy/super-Earth companions, with final architectures only weakly sensitive to initial planetesimal distributions. Population synthesis (Bern/NGPPS Generation III model) now undergoes detailed statistical comparison with radial velocity surveys.

## Findings

### 1. Leading N-body Codes

1. **REBOUND** (Rein & Liu 2012) — Open-source, modular Python/C framework. Integrators include WHFast (Wisdom-Holman), IAS15 (high-order adaptive, 15th order), MERCURIUS (hybrid symplectic, Rein et al. 2019), and now **TRACE** (Lu, Hernandez & Rein 2024). CPU-based with OpenMP parallelism. [REBOUND docs](https://rebound.hanno-rein.de/)

2. **GENGA** (Grimm & Stadel 2014; Grimm et al. 2022) — GPU-accelerated hybrid symplectic integrator. Uses CUDA for N²/2 force calculations and close-encounter handling. Capable of simulating thousands of fully-interacting bodies; extended in 2022 to handle non-Newtonian forces and high particle counts. Run times of days–weeks for ~4000 bodies over 100 Myr on Nvidia A30 GPUs. [Grimm et al. 2022, ApJ 932, 124](https://iopscience.iop.org/article/10.3847/1538-4357/ac6dd2)

3. **Mercury/Mercury6** (Chambers 1999) — Classic hybrid symplectic (Bulirsch-Stoer for encounters). Still widely used but showing age; variants include Mercury-T (Bolmont et al. 2015) with tidal forces and Posidonius (Blanco-Cuaresma & Bolmont 2017). CPU-only.

4. **SyMBA** (Duncan, Levison & Lee 1998) — Symplectic Massive Body Algorithm with recursive subdivision for close encounters. Recently parallelized as **SyMBAp** with OpenMP for multi-core CPUs (IOP 2023). Used by Matsumura et al. for pebble accretion studies. [SyMBAp, RNAAS 2023](https://iopscience.iop.org/article/10.3847/2515-5172/accc8a)

5. **Other codes** — GLISSE (Zhang & Gladman 2022, GPU), QYMSYM (Moore & Quillen 2010, GPU), PKDGRAV3 (Potter et al. 2017, CPU parallel), Hiperion (GPU).

### 2. Recent Methodological Advances (2022–2026)

6. **TRACE integrator** (Lu, Hernandez & Rein 2024, MNRAS 533, 3708) — Almost time-reversible hybrid integrator for REBOUND. Combines WHFast and Bulirsch-Stoer/IAS15 internally. Switches methods time-reversibly during close encounters following Hernandez & Dehnen prescription. In all test cases at least as accurate and fast as MERCURIUS; vastly superior for planet-planet close encounters. [arXiv:2405.03800](https://arxiv.org/abs/2405.03800)

7. **TRACE collision/fragmentation support** (Lu et al. 2025, RNAAS) — Extended TRACE to handle collisional fragmentation, adding/removing particles mid-timestep. Enables large-N protoplanetary disk simulations with various collision prescriptions. [IOPscience 2025](https://iopscience.iop.org/article/10.3847/2515-5172/add525)

8. **GENGA II** (Grimm et al. 2022) — Extended to non-Newtonian forces (general relativity, oblateness, tides, gas drag, type I/II migration) and high particle numbers (>10⁴). Enables direct simulation from planetesimal populations to planets. [ApJ 932, 124](https://iopscience.iop.org/article/10.3847/1538-4357/ac6dd2)

9. **High-order regularised symplectic integrator** (Petit et al. 2019, A&A) — Time-regularisation approach keeping symplectic stability properties while reducing step size during close encounters. Resolves near-collisions. [A&A 2019](https://www.aanda.org/articles/aa/full_html/2019/08/aa35786-19/aa35786-19.html)

10. **Coupling to 1D disk models** — Modern N-body codes now routinely include: pebble accretion prescriptions (Lambrechts et al. 2019), gas accretion (Kelvin-Helmholtz + runaway), type I/II migration with eccentricity/inclination damping (Cresswell & Nelson 2008), and pebble flux reduction by multiple accreting bodies. Implementation as CUDA kernels in GENGA demonstrated by Lorek & Lambrechts (2026).

### 3. Key Recent Results on System Architectures

11. **From planetesimals to planets in giant-planet region** (Lorek & Lambrechts 2026, A&A 707, A244) — GPU N-body (GENGA) simulations starting from streaming-instability-inspired planetesimal mass distributions (200 km to ~0.1 M⊕) between 10–50 AU. Key findings: (a) initial spatial distribution does NOT significantly affect outcomes—narrow rings scatter within ~1 Myr; (b) consistently form 1–2 gas giants + several icy planets >2 M⊕; (c) initial total mass (1 vs 10 M⊕ of planetesimals) does not change number of planets >1 M⊕ due to pebble flux self-regulation; (d) scattered disc of planetesimals is generic outcome. [A&A 2026](https://www.aanda.org/articles/aa/full_html/2026/03/aa58429-25/aa58429-25.html)

12. **Pebble accretion N-body: giant planets** (Bitsch et al. 2019, A&A 623, A88; Raorane et al. 2024, Icarus 421, 116231) — 60 embryos between 10–40 AU; gas giants form when core growth time < disk dissipation time. Raorane et al. form on average 1.7 planets >10 M⊕ with multiplicity ~2.5; ice-giant analogs are rare.

13. **Lau et al. 2024** (A&A 683, A204) — Explored giant planet formation via pebble accretion in smooth disks starting from 1000–5000 planetesimals. Without migration: 1–2 gas giants + 1–2 ice giants. With migration: massive cores lost to inner disk unless gas accretion prescription allows early transition to type-II.

14. **Four classes of planetary system architectures** (Emsenhuber, Mordasini & Burn 2023, EPJ Plus 138, 181) — Bern model population synthesis identifies four architectural classes emerging from core accretion: (i) similar-mass, tightly-packed low-mass systems; (ii) ordered architectures; (iii) anti-ordered; (iv) mixed. Compact Kepler-like "peas-in-a-pod" systems emerge naturally. [EPJ Plus 2023](https://epjplus.epj.org/articles/epjplus/abs/2023/02/13360_2023_Article_3784/13360_2023_Article_3784.html)

15. **NGPPS VII — Statistical comparison with HARPS/Coralie** (Burn et al. 2025, A&A 701, A64) — Generation 3 Bern Model compared quantitatively against RV survey completeness. Nominal population reproduces observed distributions well. Systematic differences identified at specific mass/period regimes. [A&A 2025](https://www.aanda.org/articles/aa/abs/2025/09/aa52485-24/aa52485-24.html)

16. **NGPPS VIII — Metallicity effects** (2025, A&A) — Bern Gen-III model at different metallicities: incorporates planetesimal/gas accretion, N-body interactions, giant impacts, disk evolution, long-term contraction, atmospheric loss. Analyses occurrence rates, orbital periods, eccentricities, and radius valley morphology vs host star metallicity. [A&A 2025](https://www.aanda.org/articles/aa/full_html/2025/09/aa55380-25/aa55380-25.html)

17. **Kepler compact multis from planet-planet scattering** — N-body simulations starting with 8 close-packed planets show dynamical instability sculpts final architectures to match observed Kepler spacing distributions. Minimum dynamical separations in stable systems match observations (~8 mutual Hill radii). [AJ 2024](https://iopscience.iop.org/article/10.3847/1538-3881/ad1244/pdf)

18. **Planetary population synthesis review** (Burn & Mordasini 2024, arXiv:2410.00093) — Comprehensive review of the population synthesis method, covering Bern model evolution from Generation I to III. [arXiv:2410.00093](https://arxiv.org/abs/2410.00093)

### 4. Initial Conditions Problem

19. **Streaming instability IMF** — Current standard: exponentially-tapered power law (Schäfer et al. 2017) with slope p≈0.6 below characteristic mass, exponential cutoff above. Characteristic mass scales with disk properties (Liu et al. 2020). Upper end consistent with cold classical Kuiper belt; lower end resolution-dependent (Li et al. 2019; Simon et al. 2022; Gerbig & Li 2023).

20. **Sensitivity to initial conditions** — Lorek & Lambrechts (2026) directly tested narrow-ring vs. wide-ring vs. high-mass distributions between 10–50 AU. Result: **weak sensitivity** — scattering redistributes planetesimals within ~1 Myr, and pebble-flux filtering self-regulates final planet numbers. The biggest effect is on intermediate-mass bodies (10⁻³–1 M⊕): narrow rings produce fewer due to dynamical excitation reducing pebble accretion efficiency.

21. **Starting from embryos vs. planetesimals** — Most earlier works assumed pre-existing lunar-mass embryos (Bitsch et al. 2019: 60 embryos of 0.005–0.015 M⊕; Matsumura et al. 2021: 10 cores of 0.01 M⊕). This skips the critical planetesimal→embryo transition that depends sensitively on assumed mass distribution and surface density (Lorek & Johansen 2022). GPU codes now enable starting from full planetesimal populations (~4000 bodies).

22. **Lau et al. (2024) on smooth vs. structured disks** — Starting from planetesimals rather than embryos is crucial to avoid artificial assumptions. Migration of massive cores to inner system is a persistent problem when starting too late with too-massive seeds.

### 5. Disk Models Used in N-body Simulations

23. **Steady-state α-disk models** — Most common: surface density Σ_g from steady-state accretion (Shakura-Sunyaev). Lorek & Lambrechts (2026) use Ida et al. (2016)/Chambers (2009) scaling laws with irradiation-dominated outer disk (T ∝ r^{-3/7}, Σ ∝ r^{-15/14}). Stellar accretion rate decreases as power law with time (Hartmann et al. 2016).

24. **Disk dispersal** — In N-body works, typically modeled as decreasing accretion rate until cutoff (e.g., 3.8 Myr after simulation start). More sophisticated treatments emerging:
   - **MHD disk winds** (Tabone et al. 2022; Weder et al. 2023, A&A): population-level studies show strong magnetic torques + weak winds needed to match observed accretion rates (~10⁻⁸ M⊙/yr). [A&A 2023](https://www.aanda.org/articles/aa/full_html/2023/06/aa43453-22/aa43453-22.html)
   - **Photoevaporation**: X-ray driven mass-loss rates revised downward by factor ~10 due to enhanced cooling (Sellek et al. 2024, A&A). [A&A 2024](https://www.aanda.org/articles/aa/full_html/2024/10/aa50171-24/aa50171-24.html)
   - **Wind-planet interplay** (Gárate et al. 2024, A&A 686, A53): wind-driven gas redistribution affects planet growth and gap structure. [A&A 2024](https://www.aanda.org/articles/aa/abs/2024/06/aa48596-23/aa48596-23.html)
   - **Magnetic disk wind global model** (2025, A&A): thin-disk limit wind-driven accretion model constrained by theory and observations. [A&A 2025](https://www.aanda.org/articles/aa/full_html/2025/03/aa50236-24/aa50236-24.html)

25. **Bern model disk evolution** — Generation III (Emsenhuber et al. 2021, A&A 656, A69) uses 1D viscous disk with photoevaporation, dust evolution, and ice lines. Confronted at population level against disk observations (Emsenhuber et al. 2023, A&A 673, A78). [A&A 2021](https://www.aanda.org/articles/aa/abs/2021/12/aa38553-20/aa38553-20.html)

26. **Pebble flux models** — Two approaches: (a) fixed fraction of stellar accretion rate (ξ ~ 1–2%, as in Lorek & Lambrechts 2026); (b) pebble production front expanding outward (Lambrechts & Johansen 2014; Ida et al. 2016). Total pebble mass drifting through disk typically ~100–200 M⊕ consistent with observed disk dust masses.

### 6. PPVII Context

27. **Drążkowska et al. (2023, PPVII, ASP Conf. Ser. 534, 717)** — Major review "Planet Formation Theory in the Era of ALMA and Kepler: from Pebbles to Exoplanets." Covers paradigm shift toward early planetesimal formation, pebble accretion as dominant growth mechanism, streaming instability, and the connection between disk substructures and planet formation. Emphasizes that classical models revised to account for observed exoplanet diversity. [ADS](https://ui.adsabs.harvard.edu/abs/2023ASPC..534..717D)

## Sources

### Kept
- Lorek & Lambrechts 2026, A&A 707, A244 (https://www.aanda.org/articles/aa/full_html/2026/03/aa58429-25/aa58429-25.html) — Most detailed recent GENGA planet formation study; covers codes, methods, initial conditions, and results
- Lu, Hernandez & Rein 2024, MNRAS 533, 3708 (https://arxiv.org/abs/2405.03800) — TRACE integrator paper, state-of-art for REBOUND
- Lu et al. 2025, RNAAS (https://iopscience.iop.org/article/10.3847/2515-5172/add525) — TRACE fragmentation extension
- Grimm et al. 2022, ApJ 932, 124 (https://iopscience.iop.org/article/10.3847/1538-4357/ac6dd2) — GENGA II with non-Newtonian forces
- Burn et al. 2025, A&A 701, A64 (https://www.aanda.org/articles/aa/abs/2025/09/aa52485-24/aa52485-24.html) — NGPPS VII statistical comparison
- Emsenhuber, Mordasini & Burn 2023, EPJ Plus 138, 181 (https://epjplus.epj.org/articles/epjplus/abs/2023/02/13360_2023_Article_3784/13360_2023_Article_3784.html) — Four architecture classes
- Burn & Mordasini 2024 review (https://arxiv.org/abs/2410.00093) — Population synthesis review
- Drążkowska et al. 2023, PPVII (https://ui.adsabs.harvard.edu/abs/2023ASPC..534..717D) — Comprehensive review
- Lau et al. 2024, A&A 683, A204 (https://www.aanda.org/articles/aa/full_html/2024/03/aa47863-23/aa47863-23.html) — Pebble accretion in smooth disk
- Rein et al. 2019, MNRAS 485, 5490 (https://ui.adsabs.harvard.edu/abs/2019MNRAS.485.5490R/abstract) — MERCURIUS hybrid integrator
- Weder et al. 2023, A&A (https://www.aanda.org/articles/aa/full_html/2023/06/aa43453-22/aa43453-22.html) — MHD wind population study
- Sellek et al. 2024 (https://www.aanda.org/articles/aa/full_html/2024/10/aa50171-24/aa50171-24.html) — Revised X-ray photoevaporation rates

### Dropped
- SyMBAp parallelization (RNAAS 2023) — brief research note, limited new physics
- Chambers 1998 original Mercury paper — foundational but dated
- Various disk-only papers without N-body coupling
- Howe 2025 (NASA TRS) — classification framework, not N-body simulation paper

## Gaps

1. **GREP-team / Oslo / PHAB publications** — No specific publications from these groups found in search results. Would need targeted search for "GREP planet formation Oslo" or "PHAB group" publications.

2. **Genesis code** — Mentioned in brief but no recent publications found describing it. May be internal/unpublished or renamed.

3. **Direct comparison of code performance** — No systematic benchmark paper comparing GENGA vs REBOUND/TRACE vs SyMBA on identical problems found for 2023–2026 period.

4. **Disk wind coupling to N-body** — While disk wind models are advancing rapidly, direct coupling of MHD-wind-driven disk evolution into N-body planet formation simulations (rather than simple α-disk) remains at early stages. This is an active frontier.

5. **Machine learning / emulator approaches** — Not surveyed; emerging trend of using neural network emulators to replace expensive N-body integrations for population studies.

**Suggested next steps:** Search specifically for Genesis code, GREP-team publications, and any Oslo/Copenhagen group papers on N-body planet formation with disk winds.
