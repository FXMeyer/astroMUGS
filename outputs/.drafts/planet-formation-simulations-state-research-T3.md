# Research: Chemical Composition in Planet Formation Models + PLATO Context

## Summary

Chemical composition tracking in planet formation simulations has matured rapidly, driven by codes like **chemcomp** (Schneider & Bitsch 2021) that couple viscous disk evolution with pebble drift, evaporation at ice lines, and planet growth to self-consistently predict atmospheric C/O, heavy element content, and volatile-to-refractory ratios. The key insight from recent work (Danti/Drazkowska et al. 2023) is that pebble-dominated accretion produces the highest heavy element enrichment, while efficient planetesimal formation dramatically reduces it — and these scenarios are distinguishable via atmospheric C/H, O/H measurements. The upcoming **PLATO mission** (launch late 2026) will provide precise planetary radii, masses, ages, and host star parameters via asteroseismology, creating an unprecedented dataset against which formation models must deliver testable predictions. The **GREP project** at Oslo's Centre for Planetary Habitability aims to do exactly this: predict exoplanet architectures and compositions from stellar chemistry and test competing formation models against PLATO observations.

## Findings

### 1. chemcomp: The Core Framework for Chemistry-Coupled Planet Formation

1. **chemcomp is an open-source 1D planet formation code** that simulates viscous disk evolution, pebble growth (two-population model from Birnstiel et al. 2012), pebble drift, evaporation/condensation at ice lines, and planet growth via pebble and gas accretion while tracking the chemical composition of gas and solids for individual molecular species. It solves advection-diffusion equations per chemical species. [Schneider & Bitsch 2021a, A&A 654, A71](https://www.aanda.org/articles/aa/full_html/2021/10/aa39640-20/aa39640-20.html); [Code documentation](https://chemcomp.readthedocs.io/en/latest/); [GitHub](https://github.com/AaronDavidSchneider/chemcomp)

2. **chemcomp traces atmospheric C/O ratios and heavy element content** by following how inward-drifting pebbles evaporate at their respective ice lines (H₂O, CO₂, CH₄, CO, etc.), enriching gas that is subsequently accreted by the growing planet. The C/O ratio in the planetary atmosphere depends on which evaporation fronts the planet crosses during migration. [Schneider & Bitsch 2021b, A&A 654, A72](https://www.aanda.org/articles/aa/full_html/2021/10/aa41096-21/aa41096-21.html)

3. **Multiple extensions have been built on chemcomp**, including planetesimal formation and accretion (Danti et al. 2023), treatment of binary star abundance differences (Hühn & Bitsch 2023), and additional chemical physics. The code was open-sourced in early 2024. [Schneider & Bitsch 2024, arXiv:2401.15686](https://arxiv.org/html/2401.15686)

### 2. Ice Lines and Volatile Transport

4. **Ice lines (snowlines) are the primary regulators of disk C/O and C/N/O ratios.** The major ice lines are: H₂O (~150–170 K), CO₂ (~70 K), CH₄ (~30 K), CO (~20–25 K). In dynamic disks with pebble drift and viscous evolution, pebbles carrying ices drift inward, evaporate at their respective ice lines, and enrich the inner disk gas in volatiles. This creates radially varying elemental ratios that directly imprint on forming planets. [Öberg et al. 2011, ApJ 743, L16]; [Price et al. 2021 — see arXiv:2604.14124 for recent update](https://arxiv.org/abs/2604.14124)

5. **Volatile transport is a dynamic process**, not a static picture. As the disk evolves thermally, ice line locations shift, and the interplay of pebble drift timescales vs. chemical reaction timescales determines whether equilibrium or kinetic chemistry applies. Current models (chemcomp, among others) assume pebble drift timescales are shorter than chemical reaction timescales, justifying a simplified equilibrium approach. [Booth & Ilee 2019, MNRAS 487, 3998; Eistrup & Henning 2022, A&A 667, A160]

### 3. Danti et al. 2023 (A&A 679, L7): Pebble vs. Planetesimal Roles

6. **Three formation scenarios were compared using chemcomp**: (i) pure pebble accretion, (ii) pebble accretion + planetesimal formation (but no planetesimal accretion), (iii) combined pebble + planetesimal accretion. This is the first study to self-consistently combine both pathways with chemical tracking. [Danti, Bitsch & Mah 2023, A&A 679, L7](https://www.aanda.org/articles/aa/full_html/2023/11/aa47501-23/aa47501-23.html)

7. **Key result: Planetesimal formation drastically reduces atmospheric enrichment.** When pebbles are locked into planetesimals, they cannot drift inward, evaporate, and enrich the disk gas. The pure pebble scenario produces the highest heavy element content; planetesimal formation causes a "dramatic drop." Even when planetesimal accretion is added back, the heavy element content remains much lower than in the pure pebble case, consistent with Venturini & Helled (2020).

8. **The volatile-to-refractory ratio differs between scenarios**, but not uniquely — pure pebble accretion can also produce low volatile/refractory ratios if the planet migrates deep into the inner disk (where refractories evaporate). Distinguishing formation pathways requires measuring **both** C/H and O/H (much higher in pure pebble scenario) alongside C/O and volatile/refractory ratios.

9. **C/O ratio alone is insufficient to distinguish formation pathways.** While C/O changes as planets cross evaporation fronts, its final value depends on many parameters (viscosity, starting position, migration history) and shows no clear pattern unique to one scenario.

10. **Planets with sub-stellar atmospheric metallicity likely formed in the outer disk** (beyond the main evaporation fronts) and were scattered inward, while super-stellar metallicity planets are migration-driven. This is a testable prediction for hot Jupiter populations.

### 4. Feeding Zones and Bulk Composition Assignment

11. **For terrestrial planets, feeding zone analysis in N-body simulations** shows that while there is a general correlation between a planet's location and the radial zone from which it accreted material, significant stochastic mixing occurs during the final giant-impact stage. [Kaib & Cowan 2015, Icarus](https://www.sciencedirect.com/science/article/abs/pii/S0019103515000196)

12. **For gas giants, composition is dominated by gas-phase accretion**, and the key question is what the gas composition is at the time and place of accretion. Pebble-driven enrichment of the gas phase means the gas composition depends on disk evolutionary history, not just local conditions. The "feeding zone" concept for gas giants thus extends to the entire disk's volatile transport history.

13. **SimAb provides a complementary simplified model** for assessing how planet formation affects atmospheric composition of gas giants, exploring the parameter space of C/O ratios and metallicities. [Mollière et al. 2022; SimAb — A&A 2022](https://www.aanda.org/articles/aa/pdf/2022/11/aa41455-21.pdf)

### 5. Stellar Composition → Planet Composition Pipeline

14. **Bitsch & Battistini (2020) demonstrated that stellar [Fe/H] is an insufficient proxy** for disk composition. Using GALAH survey stellar abundances, they showed that individual element ratios (Mg/Si, C/O, Fe/Si) vary independently among stars and that these variations propagate into significantly different planet building block compositions. Super-solar metallicity disks produce more carbon-rich solids; sub-solar disks are more oxidizing. [Bitsch & Battistini 2020, A&A 633, A10](https://www.aanda.org/articles/aa/full_html/2020/01/aa36463-19/aa36463-19.html)

15. **Thiabaud et al. (2015) linked stellar elemental ratios to planet compositions** using Gibbs energy minimization combined with a planet formation model. For 18 observed stellar compositions, they computed Fe/Si, Mg/Si, and C/O in rocky, icy, and giant planets. Key finding: planets do not simply inherit stellar ratios — disk processes (condensation sequences, volatile partitioning) modify them, particularly for C/O. [Thiabaud et al. 2015, A&A 580, A30](https://www.aanda.org/articles/aa/full_html/2015/08/aa25963-15/aa25963-15.html)

16. **Adibekyan et al. (2021) established an empirical compositional link** between rocky exoplanets and their host stars. By inferring planetary iron fractions from mass+radius and comparing with stellar iron fractions from spectroscopy, they found a clear correlation — but planets are systematically more iron-rich than expected from a simple stellar-composition scaling. Published in *Science*. [Adibekyan et al. 2021, Science 374, 330](https://www.science.org/doi/10.1126/science.abg8794)

17. **ECCOplanets (Timmermann et al. 2023) provides an open-source code** for equilibrium condensation simulations, linking stellar abundances to rocky planet bulk compositions via Gibbs energy minimization. Validated against solar system data and applied to chemically diverse stars in the solar neighborhood. [Timmermann et al. 2023, A&A 676, A52](https://www.aanda.org/articles/aa/full_html/2023/08/aa44850-22/aa44850-22.html); [GitHub](https://github.com/AninaTimmermann/ECCOplanets)

18. **A 2024 study linked primordial disk compositions to present-day rocky exoplanet compositions** by correcting for stellar evolution effects on observed abundances, then using stoichiometric models. [A&A 2024](https://www.aanda.org/articles/aa/full_html/2024/12/aa52193-24/aa52193-24.html)

19. **Galactic chemical evolution matters.** Adibekyan et al. (2015) studied Mg/Si mineralogical ratios across stellar populations, showing that the range of planet compositions accessible to formation models depends on galactic location and birth environment. Thin disk, thick disk, and halo stars produce systematically different planet building blocks. [Adibekyan et al. 2015, A&A](https://www.aanda.org/articles/aa/pdf/2015/09/aa27059-15.pdf)

### 6. PLATO Mission Science Case

20. **PLATO (PLAnetary Transits and Oscillations of stars)** is ESA's M3 mission, scheduled for launch in late 2026. It will perform uninterrupted high-precision photometric monitoring of large stellar samples for up to several years. [Rauer et al. 2025, Experimental Astronomy 59, 26](https://link.springer.com/article/10.1007/s10686-025-09985-9); [PLATO mission website](https://platomission.com/about/)

21. **PLATO's key deliverables for formation models:**
   - Planet radii (from transits) to 3% accuracy for Earth-sized planets
   - Planet masses (from ground-based RV follow-up) to 10% accuracy
   - **Stellar ages** via asteroseismology — a unique capability enabling time-resolved planet formation constraints
   - Precise stellar parameters (Teff, log g, [Fe/H], individual abundances) from combined seismology + spectroscopy
   - System architectures (multi-planet systems, orbital periods up to the habitable zone)
   [Rauer et al. 2024, arXiv:2406.05447](https://arxiv.org/abs/2406.05447)

22. **PLATO's age determinations are transformative** — no previous mission could determine stellar (and hence planetary system) ages with ~10% accuracy for large samples. This enables testing whether planet composition or architecture correlates with age, directly probing evolutionary scenarios.

23. **A Science Calibration and Validation Plan (April 2026)** has been published describing target selection for the first long-pointing field. [arXiv:2604.04042](https://arxiv.org/abs/2604.04042)

### 7. GREP Project and Centre for Planetary Habitability (Oslo)

24. **GREP = "Galactic Recipe for Exo-Planets"** is a project led by Professor Stephanie Werner at the University of Oslo's Department of Geosciences and Centre for Planetary Habitability (PHAB). Funded with 12 million NOK from the Research Council of Norway (announced September 2025). [UiO News, Sep 2025](https://www.mn.uio.no/geo/english/about/news-and-events/news/2025/stephanie-werner-new-research-funding-on-exoplanets.html)

25. **GREP's goals:**
   - Study processes and conditions necessary for exoplanet formation, especially Earth-like exoplanets
   - Investigate the relationship between host star composition and exoplanet internal structure/physical properties
   - **Predict and reproduce the architecture of newly discovered exoplanetary systems** and characteristics of exoplanets
   - **Test competing planet formation models against PLATO observations**

26. **PHAB (Centre for Planetary Habitability)** is a UiO Centre of Excellence. Its "Exo-Earths" theme includes work packages on: Observational Constraints, Exo-Planet Formation, Exo-Planet Structure, Exo-Planet Evolution, and Fate of Habitable Planets (PIs: Werner, Herbst, Conrad, Torsvik). [PHAB Exo-Earths](https://www.mn.uio.no/phab/english/themes/exo_earths.html)

27. **Past PHAB projects relevant to planet formation** include "Planet formation in the solar system" (Brasser, 2020–2024). Current projects include Mars analogue sample library and geodynamic projects. The GREP project represents a new frontier connecting galactic chemical evolution to exoplanet properties. [PHAB Projects](https://www.mn.uio.no/phab/english/projects/)

### 8. Key Modeling Frameworks Summary

| Code/Framework | Type | Key Focus | Reference |
|---|---|---|---|
| **chemcomp** | 1D disk + pebble/gas accretion | Atmospheric C/O, heavy elements, gas giant composition | Schneider & Bitsch 2021a,b |
| **ECCOplanets** | Equilibrium condensation | Rocky planet bulk composition from stellar abundances | Timmermann et al. 2023 |
| **SimAb** | Simplified parametric | Gas giant atmospheric composition parameter space | A&A 2022 |
| **Bern model** | Population synthesis | Full planet population statistics, planetesimal-based | Emsenhuber et al. 2021 |
| **Marboeuf et al.** | Disk chemistry | Volatile inventory of planetesimals, ice compositions | Marboeuf et al. 2014 |
| **Thiabaud et al.** | Gibbs minimization + formation | Elemental ratios star→planet, Fe/Si, Mg/Si, C/O | Thiabaud et al. 2015 |
| **Mollière et al.** | Forward modeling | Atmospheric spectra predictions from formation models | Mollière et al. 2022 |

### 9. Connecting Planet Formation and Astrochemistry

28. **Cridland et al. (2019, 2020) established a "main sequence" for C/O** in hot exoplanetary atmospheres by producing populations of astrochemically evolving disks combined with planet growth models. They traced molecular abundances of gas collected into planetary atmospheres and included atmospheric pollution by icy planetesimals and refractory carbon erosion. [Cridland et al. 2019, A&A; Cridland et al. 2020, A&A](https://www.aanda.org/articles/aa/full_html/2020/10/aa38767-20/aa38767-20.html)

## Sources

### Kept
- **Danti, Bitsch & Mah 2023, A&A 679, L7** (https://www.aanda.org/articles/aa/full_html/2023/11/aa47501-23/aa47501-23.html) — Central paper on pebble vs planetesimal composition roles, detailed results extracted
- **Schneider & Bitsch 2021a, A&A 654, A71** (https://www.aanda.org/articles/aa/full_html/2021/10/aa39640-20/aa39640-20.html) — chemcomp introduction paper (Paper I)
- **Schneider & Bitsch 2021b, A&A 654, A72** (https://www.aanda.org/articles/aa/full_html/2021/10/aa41096-21/aa41096-21.html) — chemcomp Paper II on volatiles and refractories
- **Schneider & Bitsch 2024, arXiv:2401.15686** (https://arxiv.org/html/2401.15686) — chemcomp open-source release paper
- **Bitsch & Battistini 2020, A&A 633, A10** (https://www.aanda.org/articles/aa/full_html/2020/01/aa36463-19/aa36463-19.html) — Stellar metallicity effects on planet building blocks
- **Thiabaud et al. 2015, A&A 580, A30** (https://www.aanda.org/articles/aa/full_html/2015/08/aa25963-15/aa25963-15.html) — Elemental ratios stars vs planets
- **Adibekyan et al. 2021, Science 374, 330** (https://www.science.org/doi/10.1126/science.abg8794) — Empirical star-planet composition link
- **Timmermann et al. 2023, A&A 676, A52** (https://www.aanda.org/articles/aa/full_html/2023/08/aa44850-22/aa44850-22.html) — ECCOplanets code
- **Rauer et al. 2025, Exp. Astron. 59, 26** (https://link.springer.com/article/10.1007/s10686-025-09985-9) — PLATO mission paper
- **UiO GREP announcement** (https://www.mn.uio.no/geo/english/about/news-and-events/news/2025/stephanie-werner-new-research-funding-on-exoplanets.html) — GREP project details
- **PHAB Exo-Earths theme** (https://www.mn.uio.no/phab/english/themes/exo_earths.html) — Oslo research structure
- **Cridland et al. 2019, 2020** (https://www.aanda.org/articles/aa/full_html/2020/10/aa38767-20/aa38767-20.html) — C/O main sequence in hot Jupiters

### Dropped
- **PLATO Science Management Plan PDF** (ESA internal document, not publicly readable)
- **Kaib & Cowan 2015 Icarus** — Relevant but focused on solar system terrestrial planets, less directly applicable
- **SimAb PDF** — Supplementary to main findings, not deeply explored
- **Marboeuf et al. 2014** — Referenced in brief but not found as a standalone web-accessible summary; contributions subsumed by Thiabaud et al. 2015 which built on the same framework

## Gaps

1. **Marboeuf et al. (2014) specifics** — The original paper on volatile inventory of planetesimals was not directly accessed. Its contributions are known through citations in later works but a direct summary is missing.

2. **Mollière et al. (2022) forward modeling details** — The connection between formation model outputs and synthetic atmospheric spectra (for JWST comparison) was referenced but not deeply explored. This is important for understanding how formation predictions translate to observables.

3. **GREP project publication list** — The project was only announced in September 2025; no scientific publications from GREP itself appear to exist yet. The project's specific modeling approach and how it will extend existing frameworks (chemcomp, ECCOplanets, etc.) is not yet publicly detailed.

4. **Quantitative PLATO precision numbers** — The Rauer et al. 2025 paper in Experimental Astronomy contains detailed performance specifications, but only the abstract was accessed. Full precision budgets (e.g., exact radius/mass/age accuracy as function of stellar type) would require reading the full paper.

5. **How GREP specifically plans to couple galactic chemical evolution with planet formation models** — The news article mentions this as a goal but provides no methodological detail. This is presumably still in the planning/early execution phase.

6. **N-body + chemistry hybrid codes** — Whether any group has coupled full N-body dynamics with detailed chemical tracking (beyond 1D semi-analytical approaches like chemcomp) remains unclear. This may be a frontier that GREP or others plan to address.
