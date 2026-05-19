# Deep Research Plan: Current State of Planet Formation Simulations

**Slug:** planet-formation-simulations-state
**Date:** 2026-05-13
**Purpose:** Prepare a postdoc interview briefing on the state of the art in planet formation simulations, covering N-body, dust evolution, and chemistry — tailored to the GREP project at the Centre for Planetary Habitability (Oslo).

**Reference paper:** Drazkowska et al. (2023), "Composition of giant planets: The roles of pebbles and planetesimals", A&A 679, L7
**Target project:** GREP — "Galactic Recipe for Exo-Planets", focusing on chemical composition of forming exoplanets, PLATO mission context, condensation sequences, feeding zones, building blocks.

---

## Key Questions (Updated for GREP)

1. **N-body simulations of planet formation:** What are the leading codes and methods for simulating planetary system assembly from planetesimals/embryos? How do they handle disk–planet interaction, migration, and multi-body dynamics? Key codes: REBOUND, Mercury, SyMBA, GENGA, Genesis. What recent results address compact vs extended system architectures?

2. **Dust evolution & planetesimal formation:** State of the art in dust growth, drift, fragmentation → planetesimal formation. How do pebble flux and streaming instability feed into N-body initial conditions? DustPy, two-population models. Connection to disk structure.

3. **Chemical composition in planet formation models:** How is elemental/mineral composition tracked through disk evolution and planet assembly? Condensation sequences, ice lines, volatile transport. How are feeding zones used to assign bulk composition? Key frameworks: chemcomp, Bitsch & Battistini models, Thiabaud/Marboeuf approaches. How does stellar composition (PLATO context) propagate to planet composition?

4. **Coupling chemistry to N-body dynamics:** This is the GREP frontier — how do existing models connect chemical reservoirs in disks to the dynamical accretion history of individual planets? What are the gaps? What new disc models feed into N-body codes?

5. **PLATO mission context:** What observational constraints will PLATO provide? How do formation models need to prepare for PLATO-era data (mass, radius, age, host star composition)?

6. **Interview readiness:** Key open debates, unsolved problems, and where the GREP approach fits in the landscape.

---

## Evidence Needed

- Recent reviews (2022–2026) on planet formation, especially Drazkowska et al. (2023 PPVII review)
- N-body planet formation codes and recent benchmark/comparison studies
- Dust evolution codes and their coupling to planet formation
- Chemical composition models in planet formation (condensation, volatile delivery)
- PLATO mission science case and its implications for formation models
- The GREP team's previous work (Oslo/PHAB group)
- Competing models: pebble accretion vs planetesimal accretion for composition

---

## Scale Decision

**Multi-faceted survey across 3+ domains with interview preparation angle → 3 researcher subagents**

- T1: N-body planet formation simulations + system architecture
- T2: Dust evolution + planetesimal formation + pebble flux
- T3: Chemical composition tracking + PLATO context + GREP-relevant coupling

---

## Task Ledger

| ID | Owner | Description | Status |
|----|-------|-------------|--------|
| T1 | researcher | N-body simulations: codes, migration, architectures, recent results | PENDING |
| T2 | researcher | Dust evolution: growth models, streaming instability, pebble flux | PENDING |
| T3 | researcher | Chemistry in planet formation + PLATO + composition-dynamics coupling | PENDING |
| T4 | lead | Synthesize into interview-ready briefing | PENDING |
| T5 | verifier | Citation verification | PENDING |
| T6 | reviewer | Verification pass | PENDING |

---

## Verification Log

| Check | Result | Notes |
|-------|--------|-------|
| Plan written | ✓ | Updated with GREP context |
| Research complete | PENDING | |
| Draft written | PENDING | |
| Citations verified | PENDING | |
| Review pass | PENDING | |
| Final artifacts | PENDING | |

---

## Decision Log

| Decision | Rationale |
|----------|-----------|
| 3 subagents | Three distinct simulation domains warrant parallel evidence gathering |
| GREP/PLATO focus | Job ad emphasizes composition tracking + PLATO + N-body coupling |
| Include feeding zones topic | Central to GREP: tracing composition back to dynamical formation |
| Emphasize chemistry–dynamics coupling | This is the novel frontier the position targets |
