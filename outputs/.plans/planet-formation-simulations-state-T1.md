# T1 Brief: N-body Planet Formation Simulations

## Goal
Survey the current state of N-body simulations for planet formation and planetary system assembly. Write findings to `outputs/.drafts/planet-formation-simulations-state-research-T1.md`.

## Key Topics to Cover

1. **Leading N-body codes:** REBOUND (Rein & Liu 2012), Mercury/Mercury6 (Chambers 1999), SyMBA (Duncan et al. 1998), GENGA (Grimm & Stadel 2014), Genesis. What integrators do they use (WHM, IAS15, MERCURIUS, hybrid symplectic)? GPU acceleration?

2. **Recent methodological advances (2022–2026):** Hybrid integrators, close-encounter handling, GPU-based N-body (GENGA), coupling to 1D/2D disk models, treatment of gas drag and type I/II migration within N-body.

3. **Key recent results on system architectures:** What do N-body simulations predict about compact vs extended systems? Can they reproduce the observed exoplanet diversity? What about the "missing" Solar System analogs? Key papers on synthetic planet populations (Bern model, Emsenhuber et al., Burn et al.).

4. **Initial conditions problem:** How do N-body simulations set up initial planetesimal/embryo distributions? How sensitive are results to these choices? Connection to dust evolution models.

5. **Disk models used in N-body:** What 1D disk evolution models feed into N-body codes? How is disk dispersal handled? Recent work on disk winds vs photoevaporation.

## Search Strategy
- Use web search for recent papers and reviews on N-body planet formation
- Use paper search for key references (REBOUND, GENGA, Bern model, etc.)
- Look for the PPVII review by Drazkowska et al. (2023) for context
- Check for any GREP-team publications from Oslo/PHAB

## Output
Write a structured research note covering each topic above with citations (author, year, journal, URL where available). Save to `outputs/.drafts/planet-formation-simulations-state-research-T1.md`.
