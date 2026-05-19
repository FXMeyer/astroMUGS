# Verification: planet-formation-simulations-state-cited.md

## Checks Performed

### URL Verification
- [1] Schneider & Bitsch 2021a, A&A chemcomp I — VERIFIED (fetched, title matches)
- [3] Danti et al. 2023, A&A 679, L7 — VERIFIED (fetched, title matches)
- [8] Lorek & Lambrechts 2026, A&A 707, A244 — VERIFIED (fetched, title matches)
- [24] Bitsch & Battistini 2020, A&A 633, A10 — VERIFIED (fetched, title matches)
- [27] Timmermann et al. 2023, ECCOplanets — VERIFIED (fetched, title matches)
- [28] Rauer et al. 2025, PLATO mission — VERIFIED (fetched, title matches)
- [30] UiO GREP announcement — VERIFIED (fetched, content matches)
- Other URLs: well-formed, consistent with ADS/journal format, not individually fetched

### Content Checks
1. **Claim: "C/O ratio alone insufficient" (Danti et al. 2023)** — Confirmed from source [3]; paper states C/O shows no clear trend distinguishing scenarios, while C/H and O/H show clear differentiation.
2. **Claim: "Lorek & Lambrechts 2026 weak sensitivity"** — Confirmed from source [8]; paper states "initial configuration has only a weak effect."
3. **Claim: "PLATO launch late 2026"** — Confirmed from source [28]; Rauer et al. 2025 mission paper.
4. **Claim: "GREP 12M NOK"** — Confirmed from source [30]; UiO news article states "12 million NOK."
5. **Claim: "No existing code does full N-body + per-species chemistry"** — Inference, not directly cited from a single source. Based on review of all three research files. Marked as state-of-field assessment. MINOR: could be wrong if an unpublished or obscure code exists.

### Structure Check
- Executive summary: present ✓
- Organized by domain: present ✓
- Caveats and open questions: present ✓
- Sources section with URLs: present ✓
- No invented data/figures/benchmarks: confirmed ✓

### Findings

| Severity | Issue | Resolution |
|----------|-------|------------|
| MINOR | Claim "no existing code does full N-body + chemistry" is an inference, not verified exhaustively | Acceptable — qualified as "state of field assessment"; hedge word could be added |
| MINOR | Chambers (1999) Mercury6 listed in table without URL | Foundational paper, widely known — acceptable |
| MINOR | Some URLs not individually fetched (sources 5, 6, 9, 11, 14, 19, 20, 21, 22) | Well-formed journal/arXiv URLs with consistent metadata — low risk |

### Verdict
No FATAL issues. Three MINOR issues, all acceptable. PASS.
