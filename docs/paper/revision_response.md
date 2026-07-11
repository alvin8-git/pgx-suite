# Response to reviewers (round 1)

Overall: MAJOR REVISION (5.35/10). All items below were addressed in-text; items requiring **new data/analysis** (no new samples available this round) are converted to explicit Acknowledged Limitations per honest-reporting practice, with a clear route to close them.

## Critical

**C1 — Adjudication circularity (73.0%→85.8%).** FIXED (framing) + LIMITATION. Abstract and §3.1 now report **73.0% raw as the primary endpoint** with 95% CI, and 85.8% as a *secondary, mechanistic* explanation of the residual. Each reclassification is tied to evidence independent of the pipeline (Axiom probe content for `array_FN`; raw rs9923231 genotype for `array_authoritative`). We disclose that adjudication was unblinded (§3.1, Limitations). Fully closing this needs a blinded/pre-registered adjudication — acknowledged.

**C2 — Authority tier asserted, not validated.** PARTIALLY FIXED (new §3.8) + LIMITATION. §3.8 shows Tier-3 overrides agreed with truth in the checkable cases (HG002 CYP2D6 `*2/*4` = Cyrius; VKORC1 = raw array SNP). A systematic override-vs-truth audit across all Tier-3 genes/reference materials is stated as not-yet-complete (Limitations).

**C3 — Empty tier-rescue placeholder.** FIXED. §3.6 now gives concrete counts: Tier-2 synonym collapse converted 24 DPYD replicates from disagreement to concordant; Tier-3 determined 26 VKORC1 verdicts (12 array_authoritative) + CYP2D6/PharmCAT genes; Tier-1/Tier-4 did not alter verdicts on this gene set. Full 37-locus ledger deferred to Supplementary Table S1.

## Major

- **M1 — CIs.** FIXED. Wilson 95% CIs on both endpoints (Abstract, §3.1); per-gene intervals noted for Table 3.
- **M2 — 26 = 13×2 replicates.** FIXED. Reframed throughout as "13 samples (26 technical replicates)"; §3.3 notes paired structure and calls for McNemar (⟦to report⟧).
- **M3 — Single-caller baseline.** ADDRESSED via new §3.7 (structure + Supplementary Table S2 ⟦pending⟧) + explicit Limitation. Numbers not fabricated.
- **M4 — Selection bias (16/37 genes).** FIXED. Stated in §3.7 and Limitations.
- **M5 — Sensitivity/specificity/PPV.** LIMITATION. §4 and Limitations state concordance ≠ sensitivity; the completeness claim is bounded to DPYD.
- **M6 — Reproducibility vs non-redistributable image.** FIXED. Reconciled explicitly in Limitations (reproducible from recipe, not from a public image).
- **M7 — Overclaiming ("as sensitive as an array").** FIXED. Softened in Abstract, §4, and Conclusion; general sensitivity advantage no longer claimed.

## Minor (8)
Placeholders flagged; figures/tables now present (Tables 1–4, Figures 1–4 via `figures.py`); PyPGx Zenodo DOI flagged for author; `majority`-state scoring, HG001 68% tone, 10× gate rationale, and intensifiers addressed inline. See `peer_review.md` for the full list.

## Residual (author action required before submission)
Ethics/IRB reference; data-availability final wording; PyPGx v0.26.0 Zenodo DOI; fill Supplementary Tables S1 (tier ledger) and S2 (single-caller baseline) from run logs; McNemar paired test.
