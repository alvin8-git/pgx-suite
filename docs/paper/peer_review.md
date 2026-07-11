# Peer Review — "Reconciling nine star-allele callers into one clinical verdict"

**Reviewer role:** Clinical pharmacogenomics methods reviewer (J Mol Diagn / Pharmacogenomics J tier)
**Date:** 2026-07-06
**Recommendation: MAJOR REVISION**

---

## 1. Summary of the manuscript

The authors present *pgx-suite*, a containerized Snakemake pipeline that runs up to nine star-allele callers over one WGS BAM/CRAM (GRCh38) across 37 loci (all 19 CPIC Level A genes) and reconciles their disagreeing outputs into a single per-gene clinical verdict via a four-tier engine (coverage gate → synonym collapse → authoritative-caller override → phenotype consensus). Validation is against a ThermoFisher Axiom PGx array on 13 clinical samples re-sequenced on two platforms (26 pairs, 330 gene-calls) plus GIAB HG002 and GeT-RM HG001. Headline result: raw concordance 73.0% (241/330) rising to 85.8% (283/330) after adjudication.

The engineering is real and the reconciliation-layer framing is a genuine, publishable idea. But the evidence architecture has a structural problem — the headline number is produced by a self-defined reclassification — and the statistical treatment is essentially absent. The paper is well written and honestly hedged in places, which helps, but it currently overclaims relative to what n=13 and a post-hoc adjudication can support.

---

## 2. Scores (1–10)

| Dimension | Weight | Score | Justification |
|---|---|---|---|
| **Originality** | 20% | **7** | The reconciliation layer with an explicit, auditable authority hierarchy over nine heterogeneous callers is, as claimed, not previously published, and the synonym-collapse insight is a legitimate contribution. Novelty is in orchestration/adjudication engineering, not in any new calling algorithm — every caller is off-the-shelf and the "authority hierarchy" is ultimately a hand-curated routing table, so the scientific (as opposed to software) novelty is moderate. |
| **Methodological Rigor** | 25% | **4** | No confidence intervals anywhere; no per-gene denominator reasoning; the adjudication categories (`array_FN`, `array_authoritative`) are author-defined and directly inflate the headline metric; the authority hierarchy that drives the verdict is asserted, never validated against truth. The two "platforms" share the same source DNA, so 26 pairs are not 26 independent observations — this is not acknowledged in the statistics. |
| **Evidence Sufficiency** | 25% | **4** | n=13 clinical samples; reference materials are n=2 (HG001, HG002). No sensitivity/specificity/PPV framing despite a clinical-diagnostics venue. No comparison table vs a single-caller baseline, so the reader cannot see what reconciliation buys over "just run PharmCAT/PyPGx." The per-tier rescue counts — the direct evidence that the engine works — are a literal placeholder (§3.6 ⟦…⟧). 21 of 37 loci have no orthogonal truth at all. |
| **Argument Coherence** | 25%→15% | **6** | The three-contribution structure is clear and the narrative is disciplined. But the load-bearing claim ("at least as good as the array") rests on reclassifying the array's disagreements as array faults, which is coherent only if those reclassifications are independently justified — for DPYD they are, for the residual "review/unresolved" bucket they are not, and the paper does not separate these cleanly. |
| **Writing Quality** | 15% | **8** | Clear, well-organized, appropriate register; references verified with DOIs; limitations section is unusually candid. Marred only by unfilled ⟦PLACEHOLDER⟧ tokens (including one inside Results §3.6) and some rhetorical intensifiers ("decisive point," "centerpiece"). |

**Weighted overall = 7(0.20) + 4(0.25) + 4(0.25) + 6(0.15) + 8(0.15) = 5.35 / 10.**

**Recommendation: MAJOR REVISION.** The core idea merits publication; the evidence and its statistical treatment do not yet support the claims as stated.

---

## 3. Prioritized, actionable revisions

### [CRITICAL]

**C1 — The 85.8% is at risk of circularity; decouple the adjudication from the metric. (Abstract; §2.5; §3.1; §4)**
The improvement from 73.0%→85.8% is achieved by reclassifying 42 array disagreements (30 `array_FN` + 12 `array_authoritative`) as *not counting against the pipeline*. Because the authors both define the categories and assign the calls, a skeptical reader cannot distinguish "the array was genuinely blind" from "we credited ourselves the wins." Fix: (a) Report three numbers explicitly and side-by-side — raw agreement (73.0%), agreement after *only* nomenclature/synonym collapse (a neutral, non-self-serving adjustment), and agreement after clinical adjudication (85.8%) — and make 73.0% the primary endpoint. (b) For every `array_FN` and `array_authoritative` call, provide an *independent* justification not derived from the pipeline's own output: e.g., cite the Axiom panel's published probe list to prove the DPYD c.85T>C locus is absent (this is verifiable and makes the DPYD block airtight), and for `array_authoritative` VKORC1 confirm rs9923231 against the raw array intensity/genotype, not against the caller. (c) Have a blinded adjudicator (or a pre-registered rule set applied by someone other than the pipeline author) assign categories, and state this. Without (a)–(c) the headline metric is not defensible.

**C2 — The authority hierarchy (Tier 3) is asserted, not validated. (§2.4 Tier 3; §3.4; §4)**
Tier 3 lets a single "authoritative" caller *short-circuit the vote* — Cyrius overrides for CYP2D6, PharmCAT for UGT1A1/CYP2B6/CYP4F2, VCF-Check for the single-SNP genes. The paper never demonstrates that the authoritative caller is *right* when it disagrees with the majority; it assumes so by construction. This is the paper's second-biggest vulnerability: if the authority is wrong, the engine confidently emits a wrong clinical verdict and reports high "concordance" because it agrees with the array on the easy cases. Fix: (a) On the GeT-RM/GIAB truth sets, tabulate every case where the Tier-3 authority *overruled* a caller majority and report how often the authority matched consensus truth vs the majority it overrode. (b) Add a failure-mode analysis: what happens when Cyrius abstains or is wrong (the CYP2D6 hybrid cases)? (c) Explicitly bound the claim — Tier 3 is justified only where the specialist has independent published validation for that gene (cite it per gene), otherwise it is an untested heuristic.

**C3 — Fill the per-tier rescue counts; they are the only direct evidence the engine works. (§3.6)**
§3.6 currently ends in ⟦Quantify: exact count of calls rescued by each tier⟧. This is the single most important results paragraph — it is the *mechanistic* evidence for the method, as opposed to the aggregate concordance which could be driven by the 6 already-100% genes. A methods paper about a four-tier engine cannot go to review, let alone publication, with the tier-contribution table empty. Provide: n calls entering each tier, n rescued/changed by each tier, and n where the tier flipped the final verdict, with a worked example per tier.

### [MAJOR]

**M1 — Add statistical treatment: confidence intervals and honest denominators. (§3.1–§3.3; all of Results)**
Every percentage is a point estimate with no CI. With 330 calls the 73.0% CI (Wilson) is roughly ±4.8 pp, and per-gene figures rest on n≈13–26 each so a "100% concordance" gene (ABCG2, n≈26) has a lower 95% bound near 87%, not 100%. Fix: report Wilson 95% CIs for the overall and every per-gene figure in Table 3; state the per-gene n in the table; and stop describing single-digit-sample genes as "100% concordant" without the interval.

**M2 — The 26 pairs are not 26 independent samples. (§2.5; §3.3)**
Illumina and MGI libraries were made from the same 13 source DNAs, so platform pairs are correlated replicates, not independent samples. Presenting 330 as if it were 330 independent calls overstates power, and the "platform equivalence" claim (73.9% vs 72.1%) should be tested as a *paired* comparison (McNemar on the discordant pairs), not an unpaired 1.8-pp difference. Fix: describe the design as 13 biological samples × 2 technical (platform) replicates; do the paired test; and report the effective sample size honestly in the abstract ("13 samples" should lead, "26 pairs / 330 calls" second).

**M3 — Add a single-caller baseline comparison table. (New subsection in §3; Discussion)**
The entire premise is that reconciliation beats picking one caller. That is never shown. The reader needs a table: concordance-with-truth for PyPGx alone, Stargazer alone, Aldy alone, StellarPGx alone, PharmCAT alone, vs the reconciled verdict, on the same 330 calls (or the GeT-RM set). If reconciliation does not beat the best single caller, the paper's thesis is unsupported; if it does, this is the strongest possible evidence and its absence is inexplicable.

**M4 — Selection bias toward the 16 array-overlap genes must be stated as a validity threat, not a footnote. (§2.5; §4 Limitations)**
Only 16 of 37 loci have any orthogonal comparator; the headline concordance is computed only over genes chosen *because* the array covers them, which are largely the well-behaved SNP/small-panel genes. The 21 unvalidated loci include harder structural cases. Fix: state explicitly that concordance is measured on a favorable subset, report which 21 loci have no truth, and soften any pipeline-wide accuracy claim to "on the 16 validatable genes."

**M5 — Add clinical-diagnostics performance framing (sensitivity/specificity/PPV). (§2.5; §3)**
For a clinical molecular-diagnostics venue, "concordance" is insufficient. Reframe at least the actionable-variant detection (e.g., DPYD reduced-function, CYP2D6 metabolizer class) as sensitivity/specificity vs the best available truth, with CIs. This also directly answers whether the DPYD `array_FN` calls are true positives (they need orthogonal confirmation — GIAB/PharmVar — before being counted as pipeline wins, see C1).

**M6 — Reconcile the reproducibility claim with the non-redistributable image. (Abstract; §1 contribution 3; §4; Data availability)**
The abstract sells "clinical-grade reproducibility from one command," but §4 concedes the image cannot be redistributed (Stargazer/Aldy/Cyrius licenses) and reproduction requires a private rebuild. These are in tension. Fix: downgrade the abstract claim to "reproducible from a versioned recipe" and commit concretely in Data Availability to releasing what legally can be released (`genes.tsv`, `allele_synonyms.json`, the reconciliation code, summary call tables, and the exact Dockerfile/lockfile) — otherwise the reproducibility contribution is not evaluable.

**M7 — Overclaiming in Abstract/Discussion/Conclusion. (Abstract; §4; §5)**
"at least as sensitive as, and often more complete than, a targeted PGx array" is a sensitivity claim proven on one gene (DPYD) and never measured as sensitivity. Similarly "matched or exceeded the array" folds the self-adjudicated wins into the comparison. Fix: restrict these statements to what is shown ("on DPYD, the callers concordantly detected a reduced-function variant outside the array's fixed panel"), and remove pipeline-wide sensitivity superiority language until M3/M5 support it.

### [MINOR]

- **m1 — Remove all ⟦PLACEHOLDER⟧ tokens before resubmission**, especially the one inside Results §3.6 and the front-matter/author/ethics blocks; a Results placeholder will trigger a desk reject at many venues. (Throughout)
- **m2 — Figure/Table callouts reference Figures 1–4 and Tables 1–4 that are not in the provided file.** Ensure they are included and that Table 3 carries per-gene n and CIs (see M1). (§2–§3)
- **m3 — Ref [8] (PyPGx) still lacks a Zenodo DOI**, and ref [17] cites mtDNA-Server as a proxy for mutserve; add the Haplocheck cite if contamination is discussed. (References)
- **m4 — Define "verdict status" states on first use** ({concordant, majority, discordant, no_call}) and clarify how `majority` counts toward the concordance numerator — is a `majority` verdict scored as concordant with the array or not? This materially affects the 73.0%. (§2.4; §3.1)
- **m5 — HG001 "19 of 28 comparable genes" (68%) is presented without concern** while the clinical set's raw 73% is framed as good; reconcile the tone and explain the 9 residuals individually or in a supplementary table. (§3.5)
- **m6 — "up to nine callers" varies per gene** — state the mean/median number of callers actually contributing per verdict, since a "reconciliation of nine" that is usually 2–3 callers per gene should be described as such. (§2.2; §2.4)
- **m7 — Coverage gate default 10×** is low for confident diplotype/CNV calling in a clinical context; justify or cite the threshold. (§2.4 Tier 1)
- **m8 — Tone down intensifiers** ("centerpiece," "decisive point," "not decorative") for a methods venue. (§1, §3.6, §4)

---

## 4. What would move this to Accept

Fill §3.6 (C3), add the single-caller baseline (M3) and Tier-3 override validation against truth (C2), report CIs with honest denominators (M1/M2), and split the headline into raw vs nomenclature-only vs clinically-adjudicated with independently-sourced justification for each `array_FN`/`array_authoritative` call (C1). None of these require new samples — they are re-analyses of data the authors already hold — which is why MAJOR (not REJECT) is appropriate: the study can likely support a more careful version of its thesis.
