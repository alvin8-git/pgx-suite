# Chapter Plan + INSIGHT Collection — pgx-suite manuscript

**Mode:** academic-paper `plan` · **Type:** Full research/methods article (IMRaD, ~4,500–6,000 words)
**Venue class:** Clinical pharmacogenomics (candidate journals: *The Pharmacogenomics Journal*, *J Molecular Diagnostics*, *Clin Pharmacol Ther*)
**Citations:** Vancouver (numbered) · **Contribution posture:** tiered (method → validation → engineering)

---

## INSIGHT Collection

- **[INSIGHT: gap]** Individual NGS star-allele callers disagree with one another — both as pure *nomenclature* differences (`*9A` vs `c.85T>C` vs `rs1801265`) and as *substantive* algorithmic differences (ILP vs genome-graph vs Bayesian vs targeted-SV). A clinical lab running more than one caller has no principled way to collapse that disagreement into the single CPIC-actionable call a report requires. Existing tools are either single callers (Aldy, StellarPGx, Cyrius, PharmCAT) or caller *benchmarks*; none is a *reconciliation layer*.

- **[INSIGHT: thesis — DRAFT, author to confirm/edit]** A fixed-hierarchy verdict-reconciliation engine — coverage gate → synonym collapse → authoritative-caller override → phenotype consensus — converts discordant nine-caller star-allele output into one defensible, auditable clinical verdict. Validated on 26 WGS samples (two platforms) against an orthogonal ThermoFisher Axiom PGx array, reconciliation raises array-concordance from **73.0% (241/330) to 85.8% (283/330)** after adjudication, and the residual discordance is dominated not by pipeline error but by **array panel gaps that WGS correctly rescues** (DPYD). Reconciled WGS multi-calling is therefore at least as sensitive as, and frequently more complete than, a targeted genotyping array.

- **[INSIGHT: contribution_claim — tiered, per author steer "All three"]**
  1. *Method (primary):* the 4-tier reconciliation algorithm as a reusable design pattern for multi-caller adjudication.
  2. *Evidence:* orthogonal array validation quantifying the reconciliation gain and characterizing every discordance class.
  3. *Substrate:* a single containerized Snakemake pipeline unifying 9 heterogeneous callers over 37 loci (19/19 CPIC Level A) with non-fatal orchestration — clinical-grade reproducibility.

---

## Chapter-by-chapter plan

### 1. Introduction (~700 w)
- **Core argument:** PGx is clinically actionable (CPIC Level A) but star-allele calling from NGS is hard (SVs at CYP2D6, homology, indels, phasing) → multiple callers exist → they disagree → labs are stuck. Reconciliation is the missing layer.
- Evidence: CPIC guideline landscape [cite CPIC], per-caller method papers [PyPGx, Stargazer, Aldy, StellarPGx, Cyrius, PharmCAT], prior benchmark studies showing inter-caller discordance.
- End on the tiered contribution.
- **Gap (one sentence):** No tool reconciles ≥9 heterogeneous star-allele callers into one CPIC-actionable verdict with an auditable authority trail.

### 2. Materials & Methods (~1,600 w) — methodological centerpiece
- **2.1 Pipeline architecture:** Docker multi-stage build; Snakemake wildcard-`{gene}` DAG; `genes.tsv` as single source of truth; sentinel/`.status` non-fatal caller pattern (a failing caller never aborts the batch); sanity gates (truncated-BAM/wrong-reference). GRCh38-only.
- **2.2 The nine callers** → **Table 1** (caller, version, algorithm class, genes, license).
- **2.3 Gene panel** → **Table 2** (37 loci, per-caller support matrix, CPIC level; 19/19 Level A).
- **2.4 Reconciliation engine (the core):** ordered tiers with worked examples —
  (1) **coverage gate** (NO_CALL below min-depth 10);
  (2) **variant-tier synonym collapse** (`reconcile.py` + `allele_synonyms.json`; DPYD 5-rsID cluster, CYP2C19 `*38→*2`, TPMT `*1S→*1`; VKORC1/NUDT15 *deliberately excluded*);
  (3) **authoritative-caller override** (`AUTHORITATIVE`: Cyrius→CYP2D6, PharmCAT→UGT1A1/CYP2B6/CYP4F2, VCF-Check→RYR1/CACNA1S/G6PD/VKORC1/IFNL3/CYP2C-cluster);
  (4) **phenotype consensus** (rescue to concordant when ≥2 phenotype-emitting callers agree on normalized CPIC phenotype; abstainers ignored).
  → **Figure 1** = reconciliation flowchart.
- **2.5 Validation design:** 13 clinically-genotyped samples (Axiom array) re-sequenced on Illumina + MGI = 26 platform-sample pairs; 16 overlapping genes; 330 gene-calls; adjudication categories (concordant, array_FN, array_authoritative, array_FP_or_NGS_FN, review, ngs_unresolved). Provenance: all numbers script-generated (`build_axiom_concordance.py`), none hand-entered.
- **2.6 Reference materials:** HG002/NA24385 (GIAB) and HG001/NA12878 vs GeT-RM.

### 3. Results (~1,400 w)
- **3.1** Headline: 73.0% raw → 85.8% adjudicated → **Figure 2** waterfall (concordant / array_FN / array_authoritative / residual).
- **3.2** Per-gene concordance → **Table 3 / Figure 3** (6 genes at 100%: ABCG2, CYP2B6, CYP3A4, CYP3A5, NAT2, TPMT; 10 at ≥90%).
- **3.3** Platform equivalence: MGI 73.9% vs ILMN 72.1% (1.8 pp) → **Figure 4**.
- **3.4** Discordance anatomy → **Table 4**: DPYD (24/26 array_FN — array reports `*1/*1`, callers concordantly detect reduced-function variant); VKORC1 (12 array_authoritative — rs9923231 correct, Aldy H-haplotype notation noise); CYP2D6 (Cyrius authority, hybrid-allele nomenclature residual, honest abstention).
- **3.5** Reference-material accuracy: HG002 CYP2D6 `*2/*4` (AS 1.0, IM), 4/4 concordant; HG001 19/28 GeT-RM exact.
- **3.6** Reconciliation contribution: count of calls rescued by synonym collapse and by phenotype consensus (quantifies each tier's value).

### 4. Discussion (~800 w)
- Reconciled WGS ≥ array; array gaps are clinically material (DPYD before fluoropyrimidine dosing).
- The authority hierarchy as a transferable principle: a specialist caller (targeted SV/CNV, CPIC matcher, direct-variant) outranks a general star-allele caller *for its specialty gene*.
- Nomenclature as an under-appreciated pseudo-discordance; synonym collapse as the fix.
- **Limitations:** small clinical n; several genes lack orthogonal truth (RYR1); CYP2D6 hybrid nomenclature unresolved; GRCh38-only; image not publicly redistributable (Stargazer/Aldy non-commercial, Cyrius PolyForm Strict); phenotype tier depends on callers emitting phenotype.
- Prior work: reconciler vs single-caller vs benchmark distinction.

### 5. Conclusion (~200 w)
One-paragraph: reconciliation makes multi-caller WGS operationally clinical; it matches/exceeds an orthogonal array and exposes array blind spots.

### Mandatory back-matter
Data Availability (see open Q3) · Ethics (Q2) · Author Contributions/CRediT (Q4) · Conflict of Interest · Funding · AI-usage disclosure · Limitations (in Discussion).

---

## Figures & Tables (7)
| # | Type | Content |
|---|------|---------|
| Table 1 | table | Nine callers: version, algorithm, genes, license |
| Table 2 | table | 37-locus panel × caller-support matrix, CPIC level |
| Table 3 | table | Per-gene concordance (MATCH/PARTIAL/MISMATCH/n/%) |
| Table 4 | table | Discordance adjudication categories + worked examples |
| Figure 1 | diagram | 4-tier reconciliation flowchart |
| Figure 2 | waterfall | 73.0% → 85.8% concordance by category |
| Figure 3 | bar | Per-gene concordance (sorted) |
| Figure 4 | paired bar | MGI vs Illumina platform equivalence |

---

## OPEN QUESTIONS — author-owned, cannot be derived from the code
1. **Thesis wording** — confirm or edit the DRAFT thesis above.
2. **Ethics/IRB** — the 13 clinical samples need an IRB/ethics approval reference (protocol #, institution). Required for any clinical PGx journal. *What is the approval and consent basis?*
3. **Data availability** — image is non-redistributable and samples are clinical. What *can* be shared: source code? `genes.tsv`/`allele_synonyms.json`? summary concordance tables? (Patient BAMs and the built image presumably cannot.)
4. **Authorship, affiliations, funding** — author list + CRediT roles + funding sources.
5. **Final journal** — pick within the clinical-PGx set (word limits differ: Pharmacogenomics J ~5k, J Mol Diagn ~5k, CPT variable).
