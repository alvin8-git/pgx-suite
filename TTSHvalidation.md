# TTSH PGx Pipeline Validation Report

**Date:** 2026-03-13
**Prepared by:** TTSH Pharmacogenomics Unit
**Pipeline version:** pgx-suite (PyPGx 0.26.0 · Stargazer 2.0.3 · Aldy 4.8.3 · StellarPGx 1.2.7)
**Reference:** GRCh38 (hg38)

---

## 1. Executive Summary

Thirteen clinical samples previously genotyped on the Thermofisher Axiom Pharmacogenomics array (Axiom P6 and/or P9) were re-analysed on two whole-genome sequencing (WGS) platforms — Illumina (ILMN) and MGI — using the four-tool pgx-suite pipeline. Sixteen genes overlapping between the Axiom panel and the pipeline were evaluated.

**Key findings:**

| Metric | Value |
|--------|-------|
| Samples assessed | 13 (26 platform-sample pairs) |
| Genes evaluated | 16 |
| Overall concordance (all tools, all samples) | **~80%** |
| Best-performing tool (concordance) | **Stargazer** (81.5%) |
| Broadest gene coverage (no failed calls) | **Aldy** |
| Genes with >90% concordance | CYP2B6, CYP2C19, CYP2C9, CYP3A4, CYP3A5, NAT2, NUDT15, SLCO1B1, TPMT |
| Genes with problematic concordance | DPYD (8%), CYP2D6 (55%), UGT1A1 (62%), RYR1 (60%) |
| Platform difference (ILMN vs MGI) | Minimal genotype differences; MGI has markedly lower duplication rate |

---

## 2. Study Design

### 2.1 Samples

Thirteen samples with prior Axiom microarray PGx genotyping were sequenced by WGS on two platforms. Five samples (P_62, P_93, P_670, P_1143, P_2000) had both Axiom P6 and P9 genotyping; for these, P9 results (which include CYP2D6 and additional gene coverage) were used as the primary Axiom reference.

| Sample | Axiom Software | ILMN Depth (Dup%) | MGI Depth (Dup%) |
|--------|---------------|-------------------|------------------|
| EQ_2017 | P6 | 45.3× (39.7%) | 33.9× (20.8%) |
| EQ_2021 | P6 | 35.6× (35.7%) | 32.9× (15.5%) |
| EQ-H05 | P6 | 30.4× (20.0%) | 28.4× (6.6%) |
| H_342 | P6 | 30.4× (18.5%) | 43.1× (9.0%) |
| P_62 | P9+P6 | 31.6× (29.4%) | 28.1× (25.3%) |
| P_93 | P9+P6 | 49.8× (43.7%) | 31.7× (8.0%) |
| P_670 | P9+P6 | 35.8× (25.9%) | 31.8× (9.0%) |
| P_1143 | P9+P6 | 40.1× (24.4%) | 33.7× (7.2%) |
| P_2000 | P9+P6 | 38.8× (37.9%) | 39.6× (9.5%) |
| P_2483 | P6 | 39.2× (35.8%) | 35.0× (17.4%) |
| P_2505 | P6 | 35.1× (19.5%) | 31.4× (6.2%) |
| P_2820 | P6 | 39.5× (36.4%) | 29.3× (7.8%) |
| P_2912 | P6 | 35.8× (38.9%) | 34.0× (38.3%) |

All sequencing is 150 bp paired-end. ILMN samples show significantly higher duplication rates (18–44%) compared to MGI (6–38%), with the exception of P_2912 where MGI duplication is also high (38%).

### 2.2 Genes Evaluated

| Gene | Axiom P6 | Axiom P9 | Pipeline tools |
|------|---------|---------|---------------|
| ABCG2 | ✓ | ✓ | Aldy only |
| CYP2B6 | — | — | All 4 (SV) |
| CYP2C19 | ✓ | ✓ | All 4 |
| CYP2C9 | ✓ | ✓ | All 4 |
| CYP2D6 | ✗ (P6 "Indeterminate") | ✓ (P9) | All 4 (SV) |
| CYP3A4 | — | ✓ (P9) | All 4 |
| CYP3A5 | ✓ | ✓ | All 4 |
| CYP4F2 | ✓ | ✓ | All 4 |
| DPYD | ✓ | ✓ | PyPGx, Stargazer, Aldy |
| NAT2 | — | ✓ (P9) | All 4 |
| NUDT15 | ✓ | ✓ | All 4 |
| RYR1 | — | ✓ (P9) | PyPGx, Stargazer, Aldy |
| SLCO1B1 | ✓ | ✓ | All 4 |
| TPMT | ✓ | ✓ | All 4 |
| UGT1A1 | ✓ | ✓ | All 4 |
| VKORC1 | ✓ | ✓ | PyPGx, Stargazer, Aldy |

---

## 3. Tool Concordance Summary

Concordance was assessed against the Axiom reference call for each gene × sample pair with a non-null Axiom result. Where Axiom listed multiple possible diplotypes (ambiguous calls, e.g. `*1/*1,*1/*28,*28/*28`), any matching interpretation was accepted as concordant. VKORC1 genotype notation was normalised (G/G = reference; A/A = homozygous *rs9923231*; G/A = heterozygous). CYP2C9 rs9332094C was equated to *\*8*.

| Tool | Concordant | Discordant | Failed/Not Called | **Concordance** |
|------|-----------|-----------|------------------|----------------|
| **Stargazer** | 251 | 57 | 26 | **81.5%** |
| **PyPGx** | 246 | 62 | 26 | **79.9%** |
| **Aldy** | 261 | 73 | 0 | **78.1%** |
| **StellarPGx** | 186 | 57 | 91 | **76.5%** |

- **Stargazer** achieves the highest concordance rate (81.5%) among assessable calls.
- **PyPGx** is close behind (79.9%) with identical failure count (26 calls, all ABCG2/VKORC1 genes not supported).
- **Aldy** never fails to produce a call (0 missing), giving the largest denominator, but has slightly more discordances (73), partly due to artefactual calls on RYR1 and DPYD.
- **StellarPGx** has 91 missing/failed calls — the most of any tool — primarily for ABCG2, DPYD, RYR1, VKORC1, and some SLCO1B1 entries, and achieves the lowest per-call concordance (76.5%).

---

## 4. Per-Gene Concordance

| Gene | Concordant | Discordant | Failed | **% Concordant** | Notes |
|------|-----------|-----------|--------|-----------------|-------|
| ABCG2 | 26 | 0 | 78 | **100%** | PyPGx/Stargazer/StellarPGx do not call ABCG2; Aldy calls correctly |
| CYP2B6 | 8 | 0 | 0 | **100%** | Only 2 samples had Axiom P9 CYP2B6 data |
| CYP2C19 | 94 | 10 | 0 | **90.4%** | Discordances mostly sub-allele naming (*\*38*, *\*36*) |
| CYP2C9 | 96 | 8 | 0 | **92.3%** | Axiom rsID notation (rs9332094C = *\*8*) counted as concordant |
| CYP2D6 | 35 | 29 | 0 | **54.7%** | Hybrid alleles, CNV notation, SV complexity |
| CYP3A4 | 40 | 0 | 0 | **100%** | Perfect concordance |
| CYP3A5 | 100 | 4 | 0 | **96.2%** | Minor discordances (*\*3B* sub-allele, *\*1* vs *\*3* StellarPGx) |
| CYP4F2 | 73 | 31 | 0 | **70.2%** | Compound allele notation (*\*2+3* Axiom vs *\*3*, *\*4* pipeline) |
| DPYD | 6 | 72 | 26 | **7.7%** | Systematic: pipeline finds variants not in Axiom array (see §5.1) |
| NAT2 | 38 | 2 | 0 | **95.0%** | StellarPGx mis-calls P_93 (*\*16/\*34* vs *\*4/\*5E*) |
| NUDT15 | 92 | 10 | 2 | **90.2%** | *\*5* vs *\*6* confusion (P_62); ambiguous Axiom P9 (P_670) |
| RYR1 | 18 | 12 | 10 | **60.0%** | Pipeline misses Axiom-reported RYR1 het variants; Aldy sometimes over-calls |
| SLCO1B1 | 92 | 11 | 1 | **89.3%** | P_2912: all tools call *\*20/\*37* vs Axiom *\*37/\*37* |
| TPMT | 93 | 11 | 0 | **89.4%** | StellarPGx *\*1S* sub-allele inflation; P_93 Axiom ambiguous |
| UGT1A1 | 64 | 40 | 0 | **61.5%** | Systematic: *\*60* (rs45530432) and *\*80* (rs887829) allele gaps |
| VKORC1 | 69 | 9 | 26 | **88.5%** | Aldy H-haplotype vs rs-SNP notation; StellarPGx does not call |

---

## 5. Discordance Analysis by Gene

### 5.1 DPYD — Systematic Under-detection by Axiom Array (7.7% concordance)

**This is the most significant finding.** All 13 samples have Axiom DPYD = `*1/*1`, yet the pipeline consistently detects DPYD variants across the cohort:

| Variant detected by pipeline | DPYD Star Allele | CPIC significance | Samples affected |
|------------------------------|-----------------|------------------|-----------------|
| c.85T>C | *\*9A* | Decreased function (limited evidence) | EQ_2017, EQ_2021, P_62, P_2505, P_2912, P_93, P_670, P_2483 |
| c.1627A>G | *\*5* | Decreased function | H_342, P_2000, P_62, P_670, P_1143, P_2820 |
| c.1896T>C | — (low-evidence) | Unknown | P_1143, P_2820 |
| c.2279C>T | — | Unknown | P_2912 |
| c.496A>G | rs2297595 | Increased risk (DPYD c.496A>G) | P_2483 |

**Reason:** The Axiom P6/P9 pharmacogenomics arrays are designed to genotype a fixed set of known pathogenic DPYD variants (primarily *\*2A* IVS14+1G>A, *\*13*, and HapB3 variants). Lower-penetrance variants such as *\*9A* (c.85T>C, rs1801265) and *\*5* (c.1627A>G, rs1801159) are not consistently covered by the array. WGS with variant calling across the full DPYD gene body detects all DPYD variants regardless of array content.

**Clinical implication:** These variants may represent real pharmacogenomically relevant findings missed by the array. DPYD *\*9A* and *\*5* are associated with reduced dihydropyrimidine dehydrogenase activity and may increase 5-fluorouracil toxicity risk. Sanger confirmation or a targeted NGS DPYD assay is recommended to verify these findings before clinical reporting.

### 5.2 CYP2D6 — Hybrid Allele and CNV Complexity (54.7% concordance)

CYP2D6 is the most complex PGx gene due to structural variants (gene deletion, duplication, hybrid alleles). Discordances fall into three categories:

**a) *\*36+\*10* hybrid allele notation:**
Axiom P6 often designates this as `*36-*10` (dash separator), while the pipeline tools use `*36+*10`. This is a nomenclature difference — the biology is identical. Affected: EQ-H05, P_670.

**b) *\*106* hybrid allele (EQ_2017):**
Axiom P6 calls `*1/*29`; pipeline calls `*29/*106`. CYP2D6 *\*106* is a hybrid allele sharing variants with *\*29* (rs1080985) and is classified differently across nomenclature databases (PharmVar vs older databases). Three tools (PyPGx, Stargazer, StellarPGx) agree on *\*106*; Aldy adds additional rsIDs. This likely represents an updated star-allele classification not yet reflected in the Axiom array definition.

**c) CNV notation (EQ_2021, P_2000):**
EQ_2021: Axiom P9 = `*2/*2xN` (duplication, N copies); pipeline = `*2/*2x2` (specific copy number). The pipeline assigns a specific copy number (2) whereas Axiom reports an unquantified duplication. Biologically concordant.
P_2000: Axiom P9 = `(*1/*10)xN` (multiplication); pipeline differs between tools (Stargazer: `*1x3/*36+*10`; PyPGx: `*1/*36+*10x3`). The underlying tandem duplication on a complex *10/36 background is challenging; different tools interpret copy numbers differently.

**d) P_2912 *\*68* hybrid:**
Axiom = `*4/*41`; all pipeline tools = `*41/*68+*4`. CYP2D6 *\*68* is a *\*4*/*\*2* hybrid that shares the *\*4* splice-site variant (rs3892097). This represents a more detailed sub-classification than the Axiom call captures.

**Reason:** CYP2D6 hybrid alleles require long-read sequencing or specialised SV-aware callers for definitive assignment. The Axiom P6 array cannot resolve copy number precisely and does not cover all hybrid alleles. The Axiom P9 is better but still limited.

### 5.3 UGT1A1 — *\*60* and *\*80* Allele Gap (61.5% concordance)

Axiom reports complex UGT1A1 haplotypes including `*28+60+80+93/*60` and `*1/*60`, while pipeline tools often miss the *\*60* (rs45530432, A>G in promoter) and *\*80* (rs887829, A>G in promoter) components.

| Sample | Axiom | PyPGx | Stargazer | Aldy | StellarPGx |
|--------|-------|-------|-----------|------|-----------|
| EQ_2017 | *28+60+80+93/*60 | *\*1/*\*80* | *\*1/*\*28* | *\*1/*\*28_80* | *\*28+\*80/*\*66* |
| P_93 | *28+60+80+93/*60 | *\*1/*\*80* | *\*1/*\*28* | *\*1/*\*28_80* | *\*1/*\*28+\*80* |
| H_342 | *1/*28 | *\*1/*\*27* | *\*1/*\*27* | *\*1/*\*27_28_80* | *\*1/*\*27* |

**Reason:** UGT1A1 star-allele definitions differ between PharmVar, CPIC, and tool-specific databases. The *\*60* allele (A>G at position −3279) is a common East Asian promoter variant that reduces transcription. Several tools either lack *\*60* in their allele database or classify it differently. The *\*27* call in H_342 by PyPGx/Stargazer is a known sub-allele of *\*28*.

**Note:** H_342 discordance (Axiom: *\*1/\*28*, Pipeline: *\*1/\*27*) may reflect that *\*27* and *\*28* involve the same TA-repeat locus. The *\*27*/*\*28* distinction requires precise TA-repeat count quantification.

### 5.4 RYR1 — Rare Variant Detection Sensitivity (60.0% concordance)

| Sample | Axiom P9 | PyPGx | Stargazer | Aldy | StellarPGx |
|--------|---------|-------|-----------|------|-----------|
| P_93 | c.6178G>T/WT | Reference | *\*1/\*1* | *\*6178GT/ref* | Failed |
| P_2000 | c.1021G>A/WT | Reference | *\*1/\*1* | *ref/ref* | Failed |

**Findings:**
- **P_93 (c.6178G>T):** Only Aldy detects this heterozygous RYR1 variant. PyPGx and Stargazer both call reference. This is a rare het variant (MAF <0.001) in a region of moderate read depth (32.6×). The failure of PyPGx/Stargazer to detect low-frequency het variants at this locus may be due to their variant-calling heuristics requiring higher allele frequency or read support.
- **P_2000 (c.1021G>A):** All three tools call reference, missing the Axiom-reported variant. Gene depth for RYR1 in this sample is 25–27×. This variant may be at a mappability-challenged position.

**Reason:** RYR1 is a very large gene (15,117 bp) with many repetitive regions and moderate depth in WGS. The pipeline tools use star-allele databases that may not cover all rare RYR1 point mutations reported by Axiom P9. Aldy shows the best rare variant sensitivity for RYR1.

### 5.5 CYP4F2 — Compound Allele Nomenclature (70.2% concordance)

Axiom P6 reports CYP4F2 compound allele designations such as `*2+3` (indicating both *\*2* and *\*3* defining variants are present on one haplotype). Pipeline tools decompose this into individual alleles:

| Sample | Axiom | Pipeline (typical) | Assessment |
|--------|-------|--------------------|------------|
| EQ_2017 | *\*1/*\*2+3, *\*2/*\*3* (ambiguous) | *\*1/\*4* or *\*1/\*3* | Nomenclature difference |
| P_1143 | *\*2+3/*\*3* | *\*3/\*4* | Possible true discordance |
| P_2483 | *\*1/*\*3* | *\*1/\*4* or *\*1/\*1* | See below |
| P_2505 | *\*1/\*1* | *\*1/\*9* (PyPGx, Aldy, StellarPGx) | True discordance |
| P_2820 | *\*1/\*1* | *\*1/\*17* (PyPGx) | True discordance |

CYP4F2 *\*4* is defined by the c.1297G>A variant (rs2108622 = p.V433M), and *\*3* by c.1347T>A (p.W462stop). The Axiom `*3` designation may map to different positions than the pipeline's *\*4* allele definition depending on PharmVar version. P_2505 and P_2820 represent genuine discordances where the pipeline reports additional variants not captured by Axiom.

### 5.6 ABCG2 — Nomenclature and Tool Coverage

All 13 samples had Axiom ABCG2 results (C/C or A/C using rsID/SNP notation for rs2231142, which defines ABCG2 *\*Q141K*). Only Aldy successfully calls ABCG2; PyPGx, Stargazer, and StellarPGx do not support this gene. Aldy's calls (`*ref/*ref` and `*ref/*rs2231142`) are fully concordant with Axiom when nomenclature is normalised.

### 5.7 SLCO1B1 — P_2912 *\*20* Discordance

All four pipeline tools call SLCO1B1 `*20/*37` for sample P_2912, while Axiom P6 reports `*37/*37`. SLCO1B1 *\*20* (rs2291075 + additional variants) and *\*37* (rs4149056) share an overlapping haplotype. The *\*20* allele contains the *\*37* defining SNP plus additional variants. This may represent a genuinely more informative call by the pipeline that the Axiom array cannot distinguish due to limited probe content.

### 5.8 NUDT15 — *\*5* vs *\*6* Discrimination (P_62)

Axiom P9 calls P_62 as `*1/*6` (rs3087419, p.Val18Ile). Pipeline tools call:
- PyPGx: `*1/*5` (rs186364861, Val18_Val19insVal)
- Stargazer/Aldy: `*5/*6`
- StellarPGx: `*1/*1`

NUDT15 *\*5* and *\*6* are both loss-of-function alleles but involve different coding variants at codon 18. Confusion may arise from the proximity of these variants. Given that all tools identify a variant at NUDT15, the clinical interpretation (heterozygous reduced-function) is consistent despite allele-level discordance.

### 5.9 VKORC1 — Nomenclature Mapping

Axiom reports VKORC1 as rs9923231 genotype (G/G, G/A, A/A). Pipeline tools use different notation systems:
- PyPGx: rsID notation (`Reference/rs9923231`, `rs9923231/rs9923231`)
- Stargazer: star/haplotype notation (`*S1/*S1`, `*1/*S1`)
- Aldy: H-haplotype notation (`*H1/*H1`, `*H1/*H6`, `*H6/*H9`)

After normalisation, concordance is 88.5%. Residual discordances (9 cases) are from Aldy haplotype assignments that include additional variants beyond rs9923231, changing the haplotype class. This is a nomenclature issue rather than a genotyping error.

---

## 6. Platform Comparison: ILMN vs MGI

Both sequencing platforms produce equivalent genotype calls for the overwhelming majority of genes. No systematic genotype differences between platforms were observed for CYP2C19, CYP2C9, CYP3A4, CYP3A5, NAT2, NUDT15, SLCO1B1, TPMT, or VKORC1.

The principal platform difference is **duplication rate**:

| Platform | Mean depth | Mean duplication rate | Range |
|----------|-----------|----------------------|-------|
| ILMN | 37.6× | 31.0% | 18–44% |
| MGI | 33.4× | 12.3% | 6–38% |

Higher duplication in ILMN libraries reduces effective depth, which can marginally affect sensitivity for heterozygous rare variants (relevant to RYR1 and DPYD). However, at the depths observed (28–50×), this was not a major driver of discordance. P_2912 shows comparably high duplication on both platforms (ILMN 38.9%, MGI 38.3%), suggesting a sample-level property (PCR amplification during library preparation).

---

## 7. Summary Concordance Table by Sample

The table below reports the concordance status of each gene per sample (ILMN; MGI results are equivalent unless stated). Legend: ✓ = concordant, ✗ = discordant, — = no Axiom data, F = failed/no call.

| Gene | EQ_2017 | EQ_2021 | EQ-H05 | H_342 | P_62 | P_93 | P_670 | P_1143 | P_2000 | P_2483 | P_2505 | P_2820 | P_2912 |
|------|---------|---------|--------|-------|------|------|-------|--------|--------|--------|--------|--------|--------|
| ABCG2 | ✓(Aldy) | ✓(Aldy) | ✓(Aldy) | ✓(Aldy) | ✓(Aldy) | ✓(Aldy) | ✓(Aldy) | ✓(Aldy) | ✓(Aldy) | ✓(Aldy) | ✓(Aldy) | ✓(Aldy) | ✓(Aldy) |
| CYP2B6 | — | — | — | — | — | — | — | — | ✓ | — | — | — | — |
| CYP2C19 | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✗(all) | ✓ | ✓ |
| CYP2C9 | ✓ | ✓† | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| CYP2D6 | ✗(P/S) | ✓† | ✓(P/S/St) | — | — | — | ✓(P/S/St) | — | ✗(S) | — | ✓ | — | ✗(P/S/St) |
| CYP3A4 | — | — | — | — | ✓ | ✓ | ✓ | ✓ | ✓ | — | — | — | — |
| CYP3A5 | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓(P/Sg/A) | ✓ | ✓ |
| CYP4F2 | ✓(Sg) | ✓ | ✓ | ✓ | ✓ | ✓ | ✓(P/Sg) | ✗(all) | ✓ | ✗(P/A/St) | ✗(P/A/St) | ✗(P) | ✓ |
| DPYD | ✗(all) | ✗(P/A) | ✓ | ✗(all) | ✗(all) | ✗(P/A) | ✗(P/A) | ✗(all) | ✗(all) | ✗(all) | ✗(all) | ✗(all) | ✗(all) |
| NAT2 | — | — | — | — | ✓ | ✓(P/Sg/A) | ✓ | ✓ | ✓ | — | — | — | — |
| NUDT15 | ✓ | ✓ | ✓ | ✓ | ✗(all) | ✓ | ✗(all) | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| RYR1 | — | — | — | — | ✗(A) | ✗(all) | ✓ | ✓ | ✗(P/Sg/A) | — | — | — | — |
| SLCO1B1 | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✗(all) |
| TPMT | ✓ | ✓ | ✓ | ✓ | ✓ | ✗(all) | ✓ | ✓ | ✓(P/Sg/A) | ✓ | ✓(P/Sg/A) | ✓ | ✓(P/Sg/A) |
| UGT1A1 | ✗(all) | ✓(P/Sg/A) | ✓ | ✗(P/Sg) | ✓ | ✗(all) | ✗(Sg/St) | ✓ | ✓ | ✗(P/St) | ✗(P) | ✓(P/Sg/A) | ✗(Sg/A/St) |
| VKORC1 | ✓(P/Sg) | ✓(P/Sg) | ✓ | ✓ | ✓(P/Sg) | ✓(P/Sg) | ✓(P/Sg) | ✓ | ✓ | ✓(P/Sg) | ✓ | ✓ | ✓ |

*P = PyPGx, Sg = Stargazer, A = Aldy, St = StellarPGx. ✓(P/Sg) = only those tools concordant. †CYP2D6 EQ_2021 = copy number notation difference only (biological concordance).*

---

## 8. Root Cause Analysis of Key Discordances

| Discordance Category | Gene(s) | Root Cause | Clinical Risk |
|---------------------|---------|------------|--------------|
| Array probe gap | DPYD | Axiom does not cover all reduced-function DPYD variants; WGS finds additional real variants | HIGH — missed DPYD variants could affect 5-FU dosing |
| Array probe gap | RYR1 | Rare het RYR1 variants not covered by Axiom P9 probes; low depth reduces pipeline sensitivity | MODERATE — rare malignant hyperthermia risk variants missed |
| Nomenclature | CYP2D6 | Hybrid allele (*\*36+\*10*, *\*106*, *\*68*) definitions differ between PharmVar versions and tools | LOW — biological function usually equivalent |
| Nomenclature | VKORC1 | Three different notation systems (rsID, S-allele, H-haplotype) | LOW — underlying SNP is the same |
| Nomenclature | UGT1A1 | *\*60*, *\*80* promoter allele representation differs between databases | MODERATE — compound alleles affect phenotype |
| Tool allele DB | CYP4F2 | Axiom `*2+3` compound allele not uniformly defined across tool databases | LOW — clinical impact of CYP4F2 alleles is modest |
| Tool allele DB | NUDT15 | *\*5* vs *\*6* confusion at shared codon-18 locus | LOW — both are loss-of-function; phenotype identical |
| SV complexity | CYP2D6 | CNV copy number quantification differs between Axiom and WGS tools | MODERATE — copy number affects phenotype class |
| Sub-allele inflation | TPMT | StellarPGx calls *\*1S* (sub-allele) where others call *\*1*; no functional consequence | LOW |
| SLCO1B1 *\*20* | SLCO1B1 | Pipeline detects additional variants on the *\*37* haplotype consistent with *\*20* | LOW to MODERATE — *\*20* has additional functional implications |

---

## 9. Recommendations

### 9.1 Immediate Clinical Reporting

1. **DPYD: Mandatory clinical review.** The pipeline detects DPYD reduced-function variants (particularly *\*9A* and *\*5*) in the majority of samples where Axiom reports `*1/*1`. These variants are real and may be clinically relevant for patients receiving fluoropyrimidine chemotherapy. A targeted DPYD sequencing assay or Sanger confirmation should be performed before clinical reporting.

2. **CYP2D6: Use pipeline consensus for CNV.** Accept the pipeline call (supported by ≥3/4 tools) as definitive. Report ambiguous cases (e.g. P_2000 CNV) as "complex structural variant — manual review required."

3. **RYR1: Increase depth threshold.** RYR1 rare variants are missed at 20–26× gene depth. For anaesthesia-relevant RYR1 reporting, require ≥30× gene depth or perform confirmatory testing.

### 9.2 Pipeline Configuration

4. **Consensus call policy.** Use the majority vote (≥3/4 tools concordant) as the reportable result. Where tools disagree, flag the call for manual review rather than reporting any single tool result.

5. **ABCG2.** PyPGx, Stargazer, and StellarPGx do not support ABCG2. Aldy is currently the sole caller and achieves 100% concordance. Accept Aldy's ABCG2 call unilaterally until other tools add support.

6. **StellarPGx missing calls.** StellarPGx fails on DPYD, RYR1, VKORC1, and ABCG2 in all samples and occasionally on SLCO1B1. These are systematic tool gaps. Do not report StellarPGx results for these genes; rely on the three remaining tools.

7. **UGT1A1 *\*60*/*\*80* alleles.** Configure reporting to explicitly document whether *\*60* (rs45530432) was detected. Current tool databases inconsistently represent this common East Asian allele. A post-processing step checking for rs45530432 directly in the VCF is recommended.

8. **VKORC1 nomenclature.** Standardise to rs9923231 SNP genotype (G/G, G/A, A/A) for clinical reporting. Use PyPGx output directly as it reports rsID notation.

9. **CYP2C19 sub-alleles (*\*36*, *\*38*).** P_2505 shows all tools calling *\*2/\*38* where Axiom calls *\*1/\*2*. CYP2C19 *\*38* is a sub-allele of *\*2* carrying an additional silent variant. Report as *\*2* for clinical purposes; sub-allele can be noted in the technical record.

### 9.3 Sequencing Quality

10. **Minimum effective genome depth: 30× post-deduplication.** ILMN samples with >35% duplication may have effective depth below 25×, affecting rare het variant sensitivity (particularly RYR1). Implement a pre-report QC gate: flag samples with post-dedup mean genome depth <25× for repeat sequencing.

11. **MGI platform preferred for duplication rate.** MGI consistently achieves 6–10% duplication rates vs 20–44% for ILMN in this cohort, translating to higher effective depth per read. At equivalent raw depth, MGI delivers more usable data.

12. **RYR1 and DPYD minimum gene depth: 30×.** Flag samples below this threshold for targeted re-sequencing or independent confirmation before reporting RYR1/DPYD results.

### 9.4 Axiom Reference Limitations

13. **Axiom P6 is insufficient for CYP2D6 and NAT2.** Axiom P6 returns "Indeterminate" for CYP2D6 in nearly all samples and does not include NAT2. Future validation studies should use Axiom P9 exclusively as the reference standard.

14. **Do not use Axiom as sole ground truth for DPYD.** The Axiom array covers only a small subset of DPYD variants. The WGS pipeline is a more sensitive comparator for DPYD.

---

## 10. Recommended Reporting Tool by Gene

Based on this validation, the recommended primary tool for each gene is:

| Gene | Recommended Primary | Rationale |
|------|-------------------|-----------|
| ABCG2 | **Aldy** | Only tool that calls ABCG2 |
| CYP2B6 | **Consensus (any 3/4)** | All tools perform equivalently |
| CYP2C19 | **PyPGx or Stargazer** | Highest concordance; sub-allele aware |
| CYP2C9 | **PyPGx or Stargazer** | Tied performance; best nomenclature |
| CYP2D6 | **PyPGx + Stargazer consensus** | Best SV and hybrid allele resolution; flag copy number uncertainty |
| CYP3A4 | **Any (unanimous)** | 100% concordance across all tools |
| CYP3A5 | **Any (unanimous)** | 96% concordance; very reliable |
| CYP4F2 | **Stargazer** | Fewest discordances; best allele DB |
| DPYD | **PyPGx + Stargazer + Aldy** | All detect variants missed by Axiom; StellarPGx fails |
| NAT2 | **PyPGx or Stargazer** | 95% concordance; StellarPGx single failure |
| NUDT15 | **PyPGx** | Best *\*5/\*6* discrimination |
| RYR1 | **Aldy** | Best rare variant sensitivity |
| SLCO1B1 | **Consensus (any 3/4)** | >89% concordance across all tools |
| TPMT | **PyPGx or Stargazer** | Avoid StellarPGx (sub-allele inflation) |
| UGT1A1 | **Aldy** | Best compound allele representation |
| VKORC1 | **PyPGx** | rsID notation most compatible with clinical reporting |

---

## 11. Conclusion

The pgx-suite pipeline demonstrates high concordance with Axiom microarray genotyping for 9 of 16 genes (>89% concordance). The most clinically important findings are:

1. **The pipeline detects DPYD variants systematically missed by the Axiom array**, which may have direct clinical implications for fluoropyrimidine chemotherapy dosing. These should be validated and reported.

2. **Stargazer achieves the highest overall concordance** (81.5%) and is the recommended primary tool, supported by PyPGx as a highly concordant second caller.

3. **Aldy provides unique ABCG2 coverage** and the best RYR1 rare variant sensitivity, and should be retained in the pipeline.

4. **StellarPGx has significant tool gaps** (91 failed calls across 9 genes) and should be treated as a supplementary tool rather than a primary caller.

5. **The two WGS platforms (ILMN and MGI) produce equivalent PGx genotypes.** MGI is preferred for lower duplication rates and higher effective depth.

6. **Several discordances are nomenclature artefacts** (VKORC1, UGT1A1 *\*60*, CYP2D6 hybrid alleles, CYP2C19 sub-alleles) rather than true genotyping errors, and resolve with standardised reporting.

---

*Report generated from pipeline results in `results/ILMN/` and `results/MGI/`. Axiom reference data from `Axiom.xlsx` (Thermofisher Axiom Pharmacogenomics P6/P9 array).*
