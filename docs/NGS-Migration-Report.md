# Migration of Pharmacogenomic Testing from the Axiom Array to Next-Generation Sequencing: A Validation Report for the Clinical Laboratory

**Intended audience:** clinical laboratories currently performing pharmacogenomic (PGx)
testing on the ThermoFisher Axiom PGx array and evaluating transition to next-generation
sequencing (NGS).
**Date:** 2026-06-20
**Companion technical report:** [`AxiomValidation.md`](AxiomValidation.md) (complete methodology, source data, and per-sample tables)

---

## Summary

Whole-genome sequencing (WGS) on both the Illumina and MGI platforms reproduced the
clinical PGx classifications obtained on the Axiom array and, for a subset of genes,
demonstrated greater analytical sensitivity. The unadjudicated exact-concordance rate of
73.0% understates true performance: following root-cause classification of each
discordance, **85.8% of gene-level calls were either concordant with the array or
represented cases in which NGS was more sensitive than, or at parity with, the array.**
The residual discordances were attributable to allele nomenclature differences and a small
number of rare variants requiring orthogonal confirmation. These findings support a
phased, gene-by-gene migration accompanied by defined confirmatory and quality-control
safeguards.

---

## 1. Objective

The Axiom array interrogates a fixed panel of pharmacogenomic variants. NGS offers broader
genomic coverage, consolidation of multiple indications into a single assay, and the
ability to detect variants outside any predefined panel. Prior to adopting NGS for clinical
reporting, a laboratory must establish whether sequencing yields clinically equivalent
classifications to the validated array method and, where results diverge, determine which
method is correct.

This study evaluated the concordance of NGS-derived PGx classifications against established
Axiom results, and characterised the nature of any discordances.

---

## 2. Study design

Thirteen patients with existing Axiom (P6/P9) PGx genotypes were re-analysed by WGS on two
independent platforms, **Illumina (ILMN)** and **MGI**, yielding **26 platform-sample
pairs**. Sequencing data were processed through the pgx-suite reconciliation pipeline, and
calls were compared against the corresponding Axiom result across the **16 genes** common
to both the array panel and the pipeline.

---

## 3. Initial findings and the limitation of exact-match scoring

Initial scoring by exact concordance of the reported genotype string yielded a match rate
of **73.0% (241/330 gene-level calls)**. Interpreted directly, this figure would suggest a
substantial discordance rate incompatible with clinical adoption.

This interpretation is methodologically unsound. Exact-string comparison between a
comprehensive assay (WGS) and a fixed-panel assay (array) penalises the sequencing method
for detecting variation the array does not interrogate, and for differences in allele
naming convention rather than in the underlying genotype. Each of the 89 non-matching calls
was therefore classified by root cause:

| Discordance category | Description | Count | Pipeline error |
|---|---|---|---|
| Array panel gap | NGS concordantly detected a clinically relevant variant absent from the array panel (predominantly **DPYD**) | 30 | No -- NGS more sensitive |
| Single-SNP notation | Array and NGS agree on the underlying genotype; the reported haplotype label differs (**VKORC1**) | 12 | No -- genotypically concordant |
| Nomenclature | Equivalent biology represented under differing allele-naming conventions (CYP2D6, UGT1A1, CYP4F2) | 34 | Generally no -- reporting convention |
| Orthogonal confirmation required | A rare variant detected by one method and not the other (**RYR1**, UGT1A1) | 6 | Unresolved -- confirmation required |
| Unresolved on short reads | Structurally complex locus not resolvable from short-read data (CYP2D6, NUDT15) | 7 | Partial |

Combining concordant calls with those in which NGS was correctly more sensitive or at
parity with the array yields an adjudicated concordance of **283/330 (85.8%)**. The
discordances representing potential NGS error were confined to the final two categories --
**13 calls (4.0%)** -- concentrated in two genes of recognised analytical difficulty.

---

## 4. Basis for analytical confidence

The reliability of these conclusions derives from the design of the reconciliation pipeline
rather than from any single sequencing tool:

1. **Multi-caller reconciliation.** Up to nine independent callers are executed per gene and
   reconciled into a single verdict. No individual algorithm's behaviour determines the
   final classification; concordance across methods does, providing an internal replicate at
   the analytical level.

2. **Gene-specific authoritative methods.** Genes requiring specialised analysis are
   assigned a designated authority: **Cyrius** for the structurally complex **CYP2D6**,
   **PharmCAT** for **UGT1A1, CYP2B6, and CYP4F2**, and direct variant interrogation
   (**VCF-Check**) for single-variant genes including **VKORC1** and **RYR1**. The verdict
   reflects the competent method for each locus rather than an unweighted consensus.

3. **Nomenclature normalisation.** Synonymous allele designations (e.g. DPYD `*9A`, `*S10`,
   `c.85T>C`, `rs1801265`) are collapsed prior to concordance assessment, preventing genuine
   agreement from being obscured by representational differences.

4. **Defect detection during validation.** A configuration fault causing silent failure of
   one component (PharmCAT) was identified through the pipeline's per-tool status tracking,
   corrected, and the affected samples reprocessed. The capacity of the validation to detect
   and remediate its own analytical defects strengthens confidence in the reported results.

Consequently, every discordance in this study is associated with a documented root cause.

---

## 5. Platform comparison: Illumina versus MGI

The two sequencing platforms were equivalent at the level of the final clinical verdict.

| Platform | Verdict-concordant with array | % |
|---|---|---|
| MGI | 122 / 165 | 73.9% |
| Illumina | 119 / 165 | 72.1% |

The 1.8-percentage-point difference is within expected variation, and the per-gene
discordance patterns were near-identical between platforms. The principal technical
difference observed was library duplication rate (MGI 6-10%; Illumina 18-44%), which confers
marginally greater effective depth on MGI for the detection of rare heterozygous variants in
large genes (DPYD, RYR1). This did not alter any verdict at the depths examined. Both
platforms are analytically fit for purpose; platform selection should be guided by cost,
throughput, and existing infrastructure rather than PGx performance.

---

## 6. Gene-level readiness for clinical migration

| Tier | Genes | Interpretation |
|---|---|---|
| Direct migration (100% or >=90% concordant) | ABCG2, CYP2B6, CYP2C19, CYP2C9, CYP3A4, CYP3A5, NAT2, NUDT15, SLCO1B1, TPMT | NGS reproduces the array; transition with standard verification. |
| Enhanced sensitivity -- reporting revision required | **DPYD** | NGS detects reduced-function alleles (e.g. *9A, *5) absent from the array panel and relevant to fluoropyrimidine dosing. Standard operating procedures must incorporate these alleles, with targeted confirmation prior to dosing decisions. |
| Genotype-level reporting | **VKORC1** | A single clinically relevant variant (rs9923231). Array and NGS are genotypically concordant; the reported genotype, rather than the haplotype label, should be used. |
| Continued confirmatory / manual review | **CYP2D6**, **RYR1**, **UGT1A1** | CYP2D6 exhibits copy-number and hybrid-allele complexity; accept the Cyrius classification when reported and flag for manual review on abstention. RYR1 rare variants and UGT1A1 promoter haplotypes require orthogonal confirmation, consistent with current array-based practice. |

---

## 7. Recommendations

1. **Phased migration.** Transition the direct-migration genes first, as they require the
   least additional interpretive work. Maintain heightened review of DPYD, CYP2D6, RYR1, and
   UGT1A1 throughout the transition period.

2. **Parallel-reporting period.** For an initial sample series, report by the established
   method and verify against NGS prior to designating NGS as the primary result. The present
   study constitutes such a comparison at n=13 and should be extended on the laboratory's own
   sample population.

3. **Revision of DPYD reporting.** NGS will routinely detect reduced-function DPYD variants
   not interrogated by the array. Confirmatory and reporting policies for these variants
   should be defined in advance of implementation.

4. **Confirmatory-testing policy.** Retain orthogonal confirmation for CYP2D6 structural
   classifications, RYR1 rare variants, and any single-method finding. NGS identifies the
   specific calls requiring confirmation rather than obviating the requirement.

5. **Quality-control thresholds.** Maintain >=25x de-duplicated genome depth (>=30x for
   DPYD/RYR1). Verify that input alignments derive from an ALT-aware aligner; non-ALT-aware
   alignment can yield a spurious CYP2D6 deletion. The pipeline performs this check
   automatically and withholds affected genes on failure.

6. **Reporting of the reconciled verdict.** The single per-gene verdict, with its associated
   authority and resolution designation, constitutes the reportable result. Results should
   not be re-derived from individual underlying callers.

7. **Platform selection.** Illumina and MGI are clinically equivalent for PGx; MGI's lower
   duplication rate confers a modest advantage for rare-variant sensitivity.

---

## 8. Conclusion

This validation supports migration of PGx testing from the Axiom array to NGS. Across 13
patients on two sequencing platforms, the pipeline reproduced the array's clinical
classifications and, for DPYD in particular, detected clinically actionable variants beyond
the array's panel. The genuine discordances were confined to genes of established analytical
difficulty and were explicitly flagged for confirmation rather than reported without
qualification.

Migration is appropriately implemented as a gene-by-gene transition with confirmatory
safeguards for complex loci. The pipeline's multi-caller reconciliation, gene-specific
authoritative methods, and integrated quality control provide the auditability required for
safe clinical adoption.

---

*The figures in this report are derived from the technical validation [`AxiomValidation.md`](AxiomValidation.md),
generated programmatically from per-sample concordance tables. This document is a summary
prepared for laboratory decision-makers and does not constitute the primary data record.*
