# PGx Star Allele Callers — Tool Reference

This document provides a concise reference for the four pharmacogenomics (PGx) star
allele callers bundled in the **pgx-suite** Docker container. It is intended to let
users understand each tool's approach, gene coverage, and limitations without having
to read four separate documentation sites.

All four tools are configured for **GRCh38 (hg38)** in this container.

---

## Quick Comparison

| | PyPGx | Stargazer | Aldy | StellarPGx |
|---|---|---|---|---|
| **Version** | 0.26.0 | 2.0.3 | 4.8.3 | 1.2.7 |
| **License** | MIT | Non-commercial academic (UW) | Non-commercial academic (IURTC) | MIT |
| **Input** | BAM + VCF | BAM + VCF | BAM | BAM |
| **Genes** | 88 | 58 | 39 | 21 |
| **SV/CN detection** | Yes (manual preprocessing) | Partial (3 paralog genes) | Yes (automatic) | Yes (automatic) |
| **Phasing** | Beagle (statistical) | Beagle (statistical) | ILP | Graphtyper |
| **Algorithm** | Bayesian / EM | Bayesian / EM | Integer Linear Programming | Genome graph (vg) |
| **WGS** | Yes | Yes | Yes | Yes (WGS only) |
| **WES / panel** | Yes | Yes | Yes (CN unreliable) | No |
| **Long reads** | No | No | Partial (`--param sam_long_reads`) | No |
| **Multi-sample** | Yes | Yes | No | No |
| **Recommended coverage** | ≥30× | ≥30× | ≥40× | ≥30× |

---

## 1. PyPGx

### Overview

PyPGx (Python PGx) is a Python framework developed by Seung-been Lee (sbslee) that
implements a complete NGS-to-diplotype pipeline for 88 pharmacogenes. It is the
broadest caller in this suite by gene count, covering many CPIC-annotated genes
beyond the CYP450 family including transporters, transferases, and ion channels.

PyPGx models allele frequencies using a Bayesian framework and uses Beagle for
statistical haplotype phasing. It is designed for both research and clinical
interpretation workflows and outputs structured result archives (`.zip`) that include
per-sample diplotypes, activity scores, and phenotype predictions aligned with CPIC
guidelines.

### Input requirements

- **BAM/CRAM** — aligned to GRCh38, sorted and indexed (`.bai`/`.crai`)
- **VCF** — gene-region variant calls (e.g. from `bcftools mpileup | bcftools call`)
- **Depth-of-coverage archive** (SV genes only) — produced by `pypgx prepare-depth-of-coverage`
- **Control statistics archive** (SV genes only) — produced by `pypgx compute-control-statistics VDR`

For non-SV genes the VCF alone is sufficient. For SV genes the depth and control
archives are required to call copy number variation.

### Genes covered (88 target genes)

**SV genes — require depth-of-coverage + VDR control statistics (13):**

| Gene | SV type |
|------|---------|
| CYP2A6 | Deletions, duplications (paralog: CYP2A7) |
| CYP2B6 | Deletions, duplications (paralog: CYP2B7) |
| CYP2D6 | Deletions, duplications, tandems (paralog: CYP2D7) |
| CYP2E1 | Copy number variation |
| CYP4F2 | Copy number variation |
| G6PD | Copy number variation |
| GSTM1 | Whole-gene deletion polymorphism |
| GSTT1 | Whole-gene deletion polymorphism (alt contig — see limitations) |
| SLC22A2 | Copy number variation |
| SULT1A1 | Copy number variation |
| UGT1A4 | Copy number variation |
| UGT2B15 | Copy number variation |
| UGT2B17 | Copy number variation |

**Non-SV genes — VCF only (75):**
ABCB1, ABCG2, ACYP2, ADRA2A, ADRB2, ANKK1, APOE, ATM, BCHE, BDNF, CACNA1S,
CFTR, COMT, CYP1A1, CYP1A2, CYP1B1, CYP2A13, CYP2C8, CYP2C9, CYP2C19, CYP2F1,
CYP2J2, CYP2R1, CYP2S1, CYP2W1, CYP3A4, CYP3A5, CYP3A7, CYP3A43, CYP4A11,
CYP4A22, CYP4B1, CYP17A1, CYP19A1, CYP26A1, DBH, DPYD, DRD2, F2, F5, GRIK1,
GRIK4, GRIN2B, GSTP1, HTR1A, HTR2A, IFNL3, IFNL4, ITGB3, ITPA, MT-RNR1, MTHFR,
NAT1, NAT2, NUDT15, OPRK1, OPRM1, POR, PTGIS, RARG, RYR1, SLC6A4, SLC15A2,
SLC28A3, SLC47A2, SLCO1B1, SLCO1B3, SLCO2B1, TBXAS1, TPMT, UGT1A1, UGT1A6,
UGT2B7, VKORC1, XPC

Control genes (used for CN normalisation, not callable themselves): EGFR, RYR1, VDR

### Output

Results are written to `<output>/results.zip`, which contains:
- `data.tsv` — per-sample diplotype, activity score, phenotype
- `<gene>.vcf.gz` — phased variant calls
- Copy number plots (PNG) when `--do-not-plot-copy-number` is not set

### Limitations

- **GSTT1 on alt contig**: The GRCh38 locus for GSTT1 falls on `chr22_KI270879v1_alt`,
  an alternate scaffold absent from most standard GRCh38 references. bcftools mpileup
  cannot generate a VCF for it from a standard reference, making GSTT1 effectively
  uncallable in most WGS pipelines.
- **SV genes require extra preprocessing**: Unlike Aldy and StellarPGx, PyPGx does not
  infer copy number automatically from the BAM. Users must explicitly run
  `prepare-depth-of-coverage` and `compute-control-statistics` before `run-ngs-pipeline`
  for any SV gene.
- **Multi-sample workflow**: PyPGx is primarily designed to process multiple samples
  jointly. Single-sample mode works but some statistical steps (inter-sample CN
  normalisation) are more reliable with ≥10 samples.
- **Beagle phasing**: Requires the 1KGP haplotype panel (pypgx-bundle, ~500 MB) to be
  mounted. Without it, phasing falls back to a less accurate method.
- **No long-read support**: Designed for Illumina short-read WGS/WES.

---

## 2. Stargazer

### Overview

Stargazer was developed by the University of Washington (Seattl Genomics Group) and
focuses on a subset of pharmacogenes with well-characterised star allele definitions.
It is notable for being one of the earliest tools to explicitly model gene copy number
variation by comparing read depth of the target gene to a control locus.

The pipeline runs in two modes depending on whether a GDF (Genotype Depth Format)
depth-profile file is provided:

- **VCF-only mode** (default): uses only variant calls. Cannot detect deletions or
  duplications. Suitable for non-SV genes.
- **SV mode** (`-c <control_gene> -g <gdf_file>`): uses intra-sample depth
  normalisation against a control gene (VDR, EGFR, or RYR1) to call copy number,
  then detects structural variant alleles. Required for CYP2A6, CYP2B6, CYP2D6.

The GDF file is generated from the BAM using Stargazer's own `--create-gdf-file` mode.

### Input requirements

- **BAM/CRAM** — aligned to GRCh38, sorted and indexed (for GDF creation)
- **VCF** — gene-region variant calls
- **GDF file** (SV genes only) — produced by `stargazer -G <name>.gdf -c vdr -B <bam> -t <gene> -a grc38`
- **Control gene** (SV genes only) — one of `vdr`, `egfr`, `ryr1` (VDR is default in this suite)

### Genes covered (58 target genes)

**Paralog genes — GDF + control gene required for full SV detection (3):**

| Gene | Paralog | Notes |
|------|---------|-------|
| CYP2A6 | CYP2A7 | Gene deletions, hybrid alleles |
| CYP2B6 | CYP2B7 | Gene deletions, hybrid alleles |
| CYP2D6 | CYP2D7 | Deletions, duplications, gene conversions |

**All other target genes (VCF-only mode is sufficient) (55):**
ABCB1, ABCG2, CACNA1S, CFTR, CYP1A1, CYP1A2, CYP1B1, CYP2A13, CYP2C8, CYP2C9,
CYP2C19, CYP2E1, CYP2F1, CYP2J2, CYP2R1, CYP2S1, CYP2W1, CYP3A4, CYP3A5,
CYP3A7, CYP3A43, CYP4A11, CYP4A22, CYP4B1, CYP4F2, CYP17A1, CYP19A1, CYP26A1,
DPYD, G6PD, GSTM1, GSTP1, IFNL3, NAT1, NAT2, NUDT15, POR, PTGIS, RYR1, SLC15A2,
SLC22A2, SLCO1B1, SLCO1B3, SLCO2B1, SULT1A1, TBXAS1, TPMT, UGT1A1, UGT1A4,
UGT2B7, UGT2B15, UGT2B17, VKORC1, XPC, 2C_CLUSTER

Control genes (for CN normalisation, not target genes): VDR, EGFR, RYR1

### Output

Results are written to the output directory and include:
- `genotype.txt` — per-sample diplotype, SV call, phased VCF path
- `report.tsv` — detailed per-sample report
- `phased.vcf` / `finalized.vcf` — phased variant calls

### Limitations

- **Only 3 genes support paralog-based SV detection**: Unlike PyPGx which models SV
  for 13 genes, Stargazer only performs paralog-aware CN calling for CYP2A6, CYP2B6,
  and CYP2D6. Other SV-containing genes (e.g. GSTM1, GSTT1, CYP4F2) are handled
  with VCF-only mode and SV alleles may be missed.
- **Two-step SV workflow**: GDF creation (`-G` flag) must be run before the main
  genotyping call. Stargazer does not auto-generate the GDF from the BAM.
- **Single sample per run**: Stargazer processes one sample at a time (no batch mode).
  Multi-sample normalisation is not supported.
- **No WES copy number calling**: Exome/panel data lacks the uniform depth needed for
  intra-sample CN normalisation.
- **Non-commercial academic licence (University of Washington)**: Must not be published
  to any public Docker registry or used in commercial workflows.
- **Python 2-style bare imports**: Stargazer cannot be installed via `pip install` in
  the standard way. In this container a wrapper script calls `python3 __main__.py`
  directly to work around the import path issue.

---

## 3. Aldy

### Overview

Aldy (Allele Decomposition for pharmacogenomics) was developed at Indiana University
and Carnegie Mellon University. It takes a fundamentally different algorithmic approach:
rather than calling variants first and then matching to known star alleles, Aldy frames
genotyping as an **Integer Linear Programming (ILP)** problem that jointly optimises
copy number, allele decomposition, and diplotype assignment in a single step.

This means Aldy:
1. Reads the BAM directly and profiles coverage across gene-specific CN regions.
2. Builds a set of candidate CN configurations using the coverage signal.
3. Solves an ILP to find the minimum-cost explanation of the observed reads in terms
   of known star alleles and their copy numbers.
4. Assigns the highest-scoring diplotype as the final call.

Because the CN and allele steps are coupled inside the ILP, **no external preprocessing
is needed for SV genes** — Aldy detects deletions, duplications, and fusions
automatically from the BAM.

### Input requirements

- **BAM/CRAM** — aligned to GRCh38, sorted and indexed
- No VCF required
- No control gene or GDF required

### Genes covered (39 genes)

All 39 genes have CN-region definitions in Aldy's YAML databases, meaning all support
some level of copy number-aware calling. Two genes have explicit paralog models
(CYP2A6 → CYP2A7; CYP2D6 → CYP2D7) and CYP2D6 additionally models tandem
arrangements.

ABCG2, CACNA1S, CFTR, COMT, CYP1A1, CYP1A2, CYP2A6, CYP2A13, CYP2B6, CYP2C8,
CYP2C9, CYP2C19, CYP2D6, CYP2E1, CYP2F1, CYP2J2, CYP2R1, CYP2S1, CYP2W1,
CYP3A4, CYP3A43, CYP3A5, CYP3A7, CYP4F2, DPYD, G6PD, GSTM1, GSTP1, IFNL3,
MT-RNR1, NAT1, NAT2, NUDT15, RYR1, SLCO1B1, TPMT, UGT1A1, UGT2B7, VKORC1

**Paralog-aware genes (explicit hybrid allele calling):** CYP2A6, CYP2D6

**Tandem-arrangement modelling:** CYP2D6 (e.g. \*36+\*10, \*68+\*4)

### Output

Results are written to a tab-separated `.aldy` file:
- One row per solution (Aldy may report multiple solutions if the ILP has ties)
- Columns: Sample, Gene, SolutionID, Major alleles, Minor alleles, Score, Coverage, Phenotype

### Limitations

- **Coverage requirement**: Minimum 40× mean coverage recommended for reliable CN calls.
  Below 20× the CN signal becomes too noisy for accurate SV detection.
- **WES cannot call CN changes**: Whole-exome and gene panel data have non-uniform
  capture efficiency that breaks Aldy's depth-based CN model. With exome data, Aldy
  assumes exactly two gene copies and cannot detect deletions or duplications.
- **Long reads (partial support)**: PacBio and Nanopore BAMs can be used with
  `--param sam_long_reads=true`, but SV calling accuracy on long reads is not fully
  validated.
- **No multi-sample processing**: Aldy processes one BAM at a time.
- **Gene database coverage**: At 39 genes, Aldy covers fewer genes than PyPGx (88) or
  Stargazer (58). Genes such as ABCB1, APOE, DPYD (full panel), and many
  neuropharmacology genes are absent.
- **Non-commercial academic licence (IURTC)**: Must not be used in commercial workflows
  or published to public Docker registries.

---

## 4. StellarPGx

### Overview

StellarPGx was developed by David Twesigomwe and collaborators at the Sydney Brenner
Institute for Molecular Bioscience (SBIMB), University of the Witwatersrand. It is
the only tool in this suite built on a **genome graph** framework.

Instead of aligning reads to the linear reference and calling variants, StellarPGx:
1. Constructs a variation graph of the target gene using **vg** (variation graph toolkit)
   and **graphtyper**, incorporating known star allele haplotypes as alternative paths.
2. Calls variants by aligning reads to the graph, which handles highly repetitive
   regions and paralogs more robustly than linear alignment.
3. Assigns star alleles by matching phased graph-level haplotypes to the star allele
   database.

This approach makes StellarPGx particularly well-suited to **highly polymorphic genes**
with complex repeat structures (CYP2D6, CYP2A6, CYP2B6) where linear reference
alignment is error-prone. SV detection is automatic and intrinsic to the graphtyper
variant-calling step — no GDF or control gene preprocessing is required.

StellarPGx is implemented as a **Nextflow pipeline** and runs its core tools inside
an Apptainer (Singularity) container (`stellarpgx-dev.sif`). This means the pipeline
requires either a SLURM cluster or a Docker container run with `--privileged` to allow
Apptainer to unpack its overlay filesystem.

### Input requirements

- **BAM/CRAM** — aligned to GRCh38, sorted and indexed
- **GRCh38 reference FASTA** — mounted at `/pgx/ref/hg38.fa`
- **StellarPGx repository** — mounted at `/pgx/stellarpgx` (contains gene databases and scripts)
- **SIF container** — `stellarpgx-dev.sif` mounted at `/pgx/containers/stellarpgx-dev.sif`
- **`--privileged`** Docker flag (for Apptainer inside Docker)

### Genes covered (21 genes)

StellarPGx has the narrowest gene coverage of the four tools, focused on pharmacogenes
where graph-based calling provides the most benefit.

**CYP450 genes (12):**
CYP1A1, CYP1A2, CYP2A6, CYP2B6, CYP2C8, CYP2C9, CYP2C19, CYP2D6, CYP2E1,
CYP3A4, CYP3A5, CYP4F2

**Other pharmacogenes (9):**
CYPOR (POR), GSTM1, GSTT1, NAT1, NAT2, NUDT15, SLCO1B1, TPMT, UGT1A1

### Supported reference builds

StellarPGx supports **GRCh38 (hg38)**, **GRCh37 (b37)**, and **hg19**. This container
is configured for GRCh38 only.

### Output

Results are written to `<out_dir>/<gene>/alleles/` (one subdirectory per gene), with:
- `<sample>.alleles` — called star alleles (one per line)
- Summary files in `<out_dir>/<gene>/star_allele_calls/`

### Limitations

- **WGS only**: StellarPGx explicitly requires short-read high-coverage WGS. Exome,
  gene panel, and long-read inputs are not supported.
- **Fewest genes (21)**: Covers only a fraction of the genes supported by PyPGx (88)
  or Stargazer (58). Many clinically important genes (DPYD, TPMT for full panel,
  VKORC1, etc.) are absent.
- **Requires `--privileged`**: Running Apptainer inside Docker requires elevated
  permissions. This limits deployment to local workstations or HPC environments where
  `--privileged` or `SYS_ADMIN` capability is available. Managed Kubernetes or
  restricted cloud environments typically prohibit this.
- **Large runtime dependency**: Requires the `stellarpgx-dev.sif` Apptainer image
  (~31 MB compressed) and the StellarPGx repository (~100 MB) to be volume-mounted
  at runtime. The SIF contains graphtyper, bcftools, samtools, and the star allele
  caller scripts.
- **Not suitable for single-gene panels or small BAMs**: The graph construction and
  graphtyper variant calling steps assume sufficient read depth and read length for
  reliable graph-based alignment.

---

## Gene Support Matrix

The following table shows which genes are supported by each tool (in this suite's
GRCh38 configuration). Genes only in the outer tools (PyPGx/Stargazer) are omitted
for brevity; the full per-tool lists are in the sections above.

| Gene | PyPGx | Stargazer | Aldy | StellarPGx | SV? |
|------|:-----:|:---------:|:----:|:----------:|:---:|
| CYP1A1 | — | ✓ | ✓ | ✓ | — |
| CYP1A2 | — | ✓ | ✓ | ✓ | — |
| CYP2A6 | ✓ | ✓ | ✓ | ✓ | ✓ (paralog) |
| CYP2B6 | ✓ | ✓ | ✓ | ✓ | ✓ (paralog) |
| CYP2C8 | ✓ | ✓ | ✓ | ✓ | — |
| CYP2C9 | ✓ | ✓ | ✓ | ✓ | — |
| CYP2C19 | ✓ | ✓ | ✓ | ✓ | — |
| CYP2D6 | ✓ | ✓ | ✓ | ✓ | ✓ (paralog + tandem) |
| CYP2E1 | ✓ | ✓ | ✓ | ✓ | ✓ (CN) |
| CYP3A4 | ✓ | ✓ | ✓ | ✓ | — |
| CYP3A5 | ✓ | ✓ | ✓ | ✓ | — |
| CYP4F2 | ✓ | ✓ | ✓ | ✓ | ✓ (CN) |
| SLCO1B1 | ✓ | ✓ | ✓ | ✓ | — |
| NUDT15 | ✓ | — | ✓ | ✓ | — |
| TPMT | ✓ | — | ✓ | ✓ | — |
| UGT1A1 | ✓ | ✓ | ✓ | ✓ | — |
| NAT1 | ✓ | ✓ | ✓ | ✓ | — |
| NAT2 | ✓ | ✓ | ✓ | ✓ | — |
| GSTM1 | ✓ | ✓ | ✓ | ✓ | ✓ (deletion) |
| GSTT1 | ✓† | — | — | ✓ | ✓ (deletion) |
| G6PD | ✓ | ✓ | ✓ | — | ✓ (CN) |
| VKORC1 | ✓ | ✓ | ✓ | — | — |
| DPYD | ✓ | — | ✓ | — | — |
| CYPOR/POR | — | ✓ | — | ✓ | — |
| IFNL3 | ✓ | ✓ | ✓ | — | — |
| RYR1 | ✓ | ✓ | ✓ | — | — |

†GSTT1 is on an alternate GRCh38 contig — bcftools cannot generate a VCF for it from
a standard reference; PyPGx depth-of-coverage preprocessing may also fail.

---

## SV Handling Summary

Structural variant (SV) alleles — gene deletions, duplications, and hybrid alleles —
are the most challenging class to call accurately. Each tool addresses them differently:

| Tool | SV approach | Genes | User action required |
|------|------------|-------|---------------------|
| **PyPGx** | Depth-of-coverage model + VDR normalisation | 13 genes | Run `prepare-depth-of-coverage` and `compute-control-statistics VDR` before `run-ngs-pipeline` |
| **Stargazer** | Intra-sample depth ratio to control gene (VDR/EGFR/RYR1) via GDF file | 3 paralog genes | Run `stargazer -G <gene>.gdf -c vdr -B <bam>` before main genotyping call |
| **Aldy** | Integrated ILP: CN inference from BAM is part of the core solver | All 39 genes | None — automatic |
| **StellarPGx** | Graph-based variant calling with graphtyper; SV is part of variant graph | All 21 genes | None — automatic |

---

## Output Field Reference (17 × 4 Table)

The table below maps equivalent output fields across all four tools. Field names use
each tool's native terminology; `—` means the field is not reported by that tool.

**Primary output files:**
- **PyPGx** → `<output>/results.zip` → extract `data.tsv`
- **Stargazer** → `<output>/report.tsv` (summary) + `<output>/genotype-calls.tsv` (full detail)
- **Aldy** → `<output>/<GENE>.aldy`
- **StellarPGx** → `<out_dir>/<gene>/star_allele_calls/<sample>.alleles`

| # | Output Field | PyPGx (`data.tsv`) | Stargazer (`report.tsv` / `genotype-calls.tsv`) | Aldy (`.aldy`) | StellarPGx (`.alleles`) |
|---|---|---|---|---|---|
| 1 | **Sample ID** | `Sample` | `name` (genotype-calls) / `Sample` (report) | `Sample` | Filename stem (e.g. `HG03130`) |
| 2 | **Gene** | Implicit (one run per gene) | `Gene` (report) | `Gene` | Header line in file |
| 3 | **Diplotype** | `Genotype` (e.g. `*2/*4`) | `Diplotype` (report) / `hap1_main`+`hap2_main` (genotype-calls) | `Major` (e.g. `*2+*4`) | `Result` (e.g. `*17/*29`) |
| 4 | **Haplotype 1** | `Haplotype1` (e.g. `*2;`) | `hap1_main` | First allele in `Major` | First allele in `Result` |
| 5 | **Haplotype 2** | `Haplotype2` (e.g. `*4;*10;*74;`) | `hap2_main` | Second allele in `Major` | Second allele in `Result` |
| 6 | **Sub-alleles / tag variants** | Semicolon list in `Haplotype1`/`Haplotype2` (e.g. `*4;*10;*74;`) | `hap1_main_core`+`hap1_main_tag`; `hap2_main_core`+`hap2_main_tag` | `Minor` column | `Candidate alleles` (raw haplotype codes) |
| 7 | **Alternative diplotypes** | `AlternativePhase` (semicolon-separated list) | `dip_cand` / `May also be` (report) | Multiple `SolutionID` rows (each row is a valid alternative) | — |
| 8 | **Phenotype** | `Phenotype` (CPIC class, e.g. `Intermediate Metabolizer`) | `Phenotype` (e.g. `intermediate_metabolizer`) | `Status` (CPIC class) | `Metaboliser status` (e.g. `Intermediate metaboliser (IM)`) |
| 9 | **Activity score** | `ActivityScore` (float) | `Score` (report) | — | `Activity score` (float) |
| 10 | **SV / CNV type** | `CNV` (e.g. `Normal`, `Deletion`, `Duplication`) | `dip_sv` / `hap1_sv` / `hap2_sv` (e.g. `no_sv`, `dup`) | Implicit in `Major` allele name + `Copy` column | `Initially computed CN` (integer) |
| 11 | **Copy number** | Derived from `CNV` field | Derived from `dip_sv` | `Copy` column (integer copies per allele) | `Initially computed CN` |
| 12 | **Supporting variants** | `VariantData` (`allele:pos-ref-alt:AF;…`) | `hap1_main_core`/`hap2_main_core` (pos, AD, AF, consequence) | `Location` + `Coverage` columns | `Sample core variants` (`pos~ref>alt~GT`) |
| 13 | **Functional effect** | Embedded in `VariantData` | Embedded in `hap1_main_core` (e.g. `splice_acceptor_variant:high_impact`) | `Effect` + `Type` columns | — |
| 14 | **dbSNP rsID** | Embedded in `VariantData` (where available) | — | `dbSNP` column | — |
| 15 | **Allele score / confidence** | — | `dip_score`, `hap1_score`, `hap2_score` | `SolutionID` rank (lower = better fit) | — |
| 16 | **Mean allele fraction** | Per-variant AF embedded in `VariantData` | `hap1_af_mean_gene` / `hap2_af_mean_gene` / `hap1_af_mean_main` / `hap2_af_mean_main` | `Coverage` (read depth per variant position) | — |
| 17 | **Phasing method** | Beagle statistical phasing (1KGP panel) | Beagle; `BEAGLE imputed` flag (report) + `ssr` short-sequence-repeat marker | ILP joint optimisation (phasing is implicit in solution) | Graph-based variant calling (graphtyper) |

### Example values for CYP2D6

The following shows the same CYP2D6 result expressed in each tool's output format:

| Field | PyPGx | Stargazer | Aldy | StellarPGx |
|---|---|---|---|---|
| Sample | `SPL250903` | `021ab129bb` | `SPL250903` | `HG03130` |
| Diplotype | `*2/*4` | `*2/*4` | `*2+*4` | `*17/*29` |
| Haplotype 1 | `*2;` | `*2` | `*2` | `*17` |
| Haplotype 2 | `*4;*10;*74;` | `*4` | `*4` | `*29` |
| Phenotype | `Intermediate Metabolizer` | `intermediate_metabolizer` | `Intermediate Metabolizer` | `Intermediate metaboliser (IM)` |
| Activity score | `1.0` | `1.0` | — | `1.0` |
| CNV / CN | `Normal` | `no_sv,no_sv` | `Copy: 1` per allele | `2` |
| Supporting variants | `*4:22-42128945-C-T:0.462;…` | `<42128945:C>T:22/48:0.46:…>` | `42128945` + coverage | `42127608~C>T~0/1;…` |

---

## Algorithmic Deep Dive

### Why do SV genes need special handling in PyPGx and Stargazer?

CYP2D6 (and similarly CYP2A6 and CYP2B6) has a near-identical pseudogene (CYP2D7)
~5 kb away. Standard short-read aligners cannot reliably distinguish reads originating
from the gene vs the pseudogene, leading to incorrect variant calls if CN is not
explicitly modelled. Both PyPGx and Stargazer solve this by normalising the observed
read depth at the target gene against a control locus on a different chromosome (VDR
on chr12), then inferring copy number from the depth ratio. This step must be performed
before star allele assignment.

Aldy and StellarPGx integrate this step into their core algorithms so users do not
need to run it separately.

### Why does Aldy use ILP?

Standard variant-calling pipelines identify heterozygous SNVs and assign the two
alleles to haplotypes, then look up which star allele each haplotype matches. For
genes with high copy number (e.g. CYP2D6 with 3 or 4 copies), this approach
breaks because the allele fraction of each SNV depends on how many copies carry it.
Aldy avoids this problem by framing the question as: "what combination of known star
alleles, at what copy numbers, best explains the observed allele fractions across all
variants simultaneously?" — which is naturally expressed as an ILP.

### Why does StellarPGx use a genome graph?

In regions of high sequence similarity (repeats, paralogs, segmental duplications)
linear reference alignment produces systematic artefacts: reads from the pseudogene
may map to the gene, and vice versa. A variation graph represents the gene and its
known variants as a directed graph of sequences, and aligning reads to the graph
rather than a single linear sequence is more robust to this ambiguity. Graphtyper
then calls variants directly on the graph-aligned reads.

---

## References

| Tool | Publication |
|------|-------------|
| **PyPGx** | Lee S-B et al. (2022) *Pharmacogenomics J.* doi:10.1038/s41397-021-00257-5 |
| **Stargazer** | Lee S-B et al. (2019) *Genet. Med.* doi:10.1038/s41436-018-0267-1 |
| **Aldy** | Numanagić I et al. (2018) *Nat. Methods* doi:10.1038/s41592-018-0078-y; v4 update 2023 |
| **StellarPGx** | Twesigomwe D et al. (2021) *Clin. Pharmacol. Ther.* doi:10.1002/cpt.2173 |

---

## Licensing Notice

| Tool | Licence | Restriction |
|------|---------|-------------|
| PyPGx | MIT | Open use |
| Stargazer | Non-commercial academic (University of Washington) | **No commercial use; do not push to public registries** |
| Aldy | Non-commercial academic (IURTC) | **No commercial use; do not push to public registries** |
| StellarPGx | MIT | Open use |

The `pgx-suite` Docker image must not be uploaded to Docker Hub, GHCR, or any other
public container registry because it bundles Stargazer and Aldy source code.
