# PGx Star Allele Callers — Tool Reference

This document provides a concise reference for the five pharmacogenomics (PGx)
callers bundled in the **pgx-suite** Docker container. Four tools handle star allele
genotyping across the pharmacogenome; one (OptiType) performs HLA class I typing for
HLA-A and HLA-B. It is intended to let users understand each tool's approach, gene
coverage, and limitations without having to read five separate documentation sites.

All tools are configured for **GRCh38 (hg38)** in this container.

---

## Quick Comparison

| | PyPGx | Stargazer | Aldy | StellarPGx | OptiType |
|---|---|---|---|---|---|
| **Version** | 0.26.0 | 2.0.3 | 4.8.3 | 1.2.7 | 1.3.5 |
| **License** | MIT | Non-commercial academic (UW) | Non-commercial academic (IURTC) | MIT | MIT |
| **Input** | BAM + VCF | BAM + VCF | BAM | BAM | BAM (MHC reads) |
| **Genes** | 88 | 58 | 39 | 21 | HLA-A, HLA-B, HLA-C |
| **SV/CN detection** | Yes (manual preprocessing) | Partial (3 paralog genes) | Yes (automatic) | Yes (automatic) | N/A |
| **Phasing** | Beagle (statistical) | Beagle (statistical) | ILP | Graphtyper | ILP |
| **Algorithm** | Bayesian / EM | Bayesian / EM | Integer Linear Programming | Genome graph (vg) | ILP on HLA graph |
| **WGS** | Yes | Yes | Yes | Yes (WGS only) | Yes |
| **WES / panel** | Yes | Yes | Yes (CN unreliable) | No | Yes (if MHC captured) |
| **Long reads** | No | No | Partial (`--param sam_long_reads`) | No | No |
| **Multi-sample** | Yes | Yes | No | No | No |
| **Recommended coverage** | ≥30× | ≥30× | ≥40× | ≥30× | ≥30× |

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

## 5. OptiType

### Overview

OptiType is a precision HLA class I genotyper developed at the Max Planck Institute
for Informatics (Sven Rahmann group). It does not use the conventional star allele
model. Instead, it reconstructs HLA allele sequences at **4-digit (2-field) resolution**
by mapping reads against the IMGT/HLA reference database and solving an Integer Linear
Programming (ILP) problem to find the pair of HLA alleles that best explains the
observed read pileup.

In this suite OptiType is used exclusively for **HLA-A** and **HLA-B** genotyping,
which have CPIC Level A clinical guidelines (e.g. HLA-B\*57:01 → abacavir
hypersensitivity; HLA-B\*15:02 → carbamazepine SJS/TEN risk). HLA-C is also called
internally but is not reported in the comparison output.

OptiType is invoked via Apptainer (Singularity) SIF image. It runs as a self-contained
pipeline: reads are re-mapped against the IMGT/HLA reference inside the container, then
ILP genotyping is performed.

### How it is integrated in pgx-suite

`pgx-hla.sh` handles the full workflow:

1. **MHC read extraction** — `samtools view` extracts reads from the primary MHC region
   (`chr6:28,510,020-33,480,577`) plus any chr6 unmapped reads. These are written to a
   temporary FASTQ pair via `samtools fastq`.
2. **OptiType Apptainer call** — the SIF at `/pgx/containers/optitype.sif` is invoked
   with `--rna` disabled and `--dna` mode. The two FASTQ files are passed via
   `--input`.
3. **Output** — OptiType writes a TSV (`*_result.tsv`) and a PDF coverage plot to the
   output directory. The TSV has columns `A1`, `A2`, `B1`, `B2`, `C1`, `C2`,
   `Objective` where allele values are in IMGT format (e.g. `A*01:01`, `B*57:01`).
4. **Parsing in `pgx-compare.py`** — `parse_optitype()` reads `A1`/`A2` as
   `HLA-A*{value}` and `B1`/`B2` as `HLA-B*{value}`, yielding standard CPIC-compatible
   allele strings (e.g. `HLA-B*57:01`).

### Input requirements

- **BAM/CRAM** — aligned to GRCh38, sorted and indexed
- **OptiType SIF** — pre-pulled Apptainer image at `/pgx/containers/optitype.sif`
  ```
  apptainer pull --name optitype.sif \
    docker://quay.io/biocontainers/optitype:1.3.5--hdfd78af_1
  ```
- **Apptainer** — available inside the Docker container (requires `--privileged`)
- ≥30× mean depth across the MHC region is recommended; lower coverage may reduce
  phasing confidence at the 4th digit

### Genes covered

| Gene | CPIC Level | Key clinical association |
|------|-----------|--------------------------|
| HLA-A | A | HLA-A\*31:01 → carbamazepine hypersensitivity (European ancestry) |
| HLA-B | A | HLA-B\*57:01 → abacavir; HLA-B\*15:02 → carbamazepine (Asian ancestry); HLA-B\*58:01 → allopurinol |
| HLA-C | — | Reported internally; not forwarded to pgx-compare.py |

### Output

The per-gene output directory (`<output>/HLA-A/optitype/` and
`<output>/HLA-B/optitype/`) contains:

- `*_result.tsv` — two-column allele calls (`A1`, `A2` for HLA-A; `B1`, `B2` for
  HLA-B) and an ILP objective score
- `*_coverage_plot.pdf` — per-allele read coverage visualisation
- `hla_reads_1.fastq` / `hla_reads_2.fastq` — extracted MHC reads (intermediate;
  may be large)

The `pgx-compare.py` parser emits a single comparison row per gene with the diplotype
formatted as `HLA-A*01:01/HLA-A*02:01` (slash-separated ordered pair).

### Limitations

- **4-digit resolution only**: OptiType reports 2-field (4-digit) HLA alleles.
  Higher-resolution typing (6- or 8-digit) requires a different tool (e.g. HLA-HD,
  arcasHLA).
- **Class I only**: HLA-DQB1, HLA-DRB1 and other class II loci are outside OptiType's
  scope. Use HLA-HD or another class II typer if those are required.
- **No phasing across haplotypes**: OptiType calls a diplotype but does not phase HLA
  alleles onto broader chromosomal haplotypes.
- **MHC read extraction is slow on large BAMs**: `samtools view` over the full MHC
  window (~5 Mb) on a 167 GB WGS BAM takes ~30–60 sec. This is handled by `pgx-hla.sh`
  before OptiType is invoked.
- **Requires `--privileged` Docker flag**: Apptainer inside Docker requires elevated
  privileges (same constraint as StellarPGx).
- **SIF must be pulled separately**: The OptiType SIF is not baked into the Docker
  image. It must be pre-pulled to `/pgx/containers/optitype.sif`.

---

## Gene Support Matrix

The following table shows which genes are supported by each tool (in this suite's
GRCh38 configuration). Genes only in the outer tools (PyPGx/Stargazer) are omitted
for brevity; the full per-tool lists are in the sections above.

| Gene | PyPGx | Stargazer | Aldy | StellarPGx | OptiType | SV? |
|------|:-----:|:---------:|:----:|:----------:|:--------:|:---:|
| ABCG2 | ✓ | — | ✓ | ✓ | — | — |
| CYP1A1 | ✓ | ✓ | ✓ | ✓ | — | — |
| CYP1A2 | ✓ | ✓ | ✓ | ✓ | — | — |
| CYP2A6 | ✓ | ✓ | ✓ | ✓ | — | ✓ (paralog) |
| CYP2B6 | ✓ | ✓ | ✓ | ✓ | — | ✓ (paralog) |
| CYP2C8 | ✓ | ✓ | ✓ | ✓ | — | — |
| CYP2C9 | ✓ | ✓ | ✓ | ✓ | — | — |
| CYP2C19 | ✓ | ✓ | ✓ | ✓ | — | — |
| CYP2D6 | ✓ | ✓ | ✓ | ✓ | — | ✓ (paralog + tandem) |
| CYP2E1 | ✓ | ✓ | ✓ | ✓ | — | ✓ (CN) |
| CYP3A4 | ✓ | ✓ | ✓ | ✓ | — | — |
| CYP3A5 | ✓ | ✓ | ✓ | ✓ | — | — |
| CYP4F2 | ✓ | ✓ | ✓ | ✓ | — | ✓ (CN) |
| DPYD | ✓ | ✓ | ✓ | — | — | — |
| G6PD | ✓ | ✓ | ✓ | — | — | ✓ (CN) |
| GSTM1 | ✓ | ✓ | ✓ | ✓ | — | ✓ (deletion) |
| GSTT1 | ✓† | — | — | ✓ | — | ✓ (deletion) |
| HLA-A | — | — | — | — | ✓ | — |
| HLA-B | — | — | — | — | ✓ | — |
| IFNL3 | ✓ | ✓ | ✓ | — | — | — |
| NAT1 | ✓ | ✓ | ✓ | ✓ | — | — |
| NAT2 | ✓ | ✓ | ✓ | ✓ | — | — |
| NUDT15 | ✓ | ✓ | ✓ | ✓ | — | — |
| CYPOR/POR | ✓ | ✓ | — | ✓ | — | — |
| RYR1 | ✓ | ✓ | ✓ | — | — | — |
| SLCO1B1 | ✓ | ✓ | ✓ | ✓ | — | — |
| TPMT | ✓ | ✓ | ✓ | ✓ | — | — |
| UGT1A1 | ✓ | ✓ | ✓ | ✓ | — | — |
| VKORC1 | ✓ | ✓ | ✓ | — | — | — |

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

### Example values for CYP2D6 — NA24385 (SPL250903104304-2), actual run

All values taken directly from each tool's native output file. `—` = field not present in output.

| # | Field | PyPGx (`data.tsv`) | Stargazer (`genotype-calls.tsv`) | Aldy (`CYP2D6.aldy`) | StellarPGx (`*.alleles`) |
|---|---|---|---|---|---|
| 1 | **Sample ID** | `SPL250903104304-2` | `SPL250903104304-2` | `T7_NA24385` | `T7_NA24385.bwa.sortdup.bqsr` (filename stem) |
| 2 | **Gene** | implicit — one run per gene | implicit | `CYP2D6` (`Gene` col) | `CYP2D6` (file header) |
| 3 | **Diplotype** | `*2/*4` (`Genotype` col) | `*2/*4` (`hap1_main` + `hap2_main`) | `*2/*4` (`Major` col) | `*2/*4` (`Result:` section) |
| 4 | **Haplotype 1** | `*2;` (`Haplotype1` col) | `*2` (`hap1_main`) | `*2` (first token of `Major`) | `*2` (first allele of `Result`) |
| 5 | **Haplotype 2** | `*4;*10;*74;` (`Haplotype2` — all sub-alleles on hap2) | `*4` (`hap2_main`) | `*4` (second token of `Major`) | `*4` (second allele of `Result`) |
| 6 | **Sub-alleles / tag variants** | `Haplotype2`: `*4;*10;*74;` (hap2 phased sub-alleles) | hap1 core: `42126611:C>G (S486T), 42127941:G>A (R296C)` · hap2 core: `42128945:C>T (splice)` · hap2 tag: `42130692:G>A (P34S)` | `Minor`: `2.015;4.015` (suballele with tag SNVs) | `Candidate alleles: ['2.v1_4.v1']` |
| 7 | **Alternative diplotypes** | `*65;` (`AlternativePhase` col) | `dip_cand`: `*4,*10,*65,*2,*160,*34,*39` | only 1 solution (`SolutionID=1`) | — |
| 8 | **Phenotype** | `Intermediate Metabolizer` (`Phenotype` col) | `intermediate_metabolizer` (`phenotype` col) | `intermediate` (comment: `cpic=intermediate`) | `Intermediate metaboliser (IM)` (`Metaboliser status:`) |
| 9 | **Activity score** | — (not in `data.tsv`) | `1.0` (`dip_score` col) | `1.0` (comment: `cpic_score=1.0`) | `1.0` (`Activity score:`) |
| 10 | **SV / CNV type** | `Normal` (`CNV` col) | `no_sv,no_sv` (`dip_sv` / `hap1_sv` / `hap2_sv`) | no SV (2 copies total across `Copy` col: 0 and 1) | CN=2, no SV (`Initially computed CN = 2`) |
| 11 | **Copy number** | 2 (inferred from `CNV=Normal`) | 2 (inferred from `dip_sv=no_sv`) | 2 (allele 0: `Copy=0`; allele 1: `Copy=1`) | `Initially computed CN = 2` |
| 12 | **Supporting variants** | `VariantData`: `*4:22-42128945-C-T:0.462; *2:22-42126611-C-G,22-42127941-G-A:1.0,0.553; *10:22-42130692-G-A,22-42126611-C-G:0.757,1.0; *74:22-42129819-G-T:0.455` | hap1 core: `42126611 C>G (37/37 reads)`, `42127941 G>A (21/38)` · hap2 core: `42128945 C>T (18/39)` · hap2 tag: `42130692 G>A (28/37)` | 30 variant rows: positions 42126309→42132216, with read depths 17–57× | `42126611~C>G~1/1; 42127941~G>A~0/1; 42128945~C>T~0/1; 42129809~T>C~0/1; 42129819~G>T~0/1; 42130692~G>A~0/1` |
| 13 | **Functional effect** | — (not labelled in `VariantData`) | `missense_variant:low_impact:S486T`, `missense_variant:low_impact:R296C`, `splice_acceptor_variant:high_impact:splicing_defect`, `missense_variant:high_impact:P34S` | `Effect` col: `S486T`, `R296C`, `splice defect`, `P34S`, `H94R`, `L91M`; `Type` col: `missense`/`splice` | — |
| 14 | **dbSNP rsID** | — | — | `dbSNP` col: `rs1135840` (S486T), `rs16947` (R296C), `rs3892097` (splice/`*4`), `rs1065852` (P34S), plus 26 more | — |
| 15 | **Allele score / confidence** | — | `hap1_score=1.0`, `hap2_score=0.0`, `dip_score=1.0` · `ssr=1867.5` (short-seq-repeat CN marker) | `SolutionID=1` (single solution; rank 1 = best) | — |
| 16 | **Mean allele fraction** | per-variant AF in `VariantData`: `*4:...0.462`, `*2:...1.0, 0.553`, `*10:...0.757, 1.0` | `hap1_af_mean_gene=0.46`, `hap2_af_mean_gene=0.54`; `hap1_af_mean_main=0.55`, `hap2_af_mean_main=0.46` | `Coverage` col: raw read depth per position (17–57×); no AF | GT field in core variants (`0/1` or `1/1`) |
| 17 | **Phasing method** | Beagle statistical phasing (1KGP reference panel) | Beagle; `ssr` marker for CN normalisation; `BEAGLE imputed` flag per variant | ILP: CN + allele decomposition solved jointly; phasing implicit in solution | Graph-based calling via graphtyper on vg variation graph |

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

## Concordance Results — NA24385 (HG002), 2026-03-07

**Sample:** SPL250903104304-2 | BAM: T7_NA24385.bwa.sortdup.bqsr.bam (167 GB WGS, HG002/NA24385)
**Run:** 26 genes, 4 parallel jobs | Wall time: 3m17s

### Overall concordance summary

| Gene | Concordance | Notes |
|------|-------------|-------|
| CYP1A1 | ✅ 4/4 | |
| CYP1A2 | ⚠️ 2+2 split | PyPGx/Stargazer: `*1F/*1F`; Aldy/StellarPGx: `*30/*30` — nomenclature difference |
| CYP2A6 | ⚠️ 3/4 | StellarPGx: `*46/*46` (outlier); other 3: `*1/*1` |
| CYP2B6 | ✅ 4/4 | |
| CYP2C8 | ✅ 4/4 | |
| CYP2C9 | ✅ 4/4 | |
| CYP2C19 | ✅ 4/4 | |
| CYP2D6 | ✅ 4/4 | |
| CYP2E1 | ⚠️ 2/3 | PyPGx/Stargazer: `*1/*7`; Aldy: `*1/*7B` (suballele); StellarPGx: no-call (IndexError bug) |
| CYP3A4 | ✅ 4/4 | |
| CYP3A5 | ✅ 4/4 | |
| CYP4F2 | ⚠️ 2/4 | PyPGx/StellarPGx: `*4/*5`; Stargazer: `*1/*3` (outlier); Aldy: `*4+rs4020346/*5` |
| DPYD | ✅ 3/3 | PyPGx `Reference/Reference` = Stargazer/Aldy `*1/*1` (nomenclature) |
| G6PD | ✅ 2/2 | Aldy: no diplotype output; PyPGx `MALE/B (reference)` = Stargazer `*1/*1` |
| GSTM1 | ⚠️ 2/3 | PyPGx/StellarPGx: `*0/*A` (deletion detected); Stargazer VCF-only: `*1/*1` |
| GSTT1 | ✅ 2/2 | Stargazer/Aldy not supported |
| IFNL3 | ✅ 2/2 | Aldy: no diplotype; `Reference/Reference` = `*1/*1` |
| NAT1 | ❌ 4-way discordant | PyPGx: `*4/*4`; Stargazer: `*1/*10`; Aldy/StellarPGx: `*4/*10` |
| NAT2 | ✅ 4/4 | |
| NUDT15 | ✅ 4/4 | |
| POR | ✅ 2/2 | StellarPGx: no-call (cypor/por gene name bug); Aldy: not supported |
| RYR1 | ✅ 2/2 | Aldy: no diplotype; `Reference/Reference` = `*1/*1` |
| SLCO1B1 | ✅ 4/4 | |
| TPMT | ⚠️ 3/4 | StellarPGx: `*1/*1S` (suballele); others: `*1/*1` |
| UGT1A1 | ✅ 4/4 | |
| VKORC1 | ✅ 2/2 | PyPGx: `rs9923231/rs9923231`; Stargazer: `*S1/*S1` (nomenclature); Aldy: no diplotype |

### Per-gene results (features × tools)

Each table shows the 4 output features extracted from each tool's native output file.
`—` = field not reported by this tool. `n/s` = gene not supported by this tool. `no-call` = tool ran but produced no diplotype.

---

#### CYP1A1

| Feature | PyPGx | Stargazer | Aldy | StellarPGx |
|---------|-------|-----------|------|------------|
| Diplotype | `*1/*1` | `*1/*1` | `*1/*1` | `*1/*1` |
| Activity Score | — | 2.0 | 1 | — |
| Phenotype | Indeterminate | normal_metabolizer | — | — |
| SV Mode | no SVs expected | no SVs expected | no SVs expected | no SVs expected |

✅ **4/4 concordant**

---

#### CYP1A2

| Feature | PyPGx | Stargazer | Aldy | StellarPGx |
|---------|-------|-----------|------|------------|
| Diplotype | `*1F/*1F` | `*1F/*1F` | `*30/*30` | `*30/*30` |
| Activity Score | — | 3.0 | 1 | — |
| Phenotype | Indeterminate | ultrarapid_metabolizer | — | — |
| SV Mode | no SVs expected | no SVs expected | no SVs expected | no SVs expected |

⚠️ **2+2 split.** `*1F` and `*30` share the key CYP1A2 induction variant (rs762551 A>C) but differ in their additional tag SNV sets — a star allele database definition difference between PharmVar (used by PyPGx/Stargazer) and the CPIC/Aldy/StellarPGx databases.

---

#### CYP2A6

| Feature | PyPGx | Stargazer | Aldy | StellarPGx |
|---------|-------|-----------|------|------------|
| Diplotype | `*1/*1` | `*1/*1` | `*1/*1` | `*46/*46` |
| Activity Score | — | 2.0 | 1 | — |
| Phenotype | Indeterminate | normal_metabolizer | — | — |
| SV Mode | SV (depth+VDR) | SV (GDF/BAM) | SV (auto) | SV (auto) |

⚠️ **3/4 — StellarPGx outlier.** StellarPGx calls `*46/*46`, which is a CYP2A6 allele carrying a CYP2A7-derived sequence conversion. The graph-based aligner may be assigning paralog-origin reads to this allele. The 3-tool consensus `*1/*1` is the expected HG002 call.

---

#### CYP2B6

| Feature | PyPGx | Stargazer | Aldy | StellarPGx |
|---------|-------|-----------|------|------------|
| Diplotype | `*2/*5` | `*2/*5` | `*2/*5` | `*2/*5` |
| Activity Score | — | 2.0 | 1 | — |
| Phenotype | Normal Metabolizer | normal_metabolizer | — | — |
| SV Mode | SV (depth+VDR) | SV (GDF/BAM) | SV (auto) | SV (auto) |

✅ **4/4 concordant**

---

#### CYP2C8

| Feature | PyPGx | Stargazer | Aldy | StellarPGx |
|---------|-------|-----------|------|------------|
| Diplotype | `*1/*1` | `*1/*1` | `*1/*1` | `*1/*1` |
| Activity Score | — | 2.0 | 1 | — |
| Phenotype | Indeterminate | normal_metabolizer | — | — |
| SV Mode | no SVs expected | no SVs expected | no SVs expected | no SVs expected |

✅ **4/4 concordant**

---

#### CYP2C9

| Feature | PyPGx | Stargazer | Aldy | StellarPGx |
|---------|-------|-----------|------|------------|
| Diplotype | `*1/*1` | `*1/*1` | `*1/*1` | `*1/*1` |
| Activity Score | — | 2.0 | 1 | — |
| Phenotype | Normal Metabolizer | normal_metabolizer | — | — |
| SV Mode | no SVs expected | no SVs expected | no SVs expected | no SVs expected |

✅ **4/4 concordant**

---

#### CYP2C19

| Feature | PyPGx | Stargazer | Aldy | StellarPGx |
|---------|-------|-----------|------|------------|
| Diplotype | `*1/*1` | `*1/*1` | `*1/*1` | `*1/*1` |
| Activity Score | — | 2.0 | 1 | — |
| Phenotype | Normal Metabolizer | normal_metabolizer | — | — |
| SV Mode | no SVs expected | no SVs expected | no SVs expected | no SVs expected |

✅ **4/4 concordant**

---

#### CYP2D6

| Feature | PyPGx | Stargazer | Aldy | StellarPGx |
|---------|-------|-----------|------|------------|
| Diplotype | `*2/*4` | `*2/*4` | `*2/*4` | `*2/*4` |
| Activity Score | — | 1.0 | 1 | 1.0 |
| Phenotype | Intermediate Metabolizer | intermediate_metabolizer | — | Intermediate metaboliser (IM) |
| SV Mode | SV (depth+VDR) | SV (GDF/BAM) | SV (auto) | SV (auto) |

✅ **4/4 concordant**

---

#### CYP2E1

| Feature | PyPGx | Stargazer | Aldy | StellarPGx |
|---------|-------|-----------|------|------------|
| Diplotype | `*1/*7` | `*1/*7` | `*1/*7B` | no-call |
| Activity Score | — | — | 1 | — |
| Phenotype | Indeterminate | — | — | — |
| SV Mode | SV (depth+VDR) | SV (VCF-only) | SV (auto) | SV (auto) |

⚠️ **2/3 active tools.** PyPGx and Stargazer agree on `*1/*7`. Aldy calls the suballele `*1/*7B` (which carries the same key variants as `*7` plus an additional minor SNV). StellarPGx crashes with an `IndexError: list index out of range` in `sv_modules.py` (general bug for any sample — the `all_reg` list is empty for CYP2E1 in this tool version).

---

#### CYP3A4

| Feature | PyPGx | Stargazer | Aldy | StellarPGx |
|---------|-------|-----------|------|------------|
| Diplotype | `*1/*1` | `*1/*1` | `*1/*1` | `*1/*1` |
| Activity Score | — | 2.0 | 1 | — |
| Phenotype | Indeterminate | normal_metabolizer | — | — |
| SV Mode | no SVs expected | no SVs expected | no SVs expected | no SVs expected |

✅ **4/4 concordant**

---

#### CYP3A5

| Feature | PyPGx | Stargazer | Aldy | StellarPGx |
|---------|-------|-----------|------|------------|
| Diplotype | `*3/*3` | `*3/*3` | `*3/*3` | `*3/*3` |
| Activity Score | — | 0.0 | 1 | — |
| Phenotype | Poor Metabolizer | poor_metabolizer | — | — |
| SV Mode | no SVs expected | no SVs expected | no SVs expected | no SVs expected |

✅ **4/4 concordant**

---

#### CYP4F2

| Feature | PyPGx | Stargazer | Aldy | StellarPGx |
|---------|-------|-----------|------|------------|
| Diplotype | `*4/*5` | `*1/*3` | `*4+rs4020346/*5` | `*4/*5` |
| Activity Score | — | 1.5 | 1 | — |
| Phenotype | Indeterminate | normal_metabolizer | — | — |
| SV Mode | SV (depth+VDR) | SV (VCF-only) | SV (auto) | SV (auto) |

⚠️ **2/4 — Stargazer outlier.** PyPGx and StellarPGx agree on `*4/*5`. Stargazer calls `*1/*3`, which suggests it is not detecting the `*4` or `*5` defining variants (likely a VCF-only limitation for this SV gene — Stargazer does not run GDF/CN mode for CYP4F2). Aldy calls `*4+rs4020346/*5`, appending an additional tagged variant to `*4` — the core diplotype is equivalent.

---

#### DPYD

| Feature | PyPGx | Stargazer | Aldy | StellarPGx |
|---------|-------|-----------|------|------------|
| Diplotype | `Reference/Reference` | `*1/*1` | `*1/*1` | n/s |
| Activity Score | — | 2.0 | 1 | — |
| Phenotype | Normal Metabolizer | normal_metabolizer | — | — |
| SV Mode | no SVs expected | no SVs expected | no SVs expected | — |

✅ **3/3 concordant** (PyPGx uses `Reference/Reference` nomenclature where others use `*1/*1`; functionally identical)

---

#### G6PD

| Feature | PyPGx | Stargazer | Aldy | StellarPGx |
|---------|-------|-----------|------|------------|
| Diplotype | `MALE/B (reference)` | `*1/*1` | no-call | n/s |
| Activity Score | — | 2.0 | 1 | — |
| Phenotype | G6PD Normal | normal_function | — | — |
| SV Mode | SV (depth+VDR) | SV (VCF-only) | SV (auto) | — |

✅ **2/2 active diplotype callers** (PyPGx and Stargazer both call reference/normal). Aldy runs but outputs no diplotype string for G6PD (hemizygous X-linked model may not produce a standard diplotype field in this version).

---

#### GSTM1

| Feature | PyPGx | Stargazer | Aldy | StellarPGx |
|---------|-------|-----------|------|------------|
| Diplotype | `*0/*A` | `*1/*1` | no-call | `*0/*A` |
| Activity Score | — | 2.0 | 1 | — |
| Phenotype | Indeterminate | normal_function | — | — |
| SV Mode | SV (depth+VDR) | SV (VCF-only) | SV (auto) | SV (auto) |

⚠️ **2/3 active diplotype callers.** PyPGx and StellarPGx detect the GSTM1 whole-gene deletion (`*0` = null allele) heterozygous with a functional copy (`*A`). Stargazer operates in VCF-only mode for GSTM1 (it does not produce a GDF for this gene) and therefore misses the deletion, defaulting to `*1/*1`. Aldy outputs no diplotype string. The PyPGx/StellarPGx `*0/*A` call is expected for HG002.

---

#### GSTT1

| Feature | PyPGx | Stargazer | Aldy | StellarPGx |
|---------|-------|-----------|------|------------|
| Diplotype | `*0/*A` | n/s | n/s | `*0/*A` |
| Activity Score | — | — | — | — |
| Phenotype | Indeterminate | — | — | — |
| SV Mode | SV (depth+VDR) | — | — | SV (auto) |

✅ **2/2 concordant** (only supported tools)

---

#### IFNL3

| Feature | PyPGx | Stargazer | Aldy | StellarPGx |
|---------|-------|-----------|------|------------|
| Diplotype | `Reference/Reference` | `*1/*1` | no-call | n/s |
| Activity Score | — | 2.0 | 1 | — |
| Phenotype | Favorable Response | normal_function | — | — |
| SV Mode | no SVs expected | no SVs expected | no SVs expected | — |

✅ **2/2 concordant** (Aldy: no diplotype output; Reference = `*1/*1`)

---

#### NAT1

| Feature | PyPGx | Stargazer | Aldy | StellarPGx |
|---------|-------|-----------|------|------------|
| Diplotype | `*4/*4` | `*1/*10` | `*4/*10` | `*10/*4` |
| Activity Score | — | — | 1 | — |
| Phenotype | Indeterminate | — | — | — |
| SV Mode | no SVs expected | no SVs expected | no SVs expected | no SVs expected |

❌ **Discordant across all 4 tools.** PyPGx calls `*4/*4` (homozygous reference-equivalent). Aldy and StellarPGx agree on `*4/*10` (heterozygous). Stargazer calls `*1/*10` — `*1` and `*4` are functionally similar for NAT1 but differ in their tag SNV definitions. The Aldy/StellarPGx consensus `*4/*10` is likely correct; PyPGx may be phasing both haplotypes to the reference `*4` due to insufficient discriminating variants being included in the VCF.

---

#### NAT2

| Feature | PyPGx | Stargazer | Aldy | StellarPGx |
|---------|-------|-----------|------|------------|
| Diplotype | `*4/*6` | `*4/*6` | `*4/*6` | `*4/*6` |
| Activity Score | — | 2.0 | 1 | — |
| Phenotype | Indeterminate | normal_function | — | — |
| SV Mode | no SVs expected | no SVs expected | no SVs expected | no SVs expected |

✅ **4/4 concordant**

---

#### NUDT15

| Feature | PyPGx | Stargazer | Aldy | StellarPGx |
|---------|-------|-----------|------|------------|
| Diplotype | `*1/*1` | `*1/*1` | `*1/*1` | `*1/*1` |
| Activity Score | — | 2.0 | 1 | — |
| Phenotype | Normal Metabolizer | normal_metabolizer | — | — |
| SV Mode | no SVs expected | no SVs expected | no SVs expected | no SVs expected |

✅ **4/4 concordant**

---

#### POR

| Feature | PyPGx | Stargazer | Aldy | StellarPGx |
|---------|-------|-----------|------|------------|
| Diplotype | `*1/*28` | `*1/*28` | n/s | no-call |
| Activity Score | — | 1.5 | — | — |
| Phenotype | Indeterminate | decreased_function | — | — |
| SV Mode | no SVs expected | no SVs expected | — | no SVs expected |

✅ **2/2 active tools concordant.** StellarPGx produces no-call due to a gene name mismatch bug: StellarPGx's `main.nf` only defines gene coordinates for `cypor` (not `por`). The suite remaps `por → cypor` but StellarPGx still fails to produce output for this sample. Aldy does not support POR.

---

#### RYR1

| Feature | PyPGx | Stargazer | Aldy | StellarPGx |
|---------|-------|-----------|------|------------|
| Diplotype | `Reference/Reference` | `*1/*1` | no-call | n/s |
| Activity Score | — | 2.0 | 1 | — |
| Phenotype | Uncertain Susceptibility | normal_function | — | — |
| SV Mode | no SVs expected | no SVs expected | no SVs expected | — |

✅ **2/2 concordant** (Aldy: no diplotype; Reference = `*1/*1`)

---

#### SLCO1B1

| Feature | PyPGx | Stargazer | Aldy | StellarPGx |
|---------|-------|-----------|------|------------|
| Diplotype | `*1/*1` | `*1/*1` | `*1/*1` | `*1/*1` |
| Activity Score | — | 2.0 | 1 | — |
| Phenotype | Normal Function | normal_function | — | — |
| SV Mode | no SVs expected | no SVs expected | no SVs expected | no SVs expected |

✅ **4/4 concordant**

---

#### TPMT

| Feature | PyPGx | Stargazer | Aldy | StellarPGx |
|---------|-------|-----------|------|------------|
| Diplotype | `*1/*1` | `*1/*1` | `*1/*1` | `*1/*1S` |
| Activity Score | — | 2.0 | 1 | — |
| Phenotype | Normal Metabolizer | normal_metabolizer | — | — |
| SV Mode | no SVs expected | no SVs expected | no SVs expected | no SVs expected |

⚠️ **3/4.** StellarPGx calls `*1/*1S`. The `*1S` designation is a StellarPGx-internal suballele label that maps to the standard `*1` (wild-type); it indicates StellarPGx resolved a phasing ambiguity by assigning a secondary tag. Functionally equivalent to `*1/*1`.

---

#### UGT1A1

| Feature | PyPGx | Stargazer | Aldy | StellarPGx |
|---------|-------|-----------|------|------------|
| Diplotype | `*1/*1` | `*1/*1` | `*1/*1` | `*1/*1` |
| Activity Score | — | 2.0 | 1 | — |
| Phenotype | Normal Metabolizer | normal_metabolizer | — | — |
| SV Mode | no SVs expected | no SVs expected | no SVs expected | no SVs expected |

✅ **4/4 concordant**

---

#### VKORC1

| Feature | PyPGx | Stargazer | Aldy | StellarPGx |
|---------|-------|-----------|------|------------|
| Diplotype | `rs9923231/rs9923231` | `*S1/*S1` | no-call | n/s |
| Activity Score | — | 3.0 | 1 | — |
| Phenotype | Indeterminate | increased_function | — | — |
| SV Mode | no SVs expected | no SVs expected | no SVs expected | — |

✅ **2/2 concordant** (nomenclature differs: PyPGx reports the causal rsID; Stargazer uses the `*S1` star allele label for the same haplotype). Aldy runs but outputs no diplotype string for VKORC1.

---

## Genes with Partial Tool Support

Some genes in the 29-gene suite are not supported by all five callers. Key cases:

| Gene | Notes |
|------|-------|
| **ABCG2** | Supported by Aldy + StellarPGx only. PyPGx installed database does not recognise ABCG2 in this version; Stargazer similarly produces no call. CPIC Level A/B (rosuvastatin). |
| **HLA-A / HLA-B** | Supported by OptiType only. The four star allele callers (PyPGx, Stargazer, Aldy, StellarPGx) do not perform HLA typing. |
| **GSTT1** | Supported by PyPGx + StellarPGx only. Located on `chr22_KI270879v1_alt`; Stargazer and Aldy cannot call it from a standard GRCh38 reference. |
| **G6PD** | Not supported by StellarPGx. X-linked hemizygous model not in `main.nf`. |
| **DPYD / IFNL3 / RYR1 / VKORC1** | Not in StellarPGx gene list. |
| **ABCG2 / G6PD / VKORC1 / DPYD / IFNL3 / RYR1 / POR** | Aldy support varies; POR not in Aldy database. |

---

## References

| Tool | Publication |
|------|-------------|
| **PyPGx** | Lee S-B et al. (2022) *Pharmacogenomics J.* doi:10.1038/s41397-021-00257-5 |
| **Stargazer** | Lee S-B et al. (2019) *Genet. Med.* doi:10.1038/s41436-018-0267-1 |
| **Aldy** | Numanagić I et al. (2018) *Nat. Methods* doi:10.1038/s41592-018-0078-y; v4 update 2023 |
| **StellarPGx** | Twesigomwe D et al. (2021) *Clin. Pharmacol. Ther.* doi:10.1002/cpt.2173 |
| **OptiType** | Szolek A et al. (2014) *Bioinformatics* doi:10.1093/bioinformatics/btu548 |

---

## Licensing Notice

| Tool | Licence | Restriction |
|------|---------|-------------|
| PyPGx | MIT | Open use |
| Stargazer | Non-commercial academic (University of Washington) | **No commercial use; do not push to public registries** |
| Aldy | Non-commercial academic (IURTC) | **No commercial use; do not push to public registries** |
| StellarPGx | MIT | Open use |
| OptiType | MIT | Open use (SIF image via Bioconda/BioContainers) |

The `pgx-suite` Docker image must not be uploaded to Docker Hub, GHCR, or any other
public container registry because it bundles Stargazer and Aldy source code.
