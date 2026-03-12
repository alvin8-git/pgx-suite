# PGx Suite

A Docker container that bundles six pharmacogenomics (PGx) callers — **PyPGx**, **Stargazer**, **Aldy**, **StellarPGx**, **OptiType**, and **mutserve** — pre-configured for **GRCh38 (hg38)**. A single command runs all applicable tools for a gene and produces a side-by-side concordance table plus a self-contained HTML clinical report. Covers **all 19 CPIC Level A genes** across 31 genes total.

> **License notice:** Stargazer and Aldy are restricted to non-commercial academic use.
> This image **must not be published to any public registry** (Docker Hub, GHCR, etc.).

---

## Table of Contents

- [Bundled Tools](#bundled-tools)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
  - [Host launcher — recommended](#host-launcher--recommended)
  - [Run all genes (batch, manual docker)](#run-all-genes-batch-manual-docker)
  - [Run a single gene](#run-a-single-gene)
  - [Volume mounts reference](#volume-mounts-reference)
  - [Expected output](#expected-output)
- [Gene Coverage](#gene-coverage)
- [Output Fields](#output-fields)
- [Architecture](#architecture)
- [Troubleshooting](#troubleshooting)
- [Repository Layout](#repository-layout)
- [Documentation](#documentation)
- [License](#license)

---

## Bundled Tools

| Tool | Version | Algorithm | Genes | License |
|------|---------|-----------|-------|---------|
| [PyPGx](https://github.com/sbslee/pypgx) | 0.26.0 | Bayesian / Beagle phasing | 88 | MIT |
| [Stargazer](https://stargazer.gs.washington.edu/stargazerweb/) | 2.0.3 | Bayesian / Beagle phasing | 58 | Non-commercial academic (UW) |
| [Aldy](https://github.com/0xTCG/aldy) | 4.8.3 | Integer Linear Programming | 39 | Non-commercial academic (IURTC) |
| [StellarPGx](https://github.com/SBIMB/StellarPGx) | 1.2.7 | Genome graph (graphtyper) | 21 | MIT |
| [OptiType](https://github.com/FRED-2/OptiType) | 1.3.5 | ILP on HLA-ref-filtered reads | HLA-A/B/C | MIT |
| [mutserve](https://github.com/seppinho/mutserve) | 2.0.0 | Allele-fraction pileup (chrM) | MT-RNR1 | MIT |

See [`ToolsDocumentation.md`](ToolsDocumentation.md) for a detailed comparison of each tool's capabilities, gene lists, SV handling, and limitations.

---

## Requirements

### Software

- **Docker Engine** (Linux) or **Docker Desktop** (macOS / Windows) with `--privileged` support
- ~3 GB free disk space for the built image

### Data (user-supplied)

| File | Notes |
|------|-------|
| GRCh38 reference FASTA + `.fai` index | 3–50 GB; not included |
| Input BAM or CRAM + `.bai` / `.crai` index | Sorted, GRCh38-aligned |

### Repo contents (already present after cloning)

| Path | Contents |
|------|----------|
| `pypgx/` | PyPGx source (pip-installable) |
| `pypgx/pypgx-bundle/` | 1KGP VCF panel + CNV data (~500 MB, volume-mounted at runtime) |
| `stargazer-grc38-2.0.3/` | Stargazer GRCh38 source |
| `StellarPGx/` | Nextflow pipeline, gene databases, resources |
| `StellarPGx/containers/stellarpgx-dev.sif` | StellarPGx Apptainer SIF (31 MB) |
| `StellarPGx/containers/optitype.sif` | OptiType Apptainer SIF (~500 MB; pull separately — see below) |
| `nextflow` | Nextflow binary (pre-downloaded) |

**Pulling the OptiType SIF** (one-time, requires Docker + outbound internet):

```bash
docker run --privileged --rm \
  -v "$(pwd)/StellarPGx/containers:/pgx/containers" \
  pgx-suite:latest \
  apptainer pull --name /pgx/containers/optitype.sif \
    docker://quay.io/biocontainers/optitype:1.3.5--hdfd78af_1
```

---

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/<org>/pgx-suite.git
cd pgx-suite
```

### 2. Build the Docker image

```bash
docker build -t pgx-suite:latest .
```

First build takes 15–25 minutes. Subsequent builds are fast due to Docker layer caching.

### 3. Verify with smoke tests

```bash
# Fast check — no data volumes required:
docker run --rm pgx-suite:latest bash /opt/pgx/test.sh

# Full check — includes StellarPGx (volumes required):
./docker/docker-run.sh bash /opt/pgx/test.sh
```

All 12 checks should report `PASS`:

```
============================================================
 PGx Suite Docker Container — Smoke Tests
 Reference build: GRCh38
============================================================

[1/4] PyPGx
  pypgx --version                                    PASS
  pypgx import                                       PASS
  pypgx run-ngs-pipeline help                        PASS

[2/4] Stargazer
  stargazer --version                                PASS
  stargazer wrapper exists                           PASS
  stargazer CYP2D6 VCF-only (grc38)                 PASS

[3/4] Aldy
  aldy --version                                     PASS
  aldy import                                        PASS
  aldy built-in test suite (78 tests)                PASS

[4/4] StellarPGx
  nextflow --version                                 PASS
  apptainer --version                                PASS
  StellarPGx test pipeline (CYP2D6/hg38)            PASS

============================================================
 Results: 12 PASSED  |  0 FAILED
============================================================
```

---

## Usage

### Host launcher — recommended

**`run_pgx_suite.sh`** is the single-command entry point for new samples. It handles
all Docker volume mounts automatically, accepts BAM or CRAM, and delegates to
`pgx-all-genes.sh` inside the container.

```bash
# Minimal — output defaults to results/<SAMPLE>/ under the repo root
./run_pgx_suite.sh /path/to/sample.bam

# With explicit output directory
./run_pgx_suite.sh /path/to/sample.bam --output /data/pgx_out

# CRAM input, 8 parallel gene workers
./run_pgx_suite.sh /path/to/sample.cram --output /data/pgx_out --jobs 8

# Override Docker image (e.g. a tagged release)
./run_pgx_suite.sh /path/to/sample.bam --image pgx-suite:v1.0
```

**Options:**

| Option | Default | Description |
|--------|---------|-------------|
| `--output DIR` | `<repo>/results/<SAMPLE>` | Output directory on host |
| `--jobs N` | `4` | Genes to run in parallel inside the container |
| `--image NAME` | `pgx-suite:latest` | Docker image to use |

**Pre-flight checks performed before launching Docker:**
- BAM/CRAM exists and has a matching `.bai`/`.crai` index
- `StellarPGx/`, `StellarPGx/containers/`, `pypgx/pypgx-bundle/`, `GRCh38/` are present
- `pgx-suite:latest` image is built

Sample name is derived from the filename before the first `.`
(e.g. `NA12878.bwa.bam` → `NA12878`).

The container is run with `--privileged` (required for Apptainer used by StellarPGx
and OptiType HLA typing).

---

### Run all genes (batch, manual docker)

The manual `docker run` equivalent of `run_pgx_suite.sh` — useful when you need to
customise volume mounts or run from inside another script:

```bash
docker run --privileged --rm \
  -v "$(pwd)/StellarPGx:/pgx/stellarpgx" \
  -v "$(pwd)/StellarPGx/containers:/pgx/containers" \
  -v "$(pwd)/pypgx/pypgx-bundle:/pgx/bundle" \
  -v "/path/to/ref:/pgx/ref:ro" \
  -v "/path/to/data:/pgx/data:ro" \
  -v "/path/to/results:/pgx/results" \
  pgx-suite:latest \
  pgx-all-genes.sh /pgx/data/sample.bam --output /pgx/results
```

**Supported options:**

```
pgx-all-genes.sh <BAM|CRAM> [--ref PATH] [--output PATH] [--jobs N]
```

`pgx-all-genes.sh` runs every gene in the 31-gene support matrix in parallel batches,
collects all results into `log/all_genes_summary.tsv`, and generates a standalone HTML
report at `<output>/<SAMPLE>_pgx_report.html`.

| Option | Default | Description |
|--------|---------|-------------|
| `--ref PATH` | `/pgx/ref/hg38.fa` | GRCh38 reference FASTA |
| `--output PATH` | `/pgx/results` | Root output directory |
| `--jobs N` | `4` | Genes to run concurrently |

Output layout:

```
<output>/
├── <SAMPLE>_pgx_report.html          # standalone single-file HTML report (all 31 gene panels embedded)
├── Genes/
│   └── <GENE>/                       # per-tool outputs for each gene
│       ├── <GENE>_<SAMPLE>_comparison.tsv
│       ├── <GENE>_<SAMPLE>_detail.json
│       ├── pypgx/results.zip
│       ├── stargazer/genotype-calls.tsv
│       ├── aldy/<GENE>.aldy
│       ├── stellarpgx/<gene>/alleles/*.alleles
│       ├── optitype/<SAMPLE>_result.tsv   # HLA-A/HLA-B only
│       └── mt-rnr1/<SAMPLE>_mtrna1_result.json  # MT-RNR1 only
└── log/
    ├── all_genes_summary.tsv          # concordance table across all genes and tools
    ├── bam_stats.json                 # whole-BAM QC metrics (incl. per-gene depth)
    └── bamstats.log                   # pgx-bamstats.sh stdout/stderr
```

### Run a single gene

```bash
docker run --privileged --rm \
  -v "$(pwd)/StellarPGx:/pgx/stellarpgx" \
  -v "$(pwd)/StellarPGx/containers:/pgx/containers" \
  -v "$(pwd)/pypgx/pypgx-bundle:/pgx/bundle" \
  -v "/path/to/ref:/pgx/ref" \
  -v "/path/to/data:/pgx/data" \
  -v "/path/to/results:/pgx/results" \
  pgx-suite:latest \
  pgx-run.sh CYP2D6 /pgx/data/sample.bam
```

Or use the convenience wrapper (sets all standard volume mounts):

```bash
./docker/docker-run.sh pgx-run.sh CYP2D6 /pgx/data/sample.bam
```

**Supported options:**

```
pgx-run.sh <GENE> <BAM> [--ref /pgx/ref/hg38.fa] [--output /pgx/results]
```

`pgx-run.sh` automatically:
- Validates the BAM/CRAM index, reference FASTA, and gene support
- Derives the sample name from the filename before the first `.` (e.g. `NA12878.bwa.bam` → `NA12878`)
- Generates a gene-region VCF via `bcftools mpileup | bcftools call | tabix`
- Runs SV-aware preprocessing (depth-of-coverage + VDR control stats) for PyPGx SV genes
- Runs GDF depth-profile creation for Stargazer's three paralog genes (CYP2A6, CYP2B6, CYP2D6)
- Calls all applicable tools sequentially; individual failures are logged without aborting
- Invokes `pgx-compare.py` to produce the comparison table and TSV

### Volume mounts reference

| Host path | Container path | Purpose |
|-----------|---------------|---------|
| `./StellarPGx` | `/pgx/stellarpgx` | StellarPGx pipeline scripts, databases, resources |
| `./StellarPGx/containers` | `/pgx/containers` | `stellarpgx-dev.sif` Apptainer image |
| `./pypgx/pypgx-bundle` | `/pgx/bundle` | 1KGP VCFs + CNV data for PyPGx phasing |
| `/path/to/ref` | `/pgx/ref` | GRCh38 reference FASTA + `.fai` |
| `/path/to/data` | `/pgx/data` | Input BAM/CRAM + index |
| `/path/to/results` | `/pgx/results` | Analysis output directory |

### Expected output

```
========================================================================
 PGx Star Allele Results
 Gene:   CYP2D6
 Sample: HG002
 Build:  GRCh38
 SV:     PyPGx: depth-of-coverage + VDR control stats
         Stargazer: GDF from BAM (paralog CN normalisation)
========================================================================
Tool          Diplotype         Activity Score    Phenotype
────────────────────────────────────────────────────────────────────────
PyPGx         *1/*4             1.0               Normal Metabolizer
Stargazer     *1/*4             1.0               Normal Metabolizer
Aldy          *2/*4             1.0               Normal Metabolizer
StellarPGx    *1/*4             -                 -
────────────────────────────────────────────────────────────────────────
Concordance: 3/4 tools agree on *1/*4
========================================================================

Full results saved to: /pgx/results/CYP2D6_HG002_comparison.tsv
```

The TSV (`<GENE>_<SAMPLE>_comparison.tsv`) contains one row per tool with columns:
`Tool`, `Diplotype`, `ActivityScore`, `Phenotype`, `SVMode`.

A rich `<GENE>_<SAMPLE>_detail.json` is also written with all 17 fields per tool —
this feeds the HTML per-gene detail pages.

---

## Gene Coverage

**31 genes** supported across **six tools**. Covers **19/19 CPIC Level A genes**. For full per-tool gene lists and SV details see [`ToolsDocumentation.md`](ToolsDocumentation.md).

| Gene | PyPGx | Stargazer | Aldy | StellarPGx | OptiType | mutserve | CPIC Level | SV? |
|------|:-----:|:---------:|:----:|:----------:|:--------:|:--------:|:----------:|:---:|
| ABCG2 | — | — | ✓ | ✓ | — | — | B | |
| CACNA1S | ✓ | — | — | ✓ | — | — | **A** | |
| CYP1A1 | ✓ | ✓ | ✓ | ✓ | — | — | — | |
| CYP1A2 | ✓ | ✓ | ✓ | ✓ | — | — | B | |
| CYP2A6 | ✓ | ✓ | ✓ | ✓ | — | — | B | ✓ paralog |
| CYP2B6 | ✓ | ✓ | ✓ | ✓ | — | — | **A** | ✓ paralog |
| CYP2C8 | ✓ | ✓ | ✓ | ✓ | — | — | B | |
| CYP2C9 | ✓ | ✓ | ✓ | ✓ | — | — | **A** | |
| CYP2C19 | ✓ | ✓ | ✓ | ✓ | — | — | **A** | |
| CYP2D6 | ✓ | ✓ | ✓ | ✓ | — | — | **A** | ✓ paralog + tandem |
| CYP2E1 | ✓ | ✓ | ✓ | ✓ | — | — | — | ✓ CN |
| CYP3A4 | ✓ | ✓ | ✓ | ✓ | — | — | A† | |
| CYP3A5 | ✓ | ✓ | ✓ | ✓ | — | — | **A** | |
| CYP4F2 | ✓ | ✓ | ✓ | ✓ | — | — | B | ✓ CN |
| DPYD | ✓ | ✓ | ✓ | — | — | — | **A** | |
| G6PD | ✓ | ✓ | ✓ | — | — | — | **A** | ✓ CN |
| GSTM1 | ✓ | ✓ | ✓ | ✓ | — | — | — | ✓ deletion |
| GSTT1 | ✓‡ | — | — | ✓ | — | — | — | ✓ deletion |
| HLA-A | — | — | — | — | ✓ | — | **A** | |
| HLA-B | — | — | — | — | ✓ | — | **A** | |
| IFNL3 | ✓ | ✓ | ✓ | — | — | — | **A** | |
| MT-RNR1 | — | — | — | — | — | ✓ | **A** | |
| NAT1 | ✓ | ✓ | ✓ | ✓ | — | — | B | |
| NAT2 | ✓ | ✓ | ✓ | ✓ | — | — | **A** | |
| NUDT15 | ✓ | ✓ | ✓ | ✓ | — | — | **A** | |
| POR | ✓ | ✓ | — | ✓ | — | — | — | |
| RYR1 | ✓ | ✓ | ✓ | — | — | — | **A** | |
| SLCO1B1 | ✓ | ✓ | ✓ | ✓ | — | — | **A** | |
| TPMT | ✓ | ✓ | ✓ | ✓ | — | — | **A** | |
| UGT1A1 | ✓ | ✓ | ✓ | ✓ | — | — | **A** | |
| VKORC1 | ✓ | ✓ | ✓ | — | — | — | **A** | |

† CYP3A4 appears in the tacrolimus guideline alongside CYP3A5; no standalone CPIC Level A guideline for CYP3A4 genotyping alone.
‡ GSTT1 is on `chr22_KI270879v1_alt` (alternate contig) — bcftools mpileup is skipped; PyPGx depth preprocessing may also fail on standard GRCh38 references.

---

## Output Fields

Each tool reports different field names for equivalent concepts. The table below provides the mapping (N = 17 fields × 4 tools).

| # | Field | PyPGx (`data.tsv`) | Stargazer (`report.tsv` / `genotype-calls.tsv`) | Aldy (`.aldy`) | StellarPGx (`.alleles`) |
|---|---|---|---|---|---|
| 1 | **Sample ID** | `Sample` | `name` / `Sample` (report) | `Sample` | Filename stem |
| 2 | **Gene** | Implicit (one run per gene) | `Gene` (report) | `Gene` | Header line |
| 3 | **Diplotype** | `Genotype` (e.g. `*2/*4`) | `Diplotype` / `hap1_main`+`hap2_main` | `Major` (e.g. `*2+*4`) | `Result` (e.g. `*17/*29`) |
| 4 | **Haplotype 1** | `Haplotype1` (e.g. `*2;`) | `hap1_main` | First allele in `Major` | First allele in `Result` |
| 5 | **Haplotype 2** | `Haplotype2` (e.g. `*4;*10;*74;`) | `hap2_main` | Second allele in `Major` | Second allele in `Result` |
| 6 | **Sub-alleles / tag variants** | Semicolon list in `Haplotype1`/`Haplotype2` | `hap1_main_core`+`hap1_main_tag`; `hap2_main_core`+`hap2_main_tag` | `Minor` column | `Candidate alleles` |
| 7 | **Alternative diplotypes** | `AlternativePhase` (semicolon list) | `dip_cand` / `May also be` | Multiple `SolutionID` rows | — |
| 8 | **Phenotype** | `Phenotype` (e.g. `Intermediate Metabolizer`) | `Phenotype` (e.g. `intermediate_metabolizer`) | `Status` | `Metaboliser status` |
| 9 | **Activity score** | `ActivityScore` | `Score` (report) | — | `Activity score` |
| 10 | **SV / CNV type** | `CNV` (e.g. `Normal`, `Deletion`) | `dip_sv` / `hap1_sv` / `hap2_sv` | Implicit in `Major` + `Copy` | `Initially computed CN` |
| 11 | **Copy number** | Derived from `CNV` | Derived from `dip_sv` | `Copy` (integer per allele) | `Initially computed CN` |
| 12 | **Supporting variants** | `VariantData` (`allele:pos-ref-alt:AF;…`) | `hap1_main_core`/`hap2_main_core` | `Location` + `Coverage` | `Sample core variants` (`pos~ref>alt~GT`) |
| 13 | **Functional effect** | Embedded in `VariantData` | Embedded in `hap1_main_core` | `Effect` + `Type` columns | — |
| 14 | **dbSNP rsID** | Embedded in `VariantData` | — | `dbSNP` column | — |
| 15 | **Allele score / confidence** | — | `dip_score`, `hap1_score`, `hap2_score` | `SolutionID` rank | — |
| 16 | **Mean allele fraction** | Per-variant AF in `VariantData` | `hap1_af_mean_gene` / `hap2_af_mean_gene` / `hap1_af_mean_main` / `hap2_af_mean_main` | `Coverage` per variant | — |
| 17 | **Phasing method** | Beagle statistical (1KGP panel) | Beagle; `BEAGLE imputed` flag + `ssr` marker | ILP joint optimisation | Graph-based (graphtyper) |

---

## Architecture

```
docker run --privileged pgx-suite:latest
┌───────────────────────────────────────────────────────────────┐
│  Ubuntu 22.04                                                  │
│                                                                │
│  Python 3.11                                                   │
│  ├── pypgx 0.26.0      (pip from source)                      │
│  ├── aldy 4.8.3        (pip install)                          │
│  └── shared deps: pysam, pandas, numpy, scipy, ortools, …    │
│                                                                │
│  Stargazer 2.0.3       (/usr/local/bin/stargazer wrapper)     │
│                                                                │
│  Java 21 JRE           (Beagle phasing for PyPGx + Stargazer; │
│                         mutserve.jar for MT-RNR1)             │
│                                                                │
│  mutserve 2.0.0        (JAR; chrM allele-fraction calling)    │
│                                                                │
│  Nextflow              (copied from host binary)              │
│                                                                │
│  Apptainer             (Singularity fork, from PPA)           │
│  ├── runs stellarpgx-dev.sif  [volume-mounted]                │
│  │   ├── graphtyper (graph-based variant calling)             │
│  │   ├── bcftools, samtools, tabix                            │
│  │   └── stellarpgx.py (star allele caller)                   │
│  └── runs optitype.sif  [volume-mounted]                      │
│      ├── razers3 (HLA-reference read filtering)               │
│      └── OptiTypePipeline.py (ILP HLA Class I typing)        │
│                                                                │
│  bcftools + samtools + tabix + mosdepth  (system packages)    │
└───────────────────────────────────────────────────────────────┘
```

**Why `--privileged`?** StellarPGx runs its core tools inside `stellarpgx-dev.sif` via Apptainer. Apptainer inside Docker requires `SYS_ADMIN` capability to unpack SIF overlay filesystems; `--privileged` is the standard approach for local workstations.

---

## Troubleshooting

| Symptom | Cause | Fix |
|---------|-------|-----|
| `FATAL: could not open image` | SIF not mounted | Add `-v $(pwd)/StellarPGx/containers:/pgx/containers` |
| `FATAL: kernel too old` | Missing privilege | Run with `--privileged` |
| `ModuleNotFoundError: apt_pkg` during build | Build order issue | Apptainer PPA step must run before Python 3.11 switch in Dockerfile |
| `pypgx` not found in container | Python path | Verify `python3 --version` → 3.11 inside container |
| Nextflow hangs on first start | JAR download | Ensure outbound internet access on first `nextflow` run |
| `beagle.jar` error | Java not in PATH | `JAVA_HOME=/usr/lib/jvm/java-21-openjdk-amd64`; verify with `java -version` |
| Stargazer: `no data for target gene` with GDF | hg19 GDF used in grc38 mode | Regenerate GDF with `-a grc38`, or use VCF-only mode (omit `-c` and `-g`) |
| bcftools VCF empty for GSTT1 | Alt contig `chr22_KI270879v1_alt` | Expected — `pgx-run.sh` skips bcftools for GSTT1 automatically |

---

## Repository Layout

```
pgx-suite/
├── Dockerfile                          # Ubuntu 22.04 image with all 5 tools
├── .dockerignore
├── run_pgx_suite.sh                    # Host launcher — single command for new samples
├── nextflow                            # Nextflow binary (pre-downloaded)
├── docker/
│   ├── pgx-all-genes.sh                # Batch orchestrator: all 31 genes in parallel
│   ├── pgx-run.sh                      # Single-gene orchestration entry point
│   ├── pgx-hla.sh                      # HLA typing via OptiType Apptainer SIF
│   ├── pgx-mt.sh                       # MT-RNR1 calling via mutserve JAR
│   ├── pgx-compare.py                  # Result parser → comparison TSV + detail JSON
│   ├── pgx-report.py                   # HTML report generator (standalone single-file)
│   ├── pgx-bamstats.sh                 # Whole-BAM QC → bam_stats.json (29 genes + chrM)
│   ├── test.sh                         # Smoke tests (no BAM required)
│   ├── docker-run.sh                   # Convenience docker run wrapper
│   └── README.md                       # Container-specific notes
├── pypgx/                              # PyPGx source
│   └── pypgx-bundle/                   # 1KGP VCF panel (volume-mounted)
├── stargazer-grc38-2.0.3/              # Stargazer GRCh38 source
├── StellarPGx/                         # Nextflow pipeline (volume-mounted)
│   ├── containers/stellarpgx-dev.sif   # Apptainer image (31 MB)
│   ├── database/                       # Star allele databases
│   └── resources/                      # Reference sequences, panel data
├── TestData/                           # Test BAM files
├── docs/plans/                         # Design and implementation documents
├── CHANGES.md                          # Changelog
├── ToolsDocumentation.md               # Detailed tool reference
├── README.md                           # This file
└── TODO.md                             # Roadmap
```

---

## Documentation

| Document | Description |
|----------|-------------|
| [`ToolsDocumentation.md`](ToolsDocumentation.md) | Full reference for all 6 callers: algorithms, gene lists, SV handling, output formats, limitations |
| [`CHANGES.md`](CHANGES.md) | Reverse-chronological changelog |
| [`TODO.md`](TODO.md) | Roadmap and open tasks |
| [`docker/README.md`](docker/README.md) | Container-specific notes and examples |

---

## License

| Component | License |
|-----------|---------|
| pgx-suite (this repo) | MIT |
| PyPGx | MIT |
| StellarPGx | MIT |
| Stargazer | Non-commercial academic (University of Washington) |
| Aldy | Non-commercial academic (Indiana University Research and Technology Corporation) |
| mutserve | MIT |

Stargazer and Aldy source code is bundled in this image for academic research use only. **Do not publish this image to Docker Hub, GHCR, or any public container registry.**
