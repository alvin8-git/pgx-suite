# PGx Suite

A Docker container that bundles four open-source pharmacogenomics (PGx) star allele callers — **PyPGx**, **Stargazer**, **Aldy**, and **StellarPGx** — pre-configured for **GRCh38 (hg38)**. A single command runs all applicable tools for a gene and produces a side-by-side concordance table.

> **License notice:** Stargazer and Aldy are restricted to non-commercial academic use.
> This image **must not be published to any public registry** (Docker Hub, GHCR, etc.).

---

## Table of Contents

- [Bundled Tools](#bundled-tools)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
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
| `StellarPGx/containers/stellarpgx-dev.sif` | Apptainer SIF image (31 MB) |
| `nextflow` | Nextflow binary (pre-downloaded) |

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
- Validates the BAM index, reference FASTA, and gene support
- Extracts the sample name from the BAM `@RG SM:` tag
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

The TSV (`<GENE>_<SAMPLE>_comparison.tsv`) contains one row per tool with columns: `Tool`, `Diplotype`, `ActivityScore`, `Phenotype`, `SVMode`.

---

## Gene Coverage

Genes supported by all four tools are shown below. For full per-tool gene lists and SV details see [`ToolsDocumentation.md`](ToolsDocumentation.md).

| Gene | PyPGx | Stargazer | Aldy | StellarPGx | SV? |
|------|:-----:|:---------:|:----:|:----------:|:---:|
| CYP2A6 | ✓ | ✓ | ✓ | ✓ | ✓ paralog |
| CYP2B6 | ✓ | ✓ | ✓ | ✓ | ✓ paralog |
| CYP2C8 | ✓ | ✓ | ✓ | ✓ | |
| CYP2C9 | ✓ | ✓ | ✓ | ✓ | |
| CYP2C19 | ✓ | ✓ | ✓ | ✓ | |
| CYP2D6 | ✓ | ✓ | ✓ | ✓ | ✓ paralog + tandem |
| CYP2E1 | ✓ | ✓ | ✓ | ✓ | ✓ CN |
| CYP3A4 | ✓ | ✓ | ✓ | ✓ | |
| CYP3A5 | ✓ | ✓ | ✓ | ✓ | |
| CYP4F2 | ✓ | ✓ | ✓ | ✓ | ✓ CN |
| SLCO1B1 | ✓ | ✓ | ✓ | ✓ | |
| UGT1A1 | ✓ | ✓ | ✓ | ✓ | |
| NAT1 | ✓ | ✓ | ✓ | ✓ | |
| NAT2 | ✓ | ✓ | ✓ | ✓ | |
| GSTM1 | ✓ | ✓ | ✓ | ✓ | ✓ deletion |
| NUDT15 | ✓ | — | ✓ | ✓ | |
| TPMT | ✓ | — | ✓ | ✓ | |
| CYP1A1 | — | ✓ | ✓ | ✓ | |
| CYP1A2 | — | ✓ | ✓ | ✓ | |
| DPYD | ✓ | — | ✓ | — | |
| G6PD | ✓ | ✓ | ✓ | — | ✓ CN |
| GSTT1 | ✓† | — | — | ✓ | ✓ deletion |
| VKORC1 | ✓ | ✓ | ✓ | — | |
| CYPOR/POR | — | ✓ | — | ✓ | |
| IFNL3 | ✓ | ✓ | ✓ | — | |
| RYR1 | ✓ | ✓ | ✓ | — | |

† GSTT1 is on `chr22_KI270879v1_alt` (alternate contig) — bcftools mpileup is skipped; PyPGx depth preprocessing may also fail on standard GRCh38 references.

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
│  Java 21 JRE           (Beagle phasing for PyPGx + Stargazer)│
│                                                                │
│  Nextflow              (copied from host binary)              │
│                                                                │
│  Apptainer             (Singularity fork, from PPA)           │
│  └── runs stellarpgx-dev.sif  [volume-mounted]                │
│      ├── graphtyper (graph-based variant calling)             │
│      ├── bcftools, samtools, tabix                            │
│      └── stellarpgx.py (star allele caller)                   │
│                                                                │
│  bcftools + samtools + tabix  (system packages)               │
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
├── Dockerfile                          # Ubuntu 22.04 image with all 4 tools
├── .dockerignore
├── nextflow                            # Nextflow binary (pre-downloaded)
├── docker/
│   ├── pgx-run.sh                      # Main orchestration entry point
│   ├── pgx-compare.py                  # Result parser and comparison table
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
| [`ToolsDocumentation.md`](ToolsDocumentation.md) | Full reference for all 4 callers: algorithms, gene lists, SV handling, output formats, limitations |
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

Stargazer and Aldy source code is bundled in this image for academic research use only. **Do not publish this image to Docker Hub, GHCR, or any public container registry.**
