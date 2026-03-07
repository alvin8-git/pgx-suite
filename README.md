# PGx Suite

> A Docker container that installs and runs four open-source pharmacogenomics (PGx) star allele callers — **PyPGx**, **Stargazer**, **Aldy**, and **StellarPGx** — all configured for **GRCh38 (hg38)**.

---

## Tools

| Tool | Version | Method | License |
|------|---------|--------|---------|
| [PyPGx](https://github.com/sbslee/pypgx) | 0.26.0 | pip (local source) | MIT |
| [Stargazer](https://stargazer.gs.washington.edu/stargazerweb/) | 2.0.3 | wrapper script | **Non-commercial academic only (UW)** |
| [Aldy](https://github.com/0xTCG/aldy) | 4.8.3 | pip | **Non-commercial academic only (IURTC)** |
| [StellarPGx](https://github.com/SBIMB/StellarPGx) | 1.2.7 | Nextflow + Apptainer SIF | MIT |

> **⚠️ License Notice:** Stargazer and Aldy are restricted to non-commercial academic use.
> This image **must not be published to any public registry** (Docker Hub, GHCR, etc.).

---

## Prerequisites

- Docker Engine (Linux) or Docker Desktop (macOS/Windows) with `--privileged` support
- ~3 GB disk for the built image
- The following local directories (already present in this repo):
  - `pypgx/` — PyPGx source
  - `pypgx/pypgx-bundle/` — 500 MB 1KGP VCF panel + CNV data (volume-mounted at runtime)
  - `stargazer-grc38-2.0.3/` — Stargazer GRCh38 source
  - `StellarPGx/` — Nextflow pipeline, database, resources
  - `StellarPGx/containers/stellarpgx-dev.sif` — 31 MB Apptainer SIF image
  - `nextflow` — Nextflow binary (already downloaded)
- **User-supplied** (for Phase 2 full analysis):
  - GRCh38 reference FASTA + `.fai` index (3–50 GB, not included)
  - Input BAM/CRAM file(s) + index

---

## Quick Start

### 1. Build the image

```bash
cd /data/alvin/PGxCallers
docker build -t pgx-suite:latest .
```

First build: ~15–25 minutes. Subsequent builds use Docker layer cache.

### 2. Run smoke tests

```bash
# PyPGx, Aldy, Stargazer — fast, no data volumes required:
docker run --rm pgx-suite:latest bash /opt/pgx/test.sh

# Full suite including StellarPGx (volumes required):
./docker/docker-run.sh bash /opt/pgx/test.sh
```

Expected output (all PASS):

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

### 3. Interactive shell

```bash
./docker/docker-run.sh bash
```

Or manually:

```bash
docker run --privileged --rm -it \
  -v "$(pwd)/StellarPGx:/pgx/stellarpgx" \
  -v "$(pwd)/StellarPGx/containers:/pgx/containers" \
  -v "$(pwd)/pypgx/pypgx-bundle:/pgx/bundle" \
  pgx-suite:latest bash
```

---

## Volume Mounts

| Host Path | Container Path | Purpose |
|-----------|---------------|---------|
| `./StellarPGx` | `/pgx/stellarpgx` | StellarPGx pipeline scripts, database, resources |
| `./StellarPGx/containers` | `/pgx/containers` | `stellarpgx-dev.sif` Apptainer image |
| `./pypgx/pypgx-bundle` | `/pgx/bundle` | 1KGP VCFs + CNV data for PyPGx phasing |
| `/path/to/ref` | `/pgx/ref` | GRCh38 reference FASTA + `.fai` |
| `/path/to/data` | `/pgx/data` | Input BAM/CRAM + index |
| `/path/to/results` | `/pgx/results` | Analysis output |

---

## Gene Support Matrix

| Gene | PyPGx | Stargazer | Aldy | StellarPGx |
|------|:-----:|:---------:|:----:|:----------:|
| CYP2D6 | ✓ | ✓ | ✓ | ✓ |
| CYP2C19 | ✓ | ✓ | ✓ | ✓ |
| CYP2C9 | ✓ | ✓ | ✓ | ✓ |
| CYP2B6 | ✓ | ✓ | ✓ | ✓ |
| CYP2C8 | ✓ | ✓ | ✓ | ✓ |
| CYP3A4 | ✓ | ✓ | ✓ | ✓ |
| CYP3A5 | ✓ | ✓ | ✓ | ✓ |
| CYP4F2 | ✓ | ✓ | ✓ | ✓ |
| SLCO1B1 | ✓ | ✓ | ✓ | ✓ |
| NUDT15 | ✓ | — | ✓ | ✓ |
| TPMT | ✓ | — | ✓ | ✓ |
| UGT1A1 | ✓ | — | ✓ | ✓ |
| DPYD | ✓ | — | ✓ | — |
| NAT1 | ✓ | — | — | ✓ |
| NAT2 | ✓ | — | — | ✓ |
| G6PD | ✓ | ✓ | — | — |
| GSTM1 | ✓ | — | — | ✓ |
| GSTT1 | — | — | — | ✓ |
| POR/CYPOR | — | — | — | ✓ |
| VKORC1 | ✓ | ✓ | — | — |
| CYP1A1 | — | ✓ | ✓ | ✓ |
| CYP1A2 | — | ✓ | ✓ | ✓ |
| CYP2A6 | — | ✓ | ✓ | ✓ |
| CYP2E1 | — | ✓ | — | ✓ |
| IFNL3 | ✓ | — | — | — |
| RYR1 | ✓ | — | — | — |

---

## Output Field Reference

The table below maps equivalent output fields across all four tools (N = 17 rows of output fields, 4 columns of tools). Field names use each tool's native terminology; `—` means the field is not reported by that tool.

| # | Output Field | PyPGx (`data.tsv`) | Stargazer (`report.tsv` / `genotype-calls.tsv`) | Aldy (`.aldy`) | StellarPGx (`.alleles`) |
|---|---|---|---|---|---|
| 1 | **Sample ID** | `Sample` | `name` / `Sample` (report) | `Sample` | Filename stem (e.g. `HG03130`) |
| 2 | **Gene** | Implicit (one run per gene) | `Gene` (report) | `Gene` | Header line in file |
| 3 | **Diplotype** | `Genotype` (e.g. `*2/*4`) | `Diplotype` (report) / `hap1_main`+`hap2_main` | `Major` (e.g. `*2+*4`) | `Result` (e.g. `*17/*29`) |
| 4 | **Haplotype 1** | `Haplotype1` (e.g. `*2;`) | `hap1_main` | First allele in `Major` | First allele in `Result` |
| 5 | **Haplotype 2** | `Haplotype2` (e.g. `*4;*10;*74;`) | `hap2_main` | Second allele in `Major` | Second allele in `Result` |
| 6 | **Sub-alleles / tag variants** | Semicolon list in `Haplotype1`/`Haplotype2` (e.g. `*4;*10;*74;`) | `hap1_main_core`+`hap1_main_tag`; `hap2_main_core`+`hap2_main_tag` | `Minor` column | `Candidate alleles` (raw haplotype codes) |
| 7 | **Alternative diplotypes** | `AlternativePhase` (semicolon list) | `dip_cand` / `May also be` (report) | Multiple `SolutionID` rows (each is a valid alternative) | — |
| 8 | **Phenotype** | `Phenotype` (e.g. `Intermediate Metabolizer`) | `Phenotype` (e.g. `intermediate_metabolizer`) | `Status` | `Metaboliser status` (e.g. `Intermediate metaboliser (IM)`) |
| 9 | **Activity score** | `ActivityScore` | `Score` (report) | — | `Activity score` |
| 10 | **SV / CNV type** | `CNV` (e.g. `Normal`, `Deletion`, `Duplication`) | `dip_sv` / `hap1_sv` / `hap2_sv` (e.g. `no_sv`, `dup`) | Implicit in `Major` allele name + `Copy` | `Initially computed CN` (integer) |
| 11 | **Copy number** | Derived from `CNV` | Derived from `dip_sv` | `Copy` column (integer per allele) | `Initially computed CN` |
| 12 | **Supporting variants** | `VariantData` (`allele:pos-ref-alt:AF;…`) | `hap1_main_core`/`hap2_main_core` (pos, AD, AF, effect) | `Location` + `Coverage` columns | `Sample core variants` (`pos~ref>alt~GT`) |
| 13 | **Functional effect** | Embedded in `VariantData` | Embedded in `hap1_main_core` (e.g. `splice_acceptor_variant`) | `Effect` + `Type` columns | — |
| 14 | **dbSNP rsID** | Embedded in `VariantData` (where available) | — | `dbSNP` column | — |
| 15 | **Allele score / confidence** | — | `dip_score`, `hap1_score`, `hap2_score` | `SolutionID` rank (lower = better fit) | — |
| 16 | **Mean allele fraction** | Per-variant AF embedded in `VariantData` | `hap1_af_mean_gene` / `hap2_af_mean_gene` / `hap1_af_mean_main` / `hap2_af_mean_main` | `Coverage` (read depth per variant) | — |
| 17 | **Phasing method** | Beagle statistical (1KGP panel) | Beagle; `BEAGLE imputed` flag (report) + `ssr` marker | ILP joint optimisation (no separate phase step) | Graph-based (graphtyper) |

---

## Container Architecture

```
docker run --privileged pgx-suite:latest
┌───────────────────────────────────────────────────────────────┐
│  Ubuntu 22.04                                                  │
│                                                                │
│  Python 3.11                                                   │
│  ├── pypgx 0.26.0      (pip from copied source)               │
│  ├── aldy 4.8.3        (pip install)                          │
│  └── shared deps: pysam, pandas, numpy, scipy, ortools, etc. │
│                                                                │
│  Stargazer 2.0.3       (/usr/local/bin/stargazer wrapper)     │
│                                                                │
│  Java 21 JRE           (Beagle phasing in PyPGx + Stargazer) │
│                                                                │
│  Nextflow              (copied from host binary)              │
│                                                                │
│  Apptainer             (Singularity fork, from PPA)           │
│  └── runs stellarpgx-dev.sif [volume-mounted]                 │
│      ├── graphtyper (variant calling)                         │
│      ├── bcftools, samtools, tabix                            │
│      └── stellarpgx.py (star allele caller)                   │
│                                                                │
│  bcftools + samtools + tabix  (system packages)               │
└───────────────────────────────────────────────────────────────┘
```

### Why `--privileged`?

StellarPGx runs Nextflow processes inside `stellarpgx-dev.sif` via Apptainer (the open-source Singularity fork). Running Apptainer inside Docker requires `SYS_ADMIN` capability to unpack SIF overlay filesystems — `--privileged` is the simplest way to grant this on a local workstation.

---

## Phase 2 — `pgx-run.sh` (coming soon)

A single command will run all supported tools for a given gene and produce a side-by-side comparison table:

```bash
# Example (awaiting test BAM file):
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

Expected output:

```
===== PGx Star Allele Results =================================
Gene:   CYP2D6
Sample: sample
Build:  GRCh38
──────────────────────────────────────────────────────────────
Tool         Diplotype    Activity Score   Phenotype
──────────────────────────────────────────────────────────────
PyPGx        *1/*4        1.0              Normal Metabolizer
Stargazer    *1/*4        1.0              Normal Metabolizer
Aldy         *2/*4        1.0              Normal Metabolizer
StellarPGx   *1/*4        -                -
──────────────────────────────────────────────────────────────
Concordant: 3/4 tools agree on *1/*4
Full results saved to: /pgx/results/CYP2D6_sample_comparison.tsv
```

See `TODO.md` for the full Phase 2 roadmap.

---

## Troubleshooting

| Symptom | Cause | Fix |
|---------|-------|-----|
| `FATAL: could not open image` | SIF not mounted | Add `-v $(pwd)/StellarPGx/containers:/pgx/containers` |
| `FATAL: kernel too old` | Need privilege | Run with `--privileged` |
| `ModuleNotFoundError: apt_pkg` during build | Build order issue | Ensure Apptainer PPA step is before Python 3.11 switch in Dockerfile |
| `pypgx` not found | Python path | Verify `python3 --version` is 3.11 in container |
| Nextflow hangs at startup | JAR download | Ensure internet access on first `nextflow` invocation |
| `beagle.jar` error | Java not in PATH | Verify `java -version` in container; `JAVA_HOME=/usr/lib/jvm/java-21-openjdk-amd64` |
| Stargazer: `no data for target gene` with GDF | hg19-aligned GDF on grc38 mode | Use VCF-only mode (omit `-c` and `-g` flags) |

---

## Repository Structure

```
PGxCallers/
├── Dockerfile                         # Multi-stage Ubuntu 22.04 build
├── .dockerignore                      # Exclude large data dirs from build context
├── nextflow                           # Nextflow binary (pre-downloaded)
├── docker/
│   ├── test.sh                        # Phase 1 smoke tests
│   ├── docker-run.sh                  # Convenience run wrapper
│   └── README.md                      # Container-specific docs
├── pypgx/                             # PyPGx source (pip-installable)
│   └── pypgx-bundle/                  # 500MB 1KGP VCF panel (volume-mounted)
├── stargazer-grc38-2.0.3/             # Stargazer GRCh38 source + example data
├── StellarPGx/                        # Nextflow pipeline (volume-mounted)
│   ├── containers/
│   │   └── stellarpgx-dev.sif         # Apptainer image (31 MB)
│   ├── database/                      # CYP450 star allele database
│   └── resources/                     # Reference sequences, panel data
├── docs/
│   └── plans/
│       ├── 2026-03-06-pgx-docker-design.md
│       └── 2026-03-06-pgx-docker-implementation.md
├── README.md                          # This file
└── TODO.md                            # Phase 2 roadmap
```

---

## Development Notes

- Stargazer cannot be installed with `pip install` because its source tree lacks `__init__.py` and uses Python 2-style bare imports. The wrapper at `/usr/local/bin/stargazer` runs `python3 /opt/stargazer/stargazer/__main__.py "$@"` directly, which puts the script directory in `sys.path[0]` and resolves all internal imports correctly.
- Java 21 is required (not 17) because Nextflow 25.x dropped support for Java < 17, and Ubuntu 22.04's `openjdk-17-jre-headless` had compatibility issues with Nextflow's JAR resolution.
- The Apptainer PPA (`add-apt-repository`) must run *before* switching the Python 3 default to 3.11, because `add-apt-repository` depends on `python3/apt_pkg` which is compiled for the system Python (3.10).
