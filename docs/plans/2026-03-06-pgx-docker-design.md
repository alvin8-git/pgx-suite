# PGx Suite Docker Container — Design Document

**Date:** 2026-03-06
**Reference genome:** GRCh38 (hg38) — all tools configured for GRCh38 only
**Target deployment:** Local workstation (Docker Desktop / Docker Engine)

---

## 1. Goal

Build a single Docker container that installs and runs four open-source pharmacogenomics (PGx) star allele callers:

| Tool | Version | License |
|------|---------|---------|
| PyPGx | 0.26.0 | MIT |
| Stargazer | 2.0.3 (GRCh38 build) | **Non-commercial academic only (UW)** |
| Aldy | 4.8.3 | Non-commercial academic (IURTC) |
| StellarPGx | 1.2.7 | MIT |

**⚠️ Licensing:** Stargazer and Aldy are restricted to non-commercial/academic use. This image **must not be published to any public registry** (Docker Hub, GHCR, etc.).

**Phased delivery:**
- **Phase 1** *(immediate)*: Build container, install all 4 tools, run bundled smoke tests for each.
- **Phase 2** *(after BAM test file supplied)*: `pgx-run.sh <gene> <bam>` orchestration + comparison table output.

---

## 2. Container Architecture

### Base Image

`ubuntu:22.04` — required because Apptainer (Singularity) needs glibc and kernel-namespace tools not available on Alpine or slim images.

### Why `--privileged` is needed

StellarPGx uses **Nextflow + Singularity (Apptainer)** internally. Running Singularity inside Docker requires unpacking SIF overlay filesystems, which needs `SYS_ADMIN` capability. On a local workstation, `--privileged` is the simplest approach.

### Component Stack

```
docker run --privileged pgx-suite:latest
┌────────────────────────────────────────────────────────────────┐
│  Ubuntu 22.04                                                  │
│                                                                │
│  Python 3.11                                                   │
│  ├── pypgx 0.26.0      (pip from copied source)               │
│  ├── aldy 4.8.3        (pip install aldy)                      │
│  ├── stargazer 2.0.3   (pip from copied source)               │
│  └── shared deps: pysam, pandas, numpy, matplotlib,           │
│                   scikit-learn, scipy, ortools, mappy          │
│                                                                │
│  Java 11 JRE           (Beagle phasing in pypgx + stargazer)  │
│                                                                │
│  Nextflow              (copied from host binary)              │
│                                                                │
│  Apptainer (Singularity fork)                                  │
│  └── runs stellarpgx-dev.sif at runtime via Nextflow           │
│      ├── graphtyper (variant calling)                          │
│      ├── bcftools, samtools, tabix (VCF/BAM manipulation)     │
│      └── stellarpgx.py (star allele caller scripts)           │
│                                                                │
│  bcftools + samtools   (for variant calling shared VCF)       │
│                                                                │
│  /usr/local/bin/pgx-run.sh   (Phase 2 orchestration script)  │
│  /usr/local/bin/pgx-compare.py (Phase 2 result parser)        │
└────────────────────────────────────────────────────────────────┘
```

---

## 3. Build Context & File Layout

The Dockerfile lives at the root of `/data/alvin/PGxCallers/` and uses the following from the build context:

```
PGxCallers/
├── Dockerfile                  ← to be created
├── docker/
│   ├── pgx-run.sh              ← Phase 2 entrypoint (to be created)
│   └── pgx-compare.py          ← Phase 2 result parser (to be created)
├── nextflow                    ← binary copied into image
├── pypgx/                      ← COPIED into image at /opt/pypgx
├── stargazer-grc38-2.0.3/      ← COPIED into image at /opt/stargazer
└── StellarPGx/                 ← MOUNTED as volume at runtime
```

**Volumes mounted at `docker run` time:**

| Host Path | Container Path | Size | Contents |
|-----------|---------------|------|----------|
| `./pypgx/pypgx-bundle` | `/pgx/bundle` | 500MB | 1KGP VCFs + CNV data |
| `./StellarPGx` | `/pgx/stellarpgx` | ~100MB | DB + resources + scripts |
| `./StellarPGx/containers` | `/pgx/containers` | 31MB | stellarpgx-dev.sif |
| user-supplied | `/pgx/ref` | 3–50GB | GRCh38 reference FASTA + .fai |
| user-supplied | `/pgx/data` | varies | input BAM/CRAM + index |
| user-supplied | `/pgx/results` | varies | output directory |

---

## 4. Tool Installation Methods

### PyPGx
- Source copied from `./pypgx/` → `/opt/pypgx/`
- Installed: `pip install /opt/pypgx/`
- Bundle path configured via env var: `PYPGX_BUNDLE=/pgx/bundle`

### Stargazer
- Source copied from `./stargazer-grc38-2.0.3/` → `/opt/stargazer/`
- Installed: `pip install /opt/stargazer/`
- Invoked as: `python /opt/stargazer/stargazer`

### Aldy
- Installed: `pip install aldy==4.8.3`
- Bundled gene databases included in pip package

### StellarPGx
- Nextflow binary: copied from `./nextflow` → `/usr/local/bin/nextflow`
- Apptainer: installed from `ppa:apptainer/ppa` (Ubuntu 22.04)
- Pipeline files: mounted volume `/pgx/stellarpgx` (not baked into image)
- SIF container: `/pgx/containers/stellarpgx-dev.sif` (mounted)
- Nextflow config: `singularity.enabled=true`, `cacheDir=/pgx/containers`

---

## 5. Phase 1 — Smoke Tests

Each tool's bundled test data is used to verify the installation. No user BAM file required.

| Tool | Smoke Test | Expected Output |
|------|-----------|----------------|
| PyPGx | `pypgx --version` + `pypgx --help` | Version string |
| Aldy | `aldy test` (built-in test command) | PASS message |
| Stargazer | `python /opt/stargazer/stargazer --version` | Version string |
| StellarPGx | `nextflow run /pgx/stellarpgx/main.nf -profile standard,test` | Nextflow test pipeline success |

A `test.sh` script will run all four in sequence and report PASS/FAIL.

---

## 6. Phase 2 — `pgx-run.sh` Pipeline Design

### Usage
```bash
pgx-run.sh <GENE> <BAM_FILE> [--ref /pgx/ref/hg38.fa] [--output /pgx/results]
# Example:
pgx-run.sh CYP2D6 /pgx/data/sample.bam
```

### Gene Support Matrix (from user-supplied reference table)

```
Gene       PyPGx  Stargazer  Aldy  StellarPGx
---------  -----  ---------  ----  ----------
CYP2D6       ✓        ✓        ✓       ✓
CYP2C19      ✓        ✓        ✓       ✓
CYP2C9       ✓        ✓        ✓       ✓
CYP2B6       ✓        ✓        ✓       ✓
CYP2C8       ✓        ✓        ✓       ✓
CYP3A4       ✓        ✓        ✓       ✓
CYP3A5       ✓        ✓        ✓       ✓
CYP4F2       ✓        ✓        ✓       ✓
NUDT15       ✓        -        ✓       ✓
TPMT         ✓        -        ✓       ✓
UGT1A1       ✓        -        ✓       ✓
SLCO1B1      ✓        ✓        ✓       ✓
DPYD         ✓        -        ✓       -
NAT1         ✓        -        -       ✓
NAT2         ✓        -        -       ✓
G6PD         ✓        ✓        -       -
GSTM1        ✓        -        -       ✓
GSTT1        -        -        -       ✓
POR/CYPOR    -        -        -       ✓
VKORC1       ✓        ✓        -       -
CYP1A1       -        ✓        ✓       ✓
CYP1A2       -        ✓        ✓       ✓
CYP2A6       -        ✓        ✓       ✓
CYP2E1       -        ✓        -       ✓
IFNL3        ✓        -        -       -
RYR1         ✓        -        -       -
... (full table in pgx-compare.py)
```

### Pipeline Steps

```
pgx-run.sh CYP2D6 sample.bam
     │
     ▼
Step 1: Validate inputs (BAM indexed? Gene supported by ≥1 tool?)
     │
     ▼
Step 2: Shared preprocessing
     - bcftools mpileup | bcftools call → sample.gene.vcf.gz
       (gene-region only; feeds Stargazer and PyPGx)
     │
     ├──► Step 3a: PyPGx (if gene supported)
     │    pypgx run-ngs-pipeline --gene GENE --assembly GRCh38
     │    --variants sample.gene.vcf.gz --output results/pypgx/
     │
     ├──► Step 3b: Stargazer (if gene supported)
     │    python stargazer -t gene -d wgs -a grc38
     │    -i sample.gene.vcf.gz -o results/stargazer/ [-G gdf -B bam]
     │
     ├──► Step 3c: Aldy (if gene supported)
     │    aldy genotype -g GENE -p hg38 sample.bam
     │    -o results/aldy/
     │
     └──► Step 3d: StellarPGx (if gene supported)
          nextflow run /pgx/stellarpgx/main.nf
          --gene gene --in_bam "sample*{bam,bai}"
          --ref_file /pgx/ref/hg38.fa
          --out_dir results/stellarpgx/
     │
     ▼
Step 4: pgx-compare.py — parse all outputs, build table
     │
     ▼
Step 5: Print comparison table + save results/comparison.tsv
```

### Output Format
```
===== PGx Star Allele Results =================================
Gene:   CYP2D6
Sample: NA12878
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
──────────────────────────────────────────────────────────────
Full results saved to: results/CYP2D6_NA12878_comparison.tsv
```

---

## 7. Docker Run Command

```bash
# Phase 1 — Smoke tests (no data needed)
docker run --privileged --rm pgx-suite:latest bash /opt/test.sh

# Phase 2 — Full gene analysis
docker run --privileged --rm \
  -v $(pwd)/pypgx/pypgx-bundle:/pgx/bundle \
  -v $(pwd)/StellarPGx:/pgx/stellarpgx \
  -v $(pwd)/StellarPGx/containers:/pgx/containers \
  -v /path/to/ref:/pgx/ref \
  -v /path/to/data:/pgx/data \
  -v /path/to/results:/pgx/results \
  pgx-suite:latest \
  pgx-run.sh CYP2D6 /pgx/data/sample.bam
```

---

## 8. Estimated Image Size

| Layer | Compressed Size |
|-------|----------------|
| Ubuntu 22.04 | ~30 MB |
| Python 3.11 + pip | ~50 MB |
| Python packages (aldy + pypgx + stargazer deps) | ~350 MB |
| Java 11 JRE | ~160 MB |
| Apptainer binary | ~50 MB |
| Nextflow binary | ~50 MB |
| bcftools + samtools + tabix | ~20 MB |
| pypgx source + stargazer source (copied in) | ~5 MB |
| **Total estimate** | **~715 MB compressed** |

Reference data (volumes, not in image): ~600 MB already on disk + user reference genome.

---

## 9. Files To Be Created

```
PGxCallers/
├── Dockerfile
├── .dockerignore
├── docker/
│   ├── test.sh              (Phase 1 smoke tests)
│   ├── pgx-run.sh           (Phase 2 entrypoint)
│   └── pgx-compare.py       (Phase 2 result parser)
└── docs/
    └── plans/
        └── 2026-03-06-pgx-docker-design.md  ← this file
```
