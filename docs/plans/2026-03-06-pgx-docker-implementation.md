# PGx Suite Docker Container — Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Build and validate a single `pgx-suite` Docker container that installs PyPGx, Stargazer, Aldy, and StellarPGx (via Nextflow + Apptainer) all configured for GRCh38, and passes a smoke test for each tool using bundled test data.

**Architecture:** Ubuntu 22.04 base with Python 3.11 (deadsnakes PPA) + Java 11 JRE + Apptainer (Singularity fork, PPA) + Nextflow binary. PyPGx and Stargazer are copied from the local source tree and installed via pip. Aldy is pip-installed from PyPI. StellarPGx is NOT baked in — its repo directory is mounted as a volume at runtime, and Nextflow calls Apptainer to run `stellarpgx-dev.sif`. Requires `--privileged` for Apptainer's overlay filesystem.

**Tech Stack:** Docker, Ubuntu 22.04, Python 3.11, Java 11 JRE, pip, Apptainer 1.x, Nextflow 24.x, bcftools, samtools

**Reference:** See design doc at `docs/plans/2026-03-06-pgx-docker-design.md`

---

## Task 1: Create .dockerignore

**Files:**
- Create: `.dockerignore`

**Step 1: Write the file**

```
# Exclude large directories not needed in build context
StellarPGx/
pypgx/pypgx-bundle/
Stargazer_v2.0.3.zip
.pytest_cache/
docs/
*.log
work/
results/
.nextflow/
.nextflow.log
```

**Step 2: Verify the file exists**

```bash
cat /data/alvin/PGxCallers/.dockerignore
```
Expected: file contents printed

---

## Task 2: Create the Dockerfile

**Files:**
- Create: `Dockerfile`

**Step 1: Write the Dockerfile**

```dockerfile
FROM ubuntu:22.04

# Prevent interactive prompts during apt operations
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=UTC

# ── 1. Base system packages ─────────────────────────────────────────────────
RUN apt-get update && apt-get install -y \
        software-properties-common \
        curl wget git zip unzip \
        build-essential pkg-config \
        libssl-dev zlib1g-dev libbz2-dev libreadline-dev \
        libsqlite3-dev libffi-dev liblzma-dev libcurl4-openssl-dev \
        libncurses5-dev libhts-dev \
        fuse2fs fuse-overlayfs squashfs-tools squashfuse uidmap \
    && rm -rf /var/lib/apt/lists/*

# ── 2. Python 3.11 ──────────────────────────────────────────────────────────
RUN add-apt-repository ppa:deadsnakes/ppa \
    && apt-get update \
    && apt-get install -y \
        python3.11 python3.11-dev python3.11-venv python3.11-distutils \
    && rm -rf /var/lib/apt/lists/* \
    && curl -sS https://bootstrap.pypa.io/get-pip.py | python3.11 \
    && update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.11 1 \
    && update-alternatives --install /usr/bin/python  python  /usr/bin/python3.11 1

# ── 3. Java 11 JRE (for Beagle phasing in PyPGx and Stargazer) ──────────────
RUN apt-get update \
    && apt-get install -y openjdk-11-jre-headless \
    && rm -rf /var/lib/apt/lists/*

# ── 4. Bioinformatics CLI tools ──────────────────────────────────────────────
RUN apt-get update \
    && apt-get install -y bcftools samtools tabix \
    && rm -rf /var/lib/apt/lists/*

# ── 5. Apptainer (Singularity fork, for StellarPGx) ─────────────────────────
RUN add-apt-repository -y ppa:apptainer/ppa \
    && apt-get update \
    && apt-get install -y apptainer \
    && rm -rf /var/lib/apt/lists/*

# ── 6. Nextflow ──────────────────────────────────────────────────────────────
COPY nextflow /usr/local/bin/nextflow
RUN chmod +x /usr/local/bin/nextflow \
    && nextflow -version

# ── 7. Python pip dependencies (shared) ─────────────────────────────────────
RUN python3.11 -m pip install --no-cache-dir \
        pysam \
        pandas \
        numpy \
        matplotlib \
        scikit-learn \
        scipy \
        fuc \
        natsort \
        pyyaml \
        "importlib_resources>=1.3" \
        logbook \
        mappy \
        cython

# ── 8. ortools (for Aldy ILP solver) ────────────────────────────────────────
# ortools 9.x supports Python 3.11; pin to a known-good version
RUN python3.11 -m pip install --no-cache-dir "ortools>=9.6,<10"

# ── 9. Aldy ─────────────────────────────────────────────────────────────────
RUN python3.11 -m pip install --no-cache-dir aldy==4.8.3

# ── 10. PyPGx (from local source) ───────────────────────────────────────────
COPY pypgx/ /opt/pypgx/
RUN python3.11 -m pip install --no-cache-dir /opt/pypgx/

# ── 11. Stargazer (from local source) ───────────────────────────────────────
# NOTE: Stargazer is licensed for non-commercial academic use only (UW).
# This image must NOT be pushed to any public registry.
COPY stargazer-grc38-2.0.3/ /opt/stargazer/
RUN python3.11 -m pip install --no-cache-dir /opt/stargazer/

# ── 12. Helper scripts ───────────────────────────────────────────────────────
RUN mkdir -p /opt/pgx
COPY docker/test.sh     /opt/pgx/test.sh
RUN chmod +x /opt/pgx/test.sh

# ── 13. Runtime volume mount-points ─────────────────────────────────────────
RUN mkdir -p /pgx/bundle /pgx/stellarpgx /pgx/containers /pgx/ref /pgx/data /pgx/results

# ── 14. Environment variables ────────────────────────────────────────────────
ENV PYPGX_BUNDLE=/pgx/bundle
ENV STELLARPGX_DIR=/pgx/stellarpgx
ENV NXF_SINGULARITY_CACHEDIR=/pgx/containers
ENV JAVA_HOME=/usr/lib/jvm/java-11-openjdk-amd64
ENV PATH="${JAVA_HOME}/bin:${PATH}"

WORKDIR /pgx

CMD ["bash"]
```

**Step 2: Verify the file was created**

```bash
wc -l /data/alvin/PGxCallers/Dockerfile
```
Expected: line count printed

---

## Task 3: Create the Phase 1 smoke test script

**Files:**
- Create: `docker/test.sh`

**Step 1: Create the docker/ directory and write test.sh**

```bash
#!/usr/bin/env bash
# /opt/pgx/test.sh — Phase 1 smoke tests (no user BAM/reference needed)
# Run: docker run --privileged ... pgx-suite:latest bash /opt/pgx/test.sh
set -euo pipefail

PASS=0
FAIL=0

run_test() {
    local name="$1"
    shift
    printf "  %-50s" "$name"
    if "$@" > /tmp/_pgx_test_out 2>&1; then
        echo "PASS"
        ((PASS++)) || true
    else
        echo "FAIL"
        head -5 /tmp/_pgx_test_out | sed 's/^/    /'
        ((FAIL++)) || true
    fi
}

echo "============================================================"
echo " PGx Suite Docker Container — Smoke Tests"
echo " Reference build: GRCh38"
echo "============================================================"
echo ""

# ── PyPGx ─────────────────────────────────────────────────────────────────
echo "[1/4] PyPGx"
run_test "pypgx --version"   pypgx --version
run_test "pypgx import"      python3 -c "import pypgx; print('OK')"
run_test "pypgx list-genes"  pypgx list-genes

# Functional test: run-ngs-pipeline with bundled VCF test data (no BAM needed)
cd /tmp
unzip -qo /opt/pypgx/test-data/CYP4F2-GRCh38.zip -d pypgx_test
run_test "pypgx run-ngs-pipeline CYP4F2 (VCF-only)" \
    pypgx run-ngs-pipeline CYP4F2 /tmp/pypgx_test_out \
        --variants /tmp/pypgx_test/CYP4F2-GRCh38/data.vcf \
        --assembly GRCh38 \
        --force
cd /pgx
echo ""

# ── Stargazer ─────────────────────────────────────────────────────────────
echo "[2/4] Stargazer"
run_test "stargazer --version"   stargazer --version
run_test "stargazer import"      python3 -c "import stargazer; print('OK')" 2>/dev/null || \
    run_test "stargazer script"  python3 -c "import importlib.util; s=importlib.util.find_spec('stargazer'); print('OK' if s else 'NOTFOUND')"

# Functional test: run example WGS analysis (bundled VCF + GDF, no BAM needed)
cd /opt/stargazer/example
run_test "stargazer WGS example (CYP2D6)" \
    stargazer \
        -o /tmp/stargazer_test_out \
        -d wgs \
        -t cyp2d6 \
        -c vdr \
        -a grc38 \
        -i getrm-cyp2d6-vdr.joint.filtered.vcf \
        -g getrm-cyp2d6-vdr.gdf
cd /pgx
echo ""

# ── Aldy ──────────────────────────────────────────────────────────────────
echo "[3/4] Aldy"
run_test "aldy --version"  aldy --version
run_test "aldy import"     python3 -c "import aldy; print('OK')"
run_test "aldy built-in test suite (78 tests)" \
    aldy test
echo ""

# ── StellarPGx ────────────────────────────────────────────────────────────
echo "[4/4] StellarPGx"
run_test "nextflow --version"    nextflow -version
run_test "apptainer --version"   apptainer --version

if [[ -d "/pgx/stellarpgx" && -f "/pgx/containers/stellarpgx-dev.sif" ]]; then
    run_test "StellarPGx test pipeline (CYP2D6/hg38)" \
        nextflow run /pgx/stellarpgx/main.nf \
            -profile standard,test \
            --out_dir /tmp/stellarpgx_test_out
else
    echo "  [SKIP] StellarPGx volumes not mounted"
    echo "         Requires: -v \$(pwd)/StellarPGx:/pgx/stellarpgx"
    echo "                   -v \$(pwd)/StellarPGx/containers:/pgx/containers"
fi
echo ""

# ── Summary ───────────────────────────────────────────────────────────────
echo "============================================================"
echo " Results: ${PASS} PASSED  |  ${FAIL} FAILED"
echo "============================================================"
[[ "${FAIL}" -eq 0 ]]
```

**Step 2: Verify it's executable**

```bash
ls -la /data/alvin/PGxCallers/docker/test.sh
```
Expected: file listed with contents

---

## Task 4: Build the Docker image

**Step 1: Build from the PGxCallers directory**

```bash
cd /data/alvin/PGxCallers
docker build -t pgx-suite:latest .
```

Expected: Build completes successfully, ends with:
```
Successfully tagged pgx-suite:latest
```
Note: First build takes 10-15 minutes (downloads packages). Subsequent builds use cache.

**Step 2: Verify the image was created**

```bash
docker images pgx-suite
```
Expected: `pgx-suite   latest   <hash>   <size>`

**Step 3: Quick sanity check — all tools present**

```bash
docker run --rm pgx-suite:latest bash -c "
  echo 'Python:' && python3 --version
  echo 'Java:' && java -version 2>&1 | head -1
  echo 'PyPGx:' && pypgx --version
  echo 'Aldy:' && aldy --version 2>&1 | head -1
  echo 'Stargazer:' && stargazer --version 2>&1 | head -1
  echo 'Nextflow:' && nextflow -version 2>&1 | grep version
  echo 'Apptainer:' && apptainer --version
  echo 'bcftools:' && bcftools --version | head -1
  echo 'samtools:' && samtools --version | head -1
"
```
Expected: All version strings printed without errors.

---

## Task 5: Run Phase 1 smoke tests — PyPGx and Aldy (no volumes needed)

**Step 1: Run PyPGx + Aldy tests (no volume mounts required)**

```bash
docker run --rm pgx-suite:latest bash -c "
  set -e
  echo '=== PyPGx ==='
  pypgx --version
  python3 -c \"import pypgx; print('import OK')\"
  unzip -qo /opt/pypgx/test-data/CYP4F2-GRCh38.zip -d /tmp/pypgx_test
  pypgx run-ngs-pipeline CYP4F2 /tmp/pypgx_test_out \
      --variants /tmp/pypgx_test/CYP4F2-GRCh38/data.vcf \
      --assembly GRCh38 \
      --force
  echo 'PyPGx: PASS'

  echo '=== Aldy ==='
  aldy --version
  aldy test --gene CYP2D6
  echo 'Aldy: PASS'
"
```

Expected output includes:
- `pypgx X.X.X`
- `import OK`
- `Aldy PASS` or built-in test results showing gene genotype
- `Aldy: PASS`

**Step 2: If pypgx run-ngs-pipeline fails with 'no variants'**, check the VCF:

```bash
docker run --rm pgx-suite:latest bash -c "
  unzip -qo /opt/pypgx/test-data/CYP4F2-GRCh38.zip -d /tmp/pypgx_test
  cat /tmp/pypgx_test/CYP4F2-GRCh38/data.vcf
  cat /tmp/pypgx_test/CYP4F2-GRCh38/metadata.txt
"
```
Adjust the `run-ngs-pipeline` flags based on what the test VCF actually contains.

---

## Task 6: Run Phase 1 smoke tests — Stargazer (no volumes needed)

**Step 1: Run Stargazer example test (bundled in image)**

```bash
docker run --rm pgx-suite:latest bash -c "
  set -e
  echo '=== Stargazer ==='
  stargazer --version
  cd /opt/stargazer/example
  stargazer \
      -o /tmp/stargazer_test \
      -d wgs \
      -t cyp2d6 \
      -c vdr \
      -a grc38 \
      -i getrm-cyp2d6-vdr.joint.filtered.vcf \
      -g getrm-cyp2d6-vdr.gdf
  echo 'Stargazer: PASS'
  ls /tmp/stargazer_test/
"
```

Expected: Stargazer runs, creates output directory with genotype files.

**Step 2: Verify output was produced**

```bash
docker run --rm pgx-suite:latest bash -c "
  cd /opt/stargazer/example
  stargazer -o /tmp/sg_out -d wgs -t cyp2d6 -c vdr -a grc38 \
      -i getrm-cyp2d6-vdr.joint.filtered.vcf \
      -g getrm-cyp2d6-vdr.gdf
  cat /tmp/sg_out/genotype.txt 2>/dev/null || ls /tmp/sg_out/
"
```

---

## Task 7: Run Phase 1 smoke tests — StellarPGx (requires volume mounts)

StellarPGx needs the repository mounted (for scripts/database) AND the SIF container mounted (for Apptainer). The `standard,test` profile uses chr22-only reference + small BAM files bundled in the repo.

**Step 1: Run StellarPGx test pipeline**

```bash
cd /data/alvin/PGxCallers
docker run --privileged --rm \
  -v $(pwd)/StellarPGx:/pgx/stellarpgx \
  -v $(pwd)/StellarPGx/containers:/pgx/containers \
  pgx-suite:latest \
  bash -c "
    set -e
    echo '=== StellarPGx ==='
    nextflow -version
    apptainer --version
    cd /pgx/stellarpgx
    nextflow run main.nf \
        -profile standard,test \
        --out_dir /tmp/stellarpgx_test
    echo 'StellarPGx: PASS'
    ls /tmp/stellarpgx_test/ 2>/dev/null || true
  "
```

Expected: Nextflow runs the CYP2D6 test pipeline using the chr22 reference and HG03130 test BAM. Final output in `/tmp/stellarpgx_test/`.

**Step 2: If StellarPGx fails with Apptainer/Singularity permission errors:**

The container needs `--privileged` or `SYS_ADMIN` cap. Try:
```bash
docker run --cap-add SYS_ADMIN --security-opt seccomp=unconfined --rm \
  -v $(pwd)/StellarPGx:/pgx/stellarpgx \
  -v $(pwd)/StellarPGx/containers:/pgx/containers \
  pgx-suite:latest \
  bash -c "cd /pgx/stellarpgx && nextflow run main.nf -profile standard,test --out_dir /tmp/test_out"
```

**Step 3: If Nextflow can't find the SIF, update container path**

Check that the SIF path in nextflow.config matches `/pgx/containers/stellarpgx-dev.sif`:
```bash
docker run --rm -v $(pwd)/StellarPGx:/pgx/stellarpgx pgx-suite:latest \
  bash -c "grep 'container =' /pgx/stellarpgx/nextflow.config"
```
If path differs, override via `--with-singularity /pgx/containers/stellarpgx-dev.sif` in the nextflow run command.

---

## Task 8: Run full smoke test suite

**Step 1: Run test.sh for PyPGx + Aldy + Stargazer (no volumes)**

```bash
docker run --rm pgx-suite:latest bash /opt/pgx/test.sh 2>&1 | tee /tmp/pgx_smoke_results.txt
```

**Step 2: Run full suite including StellarPGx**

```bash
cd /data/alvin/PGxCallers
docker run --privileged --rm \
  -v $(pwd)/StellarPGx:/pgx/stellarpgx \
  -v $(pwd)/StellarPGx/containers:/pgx/containers \
  pgx-suite:latest \
  bash /opt/pgx/test.sh 2>&1 | tee /tmp/pgx_full_smoke.txt
```

Expected: Final line reads `X PASSED  |  0 FAILED`

---

## Task 9: Write a convenience docker-run.sh wrapper

**Files:**
- Create: `docker/docker-run.sh`

**Step 1: Write wrapper**

```bash
#!/usr/bin/env bash
# docker/docker-run.sh — Convenience wrapper for pgx-suite container
# Usage: ./docker/docker-run.sh [docker-run-args...] -- [container-command...]
# Example: ./docker/docker-run.sh -- bash /opt/pgx/test.sh
# Example: ./docker/docker-run.sh -v /data/ref:/pgx/ref -- bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

exec docker run --privileged --rm -it \
  -v "${PROJECT_DIR}/StellarPGx:/pgx/stellarpgx" \
  -v "${PROJECT_DIR}/StellarPGx/containers:/pgx/containers" \
  -v "${PROJECT_DIR}/pypgx/pypgx-bundle:/pgx/bundle" \
  "$@" \
  pgx-suite:latest
```

**Step 2: Make it executable**

```bash
chmod +x /data/alvin/PGxCallers/docker/docker-run.sh
```

**Step 3: Test the wrapper**

```bash
cd /data/alvin/PGxCallers
./docker/docker-run.sh -- bash -c "pypgx --version && aldy --version && stargazer --version && nextflow -version"
```

---

## Task 10: Document the image and write README

**Files:**
- Create: `docker/README.md`

**Step 1: Write README**

```markdown
# PGx Suite Docker Container

Single container for four pharmacogenomics star allele callers.
All tools configured for **GRCh38**.

## ⚠️ License Notice

Stargazer (UW) and Aldy (IURTC) are **non-commercial academic use only**.
This image must NOT be pushed to any public registry.

## Prerequisites

- Docker with `--privileged` support (local workstation)
- StellarPGx/containers/stellarpgx-dev.sif (already in repo)
- ~3 GB disk for image; ~50 GB for GRCh38 reference (user-supplied)

## Build

```bash
cd /data/alvin/PGxCallers
docker build -t pgx-suite:latest .
```

## Smoke Tests

```bash
# Without StellarPGx (fast, no volumes needed):
docker run --rm pgx-suite:latest bash /opt/pgx/test.sh

# Full suite including StellarPGx:
./docker/docker-run.sh -- bash /opt/pgx/test.sh
```

## Tools Installed

| Tool | Version | Invocation |
|------|---------|-----------|
| PyPGx | 0.26.0 | `pypgx` |
| Stargazer | 2.0.3 | `stargazer` |
| Aldy | 4.8.3 | `aldy` |
| StellarPGx | 1.2.7 | `nextflow run /pgx/stellarpgx/main.nf` |

## Volume Mounts (at runtime)

| Host | Container | Required? |
|------|-----------|----------|
| `./StellarPGx` | `/pgx/stellarpgx` | StellarPGx only |
| `./StellarPGx/containers` | `/pgx/containers` | StellarPGx only |
| `./pypgx/pypgx-bundle` | `/pgx/bundle` | PyPGx phasing |
| `/path/to/ref` | `/pgx/ref` | StellarPGx full runs |
| `/path/to/data` | `/pgx/data` | User BAM/CRAM input |
| `/path/to/results` | `/pgx/results` | Output |

## Phase 2 (coming)

`pgx-run.sh <GENE> <BAM>` — runs all supported tools and produces a
side-by-side comparison table. Awaiting user-supplied test BAM file.
```

---

## Troubleshooting Reference

| Symptom | Cause | Fix |
|---------|-------|-----|
| `FATAL: could not open image` | SIF not mounted | Add `-v $(pwd)/StellarPGx/containers:/pgx/containers` |
| `FATAL: kernel too old` | Need newer kernel / privilege | Run with `--privileged` |
| `ortools not found` | ortools pip install failed | Check Python version compatibility; pin `ortools>=9.6,<10` |
| Stargazer `ModuleNotFoundError` | Package not on path | Try `python3 -c "import stargazer"` to diagnose |
| `NXF_SINGULARITY_CACHEDIR` warning | Nextflow cache dir | Set in env or mount a writable dir |
| Beagle JAR missing | Java not found | Verify `java -version` in container |

---

## File Summary

Files created by this plan:
```
PGxCallers/
├── Dockerfile
├── .dockerignore
└── docker/
    ├── test.sh
    ├── docker-run.sh
    └── README.md
```

Phase 2 files (next plan, after BAM test file is supplied):
```
docker/
├── pgx-run.sh        (orchestration: gene + BAM → run all tools)
└── pgx-compare.py    (output parser + comparison table generator)
```
