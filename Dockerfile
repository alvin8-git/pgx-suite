# ── Stage 1: Build ───────────────────────────────────────────────────────────
# Compile all Python packages against full build toolchain.
# Nothing from this stage except /opt/venv is carried into the final image.
FROM ubuntu:22.04 AS builder

ENV DEBIAN_FRONTEND=noninteractive TZ=UTC

# Build tools + Python 3.11 dev headers (needed to compile pysam and other C extensions)
RUN apt-get update && apt-get install -y \
        software-properties-common \
        curl \
        build-essential pkg-config \
        libssl-dev zlib1g-dev libbz2-dev libreadline-dev \
        libsqlite3-dev libffi-dev liblzma-dev libcurl4-openssl-dev \
        libncurses5-dev \
    && add-apt-repository ppa:deadsnakes/ppa \
    && apt-get update \
    && apt-get install -y \
        python3.11 python3.11-dev python3.11-venv python3.11-distutils \
    && rm -rf /var/lib/apt/lists/* \
    && curl -sS https://bootstrap.pypa.io/get-pip.py | python3.11

# Create a virtualenv so the compiled packages can be copied cleanly to the runtime stage
RUN python3.11 -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# Core Python dependencies
# scikit-learn is intentionally omitted: it is only used by pypgx's train_cnv_caller()
# (a model-training helper), not by run-ngs-pipeline or any step in pgx-run.sh.
RUN pip install --no-cache-dir \
        pysam \
        pandas \
        numpy \
        matplotlib \
        scipy \
        statsmodels \
        fuc \
        natsort \
        pyyaml \
        "importlib_resources>=1.3" \
        logbook \
        mappy \
        cython

# ortools (Aldy ILP solver) — pin to 9.x which supports Python 3.11
RUN pip install --no-cache-dir "ortools>=9.6,<10"

# Aldy
RUN pip install --no-cache-dir aldy==4.8.3

# PyPGx (from local source; .dockerignore strips build/ docs/ *.rst to keep context lean)
COPY pypgx/ /opt/pypgx/
RUN pip install --no-cache-dir /opt/pypgx/

# Snakemake — workflow orchestrator (replaces the bash phase/parallelism logic in
# pgx-run.sh / pgx-all-genes.sh). Added last to preserve build-cache on the heavy
# layers above.
RUN pip install --no-cache-dir snakemake==9.23.0   # pinned: orchestrator the pipeline depends on


# ── Stage 2: Runtime ─────────────────────────────────────────────────────────
# Lean image: no compiler, no -dev headers. Build tools saved ~300–350 MB.
FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive TZ=UTC

# ── 1. System runtime libraries + PPA tooling ────────────────────────────────
# software-properties-common is needed only for add-apt-repository (two PPAs below).
# It is purged in the same layer so it does not persist in the final image.
# Runtime shared libs (no -dev variants): pysam / ortools link against these at runtime.
RUN apt-get update && apt-get install -y \
        software-properties-common \
        curl wget \
        libssl3 zlib1g libbz2-1.0 libffi8 liblzma5 libcurl4 libncurses6 \
        fuse2fs fuse-overlayfs squashfs-tools squashfuse uidmap \
    && rm -rf /var/lib/apt/lists/*

# ── 2. Apptainer PPA (must run while system Python is still 3.10) ────────────
RUN add-apt-repository -y ppa:apptainer/ppa \
    && apt-get update \
    && apt-get install -y apptainer \
    && rm -rf /var/lib/apt/lists/*

# ── 3. Python 3.11 runtime (no -dev headers needed — packages come from venv) ─
RUN add-apt-repository ppa:deadsnakes/ppa \
    && apt-get update \
    && apt-get install -y python3.11 \
    && apt-get purge -y --autoremove software-properties-common \
    && rm -rf /var/lib/apt/lists/* \
    && update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.11 110 \
    && update-alternatives --install /usr/bin/python  python  /usr/bin/python3.11 110

# ── 4. Java 21 JRE + bioinformatics CLI tools ────────────────────────────────
RUN apt-get update && apt-get install -y \
        openjdk-21-jre-headless \
        bcftools samtools tabix \
    && rm -rf /var/lib/apt/lists/*

# mosdepth — fast BAM depth (single static binary, no runtime deps)
# Used by pgx-bamstats.sh for per-gene coverage and genome-wide depth.
RUN curl -sL https://github.com/brentp/mosdepth/releases/download/v0.3.12/mosdepth \
        -o /usr/local/bin/mosdepth \
    && chmod +x /usr/local/bin/mosdepth \
    && mosdepth --version

# sambamba — parallel BAM flagstat (true multi-threaded read counting)
# Splits the BAM into per-thread chunks; 6–10× faster than samtools flagstat.
# Note: does not support CRAM natively — pgx-bamstats.sh converts CRAM→BAM first.
RUN curl -sL https://github.com/biod/sambamba/releases/download/v1.0.1/sambamba-1.0.1-linux-amd64-static.gz \
        | gunzip > /usr/local/bin/sambamba \
    && chmod +x /usr/local/bin/sambamba \
    && sambamba --version 2>&1 | head -1

# mutserve v2 — mitochondrial variant caller (SNPs + heteroplasmy) for MT-RNR1
# Java 21 JRE is already present (step 4 above); no extra runtime needed.
# v2.0.x ships a ZIP (no bare-jar release asset — the old v2.0.0/mutserve.jar URL
# 404s, leaving a 9-byte "Not Found" stub). Extract mutserve.jar and fail loudly
# if the download is not a real jar.
RUN curl -fsSL https://github.com/seppinho/mutserve/releases/download/v2.0.3/mutserve.zip \
        -o /tmp/mutserve.zip \
    && python3 -c "import zipfile; zipfile.ZipFile('/tmp/mutserve.zip').extract('mutserve.jar', '/usr/local/bin/')" \
    && rm -f /tmp/mutserve.zip \
    && java -jar /usr/local/bin/mutserve.jar --version 2>&1 | head -1

# ── 5. Copy compiled Python environment from builder ─────────────────────────
# The venv contains pypgx, aldy, pysam, ortools, pandas, numpy, matplotlib, etc.
COPY --from=builder /opt/venv /opt/venv
ENV PATH="/opt/venv/bin:/usr/local/bin:${PATH}"

# ── 6. Nextflow ──────────────────────────────────────────────────────────────
COPY nextflow /usr/local/bin/nextflow
RUN chmod +x /usr/local/bin/nextflow \
    && nextflow -version

# ── 7. Stargazer ─────────────────────────────────────────────────────────────
# NOTE: Non-commercial academic licence (University of Washington).
#       Do NOT push this image to any public registry.
# .dockerignore strips unit_test/ (163 MB), example/ (42 MB), and 1kgp_vcf/hg19/ (48 MB).
# Only the grc38 phasing panel and source code are copied.
# Wrapper calls the venv Python explicitly to guarantee the right packages are found.
COPY stargazer-grc38-2.0.3/ /opt/stargazer/
RUN printf '#!/usr/bin/env bash\nexec /opt/venv/bin/python3 /opt/stargazer/stargazer/__main__.py "$@"\n' \
        > /usr/local/bin/stargazer \
    && chmod +x /usr/local/bin/stargazer

# ── 8. Helper scripts ────────────────────────────────────────────────────────
RUN mkdir -p /opt/pgx
COPY docker/test.sh           /opt/pgx/test.sh
COPY docker/pgx-compare.py    /opt/pgx/pgx-compare.py
COPY docker/pgx_altcheck.py   /opt/pgx/pgx_altcheck.py
COPY docker/pgx-bamstats.sh   /opt/pgx/pgx-bamstats.sh
COPY docker/pgx-report.py     /opt/pgx/pgx-report.py
COPY docker/pgx-hla.sh            /opt/pgx/pgx-hla.sh
COPY docker/pgx-mt.sh             /opt/pgx/pgx-mt.sh
COPY docker/pgx_cram_regions.bed  /opt/pgx/pgx_cram_regions.bed
COPY docker/genes.tsv             /opt/pgx/genes.tsv
COPY docker/cpic.json             /opt/pgx/cpic.json
COPY docker/reconcile.py          /opt/pgx/reconcile.py
COPY docker/allele_synonyms.json  /opt/pgx/allele_synonyms.json
COPY docker/vcf_check_variants.json /opt/pgx/vcf_check_variants.json
COPY docker/pharmcat_positions.vcf.bgz     /opt/pgx/pharmcat_positions.vcf.bgz
COPY docker/pharmcat_positions.vcf.bgz.csi /opt/pgx/pharmcat_positions.vcf.bgz.csi
COPY docker/Snakefile             /opt/pgx/Snakefile
# Cyrius (vendored, PolyForm Strict 1.0.0 — use only, do NOT redistribute the image)
COPY Cyrius/                      /opt/cyrius/
# Pipeline version stamp for the provenance block (override with
# --build-arg PGX_VERSION=$(git rev-parse --short HEAD)).
ARG PGX_VERSION=dev
RUN echo "pgx-suite ${PGX_VERSION} (image built $(date -u +%Y-%m-%dT%H:%M:%SZ))" > /opt/pgx/VERSION
RUN chmod +x /opt/pgx/test.sh \
             /opt/pgx/pgx-bamstats.sh /opt/pgx/pgx-report.py /opt/pgx/pgx-hla.sh \
             /opt/pgx/pgx-mt.sh \
    && ln -s /opt/pgx/pgx-bamstats.sh  /usr/local/bin/pgx-bamstats.sh \
    && ln -s /opt/pgx/pgx-report.py    /usr/local/bin/pgx-report.py \
    && ln -s /opt/pgx/pgx-hla.sh       /usr/local/bin/pgx-hla.sh \
    && ln -s /opt/pgx/pgx-mt.sh        /usr/local/bin/pgx-mt.sh

# ── 9. Runtime volume mount-points ───────────────────────────────────────────
RUN mkdir -p /pgx/bundle /pgx/stellarpgx /pgx/containers /pgx/ref /pgx/data /pgx/results

# ── 10. Environment ──────────────────────────────────────────────────────────
ENV PYPGX_BUNDLE=/pgx/bundle
ENV STELLARPGX_DIR=/pgx/stellarpgx
ENV NXF_SINGULARITY_CACHEDIR=/pgx/containers
ENV JAVA_HOME=/usr/lib/jvm/java-21-openjdk-amd64
ENV PATH="${JAVA_HOME}/bin:${PATH}"

WORKDIR /pgx
CMD ["bash"]
