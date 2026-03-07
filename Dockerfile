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
COPY docker/pgx-run.sh        /opt/pgx/pgx-run.sh
COPY docker/pgx-all-genes.sh  /opt/pgx/pgx-all-genes.sh
COPY docker/pgx-compare.py    /opt/pgx/pgx-compare.py
RUN chmod +x /opt/pgx/test.sh /opt/pgx/pgx-run.sh /opt/pgx/pgx-all-genes.sh \
    && ln -s /opt/pgx/pgx-run.sh       /usr/local/bin/pgx-run.sh \
    && ln -s /opt/pgx/pgx-all-genes.sh /usr/local/bin/pgx-all-genes.sh

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
