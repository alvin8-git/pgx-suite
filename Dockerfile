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
        libncurses5-dev \
        fuse2fs fuse-overlayfs squashfs-tools squashfuse uidmap \
    && rm -rf /var/lib/apt/lists/*

# ── 2. Apptainer PPA + install (must run BEFORE Python 3.11 default change,
#       because add-apt-repository uses python3/apt_pkg compiled for 3.10) ───
RUN add-apt-repository -y ppa:apptainer/ppa \
    && apt-get update \
    && apt-get install -y apptainer \
    && rm -rf /var/lib/apt/lists/*

# ── 3. Python 3.11 (switch default AFTER PPA steps that need python3.10) ────
RUN add-apt-repository ppa:deadsnakes/ppa \
    && apt-get update \
    && apt-get install -y \
        python3.11 python3.11-dev python3.11-venv python3.11-distutils \
    && rm -rf /var/lib/apt/lists/* \
    && curl -sS https://bootstrap.pypa.io/get-pip.py | python3.11 \
    && update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.11 110 \
    && update-alternatives --install /usr/bin/python  python  /usr/bin/python3.11 110

# ── 4. Java 21 JRE (for Beagle phasing in PyPGx and Stargazer) ──────────────
#       + bioinformatics CLI tools (merged to reduce layers)
RUN apt-get update \
    && apt-get install -y \
        openjdk-21-jre-headless \
        bcftools samtools tabix \
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
# Stargazer uses Python 2-style bare imports that only resolve when __main__.py
# is run directly as a script (sys.path[0] = script dir). Use a wrapper instead of pip install.
COPY stargazer-grc38-2.0.3/ /opt/stargazer/
RUN printf '#!/usr/bin/env bash\nexec python3 /opt/stargazer/stargazer/__main__.py "$@"\n' \
        > /usr/local/bin/stargazer \
    && chmod +x /usr/local/bin/stargazer

# ── 12. Helper scripts ───────────────────────────────────────────────────────
RUN mkdir -p /opt/pgx
COPY docker/test.sh        /opt/pgx/test.sh
COPY docker/pgx-run.sh     /opt/pgx/pgx-run.sh
COPY docker/pgx-compare.py /opt/pgx/pgx-compare.py
RUN chmod +x /opt/pgx/test.sh /opt/pgx/pgx-run.sh \
    && ln -s /opt/pgx/pgx-run.sh /usr/local/bin/pgx-run.sh

# ── 13. Runtime volume mount-points ─────────────────────────────────────────
RUN mkdir -p /pgx/bundle /pgx/stellarpgx /pgx/containers /pgx/ref /pgx/data /pgx/results

# ── 14. Environment variables ────────────────────────────────────────────────
ENV PYPGX_BUNDLE=/pgx/bundle
ENV STELLARPGX_DIR=/pgx/stellarpgx
ENV NXF_SINGULARITY_CACHEDIR=/pgx/containers
ENV JAVA_HOME=/usr/lib/jvm/java-21-openjdk-amd64
ENV PATH="${JAVA_HOME}/bin:${PATH}"

WORKDIR /pgx

CMD ["bash"]
