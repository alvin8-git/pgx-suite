#!/usr/bin/env bash
# Build StellarPGx/containers/arcashla.sif (REQUIRED for HLA class II: DQA1/DQB1/DRB1).
#
# The arcas-hla biocontainer is not usable as-is: it ships NO IMGT/HLA reference and its
# kallisto is 0.50.1, but arcasHLA 0.6.0 parses the kallisto 0.44 `pseudo` output (with
# 0.50.1 genotype dies with "pseudoalignments.tsv not found"). This script therefore
# bakes into the SIF: (1) the IMGT/HLA reference (dat/ref + dat/IMGTHLA/hla.dat), and
# (2) kallisto 0.44.0 (+ its libs) ahead of the stock binary.
#
# Needs: Docker with the pgx-suite:latest image (apptainer, run --privileged; rootless
# Docker: do not pass -u), a built arcasHLA clone (ARCAS_DIR; see
# docs/ToolsDocumentation.md 5b), and a kallisto-0.44.0 conda env (KALLISTO_ENV).
#   ARCAS_DIR=/data/alvin/ref/omnigen/tools/arcasHLA \
#   IMGTHLA_DIR=/path/to/IMGTHLA            (hla.dat + wmda/hla_nom*.txt) \
#   KALLISTO_ENV=/path/to/kallisto-env \
#   docker/build-arcashla-sif.sh [OUT_SIF] [WORKDIR]
set -euo pipefail
: "${ARCAS_DIR:?}" "${IMGTHLA_DIR:?}" "${KALLISTO_ENV:?}"
OUT="${1:-StellarPGx/containers/arcashla.sif}"; WORK="${2:-${TMPDIR:-/data/alvin/tmp}/arcas_build}"
OUT=$(realpath -m "$OUT"); mkdir -p "$WORK/tmp" "$(dirname "$OUT")"; WORK=$(realpath "$WORK")
APP="docker run --privileged --rm -v $WORK:/work -v $(dirname "$OUT"):/out -e APPTAINER_TMPDIR=/work/tmp pgx-suite:latest"
rm -rf "$WORK/sbx"
$APP apptainer build --sandbox /work/sbx docker://quay.io/biocontainers/arcas-hla:0.6.0--hdfd78af_1
S="$WORK/sbx/usr/local/share/arcas-hla-0.6.0-1/dat"
# reference (skip hla_partial*: 1.9 GB, only for --partial)
for f in allele_groups.json cDNA.json cDNA.single.json hla.convert.json hla.fasta hla.idx hla.p.json hla_transcripts.json; do
    cp "$ARCAS_DIR/dat/ref/$f" "$S/ref/$f"; done
cp "$ARCAS_DIR"/dat/info/* "$S/info/"
mkdir -p "$S/IMGTHLA/wmda"; cp "$IMGTHLA_DIR/hla.dat" "$S/IMGTHLA/"; cp "$IMGTHLA_DIR"/wmda/hla_nom*.txt "$S/IMGTHLA/wmda/"
# kallisto 0.44.0 (hla.idx must have been built with this version)
K="$WORK/sbx/opt/kallisto-0.44.0"; mkdir -p "$K/lib"; cp "$KALLISTO_ENV/bin/kallisto" "$K/"
for l in libhdf5.so.101 libz.so.1 libstdc++.so.6 libgcc_s.so.1; do cp -L "$KALLISTO_ENV/lib/$l" "$K/lib/"; done
printf '#!/bin/sh\nexport LD_LIBRARY_PATH=/opt/kallisto-0.44.0/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}\nexec /opt/kallisto-0.44.0/kallisto "$@"\n' > "$WORK/sbx/usr/local/bin/kallisto"
chmod +x "$WORK/sbx/usr/local/bin/kallisto"
rm -f "$OUT"; $APP apptainer build "/out/$(basename "$OUT")" /work/sbx
$APP apptainer exec "/out/$(basename "$OUT")" kallisto version   # expect 0.44.0
