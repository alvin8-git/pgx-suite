#!/usr/bin/env bash
# Build the T1K IPD-IMGT/HLA index (REQUIRED for the second HLA class II typer).
#
# Offline: uses the SAME hla.dat the arcasHLA SIF was built from, so both typers see one
# IPD-IMGT/HLA release (3.64.0 at the time of writing). Needs T1K's t1k-build.pl (in the
# pgx-suite image at /opt/t1k) and perl. Output dir holds hlaidx_dna_seq.fa (+ rna);
# mount it at /pgx/containers/t1k_hlaidx (or pass `pgx-hla2.sh --t1k-idx DIR`).
#   HLA_DAT=<path>/IMGTHLA/hla.dat docker/build-t1k-index.sh <OUT_DIR> [IMAGE]
set -euo pipefail
: "${HLA_DAT:?set HLA_DAT to the IPD-IMGT/HLA hla.dat used for arcasHLA}"
OUT=$(realpath -m "${1:?usage: build-t1k-index.sh OUT_DIR [IMAGE]}"); IMG="${2:-pgx-suite:latest}"
mkdir -p "$OUT"
docker run --rm --entrypoint bash -v "$OUT:/out" -v "$(realpath "$HLA_DAT"):/in/hla.dat:ro" "$IMG" \
    -c 'perl /opt/t1k/t1k-build.pl -o /out -d /in/hla.dat && ls -l /out/hlaidx_dna_seq.fa'
