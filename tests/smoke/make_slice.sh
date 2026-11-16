#!/usr/bin/env bash
# Rebuild tests/smoke/data/HG002.smoke.bam from a whole-genome GRCh38 BAM of GIAB HG002.
#
#   tests/smoke/make_slice.sh <HG002_GRCh38.bam>
#
# Keeps reads in regions.bed (4 genes + the copy-number control regions Aldy and
# StellarPGx need), drops @PG/@CO lines and BQSR tags (BI/BD), renames reads r1, r2, ...
# (mates keep a shared name) and records the aligner in a path-free @PG line.
# Edit the @PG line below if your BAM came from a different aligner.
set -euo pipefail
SRC="${1:?usage: make_slice.sh <HG002_GRCh38.bam>}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT="$HERE/data/HG002.smoke.bam"
HDR="$(mktemp -p "$HERE/data" hdr.XXXX.sam)"; trap 'rm -f "$HDR"' EXIT

samtools view -H --no-PG "$SRC" | grep -E '^@(HD|SQ)' > "$HDR"
printf '@RG\tID:rg\tSM:HG002\tPL:DNBSEQ\tLB:lb\n@PG\tID:bwa\tPN:bwa\tVN:BWA_2.0.2_acc_2.0\tCL:bwa mem -t 40 hg38.fa\n' >> "$HDR"

samtools view -b --no-PG -M -L "$HERE/regions.bed" "$SRC" \
  | samtools reheader --no-PG "$HDR" - \
  | samtools view -h --no-PG -x BI -x BD - \
  | awk 'BEGIN{OFS="\t"} /^@/{print;next} {if(!($1 in m)) m[$1]="r"(++n); $1=m[$1]; print}' \
  | samtools view -b --no-PG -o "$OUT" -
samtools index "$OUT"
ls -l "$OUT" "$OUT.bai"
