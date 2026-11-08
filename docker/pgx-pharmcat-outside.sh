#!/usr/bin/env bash
# docker/pgx-pharmcat-outside.sh <OUT> — feed pgx-suite's HLA / CYP2D6 calls to PharmCAT.
# `pharmcat_pipeline` (3.2.0) has no outside-call option, so after the normal pipeline run
# this re-runs the core tool on the already-preprocessed VCF with `-po`
# (https://pharmcat.clinpgx.org/using/Outside-Call-Format/). Non-fatal: any failure leaves
# the stage-1 pgx.report.json untouched (saved copy: pgx.report.no_outside.json on success).
# Re-run alone for an existing sample (inside the image):
#   pgx-pharmcat-outside.sh /pgx/results/<sample>
# Env: PHARMCAT_CMD overrides the command prefix (default: apptainer exec in pharmcat.sif).
set -u
OUT=$1; PC=$OUT/pharmcat
SIF=${PHARMCAT_SIF:-/pgx/containers/pharmcat.sif}
RUN=${PHARMCAT_CMD:-"apptainer exec --cleanenv --bind $PC:$PC $SIF java -jar /pharmcat/pharmcat.jar"}
HERE=$(dirname "$(readlink -f "$0")")
n=$(python3 "$HERE/pgx_outside_calls.py" --out "$OUT" --update-provenance) || { echo "[outside-calls] builder failed; stage-1 report kept" >&2; exit 0; }
[ "$n" -gt 0 ] || { echo "[outside-calls] no calls to pass; stage-1 report kept" >&2; exit 0; }
[ -s "$PC/pgx.preprocessed.vcf.bgz" ] || { echo "[outside-calls] no preprocessed VCF; stage-1 report kept" >&2; exit 0; }
T=$PC/outside_run; rm -rf "$T"; mkdir -p "$T"
if $RUN -vcf "$PC/pgx.preprocessed.vcf.bgz" -po "$PC/outside_calls.pharmcat.tsv" \
        -reporterJson -o "$T" -bf pgx > "$PC/outside_run.log" 2>&1 && [ -s "$T/pgx.report.json" ]; then
  [ -e "$PC/pgx.report.no_outside.json" ] || cp "$PC/pgx.report.json" "$PC/pgx.report.no_outside.json"
  cp "$T/pgx.report.json" "$PC/pgx.report.json.tmp" && mv "$PC/pgx.report.json.tmp" "$PC/pgx.report.json"
  [ -s "$T/pgx.phenotype.json" ] && cp "$T/pgx.phenotype.json" "$PC/pgx.phenotype.json"
  echo "[outside-calls] PharmCAT re-run with $n outside call(s)" >&2
else
  echo "[outside-calls] PharmCAT -po run failed (see $PC/outside_run.log); stage-1 report kept" >&2
fi
exit 0
