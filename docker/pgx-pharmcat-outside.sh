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
python3 - "$OUT/provenance.json" <<'PY' || true
import json, os, sys
p = sys.argv[1]   # a stale "outside call not run (killed)" note from an earlier run must not outlive a good run
if os.path.exists(p):
    d = json.load(open(p))
    if d.pop("pharmcat_outside_run", None) is not None:
        json.dump(d, open(p + ".tmp", "w"), indent=2); os.replace(p + ".tmp", p)
PY
n=$(python3 "$HERE/pgx_outside_calls.py" --out "$OUT" --update-provenance) || { echo "[outside-calls] builder failed; stage-1 report kept" >&2; exit 0; }
[ "$n" -gt 0 ] || { echo "[outside-calls] no calls to pass; stage-1 report kept" >&2; exit 0; }
[ -s "$PC/pgx.preprocessed.vcf.bgz" ] || { echo "[outside-calls] no preprocessed VCF; stage-1 report kept" >&2; exit 0; }
T=$PC/outside_run; rm -rf "$T" "$PC/outside.killed"; mkdir -p "$T"
rc=0
$RUN -vcf "$PC/pgx.preprocessed.vcf.bgz" -po "$PC/outside_calls.pharmcat.tsv" \
        -reporterJson -o "$T" -bf pgx > "$PC/outside_run.log" 2>&1 || rc=$?
if [ "$rc" -eq 0 ] && [ -s "$T/pgx.report.json" ]; then
  [ -e "$PC/pgx.report.no_outside.json" ] || cp "$PC/pgx.report.json" "$PC/pgx.report.no_outside.json"
  cp "$T/pgx.report.json" "$PC/pgx.report.json.tmp" && mv "$PC/pgx.report.json.tmp" "$PC/pgx.report.json"
  [ -s "$T/pgx.phenotype.json" ] && cp "$T/pgx.phenotype.json" "$PC/pgx.phenotype.json"
  echo "[outside-calls] PharmCAT re-run with $n outside call(s)" >&2
elif . "$HERE/kill_codes.sh" && is_kill_code "$rc"; then
  # Killed (137 = OOM or SIGKILL): still non-fatal, but never silent. The marker is read by the report (provenance.json).
  echo "outside call not run (killed, exit $rc)" > "$PC/outside.killed"
  python3 - "$OUT/provenance.json" "$rc" <<'PY' || true
import json, os, sys
p, rc = sys.argv[1], sys.argv[2]
if os.path.exists(p):
    d = json.load(open(p)); d["pharmcat_outside_run"] = f"outside call not run (killed, exit {rc})"
    json.dump(d, open(p + ".tmp", "w"), indent=2); os.replace(p + ".tmp", p)
PY
  echo "[outside-calls] outside call not run (killed, exit $rc); stage-1 report kept" >&2
  exit "$rc"   # the caller (outside.status; the Nextflow step retries 134/137) records it; still non-fatal for the run
else
  echo "[outside-calls] PharmCAT -po run failed (see $PC/outside_run.log); stage-1 report kept" >&2
fi
exit 0
