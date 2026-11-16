#!/usr/bin/env bash
# Nextflow task helper: merge the upstream result trees staged as up*/res into ./res
# (later inputs overwrite earlier ones), then point /pgx/results at ./res so every
# pgx script sees the same paths it sees under Snakemake (OUT=/pgx/results).
# ponytail: copies the trees (cp -rL, as the port plan 2.3 asks); fine at per-gene sizes.
set -euo pipefail
mkdir -p res
for d in $(ls -d up*/res 2>/dev/null | sort -V); do
    cp -rL --preserve=mode,timestamps "$d/." res/
done
rm -rf /pgx/results
ln -s "$PWD/res" /pgx/results
