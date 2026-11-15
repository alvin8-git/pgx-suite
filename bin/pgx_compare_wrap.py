#!/usr/bin/env python3
"""Nextflow COMPARE step: the Snakefile `compare` rule's run: block, moved verbatim.

usage: pgx_compare_wrap.py GENE REGION IS_ALT(0|1) SAMPLE MINDEP BAM REF OUT [SEX]
"""
import json
import os
import subprocess
import sys

gene, region, is_alt, SAMPLE, MINDEP, BAM, REF, OUT = sys.argv[1:9]
sex = sys.argv[9] if len(sys.argv) > 9 else ""
outdir = os.path.join(OUT, "genes", gene)
logdir = os.path.join(outdir, "logs")

names = {"pypgx": "PyPGx", "stargazer": "Stargazer", "aldy": "Aldy",
         "stellarpgx": "StellarPGx", "optitype": "OptiType", "mutserve": "mutserve",
         "cyrius": "Cyrius"}
# OptiType is sample-level (one run for A/B/C): its status lives at
# OUT/optitype/optitype.status, not in the per-gene logs dir.
sample_level_status = {"optitype": os.path.join(OUT, "optitype", "optitype.status")}
failed = []
for key, name in names.items():
    sf = sample_level_status.get(key, os.path.join(logdir, f"{key}.status"))
    if os.path.exists(sf) and open(sf).read().strip() not in ("0", ""):
        failed.append(name)
args = (f"--gene {gene} --sample {SAMPLE} "
        f"--output-dir {outdir} --min-depth {MINDEP}")
if is_alt != "1":   # alt-contig genes (GSTT1) have no measurable region depth
    depth = subprocess.run(
        f"samtools coverage -r {region} --reference {REF} {BAM} "
        f"| awk 'NR==2 {{print $7}}'",
        shell=True, capture_output=True, text=True).stdout.strip()
    if depth:
        args += f" --region-depth {depth}"
if failed:
    args += " --failed-tools " + ",".join(failed)
alt_aware = "yes"
try:
    alt_aware = json.load(open(os.path.join(OUT, "alt_aware.json"))).get("alt_aware", "yes")
except Exception:
    pass
args += f" --alt-aware {alt_aware}"
args += f" --bam {BAM} --ref {REF}"
if sex:
    args += f" --sex {sex}"
sys.exit(subprocess.run(f"python3 /opt/pgx/pgx-compare.py {args}", shell=True).returncode)
