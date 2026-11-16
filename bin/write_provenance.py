#!/usr/bin/env python3
"""Nextflow PROVENANCE step: the Snakefile `provenance` run: block, moved verbatim.
Only change: the engine key is `nextflow` (passed in) instead of `snakemake`.

usage: write_provenance.py OUT SAMPLE REF N_GENES NEXTFLOW_VERSION
"""
import datetime
import json
import os
import subprocess
import sys

OUT, SAMPLE, REF, n_genes, nf_version = sys.argv[1:6]


def first_line(cmd):
    # encoding/errors guard against tools that emit non-UTF-8 bytes;
    # skip warning/traceback/source-path lines so we land on the real version.
    try:
        r = subprocess.run(cmd, shell=True, capture_output=True,
                           encoding="utf-8", errors="replace", timeout=60)
        out = (r.stdout or r.stderr or "").strip()
        for ln in out.splitlines():
            s = ln.strip()
            low = s.lower()
            if not s or "warning" in low or "traceback" in low \
               or s.startswith('File "') or '.py:' in s:
                continue
            return s
        return out.splitlines()[0] if out else "unknown"
    except Exception as e:
        return f"unavailable ({e.__class__.__name__})"


ref_fai = REF + ".fai"
prov = {
    "pipeline": (open("/opt/pgx/VERSION").read().strip()
                 if os.path.exists("/opt/pgx/VERSION") else "unknown"),
    "generated_utc": datetime.datetime.utcnow().isoformat(timespec="seconds") + "Z",
    "sample": SAMPLE,
    "reference": {
        "path": REF,
        "bytes": os.path.getsize(REF) if os.path.exists(REF) else None,
        "n_contigs": (sum(1 for _ in open(ref_fai)) if os.path.exists(ref_fai) else None),
    },
    "genes_assessed": int(n_genes),
    "tools": {
        "pypgx": first_line("python3 -c \"import importlib.metadata as m; print(m.version('pypgx'))\""),
        "aldy": first_line("python3 -c \"import importlib.metadata as m; print(m.version('aldy'))\""),
        "stargazer": first_line("PYTHONWARNINGS=ignore stargazer --version 2>&1"),
        "stellarpgx": "1.2.7 (image)",
        "optitype": "1.3.5 (image)",
        "mutserve": "2.0.0 (image)",
        "samtools": first_line("samtools --version 2>&1"),
        "bcftools": first_line("bcftools --version 2>&1"),
        "nextflow": nf_version,
    },
}
with open(os.path.join(OUT, "provenance.json"), "w") as fh:
    json.dump(prov, fh, indent=2)
