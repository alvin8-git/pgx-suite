#!/usr/bin/env python3
"""Nextflow SANITY_GATE step: the Snakefile `sanity_gate` run: block, moved verbatim.

usage: sanity_gate.py SUMMARY_TSV OUTPUT
"""
import sys

summary, output = sys.argv[1], sys.argv[2]
rows = [ln.split("\t") for ln in open(summary).read().splitlines()[1:] if ln.strip()]
total = len(rows)
bad = sum(1 for r in rows if len(r) > 1 and r[1] in ("no_call", "no_data"))
frac = (bad / total) if total else 1.0
if total == 0 or frac > 0.6:
    sys.stderr.write(
        f"SANITY FAIL: {bad}/{total} genes are no_call/no_data ({frac:.0%}). "
        "Likely a misconfigured run (wrong reference, unreadable BAM, or wholesale "
        "low coverage). The HTML report was still written for inspection.\n")
    sys.exit(1)
with open(output, "w") as fh:
    fh.write(f"ok: {total - bad}/{total} genes called\n")
