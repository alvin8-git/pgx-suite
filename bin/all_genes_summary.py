#!/usr/bin/env python3
"""Nextflow ALL_GENES_SUMMARY step: the Snakefile `all_genes_summary` run: block, moved verbatim.

usage: all_genes_summary.py OUT SAMPLE GENE [GENE ...]
"""
import json
import os
import sys

OUT, SAMPLE, TARGET = sys.argv[1], sys.argv[2], sys.argv[3:]
inputs = [os.path.join(OUT, "genes", g, f"{g}_{SAMPLE}_detail.json") for g in TARGET]
output = os.path.join(OUT, "all_genes_summary.tsv")

# HLA class II two-typer consensus (hla2/<S>.hla2_consensus.json), appended as
# NEW trailing columns (HLA2*, then DecidedBy); the first seven columns are unchanged.
h2 = {}   # "-" columns only when class II was not run
cj = os.path.join(OUT, "hla2", f"{SAMPLE}.hla2_consensus.json")
sp = os.path.join(OUT, "hla2", "hla2.status")
hla2_ran = os.path.exists(sp) and open(sp).read().strip() == "0"  # rc!=0: typer failed, reported as failed per gene
if os.path.exists(cj):
    try:
        h2 = json.load(open(cj))["genes"]
        assert isinstance(h2, dict)
    except Exception as exc:
        raise Exception(f"malformed HLA class II consensus {cj}: {exc!r}")
elif hla2_ran:
    raise Exception(f"HLA class II ran but {cj} is missing")
with open(output, "w") as out:
    out.write("Gene\tStatus\tConsensusDiplotype\tConsensusPhenotype\t"
              "nAgree\tnCalled\tMeanDepth\t"
              "HLA2Agreement\tHLA2Consensus\tHLA2Flags\tDecidedBy\n")
    for f in sorted(inputs):
        d = json.load(open(f))
        v = d.get("verdict", {})
        cov = d.get("coverage", {})
        hc = h2.get(d["gene"].replace("HLA-", ""), {}) if d["gene"].startswith("HLA-D") else {}
        out.write("\t".join(str(x) for x in [
            d["gene"], v.get("status", "?") + (" (caller killed: " + ", ".join(d["killed_tools"]) + ")" if d.get("killed_tools") else ""), v.get("consensus_diplotype", "-"),
            v.get("consensus_phenotype", "-"), v.get("n_agree", 0),
            v.get("n_called", 0), cov.get("mean_depth", "-"),
            hc.get("agreement", "-"),
            "/".join(hc["consensus"]) if hc.get("consensus") else "-",
            "; ".join(hc.get("flags", [])) or "-",
            "-" if v.get("authority") in (None, "", "Phenotype") else v["authority"],   # tool whose call is the verdict (an authority gene), else "-"
        ]) + "\n")
