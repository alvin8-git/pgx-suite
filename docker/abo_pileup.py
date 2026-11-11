#!/usr/bin/env python3
"""abo_pileup.py -- per-site allele evidence for the ABO blood-group caller.

Reads the BAM/CRAM directly (samtools mpileup) at the ABO defining sites listed in
abo_allele_defs.json, so "covered, matches the reference" is distinguishable from
"no data" (a variants-only VCF cannot tell the two apart).

Usage: abo_pileup.py <BAM|CRAM> <REF.fa> <OUT.json> [--defs abo_allele_defs.json]
Output JSON: {"sites": {rsid: {chrom,pos,depth,ref,ref_count,alt_counts,ins_count,
              del_count,other_count}}, "min_mapq":..., "min_baseq":...}
"""
import argparse
import json
import os
import re
import subprocess
import sys

MIN_MAPQ = 20
MIN_BASEQ = 20

_TOKEN = re.compile(r"\^.|\$|[.,ACGTNacgtn*#<>]|[+-]\d+")


def parse_pileup_bases(bases: str, ref_base: str) -> dict:
    """Parse one mpileup base column into allele counts.

    ins_seqs / del_seqs list the inserted / deleted sequence that FOLLOWS a read's
    base at this position (mpileup attaches an indel to the preceding position)."""
    out = {"ref": 0, "alt": {}, "ins": [], "dele": [], "star": 0}
    i = 0
    n = len(bases)
    while i < n:
        c = bases[i]
        if c == "^":
            i += 2
            continue
        if c == "$":
            i += 1
            continue
        if c in ".,":
            out["ref"] += 1
            i += 1
        elif c in "ACGTNacgtn":
            b = c.upper()
            out["alt"][b] = out["alt"].get(b, 0) + 1
            i += 1
        elif c in "*#":
            out["star"] += 1
            i += 1
        elif c in "+-":
            m = re.match(r"[+-](\d+)", bases[i:])
            ln = int(m.group(1))
            start = i + m.end()
            seq = bases[start:start + ln].upper()
            (out["ins"] if c == "+" else out["dele"]).append(seq)
            i = start + ln
        else:  # '<' '>' reference skip
            i += 1
    return out


def pileup_site(bam: str, ref: str, chrom: str, pos: int, samtools: str = "samtools") -> dict:
    cmd = [samtools, "mpileup", "-f", ref, "-r", f"{chrom}:{pos}-{pos}",
           "-q", str(MIN_MAPQ), "-Q", str(MIN_BASEQ), "-d", "2000", "-A", bam]
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    if p.returncode != 0:
        raise RuntimeError(f"samtools mpileup failed at {chrom}:{pos}: {p.stderr.strip()[:200]}")
    line = next((l for l in p.stdout.splitlines() if l.strip()), None)
    if line is None:
        return {"chrom": chrom, "pos": pos, "depth": 0, "ref": None, "ref_count": 0,
                "alt_counts": {}, "ins_count": 0, "ins_seqs": {}, "del_count": 0}
    f = line.split("\t")
    parsed = parse_pileup_bases(f[4], f[2].upper())
    ins_seqs: dict = {}
    for s in parsed["ins"]:
        ins_seqs[s] = ins_seqs.get(s, 0) + 1
    return {"chrom": chrom, "pos": pos, "depth": int(f[3]), "ref": f[2].upper(),
            "ref_count": parsed["ref"], "alt_counts": parsed["alt"],
            "ins_count": len(parsed["ins"]), "ins_seqs": ins_seqs,
            "del_count": parsed["star"] + len(parsed["dele"])}


def run(bam: str, ref: str, defs_path: str) -> dict:
    with open(defs_path, encoding="utf-8") as fh:
        defs = json.load(fh)
    sites = {}
    for v in defs["variants"]:
        sites[v["rsid"]] = pileup_site(bam, ref, v["chrom"], v["pos_grch38"])
    return {"sites": sites, "min_mapq": MIN_MAPQ, "min_baseq": MIN_BASEQ,
            "bam": os.path.basename(bam)}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("bam")
    ap.add_argument("ref")
    ap.add_argument("out")
    ap.add_argument("--defs", default=os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                   "abo_allele_defs.json"))
    a = ap.parse_args()
    res = run(a.bam, a.ref, a.defs)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with open(a.out, "w") as fh:
        json.dump(res, fh, indent=2)
    return 0


if __name__ == "__main__":
    sys.exit(main())
