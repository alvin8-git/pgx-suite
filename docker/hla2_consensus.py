#!/usr/bin/env python3
"""hla2_consensus.py -- two-typer consensus for HLA class II (DRB1, DQA1, DQB1).

Typer 1: arcasHLA (kallisto pseudoalignment to IPD-IMGT/HLA cDNA, EM).
Typer 2: T1K     (read-to-allele alignment on genomic+cDNA, EM; MIT licence).
Both run on the same MHC reads and the same IPD-IMGT/HLA release, but use different
algorithms, so a disagreement is informative and an agreement is not a shared bug.

Per gene the output reports:
  arcashla / t1k      the raw calls (full names, as the typers wrote them)
  agreement           agree_2field | agree_1field | disagree | single_typer | not_typed
  consensus           per-allele names at the resolution BOTH typers support
                      (2-field "DRB1*10:01", 1-field group "DRB1*04", or "DRB1*?")
  consensus_resolution 2-field | 1-field | partial | none
  collapsed_rare      rare tie-break names (2nd field >= 100) replaced by their group
  support             arcasHLA reads/explained fraction, T1K abundance/quality
  flags               low_read_support, rare_name, zygosity notes ...
Pure stdlib; testable without a BAM (test_hla2_consensus.py).

CLI: hla2_consensus.py --sample S --arcas-json X.genotype.json
        [--arcas-verbose arcashla_verbose.log] [--t1k t1k_genotype.tsv] --out OUT.json
"""
from __future__ import annotations

import argparse
import itertools
import json
import re
import sys

CLASS_II_GENES = ("DRB1", "DQA1", "DQB1")
RARE_SECOND_FIELD = 100       # IMGT numbering: 2nd field >= 100 = rare / late-named
MIN_ARCAS_READS = 50          # reads aligned to the gene (30X WGS typically >100)
MIN_ARCAS_EXPLAINED = 70.0    # % of aligned reads explained by the chosen pair
MIN_T1K_QUALITY = 1           # T1K advises ignoring quality <= 0
MIN_T1K_ABUNDANCE = 5.0


def _strip(allele: str) -> str:
    a = allele.strip()
    return a[4:] if a.startswith("HLA-") else a


def fields(allele: str, n: int) -> str:
    """'DRB1*10:01:01' -> n=2 'DRB1*10:01', n=1 'DRB1*10'. Suffix letters (N,Q,..) dropped."""
    a = re.sub(r"[A-Za-z]$", "", _strip(allele))
    gene, _, rest = a.partition("*")
    parts = rest.split(":")[:n]
    return f"{gene}*" + ":".join(parts)


def is_rare_name(allele: str) -> bool:
    parts = _strip(allele).partition("*")[2].split(":")
    try:
        return len(parts) >= 2 and int(re.sub(r"\D", "", parts[1])) >= RARE_SECOND_FIELD
    except ValueError:
        return False


# ── input parsers ───────────────────────────────────────────────────────────
def parse_t1k_tsv(path: str) -> dict:
    """gene -> {"alleles": [opts1, opts2], "abundance": [..], "quality": [..], "secondary": [..]}.
    opts = list of tied allele names for that slot. A '.' second allele = homozygous."""
    out = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 5 or not f[0].startswith("HLA-"):
                continue
            gene = f[0][4:]
            slots, ab, q = [], [], []
            for i in (2, 5):
                if i + 2 < len(f):
                    if f[i] == ".":
                        continue
                    slots.append([_strip(a) for a in f[i].split(",") if a])
                    ab.append(float(f[i + 1]))
                    q.append(int(float(f[i + 2])))
            sec = []
            if len(f) > 8 and f[8]:
                for item in f[8].split("|"):
                    p = item.split(";")
                    if p and p[0]:
                        sec.append({"allele": _strip(p[0]),
                                    "abundance": float(p[1]) if len(p) > 1 else None,
                                    "quality": int(float(p[2])) if len(p) > 2 else None})
            if len(slots) == 1:      # homozygous: duplicate, halve nothing, keep evidence
                slots.append(list(slots[0])); ab.append(ab[0]); q.append(q[0])
            out[gene] = {"alleles": slots, "abundance": ab, "quality": q, "secondary": sec}
    return out


def parse_arcas_verbose(path: str) -> dict:
    """Per-gene evidence from `arcasHLA genotype --verbose` output (all typed genes)."""
    ev: dict = {}
    cur = None
    mode = None
    with open(path, errors="replace") as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            m = re.match(r"\[genotype\] Genotyping HLA-(\w+)", line)
            if m:
                cur = ev.setdefault(m.group(1), {"reads_aligned": None, "classes": None,
                                                 "abundance": [], "pairs": [], "zygosity": None})
                mode = None
                continue
            if cur is None:
                continue
            m = re.match(r"\[genotype\] (\d+) reads aligned to HLA-\w+ in (\d+) classes", line)
            if m:
                cur["reads_aligned"], cur["classes"] = int(m.group(1)), int(m.group(2))
            elif "Top alleles by abundance" in line:
                mode = "ab"
            elif "Pairs by % explained reads" in line:
                mode = "pairs"
            elif line.startswith("[genotype] Likely") or "zygous" in line and "[genotype]" in line:
                cur["zygosity"] = line.replace("[genotype] ", "").strip()
                mode = None
            elif line.startswith("[genotype]"):
                mode = None
            elif mode == "ab":
                m = re.match(r"\s+(\S+\*[\d:A-Za-z]+)\s+([\d.]+)%", line)
                if m:
                    cur["abundance"].append({"allele": m.group(1), "pct": float(m.group(2))})
            elif mode == "pairs":
                m = re.match(r"\s+(\S+),\s+(\S+)\s+([\d.]+)%", line)
                if m:
                    cur["pairs"].append({"alleles": [m.group(1), m.group(2)],
                                         "explained_pct": float(m.group(3))})
    return ev


# ── consensus ───────────────────────────────────────────────────────────────
def _slot_level(arcas_allele: str, t1k_opts: list) -> int:
    """2 if a T1K option matches at 2-field, 1 if at 1-field, else 0."""
    a2, a1 = fields(arcas_allele, 2), fields(arcas_allele, 1)
    if any(fields(o, 2) == a2 for o in t1k_opts):
        return 2
    if any(fields(o, 1) == a1 for o in t1k_opts):
        return 1
    return 0


def consensus_gene(gene: str, arcas: list | None, t1k: dict | None,
                   arcas_ev: dict | None = None) -> dict:
    flags: list = []
    res = {"gene": gene, "arcashla": arcas, "t1k": None, "agreement": "not_typed",
           "consensus": [], "consensus_resolution": "none", "collapsed_rare": [],
           "support": {}, "flags": flags}
    if t1k:
        res["t1k"] = ["|".join(o) if len(o) > 1 else o[0] for o in t1k["alleles"]]
        res["support"]["t1k"] = {"abundance": t1k["abundance"], "quality": t1k["quality"],
                                 "secondary": t1k.get("secondary", [])}
        for a, q, ab in zip(t1k["alleles"], t1k["quality"], t1k["abundance"]):
            if q < MIN_T1K_QUALITY or ab < MIN_T1K_ABUNDANCE:
                flags.append(f"low_read_support: T1K {a[0]} (quality {q}, abundance {ab:.1f})")
    if arcas_ev:
        pair = arcas_ev["pairs"][0]["explained_pct"] if arcas_ev.get("pairs") else None
        res["support"]["arcashla"] = {"reads_aligned": arcas_ev.get("reads_aligned"),
                                      "explained_pct": pair, "zygosity": arcas_ev.get("zygosity"),
                                      "abundance": arcas_ev.get("abundance", [])[:6],
                                      "pairs": arcas_ev.get("pairs", [])[:6]}
        ra = arcas_ev.get("reads_aligned")
        if ra is not None and ra < MIN_ARCAS_READS:
            flags.append(f"low_read_support: arcasHLA {ra} reads aligned (< {MIN_ARCAS_READS})")
        if pair is not None and pair < MIN_ARCAS_EXPLAINED:
            flags.append(f"low_read_support: arcasHLA pair explains {pair:.0f}% of reads")
    if not arcas and not t1k:
        flags.append("not typed by either typer")
        return res

    def _single(names, typer):
        out, collapsed = [], []
        for n in names:
            if is_rare_name(n):
                out.append(fields(n, 1)); collapsed.append(_strip(n))
            else:
                out.append(fields(n, 2))
        res["consensus"] = out
        res["collapsed_rare"] = collapsed
        res["consensus_resolution"] = "1-field" if collapsed else "2-field"
        res["agreement"] = "single_typer"
        flags.append(f"only {typer} typed this gene; no independent confirmation")

    if arcas and not t1k:
        _single(arcas, "arcasHLA"); return res
    if t1k and not arcas:
        _single([o[0] for o in t1k["alleles"]], "T1K"); return res

    t_slots = t1k["alleles"]
    best = None
    for perm in itertools.permutations(range(len(t_slots))):
        lv = [_slot_level(a, t_slots[j]) for a, j in zip(arcas, perm)]
        if best is None or sum(lv) > sum(best[0]):
            best = (lv, perm)
    levels, perm = best
    cons, collapsed = [], []
    for a, lv, j in zip(arcas, levels, perm):
        if lv == 2:
            cons.append(fields(a, 2))
            if is_rare_name(a):
                flags.append(f"rare_name: {fields(a, 2)} confirmed at 2-field by both typers")
        elif lv == 1:
            cons.append(fields(a, 1))
            for nm in [a] + t_slots[j]:
                if is_rare_name(nm) and _strip(nm) not in collapsed:
                    collapsed.append(_strip(nm))
        else:
            cons.append(f"{gene}*?")
    res["consensus"] = cons
    res["collapsed_rare"] = collapsed
    if all(lv == 2 for lv in levels):
        res["agreement"], res["consensus_resolution"] = "agree_2field", "2-field"
    elif all(lv >= 1 for lv in levels):
        res["agreement"], res["consensus_resolution"] = "agree_1field", "1-field"
    else:
        res["agreement"] = "disagree"
        res["consensus_resolution"] = "partial" if any(lv for lv in levels) else "none"
        flags.append("typers disagree at the allele-group level; unresolved allele(s) not reported")
    return res


def build(sample: str, arcas_json: dict | None, t1k: dict | None, arcas_ev: dict | None) -> dict:
    genes = {}
    for g in CLASS_II_GENES:
        a = arcas_json.get(g) if arcas_json else None
        genes[g] = consensus_gene(g, a, (t1k or {}).get(g), (arcas_ev or {}).get(g))
    return {"sample": sample, "schema": "hla2_consensus/1",
            "typers": {"arcashla": bool(arcas_json), "t1k": bool(t1k)},
            "genes": genes,
            "arcashla_verbose": arcas_ev or None}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--arcas-json")
    ap.add_argument("--arcas-verbose")
    ap.add_argument("--t1k")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    aj = json.load(open(a.arcas_json)) if a.arcas_json else None
    ev = parse_arcas_verbose(a.arcas_verbose) if a.arcas_verbose else None
    t1k = parse_t1k_tsv(a.t1k) if a.t1k else None
    out = build(a.sample, aj, t1k, ev)
    with open(a.out, "w") as fh:
        json.dump(out, fh, indent=2)
    return 0


if __name__ == "__main__":
    sys.exit(main())
