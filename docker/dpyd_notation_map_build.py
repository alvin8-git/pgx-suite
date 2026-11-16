#!/usr/bin/env python3
"""Rebuild docker/dpyd_notation_map.json: caller-specific DPYD allele names -> CPIC DPYD allele names.

    python3 docker/dpyd_notation_map_build.py [--retrieved YYYY-MM-DD] [--out docker/dpyd_notation_map.json]

Why: Stargazer names DPYD alleles `*S<n>` and Aldy names a lone variant `*rs<NNN>`; CPIC's diplotype table
(cpic_diplotype.json) only knows CPIC's own allele names, so a reading in those notations has no CPIC row and the
split-phenotype rule (pgx-compare.py) could not score it. The map is derived, not written from memory:

* CPIC alleles: the CPIC API (`allele_definition`, `allele_location_value`, `sequence_location`, GRCh38 positions and
  dbSNP ids of each allele's variants), embedded in the output under `cpic_alleles`.
* Stargazer: the allele definition table shipped with the tool in this repo (stargazer-grc38-2.0.3/stargazer/
  star_table.tsv, column grc38_core, and snp_table.tsv for rs ids); its path, version and sha256 are in `source`.
* A `*S` code maps only when its GRCh38 variant set equals the variant set of exactly one CPIC allele: same
  position and ref>alt for SNVs; for indels (CPIC writes HGVS, Stargazer VCF style) the dbSNP id must match a
  single-variant CPIC allele and the Stargazer allele must be that one variant. Anything else is `unscorable`
  with the reason; nothing is guessed.
* CPIC star names: a bare `*4` / `*6` (Aldy, Stargazer other_names) maps to the CPIC allele whose own name carries
  that bracketed legacy name (`c.1601G>A (*4)`), when exactly one does.
* Aldy `*rs<NNN>` (one variant): maps to the CPIC allele defined by exactly that one variant, when CPIC has exactly
  one such allele for the rs id. Composite Aldy names (`*9A+rs...`, `*rsA+rsB`) are not in the map.

PGX_CPIC_SNAPSHOT=<dir> replays <dir>/DPYD_{allele_definition,allele_location_value,sequence_location}.json instead
of the network (same convention as cpic_cyp2d6_build.py).
"""
import argparse, csv, datetime, hashlib, json, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from cpic_cyp2d6_build import _get  # noqa: E402

SG_DIR = os.path.join(HERE, "..", "stargazer-grc38-2.0.3", "stargazer")
SG_VERSION_NOTE = "Stargazer GRCh38 build, directory stargazer-grc38-2.0.3 (stargazer/version.py __version__ 2.0.0)"
API = "https://api.cpicpgx.org/v1/"
URLS = {
    "allele_definition": API + "allele_definition?genesymbol=eq.DPYD&select=*&limit=1000",
    "sequence_location": API + "sequence_location?genesymbol=eq.DPYD&select=*&limit=1000",
    "allele_location_value": API + "allele_location_value?alleledefinitionid=in.({ids})&select=*&limit=5000",
}
OUT = os.path.join(HERE, "dpyd_notation_map.json")


def sha256(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def cpic_alleles():
    """{CPIC allele name: [{pos, change, rsid}, ...]} from the CPIC API (SNV change as plus-strand ref>alt, other
    variants as CPIC's HGVS text)."""
    defs = _get(URLS["allele_definition"], "DPYD_allele_definition")
    locs = {x["id"]: x for x in _get(URLS["sequence_location"], "DPYD_sequence_location")}
    ids = ",".join(str(d["id"]) for d in defs)
    alv = _get(URLS["allele_location_value"].format(ids=ids), "DPYD_allele_location_value")
    byid = {d["id"]: d["name"] for d in defs}
    out = {d["name"]: [] for d in defs}
    for a in alv:
        l = locs[a["locationid"]]
        m = re.match(r"^g\.(\d+)([ACGT])>([ACGT])$", l["chromosomelocation"])
        out[byid[a["alleledefinitionid"]]].append(
            {"pos": l["position"], "change": f"{m.group(2)}>{m.group(3)}" if m else l["chromosomelocation"],
             "rsid": l["dbsnpid"]})
    for v in out.values():
        v.sort(key=lambda x: (x["pos"], x["change"]))
    return out


def stargazer_codes():
    """{`*S<n>`: {core: [(pos, ref>alt)], rsids: [rs]}} for DPYD from the shipped Stargazer tables."""
    rs = {}
    with open(os.path.join(SG_DIR, "snp_table.tsv"), encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["gene"] == "dpyd":
                rs[(int(r["grc38_pos"]), f'{r["wt_allele"]}>{r["var_allele"]}')] = r["rs_id"]
    codes = {}
    with open(os.path.join(SG_DIR, "star_table.tsv"), encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["gene"] != "dpyd" or not r["name"].startswith("*S"):
                continue
            core = []
            for v in r["grc38_core"].split(","):
                p, c = v.split(":")
                core.append((int(p), c))
            codes[r["name"]] = {"core": core, "stargazer_name": r["other_names"],
                                "rsids": [rs.get(c) for c in core]}
    return codes


def build(retrieved):
    cp = cpic_alleles()
    single_by_rs = {}
    for n, v in cp.items():
        if len(v) == 1:
            single_by_rs.setdefault(v[0]["rsid"], []).append(n)
    by_set = {}
    for n, v in cp.items():
        if v:
            by_set.setdefault(tuple((x["pos"], x["change"]) for x in v), []).append(n)
    stargazer, unscorable = {}, {}
    for code, e in sorted(stargazer_codes().items(), key=lambda kv: int(kv[0][2:])):
        key = tuple(sorted(e["core"]))
        hit = by_set.get(key, [])
        method = "position+alt"
        if not hit and len(e["core"]) == 1 and e["rsids"][0]:
            # indel: CPIC writes HGVS, Stargazer VCF style -> the dbSNP id of the one variant decides
            hit = single_by_rs.get(e["rsids"][0], [])
            method = "rsid"
        if len(hit) == 1:
            stargazer[code.lower()] = {"cpic": hit[0], "method": method, "stargazer_core": [f"{p}:{c}" for p, c in e["core"]],
                                       "rsid": [r for r in e["rsids"] if r]}
        else:
            unscorable[code.lower()] = ("CPIC has several alleles for this variant set" if hit else
                                        "CPIC defines no allele with exactly this variant set")
    aldy = {}
    for rs, names in sorted(single_by_rs.items()):
        if len(names) == 1:
            aldy[f"*{rs}".lower()] = {"cpic": names[0], "rsid": rs}
    # CPIC's own allele names carry the legacy star name in brackets (`c.1601G>A (*4)`): Aldy and Stargazer's
    # other_names spell it bare (`*4`). Mapped only when exactly one CPIC allele carries that bracketed name.
    star = {}
    for n in cp:
        m = re.search(r"\((\*\d+[A-Z]?)\)$", n)
        if m:
            star.setdefault(m.group(1).lower(), []).append(n)
    star = {k: {"cpic": v[0]} for k, v in sorted(star.items()) if len(v) == 1}
    return {"retrieved": retrieved,
            "source": {"cpic": URLS,
                       "stargazer": {"tool": SG_VERSION_NOTE,
                                     "star_table": "stargazer-grc38-2.0.3/stargazer/star_table.tsv",
                                     "star_table_sha256": sha256(os.path.join(SG_DIR, "star_table.tsv")),
                                     "snp_table": "stargazer-grc38-2.0.3/stargazer/snp_table.tsv",
                                     "snp_table_sha256": sha256(os.path.join(SG_DIR, "snp_table.tsv"))},
                       "aldy": "Aldy names a lone variant `*rs<NNN>`; mapped by dbSNP id to the CPIC allele defined by "
                               "exactly that variant"},
            "DPYD": {"stargazer": stargazer, "aldy": aldy, "cpic_star_names": star, "unscorable_stargazer": unscorable},
            "cpic_alleles": cp}


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--retrieved", default=datetime.date.today().isoformat())
    ap.add_argument("--out", default=OUT)
    a = ap.parse_args()
    data = build(a.retrieved)
    with open(a.out, "w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=1, sort_keys=True)
        fh.write("\n")
    d = data["DPYD"]
    print("wrote", a.out, "stargazer", len(d["stargazer"]), "unscorable", len(d["unscorable_stargazer"]), "star names", len(d["cpic_star_names"]), "aldy", len(d["aldy"]))
