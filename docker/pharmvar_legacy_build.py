#!/usr/bin/env python3
"""Build pharmvar_legacy_snapshot.json and pharmvar_legacy_suballeles.json from rendered PharmVar gene pages.

PharmVar's Legacy Label column is only in the rendered page (https://www.pharmvar.org/gene/<GENE>), not in the
download files. Procedure: render each page with a browser, save the page text to <dir>/<GENE>.txt (one line), run
  python3 pharmvar_legacy_build.py <dir> <checked-date YYYY-MM-DD>
Page text lists, per sub-allele, `<GENE>*<core>.<NNN>[<GENE>*<legacy label>]<PVid>`. The snapshot stores every
sub-allele that has a legacy label. The fold table is derived from it by one rule, nothing else:
  fold          a lettered legacy label `*<n><L>` whose number n equals the core allele it sits under (a PharmVar
                sub-allele of that numbered allele): reconcile folds it to `*<n>`;
  redesignated  a lettered label whose number differs from its core allele (CYP2C19 `*1A` is now `*38.001`; the old
                meaning was the reference allele): never folded, the new star is recorded for the reader.
Genes PharmVar does not maintain (TPMT, UGT1A1: page says so) have no entry, so their letters are never folded."""
import json, re, sys, os

GENES = "CYP2C19 CYP2C9 CYP2B6 CYP2D6 CYP3A4 CYP3A5 CYP4F2 NUDT15 SLCO1B1 NAT2".split()
NOT_MAINTAINED = ["TPMT", "UGT1A1"]       # page text: "This gene is not maintained by the PharmVar Database"
NOT_BY_LEGACY_LABEL = ["DPYD"]            # PharmVar lists DPYD by rsID / HGVS, no star legacy labels
LETTER = re.compile(r"^(\d+)([A-Za-z])$")


def parse(gene, text):
    v = re.search(r"Version ([\d.]+) Last Updated (\w+ \d+ \d+)", text)
    rows = [[f"{gene}*{m.group(1)}.{m.group(2)}", f"*{m.group(3)}", m.group(4)]
            for m in re.finditer(re.escape(gene) + r"\*([^\s.]+?)\.(\d{3})(?:" + re.escape(gene) + r"\*(\S+?))?(PV\d+)", text)
            if m.group(3)]
    return {"pharmvar_version": v.group(1), "page_last_updated": v.group(2), "sub_alleles_with_legacy_label": rows}


def derive(snap, checked):
    out = {"source": "docker/pharmvar_legacy_snapshot.json (PharmVar gene pages, Legacy Label column), derived by pharmvar_legacy_build.py",
           "checked": checked, "fold": {}, "redesignated": {}, "not_folded_genes": NOT_MAINTAINED + NOT_BY_LEGACY_LABEL}
    for g, d in snap["genes"].items():
        for sub, legacy, pv in d["sub_alleles_with_legacy_label"]:
            m = LETTER.match(legacy[1:])
            if not m:
                continue
            core = re.match(rf"{g}\*([^.]+)\.", sub).group(1)
            ent = {"pharmvar_suballele": sub, "pharmvar_id": pv, "core_allele": f"*{core}", "pharmvar_version": d["pharmvar_version"]}
            (out["fold"] if m.group(1) == core else out["redesignated"]).setdefault(g, {})[legacy] = ent
    return out


if __name__ == "__main__":
    d, checked = sys.argv[1], sys.argv[2]
    snap = {"source": "PharmVar gene pages https://www.pharmvar.org/gene/<GENE>, rendered with a browser", "checked": checked, "genes": {}}
    for g in GENES:
        snap["genes"][g] = parse(g, open(os.path.join(d, g + ".txt"), encoding="utf-8").read())
    here = os.path.dirname(os.path.abspath(__file__))
    json.dump(snap, open(os.path.join(here, "pharmvar_legacy_snapshot.json"), "w"), indent=1)
    json.dump(derive(snap, checked), open(os.path.join(here, "pharmvar_legacy_suballeles.json"), "w"), indent=1)
