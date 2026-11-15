#!/usr/bin/env python3
"""Compare two pgx-suite result trees (e.g. Snakemake vs Nextflow) file by file.

    tests/compare_runs.py A_DIR B_DIR      # each is results/<SAMPLE>/
    tests/compare_runs.py --selftest

Rules (docs/nextflow-port-plan.md section 5.2): a file on one side only is a failure unless
allow-listed; TSV/status/text byte-identical; JSON field-identical minus run-time keys;
VCF bodies identical (header command/date lines ignored); PyPGx zips by member content
(the random tmpXXXX/ prefix and zip mtimes ignored); HTML minus the generation timestamp.
Exit 0 = equivalent, 1 = differences (listed), 2 = usage error.
"""
import gzip
import json
import os
import re
import sys
import tempfile
import zipfile

# Present on one side only, by design.
ONE_SIDED = [r"^pipeline_info/", r"(^|/)\.snakemake_timestamp$", r"^smoke\.log$"]
# Presence only: binary indexes derived from a file that is compared; PDFs embed dates.
PRESENCE_ONLY = [r"\.(tbi|csi|bai)$", r"\.pdf$"]
# JSON keys ignored everywhere (run time, not results): our provenance stamp; PharmCAT's
# report/matcher timestamps.
JSON_IGNORE_KEYS = {"generated_utc", "timestamp"}
# Text log lines that carry wall-clock times or durations (Stargazer, Beagle).
LOG_LINE_IGNORE = {
    r"/stargazer(_gdf)?/[^/]+\.log$": re.compile(
        rb"^(\[INFO\] (Start|Elapsed) time:|(Start|End) time:|.*phasing time:|Total time:)"),
}
# provenance.json: engine and image tag differ by design; reference.path is the same string today.
PROV_IGNORE = {"pipeline", ("tools", "snakemake"), ("tools", "nextflow")}
VCF_HEADER_IGNORE = re.compile(rb"^##(bcftools_\w+Command|fileDate)=")
HTML_IGNORE = [re.compile(rb"\d{4}-\d{2}-\d{2}[ T]\d{2}:\d{2}(:\d{2})?(Z|\s*UTC)?")]
TMP_PREFIX = re.compile(r"^tmp[^/]*/")

# Tool nondeterminism seen between two Snakemake runs of the same input (plan 5.3). A file that
# differs only in these fields is reported as "allowed", never silently accepted. The calls
# (diplotype, phenotype, status) are always compared strictly.
#  - PyPGx: the per-allele evidence column VariantData (and what pgx-compare.py derives from it)
#    changes run to run; Genotype/Phenotype do not.
#  - Stargazer: candidate-allele lists in report.tsv come out in a different order.
#  - pgx-compare.py adds a PharmCAT cross-check to non-PharmCAT genes only if PharmCAT finished
#    first under Snakemake (undeclared dependency, a race); the Nextflow port never adds it.
#  - The HTML report is rendered from the files above; it is "allowed" only when one of them was.
NONDET_DETAIL = [("tools", "PyPGx", "supporting_variants"), ("tools", "PyPGx", "mean_af"),
                 ("cross_check", "PharmCAT")]


def files(root):
    out = set()
    for d, _, fs in os.walk(root):
        for f in fs:
            out.add(os.path.relpath(os.path.join(d, f), root))
    return out


def matches(rel, pats):
    return any(re.search(p, rel) for p in pats)


def strip_json(o, path=()):
    if isinstance(o, dict):
        return {k: strip_json(v, path + (k,)) for k, v in o.items()
                if k not in JSON_IGNORE_KEYS}
    if isinstance(o, list):
        return [strip_json(v, path) for v in o]
    return o


def load_json(p, rel):
    d = strip_json(json.load(open(p)))
    if rel == "provenance.json":
        for k in PROV_IGNORE:
            if isinstance(k, tuple):
                d.get(k[0], {}).pop(k[1], None)
            else:
                d.pop(k, None)
    if rel == "bam_stats.json":
        d.pop("bam", None)
    return d


def vcf_lines(p):
    with gzip.open(p, "rb") as fh:
        return [ln for ln in fh if not VCF_HEADER_IGNORE.match(ln)]


def zip_members(p):
    with zipfile.ZipFile(p) as z:
        return {TMP_PREFIX.sub("", n): z.read(n) for n in z.namelist()}


def html_text(p):
    b = open(p, "rb").read()
    for rx in HTML_IGNORE:
        b = rx.sub(b"<TIME>", b)
    return b


def relaxed(a, b, rel):
    """True if a and b are equal once the tool-nondeterminism rules above are applied."""
    if rel.endswith("_detail.json"):
        da, db = load_json(a, rel), load_json(b, rel)
        for d in (da, db):
            for path in NONDET_DETAIL:
                node = d
                for k in path[:-1]:
                    node = node.get(k, {}) if isinstance(node, dict) else {}
                if isinstance(node, dict):
                    node.pop(path[-1], None)
            if d.get("cross_check") == {}:
                d.pop("cross_check")
        return da == db
    if rel.endswith(".zip") and "/pypgx/" in "/" + rel:
        def rows(p):
            out = {}
            for n, data in zip_members(p).items():
                lines = data.decode(errors="replace").splitlines()
                if n.endswith(".tsv") and lines:
                    hdr = lines[0].split("\t")
                    drop = hdr.index("VariantData") if "VariantData" in hdr else None
                    lines = ["\t".join(c for i, c in enumerate(ln.split("\t")) if i != drop) for ln in lines]
                out[n] = lines
            return out
        return rows(a) == rows(b)
    if re.search(r"/stargazer/report\.tsv$", "/" + rel):
        norm = lambda p: [[",".join(sorted(c.split(","))) for c in ln.rstrip("\n").split("\t")]
                          for ln in open(p)]
        return norm(a) == norm(b)
    return False


def same(a, b, rel):
    if matches(rel, PRESENCE_ONLY):
        return True
    if rel.endswith(".json"):
        try:
            return load_json(a, rel) == load_json(b, rel)
        except ValueError:
            pass   # not valid JSON: fall through to bytes
    if rel.endswith((".vcf.gz", ".vcf.bgz")):
        return vcf_lines(a) == vcf_lines(b)
    if rel.endswith(".zip"):
        return zip_members(a) == zip_members(b)
    if rel.endswith(".html"):
        return html_text(a) == html_text(b)
    for pat, rx in LOG_LINE_IGNORE.items():
        if re.search(pat, "/" + rel):
            keep = lambda p: [ln for ln in open(p, "rb") if not rx.match(ln)]
            return keep(a) == keep(b)
    return open(a, "rb").read() == open(b, "rb").read()


def compare(a_root, b_root):
    fa, fb = files(a_root), files(b_root)
    diffs, allowed = [], []
    for rel in sorted(fa ^ fb):
        if not matches(rel, ONE_SIDED):
            diffs.append(f"only in {'A' if rel in fa else 'B'}: {rel}")
    html = []
    for rel in sorted(fa & fb):
        if matches(rel, ONE_SIDED):   # run records (trace, timeline): not results
            continue
        a, b = os.path.join(a_root, rel), os.path.join(b_root, rel)
        if same(a, b, rel):
            continue
        if rel.endswith(".html"):
            html.append(rel)
        elif relaxed(a, b, rel):
            allowed.append(rel)
        else:
            diffs.append(f"differs: {rel}")
    for rel in html:   # rendered from compared files: allowed only if an upstream file was
        (allowed if allowed and not diffs else diffs).append(rel if allowed and not diffs else f"differs: {rel}")
    return diffs, len(fa & fb), allowed


def selftest():
    with tempfile.TemporaryDirectory() as t:
        a, b = os.path.join(t, "a"), os.path.join(t, "b")
        for root, tag, gen in ((a, "x", "2026-01-01"), (b, "x", "2026-02-02")):
            os.makedirs(os.path.join(root, "genes/G/pypgx"))
            open(os.path.join(root, "all_genes_summary.tsv"), "w").write("G\tconcordant\n")
            json.dump({"generated_utc": gen, "v": 1}, open(os.path.join(root, "provenance.json"), "w"))
            with zipfile.ZipFile(os.path.join(root, "genes/G/pypgx/r.zip"), "w") as z:
                z.writestr(f"tmp{gen}/data.tsv", "same")
        assert compare(a, b)[0] == [], compare(a, b)[0]
        open(os.path.join(b, "all_genes_summary.tsv"), "w").write("G\tdiscordant\n")
        open(os.path.join(b, "extra.txt"), "w").write("x")
        d = compare(a, b)[0]
        assert d == ["only in B: extra.txt", "differs: all_genes_summary.tsv"], d
    print("selftest ok")


if __name__ == "__main__":
    if sys.argv[1:] == ["--selftest"]:
        selftest()
        sys.exit(0)
    if len(sys.argv) != 3:
        print(__doc__, file=sys.stderr)
        sys.exit(2)
    diffs, n, allowed = compare(sys.argv[1], sys.argv[2])
    for d in diffs:
        print(d)
    for d in allowed:
        print(f"allowed (tool nondeterminism): {d}")
    print(f"{n} common files compared, {len(diffs)} differences, {len(allowed)} allowed")
    sys.exit(1 if diffs else 0)
