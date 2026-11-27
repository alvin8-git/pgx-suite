#!/usr/bin/env python3
"""A caller whose sentinel is a kill code (134/137/140/143) must show as "caller killed", never as an
ordinary missing caller, in comparison.tsv, detail.json and the all-genes summary. An ordinary failure
(rc 1) must not. No framework. Run: python3 docker/test_killed_caller.py"""
import json, os, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
SUMMARY = os.path.join(HERE, "..", "bin", "all_genes_summary.py")


def run(code, gene="CYP2C19", sentinel="genes/CYP2C19/logs/stargazer.status"):
    out = tempfile.mkdtemp()
    gdir = os.path.join(out, "genes", gene); os.makedirs(os.path.join(gdir, "logs"))
    sf = os.path.join(out, sentinel); os.makedirs(os.path.dirname(sf), exist_ok=True)
    open(sf, "w").write(code)
    subprocess.run([sys.executable, os.path.join(HERE, "pgx-compare.py"), "--gene", gene, "--sample", "S1",
                    "--output-dir", gdir, "--min-depth", "10"], check=True, capture_output=True)
    subprocess.run([sys.executable, SUMMARY, out, "S1", gene], check=True, capture_output=True)
    d = json.load(open(os.path.join(gdir, f"{gene}_S1_detail.json")))
    tsv = open(os.path.join(gdir, f"{gene}_S1_comparison.tsv")).read()
    summ = open(os.path.join(out, "all_genes_summary.tsv")).read().splitlines()[1]
    return d, tsv, summ


def main():
    for code in ("134", "137", "140", "143"):
        d, tsv, summ = run(code)
        assert d.get("killed_tools") == {"Stargazer": code}, (code, d.get("killed_tools"))
        assert f"caller killed (OOM/signal, exit {code})" in tsv, code
        assert "(caller killed: Stargazer (" in summ, summ
    # sample-level callers: their sentinel lives under OUT, not in the gene's logs dir
    for tool, gene, sf in (("OptiType", "HLA-A", "optitype/optitype.status"),
                           ("arcasHLA", "HLA-DQB1", "hla2/hla2.status"),
                           ("PharmCAT", "UGT1A1", "pharmcat/pharmcat.status")):
        d, tsv, summ = run("137", gene, sf)
        assert d.get("killed_tools") == {tool: "137"}, (tool, d.get("killed_tools"))
        assert "caller killed (OOM/signal, exit 137)" in tsv and f"caller killed: {tool}" in summ, (tool, summ)
        d, tsv, summ = run("1", gene, sf)
        assert "killed_tools" not in d, tool
    # the denominator includes the killed caller, in the verdict, the summary and the same text everywhere
    d, tsv, summ = run("137")
    note = d["verdict"]["killed_note"]
    n = d["verdict"]["n_called"]
    assert d["verdict"]["n_killed"] == 1 and f"caller killed: Stargazer ({d['verdict']['n_agree']} of {n + 1} callers agree, 1 killed)" == note, note
    assert note in summ, (note, summ)
    d, tsv, summ = run("1")
    assert "killed_tools" not in d and "caller killed" not in tsv and "caller killed" not in summ
    d, tsv, summ = run("0")
    assert "killed_tools" not in d and "caller killed" not in summ
    print("test_killed_caller: PASS")


if __name__ == "__main__":
    main()
