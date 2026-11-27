#!/usr/bin/env python3
"""A caller whose sentinel is a kill code (134/137/140/143) must show as "caller killed", never as an
ordinary missing caller, in comparison.tsv, detail.json and the all-genes summary. An ordinary failure
(rc 1) must not. No framework. Run: python3 docker/test_killed_caller.py"""
import json, os, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
SUMMARY = os.path.join(HERE, "..", "bin", "all_genes_summary.py")


def run(code):
    out = tempfile.mkdtemp()
    gdir = os.path.join(out, "genes", "CYP2C19"); os.makedirs(os.path.join(gdir, "logs"))
    open(os.path.join(gdir, "logs", "stargazer.status"), "w").write(code)
    subprocess.run([sys.executable, os.path.join(HERE, "pgx-compare.py"), "--gene", "CYP2C19", "--sample", "S1",
                    "--output-dir", gdir, "--min-depth", "10"], check=True, capture_output=True)
    subprocess.run([sys.executable, SUMMARY, out, "S1", "CYP2C19"], check=True, capture_output=True)
    d = json.load(open(os.path.join(gdir, "CYP2C19_S1_detail.json")))
    tsv = open(os.path.join(gdir, "CYP2C19_S1_comparison.tsv")).read()
    summ = open(os.path.join(out, "all_genes_summary.tsv")).read().splitlines()[1]
    return d, tsv, summ


def main():
    for code in ("134", "137", "140", "143"):
        d, tsv, summ = run(code)
        assert d.get("killed_tools") == {"Stargazer": code}, (code, d.get("killed_tools"))
        assert f"caller killed (OOM/signal, exit {code})" in tsv, code
        assert "(caller killed: Stargazer)" in summ, summ
    d, tsv, summ = run("1")
    assert "killed_tools" not in d and "caller killed" not in tsv and "caller killed" not in summ
    d, tsv, summ = run("0")
    assert "killed_tools" not in d and "caller killed" not in summ
    print("test_killed_caller: PASS")


if __name__ == "__main__":
    main()
