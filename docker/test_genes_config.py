#!/usr/bin/env python3
"""Regression guard for the single-source genes.tsv (C1).

Freezes the gene matrices as they existed BEFORE genes.tsv (hardcoded in
pgx-compare.py and pgx-run.sh) and asserts that genes.tsv + the loader reproduce
them exactly. This is the safety net that substitutes for a live pipeline run:
if genes.tsv ever drifts from the validated matrix, this fails.

Run: python3 docker/test_genes_config.py
"""
import csv
import importlib.util
import os

HERE = os.path.dirname(os.path.abspath(__file__))
GENES_TSV = os.path.join(HERE, "genes.tsv")

# ── Frozen ground truth (transcribed from the pre-refactor source) ────────────
ORIG_SUPPORT = {
    "CYP2D6": {"pypgx", "stargazer", "aldy", "stellarpgx", "cyrius"},  # C: SV-aware CYP2D6 caller
    "CYP2C19": {"pypgx", "stargazer", "aldy", "stellarpgx"},
    "CYP2C9": {"pypgx", "stargazer", "aldy", "stellarpgx"},
    "CYP2B6": {"pypgx", "stargazer", "aldy", "stellarpgx", "pharmcat"},  # PharmCAT authority
    "CYP2C8": {"pypgx", "stargazer", "aldy", "stellarpgx"},
    "CYP3A4": {"pypgx", "stargazer", "aldy", "stellarpgx"},
    "CYP3A5": {"pypgx", "stargazer", "aldy", "stellarpgx"},
    "CYP4F2": {"pypgx", "stargazer", "aldy", "stellarpgx", "pharmcat"},  # PharmCAT authority
    "NUDT15": {"pypgx", "stargazer", "aldy", "stellarpgx"},
    "TPMT": {"pypgx", "stargazer", "aldy", "stellarpgx"},
    "UGT1A1": {"pypgx", "stargazer", "aldy", "stellarpgx", "vcf_check", "pharmcat"},  # PharmCAT authority
    "SLCO1B1": {"pypgx", "stargazer", "aldy", "stellarpgx"},
    "DPYD": {"pypgx", "stargazer", "aldy"},
    "NAT1": {"pypgx", "stargazer", "aldy", "stellarpgx"},
    "NAT2": {"pypgx", "stargazer", "aldy", "stellarpgx"},
    "G6PD": {"pypgx", "stargazer", "aldy", "vcf_check"},
    "GSTM1": {"pypgx", "stargazer", "aldy", "stellarpgx"},
    "GSTT1": {"pypgx", "stellarpgx"},
    "POR": {"pypgx", "stargazer", "stellarpgx"},
    "CYPOR": {"pypgx", "stargazer", "stellarpgx"},
    "VKORC1": {"pypgx", "stargazer", "aldy", "vcf_check"},  # path1: single-SNP VCF-Check (rs9923231)
    "CYP1A1": {"pypgx", "stargazer", "aldy", "stellarpgx"},
    "CYP1A2": {"pypgx", "stargazer", "aldy", "stellarpgx"},
    "CYP2A6": {"pypgx", "stargazer", "aldy", "stellarpgx"},
    "CYP2E1": {"pypgx", "stargazer", "aldy", "stellarpgx"},
    "IFNL3": {"pypgx", "stargazer", "aldy", "vcf_check"},  # path1: single-SNP VCF-Check (rs12979860)
    "RYR1": {"pypgx", "stargazer", "aldy", "vcf_check"},  # D: CPIC MH-variant VCF-Check added
    "ABCG2": {"aldy"},
    "HLA-A": {"optitype"},
    "HLA-B": {"optitype"},
    "HLA-C": {"optitype"},  # 2026-07: OptiType always typed C; it was parsed and dropped
    "CACNA1S": {"pypgx", "stellarpgx", "vcf_check"},
    "MT-RNR1": {"mutserve"},
    "CYP2C-CLUSTER": {"vcf_check"},  # rs12777823 single-SNP VCF-Check (CPIC warfarin)
    # GeT-RM-comparability genes (PyPGx-supported, low CPIC value; validated on NA12878)
    "GSTP1": {"pypgx", "aldy"},
    "SLC15A2": {"pypgx"},
    "SLC22A2": {"pypgx"},
    "SLCO2B1": {"pypgx"},
    "UGT2B15": {"pypgx"},
    # OmniGen add-ons: HLA class II (arcasHLA) + provisional ABO typer (VCF-Check)
    "HLA-DQA1": {"arcashla"},
    "HLA-DQB1": {"arcashla"},
    "HLA-DRB1": {"arcashla"},
    "ABO": {"vcf_check"},
}
ORIG_PYPGX_SV = {"CYP2A6", "CYP2B6", "CYP2D6", "CYP2E1", "CYP4F2", "G6PD", "GSTM1", "GSTT1",
                 "SLC22A2", "UGT2B15"}
ORIG_STARGAZER_SV = {"CYP2A6", "CYP2B6", "CYP2D6"}
ORIG_COORDS = {
    "CYP1A1": "chr15:74716541-74728528", "CYP1A2": "chr15:74745844-74759607",
    "CYP2A6": "chr19:40833540-40890447", "CYP2B6": "chr19:40921281-41028398",
    "CYP2C8": "chr10:95033771-95072497", "CYP2C9": "chr10:94935657-94993091",
    "CYP2C19": "chr10:94759680-94858547", "CYP2D6": "chr22:42116498-42155810",
    "CYP2E1": "chr10:133517362-133549123", "CYP3A4": "chr7:99753966-99787184",
    "CYP3A5": "chr7:99645193-99682996", "CYP4F2": "chr19:15863022-15913074",
    "DPYD": "chr1:97074742-97924034", "G6PD": "chrX:154528389-154550018",
    "GSTM1": "chr1:109684816-109696745", "GSTT1": "ALT_CONTIG",
    "IFNL3": "chr19:39240552-39253525", "NAT1": "chr8:18207108-18226689",
    "NAT2": "chr8:18388281-18404218", "NUDT15": "chr13:48034725-48050221",
    "POR": "chr7:75912154-75989855", "CYPOR": "chr7:75912154-75989855",
    "RYR1": "chr19:38430690-38590564", "SLCO1B1": "chr12:21128193-21242796",
    "TPMT": "chr6:18125310-18158169", "UGT1A1": "chr2:233754269-233779300",
    "VKORC1": "chr16:31087853-31097797", "ABCG2": "chr4:88085265-88236626",
    "HLA-A": "chr6:28510020-33480577", "HLA-B": "chr6:28510020-33480577",
    "HLA-C": "chr6:28510020-33480577",
    "CACNA1S": "chr1:201006956-201095000", "MT-RNR1": "chrM:648-1601",
    "CYP2C-CLUSTER": "chr10:94644745-94646745",
    "GSTP1": "chr11:67580811-67589653", "SLC15A2": "chr3:121891400-121947188",
    "SLC22A2": "chr6:160206754-160268821", "SLCO2B1": "chr11:75148106-75209549",
    "UGT2B15": "chr4:68640596-68676652",
    "HLA-DQA1": "chr6:28510020-33480577", "HLA-DQB1": "chr6:28510020-33480577",
    "HLA-DRB1": "chr6:28510020-33480577", "ABO": "chr9:133255000-133258000",
}
_VDR = {"CYP2D6", "CYP2C19", "CYP2C9", "CYP2B6", "CYP2C8", "CYP3A4", "CYP3A5",
        "CYP4F2", "CYP1A1", "CYP1A2", "CYP2A6", "CYP2E1", "SLCO1B1", "G6PD",
        "VKORC1", "NUDT15", "TPMT", "UGT1A1", "DPYD", "NAT1", "NAT2", "GSTM1",
        "IFNL3", "RYR1", "POR"}
ORIG_CONTROL = {g: ("vdr" if g in _VDR else "-") for g in ORIG_COORDS}

ALL_TOOLS = ("pypgx", "stargazer", "aldy", "stellarpgx", "optitype", "mutserve",
             "vcf_check", "cyrius", "pharmcat", "arcashla")


def _load_compare():
    spec = importlib.util.spec_from_file_location("c", os.path.join(HERE, "pgx-compare.py"))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_support_matches_frozen():
    c = _load_compare()
    support, pypgx_sv, stargazer_sv = c.load_gene_config(GENES_TSV)
    assert set(support) == set(ORIG_SUPPORT), set(support) ^ set(ORIG_SUPPORT)
    for g, orig_true in ORIG_SUPPORT.items():
        for tool in ALL_TOOLS:
            assert bool(support[g].get(tool)) == (tool in orig_true), (g, tool)
    assert set(pypgx_sv) == ORIG_PYPGX_SV
    assert set(stargazer_sv) == ORIG_STARGAZER_SV


def test_coords_and_control_match_frozen():
    with open(GENES_TSV) as fh:
        rows = {r["gene"]: r for r in csv.DictReader(fh, delimiter="\t")}
    assert set(rows) == set(ORIG_COORDS), set(rows) ^ set(ORIG_COORDS)
    for g, region in ORIG_COORDS.items():
        assert rows[g]["region"] == region, (g, rows[g]["region"], region)
        assert rows[g]["stargazer_control"] == ORIG_CONTROL[g], (g, rows[g]["stargazer_control"])


if __name__ == "__main__":
    test_support_matches_frozen()
    test_coords_and_control_match_frozen()
    print("all genes.tsv regression checks passed")
