#!/usr/bin/env python3
"""
docker/pgx-report.py — Generate HTML clinical PGx reports

Produces:
  <output>/<sample>.pgx.html          — landing page (gene cards + BAM QC)
  <output>/<sample>.<gene>.pgx.html   — per-gene 17X4 detail table

Usage:
  pgx-report.py --sample <SAMPLE_ID> --output <OUTPUT_DIR> [--bam-stats <JSON>]

The script reads:
  <output>/all_genes_summary.tsv
  <output>/<GENE>/<GENE>_<SAMPLE>_detail.json  (one per gene)
  <bam_stats>   (optional; bam_stats.json written by pgx-bamstats.sh)
"""

import argparse
import glob
import json
import os
import sys
from collections import Counter
from datetime import date

# ── Field metadata (17 rows in the detail table) ──────────────────────────────
FIELDS = [
    ("diplotype",             "Diplotype"),
    ("haplotype1",            "Haplotype 1"),
    ("haplotype2",            "Haplotype 2"),
    ("sub_alleles",           "Sub-alleles"),
    ("alternative_diplotypes","Alternative diplotypes"),
    ("phenotype",             "Phenotype"),
    ("activity_score",        "Activity score"),
    ("sv_type",               "SV type"),
    ("copy_number",           "Copy number"),
    ("supporting_variants",   "Supporting variants"),
    ("functional_effects",    "Functional effects"),
    ("dbsnp_ids",             "dbSNP IDs"),
    ("allele_score",          "Allele score"),
    ("mean_af",               "Mean allele frequency"),
    ("phasing_method",        "Phasing method"),
    ("status",                "Status"),
    ("sv_mode",               "SV mode note"),
]

TOOLS = ["PyPGx", "Stargazer", "Aldy", "StellarPGx", "OptiType", "mutserve"]

PHENOTYPE_COLORS = {
    "poor":         "#e53e3e",   # red
    "intermediate": "#ed8936",   # orange
    "normal":       "#38a169",   # green
    "ultrarapid":   "#3182ce",   # blue
    "rapid":        "#3182ce",
    "increased":    "#3182ce",
    "decreased":    "#ed8936",
    "indeterminate":"#718096",
    "uncertain":    "#718096",
    "favorable":    "#38a169",
    "unfavorable":  "#e53e3e",
}

def phenotype_color(pheno: str) -> str:
    if not pheno or pheno == "-":
        return "#718096"
    p = pheno.lower()
    for k, c in PHENOTYPE_COLORS.items():
        if k in p:
            return c
    return "#718096"


# ── CPIC / PharmVar clinical reference database ────────────────────────────────
# Static lookup for the 15 most actionable CPIC Level A/B gene-drug pairs.
# ── Gene loci (GRCh38, 1-based display coordinates) ──────────────────────────
# Used to show "chrN:start-end (GRCh38)" under the gene name on detail pages.
# HLA-A/HLA-B use the actual gene body coordinates, not the MHC extraction region.
# GSTT1 is on an alt contig.
GENE_LOCI: dict[str, str] = {
    "ABCG2":   "chr4:88,085,265-88,236,626",
    "CYP1A1":  "chr15:74,716,541-74,728,528",
    "CYP1A2":  "chr15:74,745,844-74,759,607",
    "CYP2A6":  "chr19:40,833,540-40,890,447",
    "CYP2B6":  "chr19:40,921,281-41,028,398",
    "CYP2C8":  "chr10:95,033,771-95,072,497",
    "CYP2C9":  "chr10:94,935,657-94,993,091",
    "CYP2C19": "chr10:94,759,680-94,858,547",
    "CYP2D6":  "chr22:42,116,498-42,155,810",
    "CYP2E1":  "chr10:133,517,362-133,549,123",
    "CYP3A4":  "chr7:99,753,966-99,787,184",
    "CYP3A5":  "chr7:99,645,193-99,682,996",
    "CYP4F2":  "chr19:15,863,022-15,913,074",
    "DPYD":    "chr1:97,074,742-97,924,034",
    "G6PD":    "chrX:154,528,389-154,550,018",
    "GSTM1":   "chr1:109,684,816-109,696,745",
    "GSTT1":   "chr22_KI270879v1_alt (alt contig)",
    "HLA-A":   "chr6:29,910,247-29,913,661",
    "HLA-B":   "chr6:31,321,649-31,324,666",
    "IFNL3":   "chr19:39,240,552-39,253,525",
    "NAT1":    "chr8:18,207,108-18,226,689",
    "NAT2":    "chr8:18,388,281-18,404,218",
    "NUDT15":  "chr13:48,034,725-48,050,221",
    "POR":     "chr7:75,912,154-75,989,855",
    "RYR1":    "chr19:38,430,690-38,590,564",
    "SLCO1B1": "chr12:21,128,193-21,242,796",
    "TPMT":    "chr6:18,125,310-18,158,169",
    "UGT1A1":  "chr2:233,754,269-233,779,300",
    "VKORC1":  "chr16:31,087,853-31,097,797",
    "CACNA1S": "chr1:201,006,956-201,083,927",
    "MT-RNR1": "chrM:648-1,601 (rCRS / GRCh38 chrM)",
}

# diplotype_check: lambda(diplotype_str) -> bool for genes where the diplotype
#   string must be inspected directly (VKORC1, DPYD non-reference, CYP4F2 *3,
#   RYR1 non-reference).
# Sources: cpicpgx.org, pharmvar.org (verify URLs at publication time).
CPIC_DB: dict = {
    "CYP2D6": {
        "desc": "Hepatic enzyme metabolising ~25% of commonly prescribed drugs including opioids, antidepressants, tamoxifen and atomoxetine.",
        "pharmvar_url": "https://www.pharmvar.org/gene/CYP2D6",
        "high_pheno":     ["poor", "ultrarapid"],
        "moderate_pheno": ["intermediate"],
        "diplotype_check": None,
        "landing_notes": {
            "poor":       "Codeine/tramadol contraindicated (inefficacy). Dose reduction required for TCAs, SSRIs, tamoxifen and atomoxetine.",
            "ultrarapid": "Codeine/tramadol contraindicated (rapid morphine conversion — life-threatening toxicity). Antidepressant efficacy likely reduced.",
            "intermediate": "Reduced CYP2D6 activity. Consider dose adjustment for codeine, tramadol, TCAs and SSRIs; monitor closely.",
        },
        "drugs": [
            {"name": "Codeine / Tramadol", "level": "A",
             "url": "https://cpicpgx.org/guidelines/guideline-for-codeine-and-cyp2d6/",
             "rec": "PM: contraindicated (inefficacy / toxicity risk). UM: contraindicated (ultra-rapid morphine conversion). IM: use with caution; consider alternative."},
            {"name": "Antidepressants (TCAs / SSRIs / SNRIs)", "level": "A",
             "url": "https://cpicpgx.org/guidelines/cpic-guideline-for-antidepressants-and-cyp2d6-and-cyp2c19/",
             "rec": "PM: reduce TCA dose or select CYP2D6-independent agent. UM: standard doses may be sub-therapeutic; consider alternative."},
            {"name": "Tamoxifen", "level": "A",
             "url": "https://cpicpgx.org/guidelines/guideline-for-tamoxifen-and-cyp2d6/",
             "rec": "PM: tamoxifen may be ineffective — consider aromatase inhibitor where appropriate."},
            {"name": "Atomoxetine", "level": "A",
             "url": "https://cpicpgx.org/guidelines/guideline-for-atomoxetine-and-cyp2d6/",
             "rec": "PM: initiate at 50% standard dose; titrate slowly. UM: standard dose may be sub-therapeutic."},
        ],
    },
    "CYP2C19": {
        "desc": "Activates clopidogrel (prodrug); metabolises PPIs, antidepressants and voriconazole. Both loss-of-function (PM) and gain-of-function (UM/RM) phenotypes have distinct clinical consequences.",
        "pharmvar_url": "https://www.pharmvar.org/gene/CYP2C19",
        "high_pheno":     ["poor", "ultrarapid"],
        "moderate_pheno": ["rapid", "intermediate"],
        "diplotype_check": None,
        "landing_notes": {
            "poor":         "Clopidogrel contraindicated — inadequate antiplatelet activation. Use prasugrel or ticagrelor. PPIs and SSRIs may accumulate.",
            "ultrarapid":   "Voriconazole likely sub-therapeutic. PPI efficacy reduced. Monitor antidepressant levels.",
            "rapid":        "Increased CYP2C19 activity. Consider voriconazole dose increase. Monitor PPI efficacy.",
            "intermediate": "Mildly reduced CYP2C19 activity. Monitor clopidogrel response; consider alternative antiplatelet.",
        },
        "drugs": [
            {"name": "Clopidogrel", "level": "A",
             "url": "https://cpicpgx.org/guidelines/guideline-for-clopidogrel-and-cyp2c19/",
             "rec": "PM/IM: use alternative antiplatelet (prasugrel or ticagrelor). UM/RM: standard dose; enhanced response expected."},
            {"name": "Antidepressants (SSRIs / TCAs)", "level": "A",
             "url": "https://cpicpgx.org/guidelines/cpic-guideline-for-antidepressants-and-cyp2d6-and-cyp2c19/",
             "rec": "PM: reduce dose or select alternative. UM: standard or higher dose; monitor for sub-therapeutic levels."},
            {"name": "Voriconazole", "level": "A",
             "url": "https://cpicpgx.org/guidelines/guideline-for-voriconazole-and-cyp2c19/",
             "rec": "PM: reduce dose substantially (elevated exposure, toxicity risk). UM/RM: TDM required; standard doses may be sub-therapeutic."},
            {"name": "Proton pump inhibitors (PPIs)", "level": "A",
             "url": "https://cpicpgx.org/guidelines/cpic-guideline-for-proton-pump-inhibitors-and-cyp2c19/",
             "rec": "PM: consider dose reduction (increased AUC). UM: may require higher dose for adequate acid suppression."},
        ],
    },
    "CYP2C9": {
        "desc": "Metabolises the S-enantiomer of warfarin, NSAIDs and phenytoin. Poor metabolisers require significantly reduced doses of narrow-therapeutic-index drugs.",
        "pharmvar_url": "https://www.pharmvar.org/gene/CYP2C9",
        "high_pheno":     ["poor"],
        "moderate_pheno": ["intermediate"],
        "diplotype_check": None,
        "landing_notes": {
            "poor":         "Warfarin sensitivity: significantly reduced starting dose required. NSAID and phenytoin toxicity risk elevated.",
            "intermediate": "Moderate warfarin sensitivity. Consider dose reduction for NSAIDs and phenytoin; TDM recommended.",
        },
        "drugs": [
            {"name": "Warfarin", "level": "A",
             "url": "https://cpicpgx.org/guidelines/guideline-for-warfarin-and-cyp2c9-and-vkorc1/",
             "rec": "PM: initiate at ≤25% of standard dose and titrate to INR. IM: reduce starting dose; use CPIC dosing calculator."},
            {"name": "NSAIDs", "level": "A",
             "url": "https://cpicpgx.org/guidelines/guideline-for-nonsteroidal-anti-inflammatory-drugs-and-cyp2c9/",
             "rec": "PM: use lowest effective dose or select alternative. Increased GI and renal toxicity risk."},
            {"name": "Phenytoin / Fosphenytoin", "level": "A",
             "url": "https://cpicpgx.org/guidelines/guideline-for-phenytoin-and-cyp2c9-and-hla-b/",
             "rec": "PM/IM: reduce maintenance dose by 25–50%; TDM essential to avoid toxicity."},
        ],
    },
    "DPYD": {
        "desc": "Rate-limiting enzyme in 5-FU and capecitabine catabolism. Even heterozygous loss-of-function variants confer significant toxicity risk — CPIC recommends pre-emptive genotyping.",
        "pharmvar_url": None,
        "high_pheno":     ["poor"],
        "moderate_pheno": ["intermediate", "decreased"],
        "diplotype_check": lambda d: bool(d) and d not in ("-",) and any(
            v.strip() not in ("*1", "Reference", "") for v in d.split("/")
        ),
        "landing_notes": {
            "poor":       "Fluoropyrimidines (5-FU/capecitabine) are contraindicated — life-threatening toxicity risk.",
            "intermediate": "Fluoropyrimidine dose reduction (≥50%) required before initiating therapy. Consult CPIC/DPWG guidelines.",
            "decreased":  "Reduced DPYD function. Fluoropyrimidine dose adjustment or avoidance required.",
            "_diplotype": "Non-reference DPYD variant detected. Fluoropyrimidine (5-FU/capecitabine) dose adjustment or avoidance may be required — consult CPIC guidelines before prescribing.",
        },
        "drugs": [
            {"name": "5-Fluorouracil (5-FU) / Capecitabine", "level": "A",
             "url": "https://cpicpgx.org/guidelines/guideline-for-fluoropyrimidines-and-dpyd/",
             "rec": "PM: contraindicated. IM/Decreased function: ≥50% dose reduction with therapeutic drug monitoring. Normal: standard dosing."},
        ],
    },
    "TPMT": {
        "desc": "Inactivates thiopurines (azathioprine, 6-mercaptopurine, thioguanine). Deficient patients accumulate toxic thioguanine nucleotides causing severe, potentially fatal myelosuppression.",
        "pharmvar_url": "https://www.pharmvar.org/gene/TPMT",
        "high_pheno":     ["poor"],
        "moderate_pheno": ["intermediate"],
        "diplotype_check": None,
        "landing_notes": {
            "poor":         "Thiopurines contraindicated at standard doses — life-threatening myelosuppression. Use alternative or dramatically reduce dose with TDM.",
            "intermediate": "Reduce thiopurine starting dose by 30–70% and titrate based on toxicity and response.",
        },
        "drugs": [
            {"name": "Azathioprine / 6-Mercaptopurine / Thioguanine", "level": "A",
             "url": "https://cpicpgx.org/guidelines/guideline-for-thiopurines-and-tpmt-and-nudt15/",
             "rec": "PM: avoid or use ≤10% of standard dose with TDM. IM: reduce starting dose 30–70%; monitor CBC."},
        ],
    },
    "NUDT15": {
        "desc": "Hydrolyses thioguanine nucleotide metabolites. Loss-of-function variants cause severe thiopurine-induced myelosuppression; particularly common in East Asian and Hispanic populations.",
        "pharmvar_url": None,
        "high_pheno":     ["poor"],
        "moderate_pheno": ["intermediate"],
        "diplotype_check": None,
        "landing_notes": {
            "poor":         "Thiopurines contraindicated at standard doses — severe myelosuppression risk. Consider non-thiopurine immunosuppressant.",
            "intermediate": "Reduce thiopurine starting dose; monitor CBC closely for myelosuppression.",
        },
        "drugs": [
            {"name": "Azathioprine / 6-Mercaptopurine / Thioguanine", "level": "A",
             "url": "https://cpicpgx.org/guidelines/guideline-for-thiopurines-and-tpmt-and-nudt15/",
             "rec": "PM: avoid or markedly reduce dose with TDM. IM: reduce starting dose; CBC monitoring mandatory."},
        ],
    },
    "SLCO1B1": {
        "desc": "Hepatic uptake transporter (OATP1B1) for statins. Decreased function variants reduce hepatic statin clearance, increasing plasma exposure and myopathy/rhabdomyolysis risk.",
        "pharmvar_url": None,
        "high_pheno":     ["poor", "decreased"],
        "moderate_pheno": ["intermediate"],
        "diplotype_check": None,
        "landing_notes": {
            "poor":         "High statin-associated myopathy risk. Use ≤20 mg/day simvastatin or switch to rosuvastatin/pravastatin.",
            "decreased":    "Increased statin myopathy risk. Avoid high-dose simvastatin; consider lower-risk statin.",
            "intermediate": "Moderately increased statin exposure. Prefer rosuvastatin or pravastatin over simvastatin.",
        },
        "drugs": [
            {"name": "Simvastatin / Statins", "level": "A",
             "url": "https://cpicpgx.org/guidelines/guideline-for-simvastatin-and-slco1b1/",
             "rec": "Poor/Decreased function: avoid simvastatin >20 mg/day; use rosuvastatin or pravastatin. IM: simvastatin ≤40 mg/day."},
        ],
    },
    "UGT1A1": {
        "desc": "Glucuronidates irinotecan's active metabolite (SN-38) and bilirubin. Reduced activity (*28/*28, Gilbert's syndrome) substantially increases irinotecan toxicity risk.",
        "pharmvar_url": None,
        "high_pheno":     ["poor"],
        "moderate_pheno": ["intermediate", "decreased"],
        "diplotype_check": None,
        "landing_notes": {
            "poor":         "High irinotecan toxicity risk (*28/*28 — Gilbert's). Reduce irinotecan starting dose by at least one dose level.",
            "intermediate": "Mildly elevated irinotecan exposure. Standard starting dose with close toxicity monitoring.",
            "decreased":    "Reduced UGT1A1 function — monitor for irinotecan and atazanavir toxicity.",
        },
        "drugs": [
            {"name": "Irinotecan", "level": "A",
             "url": "https://cpicpgx.org/guidelines/guideline-for-irinotecan-and-ugt1a1/",
             "rec": "PM (*28/*28): reduce starting dose by ≥1 dose level; titrate by toxicity. IM: standard dose with close monitoring."},
            {"name": "Atazanavir", "level": "A",
             "url": "https://cpicpgx.org/guidelines/guideline-for-atazanavir-and-ugt1a1/",
             "rec": "PM: increased unconjugated bilirubin and jaundice risk — consider alternative antiretroviral."},
        ],
    },
    "CYP2B6": {
        "desc": "Primary metaboliser of efavirenz and several other antiretrovirals. Poor metabolisers have ~4–5-fold elevated efavirenz plasma concentrations causing CNS toxicity.",
        "pharmvar_url": "https://www.pharmvar.org/gene/CYP2B6",
        "high_pheno":     ["poor"],
        "moderate_pheno": ["intermediate"],
        "diplotype_check": None,
        "landing_notes": {
            "poor":         "Efavirenz plasma levels ~4–5X elevated. Reduce efavirenz to 400 mg/day or switch to alternative ART regimen.",
            "intermediate": "Moderately elevated efavirenz exposure. Monitor for CNS side effects; consider dose adjustment.",
        },
        "drugs": [
            {"name": "Efavirenz", "level": "A",
             "url": "https://cpicpgx.org/guidelines/guideline-for-efavirenz-and-cyp2b6/",
             "rec": "PM: reduce dose to 400 mg/day or switch to alternative ART. IM: standard dose with CNS toxicity monitoring."},
        ],
    },
    "CYP3A5": {
        "desc": "Major metaboliser of tacrolimus and other calcineurin inhibitors. Expressers (*1 carriers) require higher tacrolimus doses; non-expressers (*3/*3) achieve target levels at lower doses.",
        "pharmvar_url": "https://www.pharmvar.org/gene/CYP3A5",
        "high_pheno":     [],
        "moderate_pheno": ["poor", "non-expresser", "intermediate"],
        "diplotype_check": None,
        "landing_notes": {
            "poor":         "CYP3A5 non-expresser (*3/*3): tacrolimus dose requirements typical for non-expressers. Use standard weight-based dosing with TDM.",
            "intermediate": "Intermediate CYP3A5 expression. TDM-guided tacrolimus dosing essential.",
        },
        "drugs": [
            {"name": "Tacrolimus", "level": "A",
             "url": "https://cpicpgx.org/guidelines/guideline-for-tacrolimus-and-cyp3a5/",
             "rec": "Non-expresser (*3/*3): standard dose + TDM. Expresser (*1 carrier): initiate at 1.5–2X standard dose to achieve target trough."},
        ],
    },
    "VKORC1": {
        "desc": "Warfarin target enzyme. rs9923231 (−1639G>A) in the VKORC1 promoter reduces enzyme expression, substantially increasing warfarin sensitivity and bleeding risk.",
        "pharmvar_url": None,
        "high_pheno":     [],
        "moderate_pheno": [],
        # PyPGx reports rs9923231; Stargazer reports *S1 (sensitive); Aldy reports *H1/*H2.
        # All notations represent the same promoter variant — check for any of them.
        "diplotype_check": lambda d: any(x in (d or "") for x in ("rs9923231", "*S1", "*H1", "*H2")),
        "landing_notes": {
            "_diplotype": "VKORC1 rs9923231 detected — increased warfarin sensitivity. Initiate warfarin at a reduced dose. A/A homozygotes require the greatest reduction (~3 mg/day). Use CPIC warfarin dosing calculator.",
        },
        "drugs": [
            {"name": "Warfarin", "level": "A",
             "url": "https://cpicpgx.org/guidelines/guideline-for-warfarin-and-cyp2c9-and-vkorc1/",
             "rec": "A/A homozygote: initiate at ~3 mg/day and titrate to INR. A/G heterozygote: intermediate sensitivity; reduce starting dose. Use CPIC calculator incorporating CYP2C9 + VKORC1 + CYP4F2."},
        ],
    },
    "G6PD": {
        "desc": "X-linked enzyme protecting erythrocytes from oxidative damage. Deficiency causes acute drug-induced haemolytic anaemia with rasburicase, primaquine and other oxidant agents.",
        "pharmvar_url": None,
        "high_pheno":     ["deficient", "hemizygous", "homozygous"],
        "moderate_pheno": ["variable", "intermediate", "heterozygous"],
        "diplotype_check": None,
        "landing_notes": {
            "deficient":    "G6PD deficiency: rasburicase is contraindicated (severe haemolytic anaemia risk). Avoid primaquine, dapsone and high-dose aspirin.",
            "variable":     "Partial G6PD activity. Avoid high-risk oxidant drugs; rasburicase requires careful risk-benefit assessment.",
            "intermediate": "Intermediate G6PD activity. Monitor for haemolytic events with oxidant drug exposure.",
        },
        "drugs": [
            {"name": "Rasburicase", "level": "A",
             "url": "https://cpicpgx.org/guidelines/guideline-for-rasburicase-and-g6pd/",
             "rec": "Deficient: contraindicated — risk of severe, potentially fatal haemolytic anaemia. Variable: use with extreme caution, monitoring haemoglobin."},
        ],
    },
    "NAT2": {
        "desc": "N-acetyltransferase 2 metabolises isoniazid and hydralazine. Slow acetylators accumulate isoniazid, increasing peripheral neuropathy and hepatotoxicity risk.",
        "pharmvar_url": None,
        "high_pheno":     [],
        "moderate_pheno": ["slow", "poor", "intermediate"],
        "diplotype_check": None,
        "landing_notes": {
            "slow":         "NAT2 slow acetylator: isoniazid accumulates — supplement with pyridoxine (vitamin B6) 25–50 mg/day and monitor for peripheral neuropathy.",
            "poor":         "NAT2 slow acetylator phenotype. Isoniazid dose adjustment or pyridoxine supplementation may be required.",
            "intermediate": "Intermediate NAT2 activity. Monitor isoniazid tolerability.",
        },
        "drugs": [
            {"name": "Isoniazid", "level": "A",
             "url": "https://cpicpgx.org/guidelines/guideline-for-isoniazid-and-nat2/",
             "rec": "Slow acetylator: supplement with pyridoxine 25–50 mg/day; monitor for neuropathy and hepatotoxicity. Consider dose adjustment if standard dose not tolerated."},
        ],
    },
    "RYR1": {
        "desc": "Skeletal muscle calcium release channel. Pathogenic variants confer susceptibility to malignant hyperthermia (MH) — a life-threatening reaction to volatile anaesthetics and succinylcholine.",
        "pharmvar_url": None,
        "high_pheno":     ["susceptible", "pathogenic"],
        "moderate_pheno": [],
        "diplotype_check": lambda d: bool(d) and d not in (None, "-", "Reference/Reference", "*1/*1", ""),
        "landing_notes": {
            "_diplotype":  "Non-reference RYR1 variant detected. Evaluate for malignant hyperthermia susceptibility — avoid volatile anaesthetics and succinylcholine if susceptible.",
            "susceptible": "Malignant hyperthermia susceptibility: avoid all volatile anaesthetics and succinylcholine. Use TIVA. Ensure dantrolene is immediately available.",
        },
        "drugs": [
            {"name": "Volatile anaesthetic agents", "level": "A",
             "url": "https://cpicpgx.org/guidelines/guideline-for-volatile-anesthetic-agents-and-succinylcholine-and-ryr1-and-cacna1s/",
             "rec": "MH-susceptible: contraindicated. Use total intravenous anaesthesia (TIVA) only; dantrolene must be immediately available."},
            {"name": "Succinylcholine", "level": "A",
             "url": "https://cpicpgx.org/guidelines/guideline-for-volatile-anesthetic-agents-and-succinylcholine-and-ryr1-and-cacna1s/",
             "rec": "MH-susceptible: contraindicated. Use non-depolarising neuromuscular blockers (e.g. rocuronium with sugammadex reversal)."},
        ],
    },
    "CACNA1S": {
        "desc": "Skeletal muscle L-type calcium channel subunit. Pathogenic variants (e.g. *2 Arg174Trp, *3 Arg1086His) confer malignant hyperthermia susceptibility (MHS) — a life-threatening pharmacogenomic reaction to volatile anaesthetics and succinylcholine. Partner gene to RYR1 under the same CPIC guideline.",
        "pharmvar_url": None,
        "high_pheno":     ["susceptible", "pathogenic", "mhs"],
        "moderate_pheno": [],
        "diplotype_check": lambda d: bool(d) and d not in (None, "-", "Reference/Reference", "*1/*1", ""),
        "landing_notes": {
            "_diplotype":  "Non-reference CACNA1S variant detected. Evaluate for malignant hyperthermia susceptibility — avoid volatile anaesthetics and succinylcholine if susceptible.",
            "susceptible": "Malignant hyperthermia susceptibility: avoid all volatile anaesthetics and succinylcholine. Use TIVA. Ensure dantrolene is immediately available.",
            "mhs":         "Malignant hyperthermia susceptibility: avoid all volatile anaesthetics and succinylcholine. Use TIVA. Ensure dantrolene is immediately available.",
        },
        "drugs": [
            {"name": "Volatile anaesthetic agents", "level": "A",
             "url": "https://cpicpgx.org/guidelines/guideline-for-volatile-anesthetic-agents-and-succinylcholine-and-ryr1-and-cacna1s/",
             "rec": "MH-susceptible: contraindicated. Use total intravenous anaesthesia (TIVA) only; dantrolene must be immediately available."},
            {"name": "Succinylcholine", "level": "A",
             "url": "https://cpicpgx.org/guidelines/guideline-for-volatile-anesthetic-agents-and-succinylcholine-and-ryr1-and-cacna1s/",
             "rec": "MH-susceptible: contraindicated. Use non-depolarising neuromuscular blockers (e.g. rocuronium with sugammadex reversal)."},
        ],
        "min_tools": 1,
    },
    "MT-RNR1": {
        "desc": "Mitochondrial 12S ribosomal RNA gene. Variants m.1555A>G (~1/500 Europeans) and m.1494C>T confer ribosomal hypersensitivity to aminoglycosides — even a single standard dose can cause permanent bilateral sensorineural hearing loss. Maternally inherited; affects all maternal-lineage family members.",
        "pharmvar_url": None,
        "high_pheno":     ["ototoxicity risk", "carrier"],
        "moderate_pheno": [],
        "diplotype_check": lambda d: bool(d) and d not in (None, "-", "Reference", ""),
        "landing_notes": {
            "_diplotype":           "CPIC Level A MT-RNR1 variant detected. Avoid aminoglycosides unless no safe alternatives exist.",
            "aminoglycoside-ototoxicity risk": "Aminoglycoside-ototoxicity risk: avoid gentamicin, tobramycin, amikacin, streptomycin and neomycin unless essential. Counsel patient on irreversible hearing-loss risk. Inform maternal-lineage family members.",
        },
        "drugs": [
            {"name": "Aminoglycosides (gentamicin, tobramycin, amikacin, streptomycin, neomycin)", "level": "A",
             "url": "https://cpicpgx.org/guidelines/guideline-for-aminoglycosides-and-mt-rnr1/",
             "rec": "m.1555A>G or m.1494C>T carrier: avoid aminoglycosides if alternatives exist. If unavoidable, use lowest effective dose with audiological monitoring. Inform maternal relatives."},
        ],
        "min_tools": 1,
    },
    "CYP4F2": {
        "desc": "Metabolises vitamin K1; V433M (*3 allele) reduces vitamin K catabolism, modestly increasing warfarin sensitivity. Incorporated in CPIC warfarin dosing algorithm.",
        "pharmvar_url": None,
        "high_pheno":     [],
        "moderate_pheno": [],
        "diplotype_check": lambda d: "*3" in (d or ""),
        "landing_notes": {
            "_diplotype": "CYP4F2 *3 variant detected — minor contribution to warfarin sensitivity. Use CPIC warfarin dosing calculator incorporating CYP2C9 + VKORC1 + CYP4F2 genotype.",
        },
        "drugs": [
            {"name": "Warfarin", "level": "B",
             "url": "https://cpicpgx.org/guidelines/guideline-for-warfarin-and-cyp2c9-and-vkorc1/",
             "rec": "*3 carriers may require slightly higher warfarin doses. Incorporate into multi-gene CPIC dosing algorithm alongside CYP2C9 and VKORC1."},
        ],
    },
    "ABCG2": {
        "desc": "Efflux transporter (BCRP) reducing intestinal and hepatic uptake of rosuvastatin and other statins. The *2 variant (rs2231142, c.421C>A, p.Gln141Lys) decreases transporter function, increasing statin plasma exposure.",
        "pharmvar_url": "https://www.pharmvar.org/gene/ABCG2",
        "high_pheno":     [],
        "moderate_pheno": ["decreased", "intermediate"],
        "diplotype_check": lambda d: "rs2231142" in (d or ""),
        "landing_notes": {
            "_diplotype": "ABCG2 rs2231142 (c.421C>A) detected — decreased rosuvastatin transport. Consider lower rosuvastatin starting dose (≤10 mg/day) per CPIC guidelines.",
            "decreased":  "ABCG2 decreased function — reduced rosuvastatin efflux increases plasma exposure. Use lower starting dose.",
        },
        "drugs": [
            {"name": "Rosuvastatin", "level": "A",
             "url": "https://cpicpgx.org/guidelines/cpic-guideline-for-statins/",
             "rec": "ABCG2 decreased function (*2/*2 or *2/ref): initiate rosuvastatin at ≤10 mg/day. Consider alternative statin if high-intensity therapy required."},
        ],
    },
    "HLA-A": {
        "desc": "Major histocompatibility complex class I gene. HLA-A*31:01 is associated with carbamazepine-induced severe cutaneous adverse reactions (SCAR) including DRESS, SJS, and TEN across all ancestries.",
        "pharmvar_url": None,
        "high_pheno":     [],
        "moderate_pheno": [],
        "min_tools": 1,   # OptiType is the sole HLA typing tool
        "diplotype_check": lambda d: "A*31:01" in (d or ""),
        "risk_alleles": {
            "A*31:01": "carbamazepine/oxcarbazepine: increased DRESS/SJS/TEN risk across all ancestries",
        },
        "no_risk_note": "No CPIC Level A high-risk HLA-A alleles (A*31:01) detected.",
        "landing_notes": {
            "_diplotype": "HLA-A*31:01 detected — carbamazepine and oxcarbazepine associated with severe cutaneous adverse reactions (DRESS, SJS, TEN) across all ancestries. Avoid if clinically feasible.",
        },
        "drugs": [
            {"name": "Carbamazepine", "level": "A",
             "url": "https://cpicpgx.org/guidelines/guideline-for-carbamazepine-and-hla-b/",
             "rec": "HLA-A*31:01 positive: avoid carbamazepine — significantly increased risk of DRESS, SJS, TEN in all ancestries. Use alternative anticonvulsant."},
            {"name": "Oxcarbazepine", "level": "B",
             "url": "https://cpicpgx.org/guidelines/guideline-for-oxcarbazepine-and-hla-b-and-hla-a/",
             "rec": "HLA-A*31:01 positive: consider alternatives to oxcarbazepine; increased SCAR risk."},
        ],
    },
    "HLA-B": {
        "desc": "Major histocompatibility complex class I gene. Multiple HLA-B alleles are associated with severe drug hypersensitivity: HLA-B*57:01 (abacavir), HLA-B*58:01 (allopurinol), HLA-B*15:02 (carbamazepine/phenytoin in Asian ancestry), HLA-B*13:01 (dapsone-induced DRESS/HSS, prevalent in Asian populations), HLA-B*57:03 (flucloxacillin-induced DILI, prevalent in European populations).",
        "pharmvar_url": None,
        "high_pheno":     [],
        "moderate_pheno": [],
        "min_tools": 1,   # OptiType is the sole HLA typing tool
        "diplotype_check": lambda d: any(x in (d or "") for x in ("B*57:01", "B*15:02", "B*58:01", "B*13:01", "B*57:03")),
        "risk_alleles": {
            "B*57:01": "abacavir: contraindicated — high hypersensitivity risk",
            "B*15:02": "carbamazepine / phenytoin: avoid — high SJS/TEN risk (esp. Asian ancestry)",
            "B*58:01": "allopurinol: avoid — high SJS/TEN risk",
            "B*13:01": "dapsone: avoid — high risk of DRESS/HSS (esp. Asian ancestry)",
            "B*57:03": "flucloxacillin: use with caution — increased DILI risk (esp. European ancestry)",
        },
        "no_risk_note": "No high-risk HLA-B alleles (B*57:01, B*15:02, B*58:01 [CPIC Level A]; B*13:01, B*57:03 [CPIC Level B]) detected.",
        "landing_notes": {
            "_diplotype": "High-risk HLA-B allele detected — see detail page for drug-specific recommendations.",
        },
        "drugs": [
            {"name": "Abacavir", "level": "A",
             "url": "https://cpicpgx.org/guidelines/guideline-for-abacavir-and-hla-b/",
             "rec": "HLA-B*57:01 positive: abacavir contraindicated — high risk of severe, potentially life-threatening hypersensitivity reaction. Use alternative ART regimen."},
            {"name": "Allopurinol", "level": "A",
             "url": "https://cpicpgx.org/guidelines/guideline-for-allopurinol-and-hla-b/",
             "rec": "HLA-B*58:01 positive: avoid allopurinol — high risk of Stevens-Johnson syndrome / toxic epidermal necrolysis. Use alternative urate-lowering therapy."},
            {"name": "Carbamazepine", "level": "A",
             "url": "https://cpicpgx.org/guidelines/guideline-for-carbamazepine-and-hla-b/",
             "rec": "HLA-B*15:02 positive: avoid carbamazepine — high SJS/TEN risk especially in Han Chinese and Southeast Asian ancestry. Use alternative anticonvulsant."},
            {"name": "Phenytoin / Fosphenytoin", "level": "A",
             "url": "https://cpicpgx.org/guidelines/guideline-for-phenytoin-and-cyp2c9-and-hla-b/",
             "rec": "HLA-B*15:02 positive: increased SJS risk. Avoid or use with extreme caution; consider alternative anticonvulsant."},
            {"name": "Dapsone", "level": "B",
             "url": "https://cpicpgx.org/guidelines/cpic-guideline-for-dapsone-and-hla-b/",
             "rec": "HLA-B*13:01 positive: dapsone associated with high risk of DRESS (drug reaction with eosinophilia and systemic symptoms) and HSS (hypersensitivity syndrome), particularly in Asian populations. Consider alternative antimicrobial (e.g., trimethoprim-sulfamethoxazole, atovaquone)."},
            {"name": "Flucloxacillin", "level": "B",
             "url": "https://www.pharmgkb.org/chemical/PA10001/clinicalAnnotation",
             "rec": "HLA-B*57:03 positive: flucloxacillin associated with increased risk of drug-induced liver injury (DILI), predominantly reported in individuals of European ancestry. Use with caution and monitor liver function; consider an alternative narrow-spectrum beta-lactam if available."},
        ],
    },
    "CYP1A1": {
        "desc": "CYP1A1 metabolises polycyclic aromatic hydrocarbons and some drugs. No CPIC Level A/B guideline exists; informational star-allele calling only.",
        "pharmvar_url": "https://www.pharmvar.org/gene/CYP1A1",
        "high_pheno":     [],
        "moderate_pheno": [],
        "diplotype_check": None,
        "landing_notes": {},
        "drugs": [],
    },
    "CYP1A2": {
        "desc": "CYP1A2 is the primary enzyme metabolising fluvoxamine and accounts for ~13% of total hepatic CYP activity. Ultra-rapid metabolisers may require higher fluvoxamine doses; poor metabolisers are at risk of toxicity.",
        "pharmvar_url": "https://www.pharmvar.org/gene/CYP1A2",
        "high_pheno":     [],
        "moderate_pheno": ["poor", "intermediate"],
        "diplotype_check": None,
        "landing_notes": {
            "poor":         "CYP1A2 poor metaboliser: fluvoxamine plasma levels may be markedly elevated — start at lowest effective dose and titrate carefully.",
            "intermediate": "CYP1A2 intermediate metaboliser: consider dose reduction for narrow-therapeutic-index CYP1A2 substrates.",
            "ultra-rapid":  "CYP1A2 ultra-rapid metaboliser: standard fluvoxamine doses may be subtherapeutic.",
        },
        "drugs": [
            {"name": "Fluvoxamine", "level": "A",
             "url": "https://cpicpgx.org/guidelines/cpic-guideline-for-fluvoxamine-and-cyp1a2/",
             "rec": "Poor metaboliser: initiate at lowest available dose; titrate with close monitoring. Ultra-rapid: may need higher doses or alternative SSRI."},
        ],
    },
    "CYP2A6": {
        "desc": "CYP2A6 is the principal enzyme catalysing nicotine C-oxidation to cotinine. Poor and intermediate metabolisers smoke fewer cigarettes but may have lower quit rates with standard nicotine replacement therapy.",
        "pharmvar_url": "https://www.pharmvar.org/gene/CYP2A6",
        "high_pheno":     [],
        "moderate_pheno": ["poor", "intermediate"],
        "diplotype_check": None,
        "landing_notes": {
            "poor":         "CYP2A6 poor metaboliser: nicotine clearance is reduced — may smoke fewer cigarettes per day. Standard NRT patch may be sufficient; monitor for over-replacement.",
            "intermediate": "CYP2A6 intermediate metaboliser: slightly reduced nicotine metabolism. Standard cessation therapy is appropriate.",
        },
        "drugs": [
            {"name": "Nicotine / Tobacco cessation", "level": "A",
             "url": "https://cpicpgx.org/guidelines/cpic-guideline-for-nicotine-and-cyp2a6/",
             "rec": "Poor metaboliser: nicotine elimination is slow. Patch-based NRT is preferred over short-acting forms. Dose reduction may be needed to prevent adverse effects."},
        ],
    },
    "CYP2C8": {
        "desc": "CYP2C8 metabolises amodiaquine, pioglitazone, paclitaxel and other drugs. Poor metabolisers are at increased risk of amodiaquine-induced agranulocytosis.",
        "pharmvar_url": "https://www.pharmvar.org/gene/CYP2C8",
        "high_pheno":     [],
        "moderate_pheno": ["poor", "intermediate"],
        "diplotype_check": None,
        "landing_notes": {
            "poor":         "CYP2C8 poor metaboliser: avoid amodiaquine — increased risk of agranulocytosis. Use alternative antimalarial.",
            "intermediate": "CYP2C8 intermediate metaboliser: monitor for amodiaquine toxicity.",
        },
        "drugs": [
            {"name": "Amodiaquine", "level": "B",
             "url": "https://cpicpgx.org/guidelines/cpic-guideline-for-amodiaquine-and-cyp2c8/",
             "rec": "Poor metaboliser: avoid amodiaquine due to increased agranulocytosis risk. Use alternative antimalarial (e.g. artemisinin-based combination therapy)."},
        ],
    },
    "CYP2E1": {
        "desc": "CYP2E1 metabolises ethanol, paracetamol (acetaminophen), and industrial solvents. No CPIC guideline with actionable recommendations currently exists.",
        "pharmvar_url": "https://www.pharmvar.org/gene/CYP2E1",
        "high_pheno":     [],
        "moderate_pheno": [],
        "diplotype_check": None,
        "landing_notes": {},
        "drugs": [],
    },
    "CYP3A4": {
        "desc": "CYP3A4 is the most abundant hepatic CYP enzyme, metabolising ~50% of drugs. The *22 allele (rs35599367) reduces CYP3A4 expression and is relevant to tacrolimus dosing. No standalone high-evidence CPIC guideline for CYP3A4 genotyping alone.",
        "pharmvar_url": "https://www.pharmvar.org/gene/CYP3A4",
        "high_pheno":     [],
        "moderate_pheno": [],
        "diplotype_check": lambda d: "*22" in (d or ""),
        "landing_notes": {
            "_diplotype": "CYP3A4 *22 allele detected — reduced CYP3A4 expression. May contribute to altered tacrolimus and other CYP3A4 substrate metabolism. Interpret alongside CYP3A5 genotype.",
        },
        "drugs": [
            {"name": "Tacrolimus (with CYP3A5)", "level": "A",
             "url": "https://cpicpgx.org/guidelines/guideline-for-tacrolimus-and-cyp3a5/",
             "rec": "CYP3A4 *22 in a CYP3A5 non-expresser: combined poor CYP3A metaboliser — consider lower tacrolimus starting dose and target lower initial trough concentrations."},
        ],
    },
    "GSTM1": {
        "desc": "GSTM1 (glutathione S-transferase mu 1) is commonly deleted (*0/*0 null genotype), reducing glutathione conjugation capacity. The null genotype is present in ~50% of Europeans and ~35% of Asians. No CPIC Level A/B guideline with dose-adjustment recommendations currently exists.",
        "pharmvar_url": None,
        "high_pheno":     [],
        "moderate_pheno": [],
        "diplotype_check": None,
        "landing_notes": {},
        "drugs": [],
    },
    "GSTT1": {
        "desc": "GSTT1 (glutathione S-transferase theta 1) null genotype (*0/*0) is present in ~15–20% of Europeans. No CPIC Level A/B guideline with dose-adjustment recommendations currently exists.",
        "pharmvar_url": None,
        "high_pheno":     [],
        "moderate_pheno": [],
        "diplotype_check": None,
        "landing_notes": {},
        "drugs": [],
    },
    "IFNL3": {
        "desc": "IFNL3 (IL28B) encodes interferon lambda-3. The rs12979860 C/C genotype strongly predicts favourable response to pegylated interferon alfa + ribavirin therapy for chronic hepatitis C (HCV). Relevant primarily for genotype 1 and 4 HCV where interferon-based therapy is still considered.",
        "pharmvar_url": None,
        "high_pheno":     [],
        "moderate_pheno": [],
        "diplotype_check": lambda d: "Reference/Reference" not in (d or "") and bool(d) and d not in ("-", ""),
        "landing_notes": {
            "_diplotype": "Non-reference IFNL3 genotype detected — reduced likelihood of sustained virological response (SVR) with pegylated interferon alfa + ribavirin for HCV.",
        },
        "drugs": [
            {"name": "Peginterferon alfa-2a", "level": "A",
             "url": "https://cpicpgx.org/guidelines/guideline-for-peg-interferon-alpha-based-regimens-and-ifnl3/",
             "rec": "Non-CC genotype (CT or TT at rs12979860): significantly reduced SVR rate with peg-IFN + RBV. Prefer direct-acting antiviral (DAA) regimens where available."},
            {"name": "Peginterferon alfa-2b + Ribavirin", "level": "A",
             "url": "https://cpicpgx.org/guidelines/guideline-for-peg-interferon-alpha-based-regimens-and-ifnl3/",
             "rec": "Non-CC genotype: lower SVR probability. DAA therapy is strongly preferred. If interferon-based therapy must be used, counsel patient on reduced efficacy."},
        ],
    },
    "NAT1": {
        "desc": "NAT1 (N-acetyltransferase 1) acetyltransferase with activity toward aromatic and heterocyclic amines. No CPIC Level A/B guideline with actionable dose recommendations currently exists.",
        "pharmvar_url": None,
        "high_pheno":     [],
        "moderate_pheno": [],
        "diplotype_check": None,
        "landing_notes": {},
        "drugs": [],
    },
    "POR": {
        "desc": "POR (cytochrome P450 oxidoreductase) is the essential electron donor for all microsomal CYP enzymes. POR *28 (A503V) is the most common variant and modestly alters the activity of multiple CYPs. No standalone CPIC guideline exists; interpret POR variants in context of co-metabolising CYP genotypes.",
        "pharmvar_url": "https://www.pharmvar.org/gene/POR",
        "high_pheno":     [],
        "moderate_pheno": [],
        "diplotype_check": None,
        "landing_notes": {},
        "drugs": [],
    },
}


def _pheno_cat(phenotype: str) -> str:
    """Normalise a phenotype string to a lowercase keyword category."""
    p = (phenotype or "").lower().replace("-", " ").replace("_", " ").strip()
    if not p or p == "-":
        return "unknown"
    # Check ultrarapid before rapid (substring order matters)
    if "ultrarapid" in p.replace(" ", "") or "ultra rapid" in p:
        return "ultrarapid"
    if "poor" in p:       return "poor"
    if "rapid" in p:      return "rapid"      # rapid but not ultrarapid (caught above)
    if "intermediate" in p: return "intermediate"
    if "slow" in p:       return "slow"
    if "deficient" in p:  return "deficient"
    if "hemizygous" in p: return "hemizygous"
    if "homozygous" in p: return "homozygous"
    if "heterozygous" in p: return "heterozygous"
    if "variable" in p:   return "variable"
    if "decreased" in p or "reduced" in p: return "decreased"
    if "increased" in p:  return "increased"
    if "susceptible" in p or "malignant" in p: return "susceptible"
    if "non expresser" in p or "nonexpresser" in p: return "poor"  # CYP3A5 *3/*3
    if "normal" in p or "reference" in p or "wild" in p: return "normal"
    return "other"


def _get_tier(gene: str, pheno_cat: str, diplotype: str,
              all_diplotypes: list[str] | None = None) -> str | None:
    """Return 'high', 'moderate', 'informational', or None (not flagged).

    For diplotype_check genes, all_diplotypes (one per tool that called the gene)
    are checked individually so that nomenclature differences between tools
    (e.g. rs9923231 vs *S1 vs *H1 for VKORC1) do not block detection.
    The finding is raised only if ≥2 tool diplotypes trigger the lambda.
    """
    entry = CPIC_DB.get(gene)
    if not entry:
        return None
    dc = entry.get("diplotype_check")
    if dc:
        diplos = all_diplotypes if all_diplotypes else ([diplotype] if diplotype else [])
        min_t = entry.get("min_tools", 2)
        n_triggered = sum(1 for d in diplos if d and d not in ("-", "") and dc(d))
        if n_triggered >= min_t:
            if gene in ("DPYD", "RYR1", "TPMT", "NUDT15", "HLA-A", "HLA-B"):
                return "high"
            elif gene in ("VKORC1", "ABCG2"):
                return "moderate"
            else:
                return "informational"
    # Phenotype-based tier
    if pheno_cat in ("normal", "unknown", "other"):
        return None
    if pheno_cat in entry.get("high_pheno", []):
        return "high"
    if pheno_cat in entry.get("moderate_pheno", []):
        return "moderate"
    return None


def _get_landing_note(gene: str, pheno_cat: str, diplotype: str,
                      all_diplotypes: list[str] | None = None) -> str:
    """Return the brief TL;DR clinical note for a finding."""
    entry = CPIC_DB.get(gene, {})
    notes = entry.get("landing_notes", {})
    dc = entry.get("diplotype_check")
    diplos = all_diplotypes if all_diplotypes else ([diplotype] if diplotype else [])

    # HLA genes: generate note from whichever risk alleles are present
    risk_alleles = entry.get("risk_alleles")
    if risk_alleles and dc:
        found = [note for allele, note in risk_alleles.items()
                 if any(allele in (d or "") for d in diplos if d not in ("-", ""))]
        if found:
            return " | ".join(found)
        return entry.get("no_risk_note", "")

    if dc and any(dc(d) for d in diplos if d and d not in ("-", "")) and "_diplotype" in notes:
        return notes["_diplotype"]
    return notes.get(pheno_cat, notes.get("_diplotype", ""))


def build_clinical_findings_section(genes_data: list, sample: str = "",
                                     genes_rel_prefix: str = "",
                                     embedded: bool = False) -> str:
    """Build the Key Clinical Findings TL;DR section for the landing page.

    embedded — True when gene detail panels are inlined in the same HTML file;
               links use onclick='pgxShowGene()' instead of separate file hrefs.
    """
    findings: dict[str, list] = {"high": [], "moderate": [], "informational": []}

    for gd in genes_data:
        gene             = gd["gene"]
        pheno            = gd["consensus_phenotype"]
        diplo            = gd["consensus_diplotype"]
        n_agree          = gd["n_agree"]
        n_called         = gd["n_called"]
        all_tool_diplos  = gd.get("all_tool_diplotypes", [])

        entry = CPIC_DB.get(gene, {})
        dc    = entry.get("diplotype_check")
        pheno_cat = _pheno_cat(pheno)

        if dc:
            # Diplotype-check gene: count how many tools trigger the lambda.
            # Different tools use different nomenclature for the same variant
            # (e.g. rs9923231 / *S1 / *H1 for VKORC1), so we do NOT use the
            # consensus diplotype n_agree filter here — we count lambda hits.
            min_t = entry.get("min_tools", 2)
            n_triggered = sum(1 for d in all_tool_diplos if d and d not in ("-", "") and dc(d))
            if n_triggered < min_t:
                continue
            tier = _get_tier(gene, pheno_cat, diplo, all_tool_diplos)
            note = _get_landing_note(gene, pheno_cat, diplo, all_tool_diplos)
            concordance_warn = n_triggered < n_called
        else:
            # Phenotype-based gene: require ≥2 tools to agree on diplotype.
            if n_agree < 2:
                continue
            tier = _get_tier(gene, pheno_cat, diplo)
            note = _get_landing_note(gene, pheno_cat, diplo)
            concordance_warn = n_agree < n_called

        if tier is None:
            continue
        findings[tier].append({
            "gene":             gene,
            "diplotype":        diplo,
            "phenotype":        pheno,
            "note":             note,
            "concordance_warn": concordance_warn,
            "n_agree":          n_agree,
            "n_called":         n_called,
        })

    all_findings = findings["high"] + findings["moderate"] + findings["informational"]
    if not all_findings:
        return """
    <div class="section">
        <h2>Key Clinical Findings</h2>
        <div class="cf-normal">
            <span style="font-size:1.4rem">&#10003;</span>
            <div><strong>No clinically actionable pharmacogenomic variants identified.</strong><br>
            All tested genes show normal/reference function based on concordant tool results.</div>
        </div>
    </div>"""

    def render_tier(tier_findings, tier_css, tier_icon, tier_label):
        if not tier_findings:
            return ""
        html = f'<div class="cf-tier-label {tier_css}-label">{tier_icon}&nbsp; {tier_label}</div>'
        for f in tier_findings:
            pill_color = phenotype_color(f["phenotype"])
            pheno_str  = f["phenotype"] if f["phenotype"] not in ("-", "") else "—"
            diplo_str  = f["diplotype"] if f["diplotype"] not in ("-", "") else "—"
            warn_html  = (
                f'<span class="cf-concordance-warn" title="{f["n_agree"]}/{f["n_called"]} tools agree">'
                f'&#9888; Discordant ({f["n_agree"]}/{f["n_called"]})</span>'
                if f["concordance_warn"] else ""
            )
            if embedded:
                detail_href  = "#"
                detail_click = f" onclick=\"pgxShowGene('{f['gene']}'); return false;\""
            elif sample and genes_rel_prefix:
                detail_href  = f"{genes_rel_prefix}/{f['gene']}/{sample}.{f['gene']}.pgx.html"
                detail_click = ""
            elif sample:
                detail_href  = f"{sample}.{f['gene']}.pgx.html"
                detail_click = ""
            else:
                detail_href  = "#"
                detail_click = ""
            html += f"""
            <a href="{detail_href}"{detail_click} class="cf-finding-link">
            <div class="cf-finding {tier_css}">
                <div class="cf-finding-header">
                    <span class="cf-gene">{f["gene"]}</span>
                    <span class="cf-diplo">{diplo_str}</span>
                    <span class="gene-pheno-pill" style="background:{pill_color};font-size:0.72rem;padding:0.15rem 0.5rem">{pheno_str}</span>
                    {warn_html}
                    <span class="cf-detail-arrow">&#8594; Detail</span>
                </div>
                <div class="cf-note">{f["note"]}</div>
            </div>
            </a>"""
        return html

    body  = render_tier(findings["high"],         "cf-high",     "&#9888;", "High Priority")
    body += render_tier(findings["moderate"],      "cf-moderate", "&#9679;", "Moderate")
    body += render_tier(findings["informational"], "cf-info",     "&#8505;", "Informational")

    return f"""
    <div class="section">
        <h2>Key Clinical Findings</h2>
        <p class="cf-disclaimer">Based on concordant results (&#8805;2/4 tools). For research and clinical decision support only — not a standalone diagnostic.
        Consult CPIC guidelines and a clinical pharmacist before prescribing.</p>
        <div class="cf-container">
{body}
        </div>
    </div>"""


def build_gene_cpic_section(gene: str, phenotype: str, diplotype: str,
                             all_tool_diplos: list[str] | None = None,
                             id_prefix: str = "") -> str:
    """Build the comprehensive CPIC reference section for a gene detail page."""
    entry = CPIC_DB.get(gene)
    if not entry:
        return ""

    pheno_cat = _pheno_cat(phenotype)
    tier      = _get_tier(gene, pheno_cat, diplotype, all_tool_diplos)
    note      = _get_landing_note(gene, pheno_cat, diplotype, all_tool_diplos)

    tier_badge = {
        "high":          '<span class="badge badge-red" style="font-size:0.85rem;padding:0.25rem 0.75rem">High priority</span>',
        "moderate":      '<span class="badge badge-amber" style="font-size:0.85rem;padding:0.25rem 0.75rem">Moderate</span>',
        "informational": '<span class="badge badge-blue" style="font-size:0.85rem;padding:0.25rem 0.75rem">Informational</span>',
    }.get(tier or "", '<span class="badge badge-grey" style="font-size:0.85rem;padding:0.25rem 0.75rem">Normal / Not flagged</span>')

    patient_note_html = (
        f'<div class="cpic-patient-note">&#128203; {note}</div>' if note else ""
    )

    drug_rows = ""
    for d in entry.get("drugs", []):
        lvl_cls = f'cpic-level-{d["level"].lower()}'
        drug_rows += f"""
            <tr>
                <td><a href="{d['url']}" target="_blank" class="cpic-link">{d['name']}</a></td>
                <td><span class="cpic-level-badge {lvl_cls}">{d['level']}</span></td>
                <td>{d['rec']}</td>
            </tr>"""

    pharmvar_url = entry.get("pharmvar_url")
    pharmvar_html = (
        f'<a href="{pharmvar_url}" target="_blank" class="cpic-link">PharmVar — {gene} allele nomenclature &#8599;</a>'
        if pharmvar_url else ""
    )
    sep = "&ensp;|&ensp;" if pharmvar_html else ""

    return f"""
    <div class="cpic-section" id="{id_prefix}cpic-reference">
        <div class="cpic-header">
            <div>
                <h2>CPIC Clinical Reference — {gene}</h2>
                <p class="cpic-desc">{entry.get("desc", "")}</p>
            </div>
            <div style="text-align:right;flex-shrink:0">
                <div class="qc-label">This patient</div>
                <div style="margin-top:0.35rem">{tier_badge}</div>
            </div>
        </div>
        {patient_note_html}
        <div class="cpic-drug-wrap">
            <table class="cpic-drug-table">
                <thead><tr>
                    <th>Drug / Drug class</th>
                    <th>CPIC Level</th>
                    <th>Recommendation</th>
                </tr></thead>
                <tbody>{drug_rows}</tbody>
            </table>
        </div>
        <div class="cpic-footer">
            <a href="https://cpicpgx.org" target="_blank" class="cpic-link">cpicpgx.org &#8599;</a>
            {sep}{pharmvar_html}
        </div>
    </div>"""


def concordance_color(n_agree: int, n_called: int) -> tuple[str, str]:
    """Return (card_class, badge_text) based on tool agreement."""
    if n_called == 0:
        return ("card-no-data",  "No data")
    frac = n_agree / n_called
    if frac == 1.0:
        return ("card-green",  f"{n_agree}/{n_called}")
    elif frac >= 0.75:
        return ("card-amber",  f"{n_agree}/{n_called}")
    elif frac >= 0.5:
        return ("card-orange", f"{n_agree}/{n_called}")
    else:
        return ("card-red",    f"{n_agree}/{n_called}")


def normalize_diplotype(d: str) -> str:
    """Sort alleles alphabetically so *2/*4 == *4/*2."""
    if not d or d == "-":
        return "-"
    sep = "/" if "/" in d else None
    if sep is None:
        return d
    parts = [p.strip() for p in d.split(sep)]
    return "/".join(sorted(parts))


def compute_concordance(tool_data: dict) -> tuple[str, str, str, int]:
    """
    Returns (consensus_diplotype, consensus_phenotype, card_class, n_agree).
    """
    diplotypes = []
    phenotypes = []
    for tool in TOOLS:
        td = tool_data.get(tool, {})
        d = normalize_diplotype(td.get("diplotype", "-"))
        p = td.get("phenotype", "-")
        if d and d != "-":
            diplotypes.append(d)
        if p and p != "-":
            phenotypes.append(p)

    n_called = len(diplotypes)
    if n_called == 0:
        return ("-", "-", "card-no-data", 0)

    most_common_dip, n_agree = Counter(diplotypes).most_common(1)[0]
    most_common_pheno = Counter(phenotypes).most_common(1)[0][0] if phenotypes else "-"

    card_class, _ = concordance_color(n_agree, n_called)
    return (most_common_dip, most_common_pheno, card_class, n_agree)


# ── CSS / JS shared across all pages ──────────────────────────────────────────
SHARED_CSS = """
:root {
    --primary: #1a365d;
    --primary-light: #2b6cb0;
    --bg: #f7fafc;
    --card-bg: #ffffff;
    --border: #e2e8f0;
    --text: #2d3748;
    --muted: #718096;
    --green: #276749;
    --green-bg: #c6f6d5;
    --amber: #744210;
    --amber-bg: #fefcbf;
    --orange: #7b341e;
    --orange-bg: #feebc8;
    --red: #742a2a;
    --red-bg: #fed7d7;
    --blue: #2a4365;
    --blue-bg: #bee3f8;
    --grey: #4a5568;
    --grey-bg: #edf2f7;
}
* { box-sizing: border-box; margin: 0; padding: 0; }
body {
    font-family: 'Segoe UI', system-ui, -apple-system, sans-serif;
    background: var(--bg);
    color: var(--text);
    line-height: 1.5;
}
header {
    background: var(--primary);
    color: white;
    padding: 1.25rem 2rem;
    display: flex;
    align-items: center;
    gap: 1rem;
    box-shadow: 0 2px 8px rgba(0,0,0,0.3);
}
header .logo {
    font-size: 1.5rem;
    font-weight: 700;
    letter-spacing: -0.5px;
}
header .subtitle {
    font-size: 0.85rem;
    opacity: 0.75;
}
header .spacer { flex: 1; }
header .report-date { font-size: 0.8rem; opacity: 0.7; text-align: right; }
.container { max-width: 1200px; margin: 0 auto; padding: 2rem 1.5rem; }
h2 { font-size: 1.2rem; font-weight: 600; color: var(--primary); margin-bottom: 1rem; }
h3 { font-size: 1rem; font-weight: 600; color: var(--primary); margin-bottom: 0.75rem; }
.section { margin-bottom: 2.5rem; }
.badge {
    display: inline-block;
    padding: 0.15rem 0.6rem;
    border-radius: 999px;
    font-size: 0.75rem;
    font-weight: 600;
}
.badge-green  { background: var(--green-bg);  color: var(--green);  }
.badge-amber  { background: var(--amber-bg);  color: var(--amber);  }
.badge-orange { background: var(--orange-bg); color: var(--orange); }
.badge-red    { background: var(--red-bg);    color: var(--red);    }
.badge-grey   { background: var(--grey-bg);   color: var(--grey);   }
.badge-blue   { background: var(--blue-bg);   color: var(--blue);   }
table { width: 100%; border-collapse: collapse; font-size: 0.875rem; }
th {
    background: var(--primary);
    color: white;
    padding: 0.6rem 0.75rem;
    text-align: left;
    font-weight: 600;
    font-size: 0.8rem;
    text-transform: uppercase;
    letter-spacing: 0.05em;
}
td { padding: 0.5rem 0.75rem; border-bottom: 1px solid var(--border); vertical-align: top; }
tr:hover td { background: #f0f4f8; }
tr:last-child td { border-bottom: none; }
.dash { color: var(--muted); }
footer {
    text-align: center;
    padding: 2rem;
    font-size: 0.78rem;
    color: var(--muted);
    border-top: 1px solid var(--border);
    margin-top: 3rem;
}
.back-link {
    display: inline-flex;
    align-items: center;
    gap: 0.4rem;
    color: var(--primary-light);
    text-decoration: none;
    font-size: 0.875rem;
    font-weight: 500;
    margin-bottom: 1.5rem;
}
.back-link:hover { text-decoration: underline; }
"""

# ── Landing page ──────────────────────────────────────────────────────────────
LANDING_EXTRA_CSS = """
.sample-banner {
    background: white;
    border: 1px solid var(--border);
    border-radius: 8px;
    padding: 1.25rem 1.5rem;
    margin-bottom: 2rem;
    display: flex;
    align-items: flex-start;
    gap: 2rem;
    flex-wrap: wrap;
}
.sample-id { font-size: 1.5rem; font-weight: 700; color: var(--primary); }
.sample-meta { font-size: 0.85rem; color: var(--muted); margin-top: 0.2rem; }
.qc-grid {
    display: grid;
    grid-template-columns: repeat(auto-fill, minmax(170px, 1fr));
    gap: 0.75rem;
    margin-bottom: 2rem;
}
.qc-card {
    background: white;
    border: 1px solid var(--border);
    border-radius: 8px;
    padding: 0.9rem 1rem;
}
.qc-label { font-size: 0.72rem; text-transform: uppercase; letter-spacing: 0.05em; color: var(--muted); }
.qc-value { font-size: 1.25rem; font-weight: 700; color: var(--primary); margin-top: 0.2rem; }
.qc-unit  { font-size: 0.78rem; color: var(--muted); }
.gene-grid {
    display: grid;
    grid-template-columns: repeat(auto-fill, minmax(200px, 1fr));
    gap: 1rem;
}
.gene-card {
    background: white;
    border-radius: 8px;
    padding: 1rem 1.1rem;
    cursor: pointer;
    text-decoration: none;
    color: var(--text);
    transition: transform 0.12s, box-shadow 0.12s;
    border-left: 5px solid #ccc;
    display: block;
}
.gene-card:hover { transform: translateY(-2px); box-shadow: 0 4px 16px rgba(0,0,0,0.12); }
.card-green  { border-left-color: #38a169; }
.card-amber  { border-left-color: #d69e2e; }
.card-orange { border-left-color: #ed8936; }
.card-red    { border-left-color: #e53e3e; }
.card-no-data{ border-left-color: #a0aec0; }
.gene-name   { font-size: 1rem; font-weight: 700; }
.gene-diplo  { font-size: 0.8rem; color: var(--muted); margin-top: 0.25rem; white-space: nowrap; overflow: hidden; text-overflow: ellipsis; }
.gene-pheno-pill {
    display: inline-block;
    margin-top: 0.5rem;
    padding: 0.2rem 0.55rem;
    border-radius: 999px;
    font-size: 0.7rem;
    font-weight: 600;
    color: white;
}
.gene-concord { font-size: 0.7rem; margin-top: 0.4rem; }
.gene-depth   { font-size: 0.7rem; margin-top: 0.25rem; color: var(--muted); }
/* ── Sequencing depth coverage flags ── */
.depth-ok       { color: #2d7d46; font-weight: 700; }   /* ≥80% at ≥30× */
.depth-caution  { color: #b38600; font-weight: 700; }   /* 50–79% */
.depth-poor     { color: #d46b08; font-weight: 700; }   /* 10–49% */
.depth-critical { color: #c0392b; font-weight: 700; }   /* <10%  */
/* ── Key Clinical Findings section ── */
.cf-disclaimer { font-size: 0.75rem; color: var(--muted); margin-bottom: 1rem; font-style: italic; }
.cf-container { display: flex; flex-direction: column; gap: 0.4rem; }
.cf-tier-label { font-size: 0.72rem; font-weight: 700; text-transform: uppercase; letter-spacing: 0.08em; padding: 0.6rem 0 0.2rem; }
.cf-high-label     { color: #742a2a; }
.cf-moderate-label { color: #744210; }
.cf-info-label     { color: #2a4365; }
.cf-finding-link { display: block; text-decoration: none; color: inherit; border-radius: 8px; transition: transform 0.12s ease, box-shadow 0.12s ease; }
.cf-finding-link:hover { transform: translateY(-2px); box-shadow: 0 4px 12px rgba(0,0,0,0.12); }
.cf-finding-link:hover .cf-detail-arrow { opacity: 1; }
.cf-detail-arrow { font-size: 0.72rem; color: var(--primary-light); opacity: 0; transition: opacity 0.12s ease; margin-left: auto; white-space: nowrap; }
.cf-finding { background: white; border-radius: 8px; padding: 0.85rem 1.1rem; border-left: 5px solid #ccc; }
.cf-high     { border-left-color: #e53e3e; }
.cf-moderate { border-left-color: #d69e2e; }
.cf-info     { border-left-color: #3182ce; }
.cf-finding-header { display: flex; align-items: center; gap: 0.65rem; flex-wrap: wrap; margin-bottom: 0.35rem; }
.cf-gene  { font-weight: 700; font-size: 0.95rem; color: var(--primary); min-width: 70px; }
.cf-diplo { font-size: 0.8rem; color: var(--muted); font-family: monospace; }
.cf-note  { font-size: 0.82rem; color: var(--text); line-height: 1.5; }
.cf-concordance-warn { font-size: 0.72rem; color: #744210; background: #fefcbf; border-radius: 4px; padding: 0.1rem 0.4rem; }
.cf-normal { background: #c6f6d5; color: #276749; border-radius: 8px; padding: 1rem 1.25rem; font-size: 0.88rem; display: flex; align-items: center; gap: 0.75rem; }
.legend-row {
    display: flex;
    gap: 1.2rem;
    flex-wrap: wrap;
    font-size: 0.78rem;
    margin-bottom: 1.25rem;
    align-items: center;
}
.legend-item { display: flex; align-items: center; gap: 0.4rem; }
.legend-swatch { width: 14px; height: 14px; border-radius: 3px; }
"""

def fmt_num(n, suffix=""):
    """Format large numbers with comma separators."""
    try:
        return f"{int(n):,}{suffix}"
    except Exception:
        return str(n)

def bam_stats_cards(bs: dict) -> str:
    items = [
        ("Total reads",    fmt_num(bs.get("total_reads", 0)),       ""),
        ("Mapped",         f"{bs.get('mapped_pct','—')}",            "%"),
        ("Duplicates",     f"{bs.get('duplicate_pct','—')}",         "%"),
        ("Mean depth",     f"{bs.get('mean_depth_genome','—')}",     "X"),
        ("Read length",    f"{bs.get('read_length','—')}",           " bp"),
        ("Insert size",    f"{bs.get('insert_size_mean','—')}",      " bp"),
        ("MAPQ≥20",        f"{bs.get('mapq20_pct','—')}",            "%"),
        ("Inferred sex",   f"{bs.get('inferred_sex','—')}",          ""),
        ("X/Y ratio",      f"{bs.get('xy_depth_ratio','—')}",        ""),
        ("Error rate",     f"{bs.get('error_rate','—')}",            ""),
    ]
    cards = ""
    for label, val, unit in items:
        cards += f"""
        <div class="qc-card">
            <div class="qc-label">{label}</div>
            <div class="qc-value">{val}<span class="qc-unit">{unit}</span></div>
        </div>"""
    return cards


def _depth_css(pct_ge_30x) -> str:
    """Return CSS class for a ≥30× coverage fraction value.

    Thresholds align with CAP/AMP PGx guidance and tool recommendations:
      ≥80%  → depth-ok       (green)  — reliable calling
      50–79% → depth-caution (amber)  — SNV OK, SV uncertain
      10–49% → depth-poor    (orange) — results may be unreliable
       <10%  → depth-critical (red)   — insufficient for calling
    """
    try:
        v = float(pct_ge_30x)
    except (TypeError, ValueError):
        return ""
    if v >= 80:
        return "depth-ok"
    if v >= 50:
        return "depth-caution"
    if v >= 10:
        return "depth-poor"
    return "depth-critical"


def gene_depth_table(bs: dict) -> str:
    gd = bs.get("gene_depth", {})
    if not gd:
        return ""
    rows = ""
    for gene in sorted(gd.keys()):
        d = gd[gene]
        if d.get("note") == "alt_contig":
            rows += f"""
        <tr>
            <td>{gene}</td>
            <td colspan="3" style="color:var(--muted);font-style:italic">Alt contig (chr22_KI270879v1_alt) — not measurable on primary assembly</td>
        </tr>"""
        else:
            mean  = d.get('mean', '—')
            p20   = d.get('pct_ge_20x', '—')
            p30   = d.get('pct_ge_30x', '—')
            p30_cls = _depth_css(p30)
            p30_cell = f'<span class="{p30_cls}">{p30}%</span>' if p30_cls else f'{p30}%'
            rows += f"""
        <tr>
            <td>{gene}</td>
            <td>{mean}X</td>
            <td>{p20}%</td>
            <td>{p30_cell}</td>
        </tr>"""
    return f"""
    <div class="section">
        <h2>Per-gene Sequencing Depth</h2>
        <table>
            <thead><tr>
                <th>Gene</th><th>Mean depth</th><th>≥20X</th><th>≥30X</th>
            </tr></thead>
            <tbody>{rows}</tbody>
        </table>
    </div>"""


def build_landing(sample: str, bam: str, genes_data: list, bs: dict | None, out_dir: str,
                  genes_rel_prefix: str = "",
                  gene_fragments: dict | None = None):
    """Build <sample>_pgx_report.html landing page with embedded gene detail panels."""

    gene_cards_html = ""
    gene_depth_map = (bs or {}).get("gene_depth", {})
    for gd in genes_data:
        gene = gd["gene"]
        dip  = gd["consensus_diplotype"]
        pheno = gd["consensus_phenotype"]
        card_class = gd["card_class"]
        n_agree = gd["n_agree"]
        n_called = gd["n_called"]
        badge_text = f"{n_agree}/{n_called} tools" if n_called else "No data"
        pill_color = phenotype_color(pheno)
        pheno_short = pheno if pheno != "-" else "—"

        card_class_map = {
            "card-green":   "card-green",
            "card-amber":   "card-amber",
            "card-orange":  "card-orange",
            "card-red":     "card-red",
            "card-no-data": "card-no-data",
        }

        badge_cls_map = {
            "card-green":   "badge-green",
            "card-amber":   "badge-amber",
            "card-orange":  "badge-orange",
            "card-red":     "badge-red",
            "card-no-data": "badge-grey",
        }
        badge_cls = badge_cls_map.get(card_class, "badge-grey")
        css_cls    = card_class_map.get(card_class, "")

        depth_info = gene_depth_map.get(gene, {})
        if depth_info:
            if depth_info.get("note") == "alt_contig":
                depth_html = '<div class="gene-depth">Depth: Alt contig</div>'
            else:
                mean  = depth_info.get("mean", "—")
                pct30 = depth_info.get("pct_ge_30x", "—")
                p30_cls = _depth_css(pct30)
                p30_str = f'<span class="{p30_cls}">{pct30}%</span>' if p30_cls else f'{pct30}%'
                depth_html = f'<div class="gene-depth">Depth: {mean}X &nbsp;|&nbsp; &#8805;30X: {p30_str}</div>'
        else:
            depth_html = ""

        if gene_fragments is not None:
            card_onclick = f'onclick="pgxShowGene(\'{gene}\'); return false;"'
            card_href = "#"
        else:
            card_onclick = ""
            if genes_rel_prefix:
                card_href = f"{genes_rel_prefix}/{gene}/{sample}.{gene}.pgx.html"
            else:
                card_href = f"{sample}.{gene}.pgx.html"

        gene_cards_html += f"""
            <a href="{card_href}" {card_onclick} class="gene-card {css_cls}">
                <div class="gene-name">{gene}</div>
                <div class="gene-diplo" title="{dip}">{dip if dip != '-' else '—'}</div>
                <div>
                    <span class="gene-pheno-pill" style="background:{pill_color}">{pheno_short}</span>
                </div>
                <div class="gene-concord">
                    Concordance: <span class="badge {badge_cls}">{badge_text}</span>
                </div>
                {depth_html}
            </a>"""

    bam_basename = os.path.basename(bam) if bam else "—"

    clinical_html = build_clinical_findings_section(
        genes_data, sample, genes_rel_prefix,
        embedded=(gene_fragments is not None))

    qc_html = ""
    gene_depth_html = ""
    if bs:
        qc_html = f"""
        <div class="section">
            <h2>BAM Quality Control</h2>
            <div class="qc-grid">{bam_stats_cards(bs)}</div>
        </div>"""
        gene_depth_html = gene_depth_table(bs)

    today = date.today().isoformat()

    # Build embedded gene panels HTML
    gene_panels_html = ""
    if gene_fragments:
        for gfragment_gene, fragment in gene_fragments.items():
            gene_panels_html += (
                f'<div id="gene-panel-{gfragment_gene}" class="gene-panel" hidden>\n'
                f'<div class="container">\n{fragment}\n</div>\n</div>\n'
            )

    # Include DETAIL_EXTRA_CSS when embedding gene panels
    extra_style = f"\n{DETAIL_EXTRA_CSS}" if gene_fragments else ""

    # JS show/hide for embedded mode
    embedded_js = ""
    if gene_fragments:
        embedded_js = """
<script>
function pgxShowGene(gene) {
    document.getElementById('main-view').hidden = true;
    document.querySelectorAll('.gene-panel').forEach(function(el) { el.hidden = true; });
    var panel = document.getElementById('gene-panel-' + gene);
    if (panel) { panel.hidden = false; }
    window.scrollTo(0, 0);
}
function pgxShowMain() {
    document.querySelectorAll('.gene-panel').forEach(function(el) { el.hidden = true; });
    document.getElementById('main-view').hidden = false;
    window.scrollTo(0, 0);
}
</script>"""

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>PGx Report &#8212; {sample}</title>
    <style>
{SHARED_CSS}
{LANDING_EXTRA_CSS}{extra_style}
    </style>
{embedded_js}
</head>
<body>
<header>
    <div>
        <div class="logo">&#9652; PGx Suite</div>
        <div class="subtitle">Pharmacogenomics Star-Allele Report</div>
    </div>
    <div class="spacer"></div>
    <div class="report-date">Report date: {today}<br>GRCh38</div>
</header>

<div id="main-view">
<div class="container">

    <div class="sample-banner">
        <div>
            <div class="qc-label">Sample ID</div>
            <div class="sample-id">{sample}</div>
            <div class="sample-meta">BAM: {bam_basename}</div>
        </div>
        <div>
            <div class="qc-label">Tools</div>
            <div style="margin-top:0.3rem; display:flex; gap:0.5rem; flex-wrap:wrap;">
                <span class="badge badge-blue">PyPGx 0.26</span>
                <span class="badge badge-blue">Stargazer 2.0.3</span>
                <span class="badge badge-blue">Aldy 4.8</span>
                <span class="badge badge-blue">StellarPGx 1.2.7</span>
                <span class="badge badge-blue">OptiType 1.3.5</span>
            </div>
        </div>
    </div>

{clinical_html}

{qc_html}

    <div class="section">
        <h2>Gene Diplotype Summary</h2>
        <div class="legend-row">
            <span style="font-weight:600;color:var(--muted)">Concordance:</span>
            <span class="legend-item">
                <span class="legend-swatch" style="background:#38a169"></span>4/4 tools agree
            </span>
            <span class="legend-item">
                <span class="legend-swatch" style="background:#d69e2e"></span>3/4 agree
            </span>
            <span class="legend-item">
                <span class="legend-swatch" style="background:#ed8936"></span>2/4 agree
            </span>
            <span class="legend-item">
                <span class="legend-swatch" style="background:#e53e3e"></span>&lt;2/4 agree
            </span>
            <span class="legend-item">
                <span class="legend-swatch" style="background:#a0aec0"></span>No data
            </span>
        </div>
        <div class="gene-grid">
{gene_cards_html}
        </div>
    </div>

{gene_depth_html}

</div>
</div>

<!-- Gene detail panels (embedded, shown/hidden by JS) -->
{gene_panels_html}

<footer>
    PGx Suite &bull; GRCh38 &bull; PyPGx 0.26 / Stargazer 2.0.3 / Aldy 4.8 / StellarPGx 1.2.7 / OptiType 1.3.5
    &bull; For research and clinical decision support only &#8212; not a standalone diagnostic.
</footer>
</body>
</html>"""

    out_path = os.path.join(out_dir, f"{sample}_pgx_report.html")
    with open(out_path, "w") as fh:
        fh.write(html)
    print(f"Landing page: {out_path}")


# ── Per-gene detail page ───────────────────────────────────────────────────────
DETAIL_EXTRA_CSS = """
.gene-header {
    background: white;
    border: 1px solid var(--border);
    border-radius: 8px;
    padding: 1.25rem 1.5rem;
    margin-bottom: 2rem;
    display: flex;
    align-items: flex-start;
    gap: 2rem;
    flex-wrap: wrap;
}
.gene-title { font-size: 1.8rem; font-weight: 700; color: var(--primary); }
.detail-table-wrap {
    background: white;
    border: 1px solid var(--border);
    border-radius: 8px;
    overflow: hidden;
}
.detail-table th:first-child,
.detail-table td:first-child {
    background: #f7fafc;
    font-weight: 600;
    color: var(--primary);
    min-width: 180px;
    border-right: 2px solid var(--border);
}
.detail-table th { font-size: 0.85rem; }
.variant-list { font-family: monospace; font-size: 0.78rem; }
.sv-note {
    background: #ebf8ff;
    border-left: 4px solid #3182ce;
    border-radius: 4px;
    padding: 0.6rem 0.9rem;
    font-size: 0.83rem;
    color: #2a4365;
    margin-top: 1.25rem;
}
/* ── Variant evidence subtable ─────────────────────────────────────── */
.variant-section {
    background: white;
    border: 1px solid var(--border);
    border-radius: 8px;
    padding: 1.25rem 1.5rem;
    margin-top: 2rem;
}
.var-meta { display: flex; align-items: center; gap: 1.5rem; flex-wrap: wrap; margin-bottom: 0.75rem; }
.var-summary { font-size: 0.82rem; color: var(--muted); }
.var-legend { display: flex; gap: 1rem; font-size: 0.75rem; align-items: center; }
.vleg-swatch { display: inline-block; width: 12px; height: 12px; border-radius: 2px; margin-right: 3px; vertical-align: middle; }
.vleg-all   { background: #c6f6d5; }
.vleg-multi { background: #fefcbf; }
.vleg-single{ background: #edf2f7; }
.var-table-wrap { overflow-x: auto; border: 1px solid var(--border); border-radius: 6px; }
.var-table { min-width: 750px; font-size: 0.78rem; border-collapse: collapse; width: 100%; }
.var-table thead th {
    background: var(--primary);
    color: white;
    padding: 0.45rem 0.65rem;
    font-size: 0.72rem;
    text-transform: uppercase;
    letter-spacing: 0.05em;
    text-align: left;
    font-weight: 600;
    white-space: nowrap;
}
.var-table thead th.vtool-col { background: #2c5282; }
.var-table td { padding: 0.38rem 0.65rem; border-bottom: 1px solid var(--border); vertical-align: top; }
.var-table tr:last-child td { border-bottom: none; }
.vrow-all   td { background: #f0fff4; }
.vrow-multi td { background: #fffff0; }
.vrow-single td { background: inherit; }
.vpos-col { font-family: 'Consolas', 'Cascadia Code', monospace; white-space: nowrap; color: #4a5568; }
.vchange-col { font-family: 'Consolas', 'Cascadia Code', monospace; font-weight: 700; white-space: nowrap; }
.vrsid-col { white-space: nowrap; }
.veff-col { font-size: 0.73rem; max-width: 200px; }
.vtool-data { font-size: 0.73rem; font-family: 'Consolas', 'Cascadia Code', monospace; line-height: 1.6; }
.vnot-detected { color: #cbd5e0; text-align: center; font-size: 0.82rem; }
.vallele-tag {
    display: inline-block;
    background: #ebf8ff; color: #2b6cb0;
    padding: 0.05rem 0.35rem;
    border-radius: 4px;
    font-size: 0.68rem;
    font-weight: 700;
    margin-bottom: 0.1rem;
    font-family: 'Segoe UI', sans-serif;
}
.vcopy-tag {
    display: inline-block;
    background: #f0e6ff; color: #553c9a;
    padding: 0.05rem 0.35rem;
    border-radius: 4px;
    font-size: 0.68rem;
    font-weight: 600;
    margin-bottom: 0.1rem;
    font-family: 'Segoe UI', sans-serif;
}
.vgt-tag {
    display: inline-block;
    background: #fff5f5; color: #742a2a;
    padding: 0.05rem 0.35rem;
    border-radius: 4px;
    font-size: 0.68rem;
    font-weight: 600;
    font-family: 'Segoe UI', sans-serif;
}
.rsid-link { color: var(--primary-light); text-decoration: none; font-size: 0.75rem; }
.rsid-link:hover { text-decoration: underline; }
.var-count-cell { font-size: 0.78rem; color: var(--muted); }
.var-anchor-link { font-size: 0.72rem; color: var(--primary-light); margin-left: 0.4rem; }
.gene-page-nav { display: flex; gap: 1rem; margin-bottom: 1.5rem; flex-wrap: wrap; }
.gene-page-nav a { font-size: 0.85rem; color: var(--primary); text-decoration: none; padding: 0.35rem 0.8rem; border: 1px solid var(--primary-light); border-radius: 6px; transition: background 0.12s ease; }
.gene-page-nav a:hover { background: var(--primary-light); color: white; }
.gene-depth-detail { font-size: 0.82rem; color: var(--muted); margin-top: 0.4rem; }
.gene-locus { font-size: 0.78rem; color: var(--muted); font-family: monospace; margin-top: 0.2rem; }
/* ── CPIC reference section (gene detail page) ── */
.cpic-section { background: white; border: 1px solid var(--border); border-radius: 8px; padding: 1.5rem; margin-top: 2rem; }
.cpic-header { display: flex; justify-content: space-between; align-items: flex-start; gap: 1.5rem; margin-bottom: 0.9rem; flex-wrap: wrap; }
.cpic-desc { font-size: 0.85rem; color: var(--muted); margin-top: 0.4rem; max-width: 720px; line-height: 1.5; }
.cpic-patient-note { background: #ebf8ff; border-left: 4px solid #3182ce; border-radius: 4px; padding: 0.6rem 0.9rem; font-size: 0.83rem; color: #2a4365; margin-bottom: 1rem; line-height: 1.5; }
.cpic-drug-wrap { overflow-x: auto; border: 1px solid var(--border); border-radius: 6px; margin-bottom: 1rem; }
.cpic-drug-table { width: 100%; border-collapse: collapse; font-size: 0.83rem; }
.cpic-drug-table th { background: var(--primary); color: white; padding: 0.5rem 0.75rem; text-align: left; font-size: 0.78rem; text-transform: uppercase; letter-spacing: 0.05em; font-weight: 600; }
.cpic-drug-table td { padding: 0.5rem 0.75rem; border-bottom: 1px solid var(--border); vertical-align: top; }
.cpic-drug-table tr:last-child td { border-bottom: none; }
.cpic-level-badge { display: inline-block; padding: 0.1rem 0.45rem; border-radius: 4px; font-weight: 700; font-size: 0.75rem; }
.cpic-level-a { background: #c6f6d5; color: #276749; }
.cpic-level-b { background: #bee3f8; color: #2a4365; }
.cpic-footer { font-size: 0.78rem; color: var(--muted); margin-top: 0.5rem; }
.cpic-link { color: var(--primary-light); text-decoration: none; }
.cpic-link:hover { text-decoration: underline; }
"""


def fmt_value(key: str, val) -> str:
    """Format a field value for HTML display."""
    if val is None or val == "" or val == "-":
        return '<span class="dash">—</span>'
    if isinstance(val, list):
        if not val:
            return '<span class="dash">—</span>'
        return "<br>".join(f'<span class="variant-list">{v}</span>' for v in val)
    s = str(val)
    if s == "-":
        return '<span class="dash">—</span>'
    return s


# ── Variant harmonization ──────────────────────────────────────────────────────

def harmonize_variants(tools_data: dict) -> list[dict]:
    """Merge supporting_variants from all tools into a position-sorted union.

    Variants within ±1 bp with the same REF and ALT are treated as the same
    event (handles VCF coordinate convention differences between tools).

    Returns a list of cluster dicts, each with:
        pos    – canonical GRCh38 position (int, lowest seen)
        ref    – reference allele (str)
        alt    – alternate allele (str)
        rsid   – best rsID across tools (str, or "-")
        effect – best functional effect label across tools (str, or "-")
        tools  – dict {tool_name: variant_dict} for tools that detected this variant
    """
    collected: list[tuple[int, str, str, str, dict]] = []
    for tool in TOOLS:
        for v in tools_data.get(tool, {}).get("supporting_variants", []):
            try:
                pos = int(v.get("pos", 0))
            except (ValueError, TypeError):
                continue
            ref = (v.get("ref") or "").upper().strip()
            alt = (v.get("alt") or "").upper().strip()
            if not ref or not alt or ref == "-" or alt == "-":
                continue
            collected.append((pos, ref, alt, tool, v))

    if not collected:
        return []

    collected.sort(key=lambda x: (x[0], x[1], x[2]))

    clusters: list[dict] = []
    for pos, ref, alt, tool, v in collected:
        # Search for a matching cluster (same ref+alt, position within ±1 bp)
        merged = False
        for cl in reversed(clusters):
            if cl["pos"] < pos - 2:
                break          # sorted — no older cluster can be within ±1
            if cl["ref"] == ref and cl["alt"] == alt and abs(cl["pos"] - pos) <= 1:
                cl["pos"] = min(cl["pos"], pos)   # keep the lower coordinate
                # Accumulate best rsid and effect from any tool
                if cl["rsid"] in ("-", "", None) and (v.get("rsid") or "-") != "-":
                    cl["rsid"] = v["rsid"]
                if cl["effect"] in ("-", "", None) and (v.get("effect") or "-") != "-":
                    cl["effect"] = v["effect"]
                # First reporter wins within a tool (avoids duplicate Aldy rows)
                cl["tools"].setdefault(tool, v)
                merged = True
                break
        if not merged:
            clusters.append({
                "pos":    pos,
                "ref":    ref,
                "alt":    alt,
                "rsid":   (v.get("rsid") or "-"),
                "effect": (v.get("effect") or "-"),
                "tools":  {tool: v},
            })

    clusters.sort(key=lambda c: c["pos"])
    return clusters


def _tool_cell(tool: str, v: dict | None) -> str:
    """Return the <td> HTML for one tool's evidence at a variant locus."""
    if v is None:
        return '<td class="vnot-detected">—</td>'
    parts: list[str] = []
    if tool == "PyPGx":
        allele = (v.get("allele") or "-")
        af     = (v.get("af")     or "-")
        if allele != "-":
            parts.append(f'<span class="vallele-tag">{allele}</span>')
        if af != "-":
            parts.append(f"AF&nbsp;{af}")
    elif tool == "Stargazer":
        af    = (v.get("af")    or "-")
        depth = (v.get("depth") or "-")
        if af    != "-": parts.append(f"AF&nbsp;{af}")
        if depth != "-": parts.append(f"Dp&nbsp;{depth}")
    elif tool == "Aldy":
        allele = (v.get("allele") or "-")
        depth  = (v.get("depth")  or "-")
        rsid   = (v.get("rsid")   or "-")
        if allele != "-":
            parts.append(f'<span class="vcopy-tag">{allele}</span>')
        if depth  != "-": parts.append(f"Dp&nbsp;{depth}")
        if rsid   != "-":
            parts.append(
                f'<a href="https://www.ncbi.nlm.nih.gov/snp/{rsid}" '
                f'target="_blank" class="rsid-link">{rsid}</a>'
            )
    elif tool == "StellarPGx":
        gt = (v.get("gt") or "-")
        if gt != "-":
            parts.append(f'<span class="vgt-tag">GT&nbsp;{gt}</span>')
    inner = "<br>".join(parts) if parts else "✓"
    return f'<td class="vtool-data">{inner}</td>'


def render_variant_subtable(clusters: list[dict], n_tools_per_gene: dict | None = None,
                            id_prefix: str = "") -> str:
    """Return the full HTML for the harmonised cross-tool variant evidence section."""
    if not clusters:
        return ""

    n_total  = len(clusters)
    n_shared = sum(1 for c in clusters if len(c["tools"]) >= 2)

    summary = (
        f'<span class="var-summary">'
        f'{n_total} unique variant{"s" if n_total != 1 else ""} '
        f'&nbsp;·&nbsp; {n_shared} detected by ≥2 tools'
        f'</span>'
    )
    legend = """
        <span class="var-legend">
            <span class="vleg-swatch vleg-all"></span>All tools&ensp;
            <span class="vleg-swatch vleg-multi"></span>≥2 tools&ensp;
            <span class="vleg-swatch vleg-single"></span>1 tool only
        </span>"""

    tool_header_cells = "".join(
        f'<th class="vtool-col">{t}</th>' for t in TOOLS
    )

    rows = ""
    for cl in clusters:
        n_det = len(cl["tools"])
        row_cls = (
            "vrow-all"    if n_det == len(TOOLS) else
            "vrow-multi"  if n_det >= 2          else
            "vrow-single"
        )
        rsid = cl.get("rsid") or "-"
        rsid_html = (
            f'<a href="https://www.ncbi.nlm.nih.gov/snp/{rsid}" '
            f'target="_blank" class="rsid-link">{rsid}</a>'
            if rsid != "-" else '<span class="dash">—</span>'
        )
        eff = cl.get("effect") or "-"
        eff_html = eff if eff != "-" else '<span class="dash">—</span>'

        tool_cells = "".join(_tool_cell(t, cl["tools"].get(t)) for t in TOOLS)

        rows += (
            f'<tr class="{row_cls}">'
            f'<td class="vpos-col">{cl["pos"]}</td>'
            f'<td class="vchange-col">{cl["ref"]}&rarr;{cl["alt"]}</td>'
            f'<td class="vrsid-col">{rsid_html}</td>'
            f'<td class="veff-col">{eff_html}</td>'
            f'{tool_cells}'
            f'</tr>\n'
        )

    return f"""
    <div class="variant-section" id="{id_prefix}variant-evidence">
        <h2>Supporting Variant Evidence</h2>
        <div class="var-meta">{summary}{legend}</div>
        <div class="var-table-wrap">
        <table class="var-table">
            <thead><tr>
                <th>Position (GRCh38)</th>
                <th>REF&rarr;ALT</th>
                <th>rsID</th>
                <th>Functional effect</th>
                {tool_header_cells}
            </tr></thead>
            <tbody>
{rows}            </tbody>
        </table>
        </div>
    </div>"""


def _build_gene_inner(sample: str, gene: str, detail: dict, gene_depth: dict | None,
                      back_href: str, id_prefix: str = "") -> str:
    """Return the inner HTML fragment for a gene detail panel (no html/head/body wrapper).

    id_prefix — prepended to all section IDs (e.g. 'CYP2D6-') to avoid collisions
                when multiple gene panels are embedded in the same document.
    back_href — href for the <- Back link; use 'javascript:void(0)' for embedded mode.
    """
    tools_data = detail.get("tools", {})
    sv_note = detail.get("sv_mode", "")

    clusters = harmonize_variants(tools_data)
    var_subtable_html = render_variant_subtable(clusters, id_prefix=id_prefix)

    from collections import Counter as _Counter
    pheno_list = [tools_data.get(t, {}).get("phenotype", "-") for t in TOOLS
                  if tools_data.get(t, {}).get("phenotype", "-") not in ("-", "")]
    consensus_pheno = _Counter(pheno_list).most_common(1)[0][0] if pheno_list else ""
    diplo_list = [normalize_diplotype(tools_data.get(t, {}).get("diplotype", "-")) for t in TOOLS
                  if normalize_diplotype(tools_data.get(t, {}).get("diplotype", "-")) != "-"]
    consensus_diplo = _Counter(diplo_list).most_common(1)[0][0] if diplo_list else ""
    all_tool_diplos_page = [tools_data.get(t, {}).get("diplotype", "-") for t in TOOLS
                            if tools_data.get(t, {}).get("diplotype", "-") not in ("-", "")]
    cpic_section_html = build_gene_cpic_section(
        gene, consensus_pheno, consensus_diplo, all_tool_diplos_page, id_prefix=id_prefix)

    header_cells = "".join(f"<th>{t}</th>" for t in TOOLS)
    rows_html = ""
    for field_key, field_label in FIELDS:
        if field_key == "sv_mode":
            val_html = f'<td colspan="{len(TOOLS)}" style="color:var(--muted);font-style:italic">{fmt_value("sv_mode", sv_note)}</td>'
            rows_html += f"<tr><td>{field_label}</td>{val_html}</tr>\n"
        elif field_key == "supporting_variants":
            cells = ""
            for tool in TOOLS:
                vlist = tools_data.get(tool, {}).get("supporting_variants", [])
                n = len(vlist) if isinstance(vlist, list) else 0
                cells += f'<td><span class="var-count-cell">{n} variant{"s" if n != 1 else ""}</span></td>'
            anchor = f'<a href="#{id_prefix}variant-evidence" class="var-anchor-link">&#8595; see table below</a>'
            rows_html += f"<tr><td>{field_label}{anchor}</td>{cells}</tr>\n"
        else:
            cells = ""
            for tool in TOOLS:
                td = tools_data.get(tool, {})
                raw = td.get(field_key, "-")
                cells += f"<td>{fmt_value(field_key, raw)}</td>"
            rows_html += f"<tr><td>{field_label}</td>{cells}</tr>\n"

    diplotypes = []
    for tool in TOOLS:
        d = normalize_diplotype(tools_data.get(tool, {}).get("diplotype", "-"))
        if d and d != "-":
            diplotypes.append(d)
    n_called = len(diplotypes)
    n_agree = Counter(diplotypes).most_common(1)[0][1] if diplotypes else 0
    card_class, badge_text = concordance_color(n_agree, n_called)
    badge_cls_map = {
        "card-green":   "badge-green",
        "card-amber":   "badge-amber",
        "card-orange":  "badge-orange",
        "card-red":     "badge-red",
        "card-no-data": "badge-grey",
    }
    badge_cls = badge_cls_map.get(card_class, "badge-grey")

    if gene_depth and gene_depth.get("note") == "alt_contig":
        depth_detail_html = '<div class="gene-depth-detail">Depth: Alt contig (chr22_KI270879v1_alt)</div>'
    elif gene_depth:
        mean  = gene_depth.get("mean", "—")
        pct30 = gene_depth.get("pct_ge_30x", "—")
        pct20 = gene_depth.get("pct_ge_20x", "—")
        p30_cls = _depth_css(pct30)
        p30_str = (f'<span class="{p30_cls}">{pct30}%</span>'
                   if p30_cls else f'<strong>{pct30}%</strong>')
        depth_detail_html = (
            f'<div class="gene-depth-detail">'
            f'Mean depth: <strong>{mean}X</strong>'
            f' &nbsp;|&nbsp; &#8805;20X: <strong>{pct20}%</strong>'
            f' &nbsp;|&nbsp; &#8805;30X: {p30_str}'
            f'</div>'
        )
    else:
        depth_detail_html = ""

    locus = GENE_LOCI.get(gene, "")
    locus_html = (
        f'<div class="gene-locus">&#128205; {locus} (GRCh38)</div>'
        if locus else ""
    )

    sv_note_html = f'<div class="sv-note">&#128202; {sv_note}</div>' if sv_note and sv_note != "-" else ""

    # back_onclick is set when back_href signals embedded mode
    if back_href.startswith("javascript:"):
        back_el = f'<a href="#" onclick="pgxShowMain(); return false;" class="back-link">&#8592; Back to sample summary</a>'
    else:
        back_el = f'<a href="{back_href}" class="back-link">&#8592; Back to sample summary</a>'

    return f"""
    {back_el}

    <nav class="gene-page-nav">
        <a href="#{id_prefix}tool-results">Tool Results</a>
        <a href="#{id_prefix}variant-evidence">Variant Evidence</a>
        <a href="#{id_prefix}cpic-reference">CPIC Reference</a>
    </nav>

    <div class="gene-header">
        <div>
            <div class="qc-label">Gene</div>
            <div class="gene-title">{gene}</div>
            {locus_html}
            <div class="sample-meta" style="margin-top:0.3rem">Sample: {sample} &bull; GRCh38</div>
        </div>
        <div>
            <div class="qc-label">Tool concordance</div>
            <div style="margin-top:0.4rem">
                <span class="badge {badge_cls}" style="font-size:0.9rem;padding:0.3rem 0.9rem">{badge_text}</span>
            </div>
            {depth_detail_html}
        </div>
    </div>

    <div class="section" id="{id_prefix}tool-results">
        <h2>Tool Results &#8212; {gene}</h2>
        <div class="detail-table-wrap">
            <table class="detail-table">
                <thead>
                    <tr>
                        <th>Field</th>
                        {header_cells}
                    </tr>
                </thead>
                <tbody>
{rows_html}                </tbody>
            </table>
        </div>
{sv_note_html}
    </div>

{var_subtable_html}

{cpic_section_html}
"""


def build_gene_page(sample: str, gene: str, detail: dict, landing_file: str, out_dir: str,
                    gene_depth: dict | None = None):
    """Build <sample>.<gene>.pgx.html detail page (standalone file, not used in default pipeline)."""
    today = date.today().isoformat()
    inner = _build_gene_inner(sample, gene, detail, gene_depth,
                              back_href=landing_file, id_prefix="")

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>PGx Report &#8212; {sample} &#8212; {gene}</title>
    <style>
{SHARED_CSS}
{DETAIL_EXTRA_CSS}
    </style>
</head>
<body>
<header>
    <div>
        <div class="logo">&#9652; PGx Suite</div>
        <div class="subtitle">Pharmacogenomics Star-Allele Report</div>
    </div>
    <div class="spacer"></div>
    <div class="report-date">Report date: {today}<br>GRCh38</div>
</header>

<div class="container">
{inner}
</div>

<footer>
    PGx Suite &bull; GRCh38 &bull; PyPGx 0.26 / Stargazer 2.0.3 / Aldy 4.8 / StellarPGx 1.2.7 / OptiType 1.3.5
    &bull; For research and clinical decision support only &#8212; not a standalone diagnostic.
</footer>
</body>
</html>"""

    out_path = os.path.join(out_dir, f"{sample}.{gene}.pgx.html")
    with open(out_path, "w") as fh:
        fh.write(html)
    print(f"  {gene}: {out_path}")



# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser(description="Generate PGx HTML reports")
    ap.add_argument("--sample",    required=True, help="Sample ID")
    ap.add_argument("--output",    required=True, help="Root output directory (landing page written here)")
    ap.add_argument("--genes-dir", default="",    help="Directory containing per-gene subdirs (default: <output>/Genes)")
    ap.add_argument("--bam",       default="",    help="BAM file path (for display)")
    ap.add_argument("--bam-stats", default="",    help="Path to bam_stats.json (default: <output>/log/bam_stats.json)")
    args = ap.parse_args()

    sample  = args.sample
    out_dir = args.output
    genes_dir = args.genes_dir or os.path.join(out_dir, "Genes")
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(genes_dir, exist_ok=True)

    # Load BAM stats
    bs = None
    bam_stats_path = args.bam_stats or os.path.join(out_dir, "log", "bam_stats.json")
    if os.path.isfile(bam_stats_path):
        with open(bam_stats_path) as fh:
            bs = json.load(fh)
    else:
        print(f"[warn] bam_stats.json not found at {bam_stats_path}; skipping BAM QC section",
              file=sys.stderr)

    bam_path = args.bam or (bs.get("bam", "") if bs else "")
    gene_depth_map = (bs or {}).get("gene_depth", {})

    # Discover genes from summary TSV (written to <output>/log/ by pgx-all-genes.sh)
    summary_tsv = os.path.join(out_dir, "log", "all_genes_summary.tsv")
    if not os.path.isfile(summary_tsv):
        # Fallback: legacy location at output root
        summary_tsv = os.path.join(out_dir, "all_genes_summary.tsv")
    if not os.path.isfile(summary_tsv):
        print(f"ERROR: all_genes_summary.tsv not found in {out_dir}/log/ or {out_dir}/", file=sys.stderr)
        sys.exit(1)

    seen_genes: list[str] = []
    with open(summary_tsv) as fh:
        for line in fh:
            if line.startswith("Gene"):
                continue
            parts = line.rstrip("\n").split("\t")
            if parts and parts[0] and parts[0] not in seen_genes:
                seen_genes.append(parts[0])

    # Build gene fragments and collect landing data
    print(f"Collecting gene detail panels for {len(seen_genes)} genes …")
    genes_data = []
    gene_fragments: dict = {}

    for gene in seen_genes:
        # Find detail JSON — new layout: genes_dir/<gene>/; fallback: out_dir/<gene>/
        detail_json = os.path.join(genes_dir, gene, f"{gene}_{sample}_detail.json")
        if not os.path.isfile(detail_json):
            matches = glob.glob(os.path.join(genes_dir, gene, "*_detail.json"))
            if not matches:
                # Legacy fallback
                matches = glob.glob(os.path.join(out_dir, gene, "*_detail.json"))
            detail_json = matches[0] if matches else None

        if detail_json and os.path.isfile(detail_json):
            with open(detail_json) as fh:
                detail = json.load(fh)
        else:
            print(f"  [warn] no detail JSON for {gene}", file=sys.stderr)
            detail = {"gene": gene, "sample": sample, "tools": {}}

        tools_data = detail.get("tools", {})
        consensus_dip, consensus_pheno, card_class, n_agree = compute_concordance(tools_data)
        n_called = sum(
            1 for t in TOOLS
            if normalize_diplotype(tools_data.get(t, {}).get("diplotype", "-")) != "-"
        )

        # Collect raw per-tool diplotypes (unnormalised) for diplotype_check genes
        all_tool_diplotypes = [
            tools_data.get(t, {}).get("diplotype", "-")
            for t in TOOLS
            if tools_data.get(t, {}).get("diplotype", "-") not in ("-", "")
        ]

        genes_data.append({
            "gene": gene,
            "consensus_diplotype": consensus_dip,
            "consensus_phenotype": consensus_pheno,
            "card_class": card_class,
            "n_agree": n_agree,
            "n_called": n_called,
            "all_tool_diplotypes": all_tool_diplotypes,
        })

        fragment = _build_gene_inner(
            sample, gene, detail, gene_depth_map.get(gene),
            back_href="javascript:void(0)", id_prefix=f"{gene}-")
        gene_fragments[gene] = fragment
        print(f"  {gene}: fragment built")

    # Build landing page at <out_dir>/<sample>_pgx_report.html
    print(f"Generating standalone HTML report …")
    build_landing(sample, bam_path, genes_data, bs, out_dir,
                  genes_rel_prefix="Genes",
                  gene_fragments=gene_fragments)
    print("Done.")


if __name__ == "__main__":
    main()
