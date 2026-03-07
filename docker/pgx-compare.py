#!/usr/bin/env python3
"""
pgx-compare.py — Parse outputs from all 4 PGx callers and produce a comparison table.

Usage: pgx-compare.py --gene GENE --sample SAMPLE --output-dir DIR
"""
import argparse
import csv
import glob
import io
import os
import re
import sys
import zipfile
from collections import Counter
from dataclasses import dataclass, field


@dataclass
class CallerResult:
    tool: str
    diplotype: str = "-"
    activity_score: str = "-"
    phenotype: str = "-"
    status: str = "not_run"  # not_run | ok | failed


# ── SV gene sets ─────────────────────────────────────────────────────────────
# Genes where structural variant / copy number detection requires extra steps.
# Used only for annotating the output; actual SV logic is in pgx-run.sh.
PYPGX_SV_GENES = frozenset(
    {"CYP2A6", "CYP2B6", "CYP2D6", "CYP2E1", "CYP4F2", "G6PD", "GSTM1", "GSTT1"}
)
# Stargazer paralog genes: GDF depth profile created from BAM for CN normalisation.
STARGAZER_SV_GENES = frozenset({"CYP2A6", "CYP2B6", "CYP2D6"})

# ── Gene support matrix ───────────────────────────────────────────────────────
GENE_SUPPORT: dict[str, dict[str, bool]] = {
    "CYP2D6":   {"pypgx": True,  "stargazer": True,  "aldy": True,  "stellarpgx": True},
    "CYP2C19":  {"pypgx": True,  "stargazer": True,  "aldy": True,  "stellarpgx": True},
    "CYP2C9":   {"pypgx": True,  "stargazer": True,  "aldy": True,  "stellarpgx": True},
    "CYP2B6":   {"pypgx": True,  "stargazer": True,  "aldy": True,  "stellarpgx": True},
    "CYP2C8":   {"pypgx": True,  "stargazer": True,  "aldy": True,  "stellarpgx": True},
    "CYP3A4":   {"pypgx": True,  "stargazer": True,  "aldy": True,  "stellarpgx": True},
    "CYP3A5":   {"pypgx": True,  "stargazer": True,  "aldy": True,  "stellarpgx": True},
    "CYP4F2":   {"pypgx": True,  "stargazer": True,  "aldy": True,  "stellarpgx": True},
    "NUDT15":   {"pypgx": True,  "stargazer": False, "aldy": True,  "stellarpgx": True},
    "TPMT":     {"pypgx": True,  "stargazer": False, "aldy": True,  "stellarpgx": True},
    "UGT1A1":   {"pypgx": True,  "stargazer": False, "aldy": True,  "stellarpgx": True},
    "SLCO1B1":  {"pypgx": True,  "stargazer": True,  "aldy": True,  "stellarpgx": True},
    "DPYD":     {"pypgx": True,  "stargazer": False, "aldy": True,  "stellarpgx": False},
    "NAT1":     {"pypgx": True,  "stargazer": False, "aldy": False, "stellarpgx": True},
    "NAT2":     {"pypgx": True,  "stargazer": False, "aldy": False, "stellarpgx": True},
    "G6PD":     {"pypgx": True,  "stargazer": True,  "aldy": False, "stellarpgx": False},
    "GSTM1":    {"pypgx": True,  "stargazer": False, "aldy": False, "stellarpgx": True},
    "GSTT1":    {"pypgx": False, "stargazer": False, "aldy": False, "stellarpgx": True},
    "POR":      {"pypgx": False, "stargazer": False, "aldy": False, "stellarpgx": True},
    "CYPOR":    {"pypgx": False, "stargazer": False, "aldy": False, "stellarpgx": True},
    "VKORC1":   {"pypgx": True,  "stargazer": True,  "aldy": False, "stellarpgx": False},
    "CYP1A1":   {"pypgx": False, "stargazer": True,  "aldy": True,  "stellarpgx": True},
    "CYP1A2":   {"pypgx": False, "stargazer": True,  "aldy": True,  "stellarpgx": True},
    "CYP2A6":   {"pypgx": False, "stargazer": True,  "aldy": True,  "stellarpgx": True},
    "CYP2E1":   {"pypgx": False, "stargazer": True,  "aldy": False, "stellarpgx": True},
    "IFNL3":    {"pypgx": True,  "stargazer": False, "aldy": False, "stellarpgx": False},
    "RYR1":     {"pypgx": True,  "stargazer": False, "aldy": False, "stellarpgx": False},
}


# ── Parser: PyPGx ─────────────────────────────────────────────────────────────
def parse_pypgx(output_dir: str, gene: str, sample: str) -> CallerResult:
    result = CallerResult(tool="PyPGx")

    # run-ngs-pipeline writes results.zip in the output directory
    candidates = [os.path.join(output_dir, "pypgx", "results.zip")]
    candidates += glob.glob(
        os.path.join(output_dir, "pypgx", "**", "results.zip"), recursive=True
    )
    zip_path = next((p for p in candidates if os.path.exists(p)), None)
    if zip_path is None:
        return result  # status stays "not_run"

    try:
        with zipfile.ZipFile(zip_path) as zf:
            tsv_name = next(
                (n for n in zf.namelist() if n.endswith("data.tsv")), None
            )
            if tsv_name is None:
                result.status = "failed"
                result.diplotype = "no data.tsv in results.zip"
                return result

            with zf.open(tsv_name) as f:
                reader = csv.DictReader(io.TextIOWrapper(f), delimiter="\t")
                for row in reader:
                    result.diplotype = (
                        row.get("Genotype") or row.get("Diplotype")
                        or row.get("Haplotype") or "-"
                    )
                    result.activity_score = str(
                        row.get("ActivityScore") or row.get("Activity_Score") or "-"
                    )
                    result.phenotype = row.get("Phenotype") or "-"
                    result.status = "ok"
                    break  # first data row
    except Exception as exc:
        result.status = "failed"
        result.diplotype = f"parse error: {exc}"

    return result


# ── Parser: Stargazer ─────────────────────────────────────────────────────────
def parse_stargazer(output_dir: str, gene: str, sample: str) -> CallerResult:
    result = CallerResult(tool="Stargazer")
    gene_lower = gene.lower()

    # Stargazer 2.0.3 output: genotype-calls.tsv (tab-separated, columns include
    # name, hap1_main, hap2_main, dip_score, phenotype).
    # Older versions wrote genotype.txt; both are searched.
    candidates = [
        os.path.join(output_dir, "stargazer", "genotype-calls.tsv"),
        os.path.join(output_dir, "stargazer", "genotype.txt"),
        os.path.join(output_dir, "stargazer", f"{gene_lower}.genotype.txt"),
        os.path.join(output_dir, "stargazer", f"{gene}.genotype.txt"),
    ]
    candidates += glob.glob(
        os.path.join(output_dir, "stargazer", "**", "genotype*.tsv"), recursive=True
    )
    candidates += glob.glob(
        os.path.join(output_dir, "stargazer", "**", "genotype*.txt"), recursive=True
    )
    geno_path = next((p for p in candidates if os.path.exists(p)), None)
    if geno_path is None:
        return result

    try:
        with open(geno_path) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                # Strip leading # from header keys
                row = {k.lstrip("#").strip(): v for k, v in row.items()}
                # genotype-calls.tsv uses hap1_main / hap2_main for the diplotype
                hap1 = row.get("hap1_main") or ""
                hap2 = row.get("hap2_main") or ""
                if hap1 and hap2:
                    result.diplotype = f"{hap1}/{hap2}"
                else:
                    result.diplotype = (
                        row.get("Genotype") or row.get("genotype")
                        or row.get("Diplotype") or "-"
                    )
                result.activity_score = str(
                    row.get("dip_score") or row.get("Activity_score")
                    or row.get("activity_score") or "-"
                )
                result.phenotype = (
                    row.get("phenotype") or row.get("Phenotype") or "-"
                )
                result.status = "ok"
                break
    except Exception as exc:
        result.status = "failed"
        result.diplotype = f"parse error: {exc}"

    return result


# ── Parser: Aldy ──────────────────────────────────────────────────────────────
def parse_aldy(output_dir: str, gene: str, sample: str) -> CallerResult:
    result = CallerResult(tool="Aldy")

    candidates = [
        os.path.join(output_dir, "aldy", f"{gene}.aldy"),
        os.path.join(output_dir, "aldy", f"{gene.lower()}.aldy"),
        os.path.join(output_dir, "aldy", f"{sample}.{gene}.aldy"),
    ]
    candidates += glob.glob(
        os.path.join(output_dir, "aldy", "**", "*.aldy"), recursive=True
    )
    aldy_path = next((p for p in candidates if os.path.exists(p)), None)
    if aldy_path is None:
        return result

    # Aldy 4.x output is tab-separated with a comment header.
    # Key columns: Sample, Gene, Solution#, Major, Minor, Score, Coverage, ...
    PHENOTYPE_KEYWORDS = (
        "metabolizer", "intermediate", "poor", "rapid", "ultra", "normal",
        "nm", "im", "pm", "rm", "um", "indeterminate"
    )
    try:
        with open(aldy_path) as f:
            header: list[str] = []
            for line in f:
                line = line.rstrip("\n")
                if not line:
                    continue
                if line.startswith("#"):
                    # May be header row
                    header = [c.lstrip("#").strip() for c in line.split("\t")]
                    continue
                parts = line.split("\t")
                if header and len(parts) == len(header):
                    row = dict(zip(header, parts))
                    # Try named columns first
                    major = row.get("Major") or row.get("Allele1") or ""
                    minor = row.get("Minor") or row.get("Allele2") or ""
                    if major and minor:
                        result.diplotype = f"{major}/{minor}"
                    elif major:
                        result.diplotype = major
                    result.activity_score = str(
                        row.get("Activity") or row.get("ActivityScore") or "-"
                    )
                    result.phenotype = (
                        row.get("Phenotype") or row.get("phenotype") or "-"
                    )
                    result.status = "ok"
                    break
                else:
                    # Fallback: heuristic scan for star alleles
                    star_cols = [p for p in parts if re.search(r"\*\d", p)]
                    if len(star_cols) >= 2:
                        result.diplotype = f"{star_cols[0]}/{star_cols[1]}"
                    elif star_cols:
                        result.diplotype = star_cols[0]

                    for p in parts:
                        if re.fullmatch(r"[\d.]+", p) and float(p) <= 5.0:
                            result.activity_score = p
                            break

                    for p in parts:
                        if any(kw in p.lower() for kw in PHENOTYPE_KEYWORDS):
                            result.phenotype = p
                            break

                    result.status = "ok"
                    break
    except Exception as exc:
        result.status = "failed"
        result.diplotype = f"parse error: {exc}"

    return result


# ── Parser: StellarPGx ────────────────────────────────────────────────────────
def parse_stellarpgx(output_dir: str, gene: str, sample: str) -> CallerResult:
    result = CallerResult(tool="StellarPGx")
    gene_lower = gene.lower()

    # StellarPGx 1.2.7 writes:  stellarpgx/<gene>/alleles/<sample>_<gene>.alleles
    # Older versions wrote:      stellarpgx/star_allele_calls/**
    candidates = (
        glob.glob(os.path.join(
            output_dir, "stellarpgx", gene_lower, "alleles", "*.alleles"
        ))
        + glob.glob(os.path.join(
            output_dir, "stellarpgx", gene, "alleles", "*.alleles"
        ))
    )
    star_dir = os.path.join(output_dir, "stellarpgx", "star_allele_calls")
    if os.path.isdir(star_dir):
        candidates += (
            glob.glob(os.path.join(star_dir, "**", f"*{gene_lower}*"), recursive=True)
            + glob.glob(os.path.join(star_dir, "**", "*.txt"), recursive=True)
        )
    call_path = next((p for p in candidates if os.path.isfile(p)), None)
    if call_path is None:
        return result

    try:
        with open(call_path) as f:
            lines = f.readlines()

        # Parse the structured .alleles format:
        # sections headed by "Result:", "Activity score:", "Metaboliser status:"
        for i, line in enumerate(lines):
            stripped = line.strip()
            if stripped == "Result:" and i + 1 < len(lines):
                dip = lines[i + 1].strip()
                if re.match(r"^\*", dip):
                    result.diplotype = dip
                    result.status = "ok"
            elif stripped == "Activity score:" and i + 1 < len(lines):
                result.activity_score = lines[i + 1].strip() or "-"
            elif stripped == "Metaboliser status:" and i + 1 < len(lines):
                result.phenotype = lines[i + 1].strip() or "-"
            # Fallback: diplotype pattern on a line (e.g. star_allele_calls format)
            elif result.status != "ok":
                m = re.search(r"(\*\w+/\*\w+)", stripped)
                if m:
                    result.diplotype = m.group(1)
                    result.status = "ok"
    except Exception as exc:
        result.status = "failed"
        result.diplotype = f"parse error: {exc}"

    return result


# ── Output ────────────────────────────────────────────────────────────────────
def _sv_note(gene: str) -> str:
    notes = []
    if gene in PYPGX_SV_GENES:
        notes.append("PyPGx: depth-of-coverage + VDR control stats")
    if gene in STARGAZER_SV_GENES:
        notes.append("Stargazer: GDF from BAM (paralog CN normalisation)")
    if notes:
        return "SV mode — " + "; ".join(notes)
    return "no SVs expected"


def print_table(
    gene: str,
    sample: str,
    results: list[CallerResult],
    output_dir: str,
) -> None:
    W = 72
    SEP = "─" * W

    print()
    print("=" * W)
    print(f" PGx Star Allele Results")
    print(f" Gene:   {gene}")
    print(f" Sample: {sample}")
    print(f" Build:  GRCh38")
    print(f" SV:     {_sv_note(gene)}")
    print("=" * W)
    print(f"{'Tool':<14}{'Diplotype':<18}{'Activity Score':<18}{'Phenotype'}")
    print(SEP)

    called_diplotypes: list[str] = []
    for r in results:
        print(f"{r.tool:<14}{r.diplotype:<18}{r.activity_score:<18}{r.phenotype}")
        if r.status == "ok" and r.diplotype not in ("-", ""):
            called_diplotypes.append(r.diplotype)

    print(SEP)

    if called_diplotypes:
        counts = Counter(called_diplotypes)
        top_dip, top_count = counts.most_common(1)[0]
        total = len(results)
        print(f"Concordance: {top_count}/{total} tools agree on {top_dip}")
    else:
        print("Concordance: no calls available")

    print("=" * W)
    print()

    tsv_path = os.path.join(output_dir, f"{gene}_{sample}_comparison.tsv")
    with open(tsv_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(
            ["Gene", "Sample", "Build", "Tool", "Diplotype", "ActivityScore", "Phenotype", "Status", "SVMode"]
        )
        for r in results:
            sv_mode = _sv_note(gene)
            writer.writerow(
                [gene, sample, "GRCh38", r.tool, r.diplotype, r.activity_score, r.phenotype, r.status, sv_mode]
            )
    print(f"Full results saved to: {tsv_path}")


# ── Main ──────────────────────────────────────────────────────────────────────
def main() -> None:
    parser = argparse.ArgumentParser(
        description="Parse PGx caller outputs and produce a comparison table."
    )
    parser.add_argument("--gene",       required=True, help="Gene name (e.g. CYP2D6)")
    parser.add_argument("--sample",     required=True, help="Sample name")
    parser.add_argument("--output-dir", required=True, help="Results directory")
    args = parser.parse_args()

    gene = args.gene.upper()
    support = GENE_SUPPORT.get(gene, {})
    if not support:
        print(f"ERROR: Gene '{gene}' not found in support matrix.", file=sys.stderr)
        sys.exit(1)

    results: list[CallerResult] = []
    if support.get("pypgx"):
        results.append(parse_pypgx(args.output_dir, gene, args.sample))
    if support.get("stargazer"):
        results.append(parse_stargazer(args.output_dir, gene, args.sample))
    if support.get("aldy"):
        results.append(parse_aldy(args.output_dir, gene, args.sample))
    if support.get("stellarpgx"):
        results.append(parse_stellarpgx(args.output_dir, gene, args.sample))

    print_table(gene, args.sample, results, args.output_dir)


if __name__ == "__main__":
    main()
