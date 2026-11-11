#!/usr/bin/env python3
"""
pgx-compare.py — Parse outputs from all 4 PGx callers and produce a comparison
table (TSV) and a full-detail JSON used by pgx-report.py for HTML reporting.

Usage: pgx-compare.py --gene GENE --sample SAMPLE --output-dir DIR
"""
import argparse
import csv
import glob
import gzip
import io
import json
import os
import re
import subprocess
import sys
import zipfile
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from typing import Any

# Variant-tier reconciliation (A.1). pgx-compare.py is a symlinked entrypoint, so
# resolve the real dir before importing the sibling module.
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
import reconcile  # noqa: E402
import omnigen_addons as _oa  # noqa: E402  (ABO + HLA class-II helpers)

# VCF-Check variant tables (design D) — curated clinical variant lists kept as
# data, not code. Currently RYR1 (CPIC malignant-hyperthermia alleles); CACNA1S
# and G6PD may migrate here later.
_VCF_CHECK_DATA = os.path.join(os.path.dirname(os.path.realpath(__file__)), "vcf_check_variants.json")
try:
    with open(_VCF_CHECK_DATA, encoding="utf-8") as _vfh:
        _VCF_CHECK_VARIANTS: dict = json.load(_vfh)
except FileNotFoundError:
    _VCF_CHECK_VARIANTS = {}

# Mean read depth (x) below which a gene cannot be called. Without this gate a
# variant caller that sees no ALT alleles in an uncovered region emits wild-type
# (*1/*1) — a confident-looking call manufactured from missing data.
# ponytail: single global floor; CNV-heavy genes may warrant a higher value —
# override per run with --min-depth if assay validation demands it.
DEFAULT_MIN_DEPTH = 10.0


@dataclass
class CallerResult:
    tool: str
    # ── Fields 3-7: allele calls ──────────────────────────────────────────────
    diplotype: str = "-"
    haplotype1: str = "-"
    haplotype2: str = "-"
    sub_alleles: str = "-"
    alternative_diplotypes: str = "-"
    # ── Fields 8-9: clinical ──────────────────────────────────────────────────
    phenotype: str = "-"
    activity_score: str = "-"
    # ── Fields 10-11: structural variants ─────────────────────────────────────
    sv_type: str = "-"
    copy_number: str = "-"
    # ── Fields 12-16: supporting evidence ────────────────────────────────────
    supporting_variants: list = field(default_factory=list)
    functional_effects: str = "-"
    dbsnp_ids: str = "-"
    allele_score: str = "-"
    mean_af: str = "-"
    # ── Field 17: method ──────────────────────────────────────────────────────
    phasing_method: str = "-"
    # ── Operational ───────────────────────────────────────────────────────────
    status: str = "not_run"   # not_run | ok | failed

    def to_dict(self, sample: str, gene: str) -> dict[str, Any]:
        """Serialise to the 17-field dict written into the detail JSON."""
        return {
            "sample_id":              sample,
            "gene":                   gene,
            "diplotype":              self.diplotype,
            "haplotype1":             self.haplotype1,
            "haplotype2":             self.haplotype2,
            "sub_alleles":            self.sub_alleles,
            "alternative_diplotypes": self.alternative_diplotypes,
            "phenotype":              self.phenotype,
            "activity_score":         self.activity_score,
            "sv_type":                self.sv_type,
            "copy_number":            self.copy_number,
            "supporting_variants":    self.supporting_variants,
            "functional_effects":     self.functional_effects,
            "dbsnp_ids":              self.dbsnp_ids,
            "allele_score":           self.allele_score,
            "mean_af":                self.mean_af,
            "phasing_method":         self.phasing_method,
            "status":                 self.status,
        }


# ── SV gene sets ──────────────────────────────────────────────────────────────
# ── Gene support matrix — single source of truth ──────────────────────────────
# Loaded from genes.tsv, the SAME file pgx-run.sh parses. Previously this matrix
# was duplicated here (Python dict) and in pgx-run.sh (bash arrays) in different
# encodings, with nothing checking they agreed.
def load_gene_config(path: str | None = None):
    """Return (support, pypgx_sv, stargazer_sv) read from genes.tsv.

    `support[gene]` holds only the True tool flags; consumers use
    `support.get(tool)` (None is falsy), so omitting False is behaviour-identical
    to the old explicit-False dict.
    """
    path = path or os.environ.get("PGX_GENES_TSV") \
        or os.path.join(os.path.dirname(os.path.abspath(__file__)), "genes.tsv")
    support: dict[str, dict[str, bool]] = {}
    pypgx_sv: set[str] = set()
    stargazer_sv: set[str] = set()
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            g = row["gene"]
            support[g] = {
                k: True
                for k in ("pypgx", "stargazer", "aldy", "stellarpgx",
                          "optitype", "mutserve", "vcf_check", "cyrius", "pharmcat",
                          "arcashla")
                if row.get(k) == "1"
            }
            if row["pypgx_sv"] == "1":
                pypgx_sv.add(g)
            if row["stargazer_sv"] == "1":
                stargazer_sv.add(g)
    return support, frozenset(pypgx_sv), frozenset(stargazer_sv)


GENE_SUPPORT, PYPGX_SV_GENES, STARGAZER_SV_GENES = load_gene_config()


# ── Helpers ───────────────────────────────────────────────────────────────────
def _dash(v: str | None) -> str:
    return (v or "").strip() or "-"


def _parse_pypgx_variant_data(vdata: str) -> tuple[list[dict], str, str]:
    """Parse PyPGx VariantData field.

    Format: *allele:chrom-pos-ref-alt[,pos-ref-alt,...]:AF[,AF,...];...
    Returns (variants_list, functional_effects, dbsnp_ids).
    PyPGx does not report rsIDs or effect labels in data.tsv.
    """
    variants: list[dict] = []
    if not vdata or vdata == "-":
        return variants, "-", "-"
    for entry in vdata.rstrip(";").split(";"):
        entry = entry.strip()
        if not entry:
            continue
        m = re.match(r"^(\S+?):(.+)$", entry)
        if not m:
            continue
        allele = m.group(1)
        rest = m.group(2)
        # rest = "chrom-pos-ref-alt[,...]:AF[,...]"
        parts = rest.rsplit(":", 1)
        loci_str = parts[0]
        afs_str = parts[1] if len(parts) > 1 else ""
        afs = [a.strip() for a in afs_str.split(",")]
        loci = [l.strip() for l in loci_str.split(",")]
        for i, locus in enumerate(loci):
            lm = re.match(r"(\w+)-(\d+)-([A-Z]+)-([A-Z<>.]+)", locus)
            if lm:
                variants.append({
                    "allele": allele,
                    "chrom":  lm.group(1),
                    "pos":    lm.group(2),
                    "ref":    lm.group(3),
                    "alt":    lm.group(4),
                    "af":     afs[i] if i < len(afs) else "-",
                    "depth":  "-",
                    "effect": "-",
                    "rsid":   "-",
                })
    return variants, "-", "-"


def _parse_stargazer_core(core_str: str) -> tuple[list[dict], list[str]]:
    """Parse Stargazer hap_main_core field.

    Format: <pos:ref>alt:reads_ref/reads_alt:AF:effect_type:impact:protein_change>,...
    Returns (variants, effect_labels).
    """
    variants: list[dict] = []
    effects: list[str] = []
    if not core_str or core_str in ("-", "."):
        return variants, effects
    for token in core_str.split(">,"):
        token = token.strip().lstrip("<").rstrip(">")
        if not token:
            continue
        # <42126611:C>G:37/37:1.00:missense_variant:low_impact:S486T>
        m = re.match(
            r"(\d+):([A-Z]+)>([A-Z.]+):(\d+)/(\d+):([\d.]+)"
            r"(?::([^:]+))?(?::([^:]+))?(?::([^:>]+))?",
            token,
        )
        if m:
            depth_ref, depth_alt = int(m.group(4)), int(m.group(5))
            total = depth_ref + depth_alt
            effect = m.group(7) or "-"
            protein = m.group(9) or "-"
            if effect != "-":
                effects.append(f"{protein} ({effect})" if protein != "-" else effect)
            variants.append({
                "allele": "-",
                "chrom":  "-",
                "pos":    m.group(1),
                "ref":    m.group(2),
                "alt":    m.group(3),
                "af":     m.group(6),
                "depth":  str(total),
                "effect": f"{protein} ({effect})" if protein != "-" else effect,
                "rsid":   "-",
            })
    return variants, effects


# ── Parser: PyPGx ─────────────────────────────────────────────────────────────
def parse_pypgx(output_dir: str, gene: str, sample: str) -> CallerResult:
    result = CallerResult(tool="PyPGx", phasing_method="Beagle (1KGP reference panel)")

    candidates = [os.path.join(output_dir, "pypgx", "results.zip")]
    candidates += glob.glob(
        os.path.join(output_dir, "pypgx", "**", "results.zip"), recursive=True
    )
    zip_path = next((p for p in candidates if os.path.exists(p)), None)
    if zip_path is None:
        return result

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
                    result.diplotype = _dash(
                        row.get("Genotype") or row.get("Diplotype") or row.get("Haplotype")
                    )
                    result.activity_score = _dash(
                        row.get("ActivityScore") or row.get("Activity_Score")
                    )
                    result.phenotype = _dash(row.get("Phenotype"))
                    result.haplotype1 = _dash(row.get("Haplotype1"))
                    result.haplotype2 = _dash(row.get("Haplotype2"))

                    # Sub-alleles: Haplotype2 often contains phased sub-alleles (;-separated)
                    h2 = row.get("Haplotype2") or ""
                    parts = [p for p in h2.split(";") if p.strip()]
                    if len(parts) > 1:
                        result.sub_alleles = "; ".join(parts)

                    result.alternative_diplotypes = _dash(row.get("AlternativePhase"))
                    result.sv_type = _dash(row.get("CNV"))
                    result.copy_number = "2" if result.sv_type == "Normal" else "-"

                    vdata = row.get("VariantData") or ""
                    result.supporting_variants, result.functional_effects, result.dbsnp_ids = \
                        _parse_pypgx_variant_data(vdata)

                    result.allele_score = "-"
                    result.mean_af = "-"
                    if result.supporting_variants:
                        afs = [v["af"] for v in result.supporting_variants if v["af"] != "-"]
                        if afs:
                            result.mean_af = f"per-variant AF: {', '.join(afs[:4])}" + \
                                             ("…" if len(afs) > 4 else "")

                    result.status = "ok"
                    break
    except Exception as exc:
        result.status = "failed"
        result.diplotype = f"parse error: {exc}"

    return result


# ── Parser: Stargazer ─────────────────────────────────────────────────────────
def parse_stargazer(output_dir: str, gene: str, sample: str) -> CallerResult:
    result = CallerResult(
        tool="Stargazer",
        phasing_method="Beagle statistical phasing + SSR CN marker",
    )
    gene_lower = gene.lower()

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
                row = {k.lstrip("#").strip(): v for k, v in row.items()}

                hap1 = _dash(row.get("hap1_main"))
                hap2 = _dash(row.get("hap2_main"))
                if hap1 != "-" and hap2 != "-":
                    result.diplotype = f"{hap1}/{hap2}"
                else:
                    result.diplotype = _dash(
                        row.get("Genotype") or row.get("genotype") or row.get("Diplotype")
                    )

                raw_score = row.get("dip_score") or ""
                try:
                    if float(raw_score) < 0:
                        result.status = "not_run"
                        break
                except (ValueError, TypeError):
                    pass

                result.haplotype1 = hap1
                result.haplotype2 = hap2
                result.activity_score = _dash(raw_score)
                result.phenotype = _dash(row.get("phenotype") or row.get("Phenotype"))

                # Sub-alleles from candidate lists
                h1c = _dash(row.get("hap1_cand"))
                h2c = _dash(row.get("hap2_cand"))
                if h1c != "-" or h2c != "-":
                    result.sub_alleles = f"Hap1 candidates: {h1c} | Hap2 candidates: {h2c}"

                result.alternative_diplotypes = _dash(row.get("dip_cand"))

                sv1 = _dash(row.get("hap1_sv"))
                sv2 = _dash(row.get("hap2_sv"))
                result.sv_type = _dash(row.get("dip_sv")) or f"{sv1} / {sv2}"
                result.copy_number = "2" if "no_sv" in result.sv_type else "-"

                # Supporting variants from core fields
                variants: list[dict] = []
                effects: list[str] = []
                for hap_field in ("hap1_main_core", "hap2_main_core",
                                  "hap1_main_tag", "hap2_main_tag"):
                    v, e = _parse_stargazer_core(row.get(hap_field) or "")
                    variants.extend(v)
                    effects.extend(e)
                result.supporting_variants = variants
                result.functional_effects = "; ".join(dict.fromkeys(effects)) or "-"

                ssr = _dash(row.get("ssr"))
                h1s = _dash(row.get("hap1_score"))
                h2s = _dash(row.get("hap2_score"))
                result.allele_score = (
                    f"dip={raw_score}, hap1={h1s}, hap2={h2s}, SSR={ssr}"
                )

                af1g = _dash(row.get("hap1_af_mean_gene"))
                af2g = _dash(row.get("hap2_af_mean_gene"))
                af1m = _dash(row.get("hap1_af_mean_main"))
                af2m = _dash(row.get("hap2_af_mean_main"))
                result.mean_af = (
                    f"Hap1 gene={af1g}, main={af1m} | Hap2 gene={af2g}, main={af2m}"
                )

                result.status = "ok"
                break
    except Exception as exc:
        result.status = "failed"
        result.diplotype = f"parse error: {exc}"

    return result


# ── Parser: Aldy ──────────────────────────────────────────────────────────────
def parse_aldy(output_dir: str, gene: str, sample: str) -> CallerResult:
    result = CallerResult(
        tool="Aldy",
        phasing_method="Integer Linear Programming (joint CN + allele decomposition)",
    )

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

    try:
        with open(aldy_path) as f:
            lines = f.readlines()

        header: list[str] = []
        cpic_phenotype = "-"
        cpic_score = "-"
        diplotype = "-"
        minor = "-"
        variants: list[dict] = []
        effects: list[str] = []
        rsids: list[str] = []
        copies_seen: set[str] = set()

        for line in lines:
            line = line.rstrip("\n")
            if not line:
                continue

            if line.startswith("#"):
                # Extract cpic metadata from solution comment
                m_cpic = re.search(r"cpic=(\S+?)(?:;|$)", line)
                m_score = re.search(r"cpic_score=([\d.]+)", line)
                if m_cpic:
                    cpic_phenotype = m_cpic.group(1)
                if m_score:
                    cpic_score = m_score.group(1)
                # Parse column header
                candidate_header = [c.lstrip("#").strip() for c in line.split("\t")]
                if "Major" in candidate_header or "Sample" in candidate_header:
                    header = candidate_header
                continue

            if not header:
                continue

            parts = line.split("\t")
            if len(parts) < 4:   # need at least Sample, Gene, SolutionID, Major
                continue
            row = dict(zip(header, parts))

            # First data row: extract diplotype / minor
            if diplotype == "-":
                major = row.get("Major") or ""
                diplotype = major.strip() or "-"
                minor = _dash(row.get("Minor"))
                result.sv_type = "-"

            # Per-variant rows
            loc = row.get("Location") or row.get("Allele") or ""
            typ = row.get("Type") or ""
            cov = row.get("Coverage") or "-"
            eff = row.get("Effect") or "-"
            rsid = row.get("dbSNP") or "-"
            copy = row.get("Copy") or ""

            if loc and typ and re.search(r"[A-Z]>[A-Z]", typ):
                ref_alt = typ.split(">")
                variant = {
                    "allele": f"Copy {copy}" if copy != "" else "-",
                    "chrom":  "-",
                    "pos":    loc,
                    "ref":    ref_alt[0] if len(ref_alt) > 0 else "-",
                    "alt":    ref_alt[1] if len(ref_alt) > 1 else "-",
                    "af":     "-",
                    "depth":  cov,
                    "effect": eff if eff != "none" else "-",
                    "rsid":   rsid,
                }
                variants.append(variant)
                if eff and eff not in ("none", "-"):
                    effects.append(eff)
                if rsid and rsid != "-":
                    rsids.append(rsid)

        result.diplotype = diplotype
        result.haplotype1 = diplotype.split("/")[0].strip() if "/" in diplotype else diplotype
        result.haplotype2 = diplotype.split("/")[1].strip() if "/" in diplotype else "-"
        result.sub_alleles = minor
        result.phenotype = cpic_phenotype.replace("_", " ").title() if cpic_phenotype != "-" else "-"
        result.activity_score = cpic_score
        result.copy_number = "2"   # Aldy always solves to diploid by default
        result.supporting_variants = variants
        result.functional_effects = "; ".join(dict.fromkeys(effects)) or "-"
        result.dbsnp_ids = ", ".join(dict.fromkeys(rsids)) or "-"
        result.allele_score = f"SolutionID=1 (single best solution)"
        result.mean_af = f"{len(variants)} variants, depth in Coverage column"
        if diplotype != "-":
            result.status = "ok"

    except Exception as exc:
        result.status = "failed"
        result.diplotype = f"parse error: {exc}"

    return result


# ── Parser: StellarPGx ────────────────────────────────────────────────────────
def parse_stellarpgx(output_dir: str, gene: str, sample: str) -> CallerResult:
    result = CallerResult(
        tool="StellarPGx",
        phasing_method="Genome graph variant calling (graphtyper / vg)",
    )
    gene_lower = gene.lower()

    candidates = (
        glob.glob(os.path.join(output_dir, "stellarpgx", gene_lower, "alleles", "*.alleles"))
        + glob.glob(os.path.join(output_dir, "stellarpgx", gene, "alleles", "*.alleles"))
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

        variants: list[dict] = []
        cn = "-"

        for i, line in enumerate(lines):
            stripped = line.strip()

            if stripped.startswith("Initially computed CN"):
                m = re.search(r"=\s*(\d+)", stripped)
                if m:
                    cn = m.group(1)
                    result.copy_number = cn

            elif stripped == "Sample core variants:" and i + 1 < len(lines):
                var_str = lines[i + 1].strip()
                for tok in var_str.rstrip(";").split(";"):
                    tok = tok.strip()
                    if not tok:
                        continue
                    # format: pos~ref>alt~GT
                    m = re.match(r"(\d+)~([A-Z]+)>([A-Z]+)~([\d/|]+)", tok)
                    if m:
                        variants.append({
                            "allele": "-",
                            "chrom":  "-",
                            "pos":    m.group(1),
                            "ref":    m.group(2),
                            "alt":    m.group(3),
                            "af":     "-",
                            "depth":  "-",
                            "effect": "-",
                            "rsid":   "-",
                            "gt":     m.group(4),
                        })

            elif stripped == "Candidate alleles:" and i + 1 < len(lines):
                result.sub_alleles = lines[i + 1].strip()

            elif stripped == "Result:" and i + 1 < len(lines):
                dip = lines[i + 1].strip()
                if re.match(r"^\*", dip):
                    result.diplotype = dip
                    alleles = dip.split("/")
                    result.haplotype1 = alleles[0].strip() if alleles else "-"
                    result.haplotype2 = alleles[1].strip() if len(alleles) > 1 else "-"
                    result.status = "ok"

            elif stripped == "Activity score:" and i + 1 < len(lines):
                result.activity_score = _dash(lines[i + 1].strip())

            elif stripped == "Metaboliser status:" and i + 1 < len(lines):
                result.phenotype = _dash(lines[i + 1].strip())

            elif result.status != "ok":
                m = re.search(r"(\*\w+/\*\w+)", stripped)
                if m:
                    result.diplotype = m.group(1)
                    result.status = "ok"

        result.supporting_variants = variants
        result.sv_type = f"CN={cn}" if cn != "-" else "-"
        result.allele_score = "-"
        result.mean_af = "GT per variant (0/1 or 1/1)"

    except Exception as exc:
        result.status = "failed"
        result.diplotype = f"parse error: {exc}"

    return result


# ── Parser: OptiType (HLA-A, HLA-B, HLA-C) ──────────────────────────────────
def parse_optitype(output_dir: str, gene: str, sample: str) -> CallerResult:
    """Parse OptiType result TSV for HLA-A, HLA-B or HLA-C.

    OptiType outputs a TSV with columns:
        (index)  A1  A2  B1  B2  C1  C2  Reads  Objective
    Values are like "A*01:01" (without "HLA-" prefix).

    OptiType types all three class-I loci in a single run. Until 2026-07 only the
    A and B columns were read and the C1/C2 columns were silently discarded, so
    HLA-C never reached the per-gene outputs — which made the downstream OmniGen
    HLA-C*06:02 psoriasis association unscreenable (it rendered as "absent"
    rather than "not typed") for every sample ever reported.
    """
    result = CallerResult(
        tool="OptiType",
        phasing_method="ILP optimisation on HLA-reference-filtered MHC reads",
    )

    # OptiType now runs ONCE per sample (a sample-level Snakemake rule) and writes
    # a single shared result TSV at <sample_root>/optitype/ — all three class-I loci
    # (A/B/C) parse that same TSV. Older per-gene runs wrote it under the gene dir
    # (<sample_root>/genes/<gene>/optitype/). Check the sample-level path first, then
    # fall back to the per-gene dir so results produced by the previous layout (and
    # the parser unit tests, which pass the gene dir directly) still parse.
    sample_root = os.path.dirname(os.path.dirname(output_dir))  # .../genes/<gene> -> sample root
    optitype_dirs = [
        os.path.join(sample_root, "optitype"),   # sample-level (current)
        os.path.join(output_dir, "optitype"),    # per-gene (legacy / unit-test layout)
    ]
    # Look for <sample>_result.tsv or any *_result.tsv, sample-level location first.
    candidates: list[str] = []
    for optitype_dir in optitype_dirs:
        candidates.append(os.path.join(optitype_dir, f"{sample}_result.tsv"))
        candidates += glob.glob(os.path.join(optitype_dir, "*_result.tsv"))
    tsv_path = next((p for p in candidates if os.path.exists(p)), None)
    if tsv_path is None:
        return result

    try:
        with open(tsv_path) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                if gene == "HLA-A":
                    a1 = (row.get("A1") or "").strip()
                    a2 = (row.get("A2") or "").strip()
                    if a1 and a2:
                        # OptiType outputs "A*01:01"; prepend "HLA-" → "HLA-A*01:01"
                        result.haplotype1 = f"HLA-{a1}" if not a1.startswith("HLA") else a1
                        result.haplotype2 = f"HLA-{a2}" if not a2.startswith("HLA") else a2
                        result.diplotype  = f"{result.haplotype1}/{result.haplotype2}"
                        result.status = "ok"
                elif gene == "HLA-B":
                    b1 = (row.get("B1") or "").strip()
                    b2 = (row.get("B2") or "").strip()
                    if b1 and b2:
                        # OptiType outputs "B*57:01"; prepend "HLA-" → "HLA-B*57:01"
                        result.haplotype1 = f"HLA-{b1}" if not b1.startswith("HLA") else b1
                        result.haplotype2 = f"HLA-{b2}" if not b2.startswith("HLA") else b2
                        result.diplotype  = f"{result.haplotype1}/{result.haplotype2}"
                        result.status = "ok"
                elif gene == "HLA-C":
                    c1 = (row.get("C1") or "").strip()
                    c2 = (row.get("C2") or "").strip()
                    if c1 and c2:
                        # OptiType outputs "C*06:02"; prepend "HLA-" → "HLA-C*06:02"
                        result.haplotype1 = f"HLA-{c1}" if not c1.startswith("HLA") else c1
                        result.haplotype2 = f"HLA-{c2}" if not c2.startswith("HLA") else c2
                        result.diplotype  = f"{result.haplotype1}/{result.haplotype2}"
                        result.status = "ok"
                reads = _dash(row.get("Reads"))
                obj   = _dash(row.get("Objective"))
                result.allele_score = f"Reads={reads}, Objective={obj}"
                break
    except Exception as exc:
        result.status = "failed"
        result.diplotype = f"parse error: {exc}"

    return result


# ── Parser: mutserve (MT-RNR1) ────────────────────────────────────────────────
def parse_mutserve(output_dir: str, gene: str, sample: str) -> CallerResult:
    """Parse pgx-mt.sh JSON result for MT-RNR1."""
    result = CallerResult(
        tool="mutserve",
        phasing_method="AF-based mitochondrial variant calling (mutserve v2)",
    )
    json_path = os.path.join(output_dir, "mt-rnr1", f"{sample}_mtrna1_result.json")
    if not os.path.exists(json_path):
        return result
    try:
        import json as _json
        data = _json.load(open(json_path))
        result.diplotype  = data.get("diplotype", "-")
        result.phenotype  = data.get("phenotype", "-")
        result.status     = "ok"
        variants = data.get("variants", [])
        if variants:
            af_strs = [f"{v['label']}(AF={v['af']:.2f},{v['type']},{v.get('cpic_class', 'CPIC')})"
                       for v in variants]
            result.allele_score = "; ".join(af_strs)
            result.supporting_variants = [
                {"allele": v["label"], "chrom": "chrM", "pos": str(v.get("pos", "-")),
                 "ref": v.get("ref", "-"), "alt": v.get("alt", "-"), "af": f"{v['af']:.4f}",
                 "depth": str(v.get("coverage", "-")),
                 "effect": f"{v.get('cpic_class', '')} ({v['type']})", "rsid": "-"}
                for v in variants]
        elif data.get("diplotype") == "not typed":
            result.allele_score = "not typed (insufficient chrM coverage)"
        else:
            result.allele_score = "No CPIC increased-risk MT-RNR1 variant detected"
    except Exception as exc:
        result.status    = "failed"
        result.diplotype = f"parse error: {exc}"
    return result


# ── Parser: Cyrius (CYP2D6 SV-aware caller) ──────────────────────────────────
# Cyrius resolves CYP2D6 from the 2D6/2D7 paralog + CNV/hybrids that the
# star-allele callers struggle with. Output is a TSV: Sample / Genotype / Filter.
# Genotype is the CYP2D6 diplotype (xN for duplications, hybrids like *68+*4).
def parse_cyrius(output_dir: str, gene: str, sample: str) -> CallerResult:
    """Parse Cyrius CYP2D6 output (genes/CYP2D6/cyrius/<prefix>.tsv)."""
    result = CallerResult(tool="Cyrius")
    tsvs = glob.glob(os.path.join(output_dir, "cyrius", "*.tsv"))
    if not tsvs:
        result.status = "failed"
        result.diplotype = "no Cyrius output"
        return result
    try:
        with open(tsvs[0]) as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        row = next((r for r in rows if r.get("Sample") == sample), rows[0] if rows else None)
        if not row:
            result.status = "failed"
            result.diplotype = "empty Cyrius output"
            return result
        gt = (row.get("Genotype") or "").strip()
        filt = (row.get("Filter") or "").strip()
        if not gt or gt.lower() in ("none", "nan"):
            result.status = "failed"
            result.diplotype = f"Cyrius no-call ({filt or 'None'})"
            return result
        result.diplotype = gt
        result.haplotype1 = gt.split("/")[0].strip() if "/" in gt else gt
        result.haplotype2 = gt.split("/")[1].strip() if "/" in gt else "-"
        if filt and filt.upper() != "PASS":
            result.alternative_diplotypes = filt
        result.status = "ok"
    except Exception as exc:
        result.status = "failed"
        result.diplotype = f"parse error: {exc}"
    return result


# ── Parser: PharmCAT (CPIC reference implementation, sample-level) ───────────
# PharmCAT runs once per sample (via its Apptainer SIF) and writes one report.json
# covering all its CPIC genes. It is the authority for the genes where the
# star-allele callers genuinely disagree (UGT1A1 *28/*80 LD, CYP2B6, CYP4F2). The
# report lives at OUT/pharmcat/<base>.report.json, two levels up from the per-gene
# output dir. Structure: genes[<gene>].sourceDiplotypes[0].{label,phenotypes}.
def parse_pharmcat(output_dir: str, gene: str, sample: str) -> CallerResult:
    """Read this gene's call from the sample-level PharmCAT report.json."""
    result = CallerResult(tool="PharmCAT")
    pc_dir = os.path.normpath(os.path.join(output_dir, "..", "..", "pharmcat"))
    reports = glob.glob(os.path.join(pc_dir, "*.report.json"))
    if not reports:
        result.status = "failed"
        result.diplotype = "no PharmCAT report"
        return result
    try:
        d = json.load(open(reports[0]))
        genes = d.get("genes", {})
        node = genes.get(gene)
        if node is None:  # tolerate source-nested layouts (CPIC/DPWG/...)
            for src in genes.values():
                if isinstance(src, dict) and gene in src:
                    node = src[gene]
                    break
        if node is None:
            result.status = "failed"
            result.diplotype = f"{gene} not in PharmCAT report"
            return result
        sd = node.get("sourceDiplotypes") or []
        labels = [x.get("label") for x in sd if x.get("label")]
        if not labels:
            result.status = "failed"
            result.diplotype = "PharmCAT no-call"
            return result
        result.diplotype = labels[0]
        result.haplotype1 = labels[0].split("/")[0].strip() if "/" in labels[0] else labels[0]
        result.haplotype2 = labels[0].split("/")[1].strip() if "/" in labels[0] else "-"
        phenos = sd[0].get("phenotypes") or []
        result.phenotype = phenos[0] if phenos else "-"
        if len(labels) > 1:
            result.alternative_diplotypes = "; ".join(labels[1:4])
        result.status = "ok"
    except Exception as exc:
        result.status = "failed"
        result.diplotype = f"parse error: {exc}"
    return result


# ── Parser: CACNA1S VCF-Check (MHS pathogenic alleles *2 and *3) ──────────────
# Only PyPGx and StellarPGx call CACNA1S; if StellarPGx fails, a single-tool
# call is the only result. The entire CPIC-relevant allele set for CACNA1S is
# just two missense SNPs — both are directly detectable from the gene VCF.
# NOTE: CACNA1S is on the negative strand. REF/ALT below are positive-strand
# (genomic) alleles as reported in VCF FORMAT.
# *2 (rs772226819) sits at chr1:201091993 — 8 kb beyond the old GENE_COORDS
# endpoint; pgx-run.sh CACNA1S coords have been extended to chr1:201006956-
# 201095000 to capture it in the bcftools VCF.
_CACNA1S_SNP_ALLELES: list[dict] = [
    {"pos": 201091993, "ref": "G", "alt": "A", "allele": "*2", "rsid": "rs772226819",
     "note": "p.Arg174Trp (c.520C>T); most common CACNA1S MHS allele in European families; CPIC Level A"},
    {"pos": 201060815, "ref": "C", "alt": "T", "allele": "*3", "rsid": "rs1800559",
     "note": "p.Arg1086His (c.3257G>A); most common overall CACNA1S MHS allele; CPIC Level A"},
]


def parse_cacna1s_vcf(output_dir: str, gene: str, sample: str) -> CallerResult:
    """Query CACNA1S.vcf.gz for *2 (Arg174Trp) and *3 (Arg1086His) MHS-associated alleles."""
    result = CallerResult(
        tool="VCF-Check",
        phasing_method="bcftools mpileup — direct variant query (unphased)",
    )

    vcf = os.path.join(output_dir, f"{gene}.vcf.gz")
    if not os.path.exists(vcf):
        result.status = "failed"
        result.diplotype = "VCF not found"
        return result

    region = "chr1:201060815-201091993"
    try:
        proc = subprocess.run(
            ["bcftools", "view", "-r", region, vcf],
            capture_output=True, text=True, timeout=30,
        )
    except FileNotFoundError:
        result.status = "failed"
        result.diplotype = "bcftools not in PATH"
        return result
    except subprocess.TimeoutExpired:
        result.status = "failed"
        result.diplotype = "bcftools timed out"
        return result

    alleles_found: dict[str, tuple] = {}
    variants_list: list[dict] = []

    for line in proc.stdout.splitlines():
        if line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 10:
            continue
        pos = int(fields[1])
        ref, alt_str = fields[3], fields[4]
        fmt_keys = fields[8].split(":")
        fmt_vals = fields[9].split(":")
        fmt = dict(zip(fmt_keys, fmt_vals))
        gt_raw = fmt.get("GT", "./.")
        gt = gt_raw.replace("|", "/")
        dp  = fmt.get("DP", "-")
        ad  = fmt.get("AD", "-")
        alt_copies = sum(1 for g in gt.split("/") if g not in ("0", "."))
        if alt_copies == 0:
            continue
        alts = alt_str.split(",")
        for adef in _CACNA1S_SNP_ALLELES:
            if adef["pos"] != pos or adef["alt"] not in alts:
                continue
            alleles_found[adef["allele"]] = (gt, dp, alt_copies, adef["rsid"], adef["note"])
            variants_list.append({
                "allele": adef["allele"], "chrom": "chr1", "pos": str(pos),
                "ref": ref, "alt": adef["alt"],
                "af": ad, "depth": dp, "effect": adef["note"],
                "rsid": adef["rsid"],
            })

    result.status = "ok"
    result.supporting_variants = variants_list

    if not alleles_found:
        result.diplotype = "*1/*1"
        result.phenotype = ("No CPIC-listed CACNA1S MHS variant found among the "
                            f"{len(_CACNA1S_SNP_ALLELES)} sites checked (c.520C>T, c.3257G>A)")
        return result

    def _sort_key(name: str) -> int:
        m = re.search(r"(\d+)", name)
        return int(m.group(1)) if m else 999

    parts = []
    for aname in sorted(alleles_found, key=_sort_key):
        gt, dp, copies, rsid, note = alleles_found[aname]
        zygosity = "hom" if copies >= 2 else "het"
        parts.append(f"{aname}({zygosity})")
    result.diplotype = "+".join(parts) + " [unphased]"
    result.phenotype = "MHS-associated variant(s) detected — evaluate for malignant hyperthermia susceptibility"
    return result


# ── Parser: G6PD VCF-Check (full CPIC/ClinPGx allele table) ───────────────────
# G6PD is typed against the complete CPIC G6PD allele definition table that
# PharmCAT 3.2.0 ships (171 allele-defining sites, 184 non-reference named
# alleles, GRCh38; stored in vcf_check_variants.json["G6PD"]). An allele is
# called when ALL of its defining sites carry the alt allele; nested calls are
# resolved to the most specific allele (A- 202A_376G beats A 376G). The CPIC
# function class of each allele (I/II/III Deficient, IV Normal, Uncertain) and
# the sample's sex give the CPIC phenotype: hemizygous male with a deficient
# allele = Deficient; female with two deficient alleles = Deficient; female with
# one = Variable (CPIC G6PD guideline, Gammal et al. 2023).
# Governing rule: ERROR is never NEGATIVE. "No allele found" is only ever worded
# "no deficiency variant among the N sites checked" - never "Normal" - and
# sites without enough depth are counted as not typed, not as reference.
# The compare CLI sets _RUN_CTX (bam/ref/sex) so site depth and sex can be read.
_RUN_CTX: dict = {"bam": None, "ref": None, "sex": None}


def _infer_sex_from_bam(bam: str) -> str | None:
    """MALE/FEMALE from the chrX:chrY idxstats read ratio (same rule as pgx-bamstats.sh)."""
    try:
        out = subprocess.run(["samtools", "idxstats", bam], capture_output=True, text=True,
                             timeout=120).stdout
    except Exception:
        return None
    x = y = xl = yl = 0
    for ln in out.splitlines():
        f = ln.split("\t")
        if len(f) < 4:
            continue
        if f[0] in ("chrX", "X"):
            xl, x = int(f[1]), int(f[2])
        elif f[0] in ("chrY", "Y"):
            yl, y = int(f[1]), int(f[2])
    if not xl or not yl or (x + y) == 0:
        return None
    xd, yd = x / xl, y / yl
    if yd <= 0:
        return None   # no chrY reads at all: stripped/extracted BAM as likely as female
    return "FEMALE" if xd / yd > 5 else "MALE"


def _site_depths(bam: str, ref: str | None, chrom: str, positions: list[int]) -> dict[int, int] | None:
    """Per-position depth (BQ>=13, as bcftools mpileup) at the table sites; None if unavailable.

    One indexed region query spanning the sites (samtools depth -b would scan the whole BAM).
    """
    if not bam or not os.path.exists(bam) or not positions:
        return None
    reg = f"{chrom}:{min(positions)}-{max(positions)}"
    if ref:   # same read filters as bcftools mpileup (BAQ, BQ>=13, orphans/dups dropped): depth = evidence used
        cmd = ["samtools", "mpileup", "-a", "-Q", "13", "-d", "500", "-f", ref, "-r", reg]
    else:
        cmd = ["samtools", "depth", "-a", "-q", "13", "-r", reg]
    try:
        proc = subprocess.run(cmd + [bam], capture_output=True, text=True, timeout=600)
    except Exception:
        return None
    if proc.returncode != 0:
        return None
    want = set(positions)
    d = {}
    for ln in proc.stdout.splitlines():
        f = ln.split("\t")
        col = 3 if ref else 2
        if len(f) > col and int(f[1]) in want:
            d[int(f[1])] = int(f[col])
    return d


# CPIC alleles whose defining sites must be covered before a G6PD negative is reported:
# A- (202A/376G), Mediterranean c.563C>T, Viangchan c.871G>A, Mahidol c.487G>A,
# Canton c.1376G>T, Kaiping c.1388G>A, Union c.1360C>T, Chinese-5 c.1024C>T.
_G6PD_KEY_ALLELES = (("A- 202A_376G", "A-"), ("Mediterranean", "c.563C>T"),
                     ("Viangchan", "c.871G>A"), ("Mahidol", "c.487G>A"),
                     ("Canton", "c.1376G>T"), ("Kaiping", "c.1388G>A"),
                     ("Union", "c.1360C>T"), ("Chinese-5", "c.1024C>T"))


def _g6pd_key_sites(table: dict) -> dict[int, str]:
    """pos -> label of the sites defining the common deficient alleles."""
    out: dict[int, str] = {}
    for a in table["alleles"]:
        for prefix, label in _G6PD_KEY_ALLELES:
            if a["name"].startswith(prefix):
                for p in a["sites"]:
                    out[int(p)] = f"{label} (pos {p})" if label == "A-" else label
    return out


G6PD_MIN_QUAL, G6PD_MIN_DP = 20, 3   # a single-read / low-QUAL record is never a Deficient call


def _trim_common_suffix(ref: str, alt: str) -> tuple[str, str]:
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    return ref, alt


def g6pd_call(records: list[tuple], table: dict, depth: dict[int, int] | None,
              sex: str | None) -> dict[str, Any]:
    """Pure G6PD typing. records: (pos, ref, [alts], gt_alleles_list) from the VCF.

    Returns dict(diplotype, phenotype, alleles, unassigned, n_sites, n_typed, untyped, sex).
    """
    # >= 5 reads per X copy: a hemizygous male needs half the depth of a female (or unknown sex),
    # whose heterozygous allele is in only ~half the reads (P(no alt read) = 0.5^n).
    min_cov = int(table.get("min_cov", 10))
    if sex == "MALE":
        min_cov = int(table.get("min_cov_male", max(1, min_cov // 2)))
    sites = {s["pos"]: s for s in table["sites"]}
    n_sites = len(sites)
    if depth is None:
        untyped, n_typed = [], None
    else:
        untyped = sorted(p for p in sites if depth.get(p, 0) < min_cov)
        n_typed = n_sites - len(untyped)
    male = sex == "MALE"

    # observed alt copies at table sites, keyed (pos, alt). A call that cannot be trusted is never
    # counted (no positive from it) and never ignored (no negative over it): it goes to `qc`.
    obs: dict[tuple, int] = {}
    qc: list[str] = []
    male_het: list[str] = []
    haploid = False   # a haploid chrX GT ("1") can only come from a hemizygous male
    for rec in records:
        pos, ref, alts, gts = rec[:4]
        qual, dp = (rec[4], rec[5]) if len(rec) >= 6 else (None, None)
        s = sites.get(pos)
        if not s:
            continue
        for i, alt in enumerate(alts, start=1):
            c = sum(1 for g in gts if g == i)
            if not c:
                continue
            r2, a2 = _trim_common_suffix(ref, alt)   # bcftools may pad REF/ALT with trailing context
            if r2 != s["ref"] or a2 not in s["alts"]:
                qc.append(f"{pos}{ref}>{alt} does not match the CPIC site {s['ref']}>{'/'.join(s['alts'])}")
                continue
            if (qual is not None and qual < G6PD_MIN_QUAL) or (dp is not None and dp < G6PD_MIN_DP):
                qc.append(f"{pos}{a2} low-quality call (QUAL {qual}, DP {dp})")
                continue
            if male and c < len(gts):
                male_het.append(f"{pos}{a2}")
                qc.append(f"{pos}{a2} heterozygous on chrX in a male (artifact, contamination or wrong sex)")
                continue
            haploid = haploid or len(gts) == 1
            obs[(pos, a2)] = 1 if male else c

    matched = []
    for a in table["alleles"]:
        keys = {(int(p), alt) for p, alt in a["sites"].items()}
        if keys <= obs.keys():
            matched.append((a, keys))
    # most specific allele first; each takes the copies still left at all its sites, so a
    # nested allele is called only on leftover copies (1/1 202A + 0/1 376G = A- + Asahi).
    left = dict(obs)
    best = []
    for a, k in sorted(matched, key=lambda m: -len(m[1])):
        n = min(left[x] for x in k)
        if n:
            best.append((a, k, n))
            for x in k:
                left[x] -= n
    unassigned = sorted(x for x, n in left.items() if n)

    def fclass(a):
        f = a["function"]
        return "CNSHA" if f.startswith("I/") else "Deficient" if "Deficient" in f else \
            "Normal" if f.endswith("Normal") else "Uncertain"

    called = []
    for a, k, copies in best:
        called.append({"name": a["name"], "function": a["function"], "copies": copies,
                       "cls": fclass(a), "sites": sorted(p for p, _ in k), "keys": sorted(k)})

    out = {"alleles": called, "unassigned": [f"{p}{a}" for p, a in unassigned], "male_het": male_het, "qc": qc,
           "n_sites": n_sites, "n_typed": n_typed, "untyped": untyped,
           "sex": sex or "unknown"}
    sex_txt = {"MALE": "male (hemizygous)", "FEMALE": "female"}.get(sex or "", "sex not determined")
    cov_txt = (f"{n_typed} of {n_sites} CPIC G6PD allele-defining sites typed at >= {min_cov}x"
               if n_typed is not None else f"{n_sites} CPIC G6PD sites queried, depth not verified")

    if n_typed == 0:
        out["diplotype"] = "not typed"
        out["phenotype"] = f"G6PD not typed - no CPIC site reached {min_cov}x"
        return out

    # ERROR is never NEGATIVE: a negative needs verified depth, incl. at the key sites.
    negative = not unassigned and not qc and all(x["cls"] in ("Normal",) for x in called)
    if negative and depth is None:
        out["diplotype"] = "not typed"
        out["phenotype"] = ("G6PD not typed - read depth at the CPIC sites could not be checked "
                            "(BAM unavailable or samtools depth failed)")
        return out
    if negative:
        keys = _g6pd_key_sites(table)
        low = [f"{lab} {depth.get(p, 0)}x" for p, lab in sorted(keys.items()) if depth.get(p, 0) < min_cov]
        if low:
            out["diplotype"] = "not typed"
            out["phenotype"] = (f"G6PD not typed (key site(s) below {min_cov}x: " + "; ".join(low) + ")")
            out["low_key_sites"] = low
            return out

    if qc and not called and not unassigned:
        out["diplotype"] = "not typed"
        out["phenotype"] = "Indeterminate - G6PD call(s) failed QC: " + "; ".join(qc)
        return out
    if not called and not unassigned:
        out["diplotype"] = "Reference (no CPIC G6PD allele-defining variant at typed sites)"
        out["phenotype"] = (f"No G6PD deficiency variant found among the sites checked ({cov_txt}; "
                            f"{sex_txt}); not a measurement of enzyme activity")
        return out

    tot = {c: sum(x["copies"] for x in called if x["cls"] == c)
           for c in ("CNSHA", "Deficient", "Normal", "Uncertain")}
    n_def = tot["CNSHA"] + tot["Deficient"]
    names = [f"{x['name']} ({x['function']})" for x in called]
    if unassigned and not n_def:
        pheno = ("Indeterminate - variant(s) at CPIC G6PD site(s) not forming a defined allele: "
                 + ", ".join(out["unassigned"]))
    elif tot["Uncertain"] and not n_def:
        pheno = "Indeterminate - allele of uncertain function: " + "; ".join(names)
    elif n_def:
        cn = "with CNSHA" if tot["CNSHA"] else ""
        if sex == "FEMALE":
            two = n_def >= 2
            pheno = (f"G6PD Deficient {cn}".strip() + " (female, two deficient alleles)") if two else \
                    ("G6PD Variable (female heterozygous for a deficient allele; "
                     "CPIC: enzyme activity test recommended)")
        elif sex == "MALE":
            pheno = f"G6PD Deficient {cn}".strip() + " (hemizygous male)"
        else:  # sex unknown
            hom = n_def >= 2 or haploid
            pheno = (f"G6PD Deficient {cn}".strip() + " (1/1 call: hemizygous male or homozygous female)") \
                if hom else ("G6PD Variable if female / Deficient if male (sex not determined, "
                             "single deficient allele copy)")
        pheno += "; " + "; ".join(names) + " - rasburicase/pegloticase contraindicated; avoid primaquine, dapsone"
    else:  # only normal-function alleles (e.g. A, class IV)
        pheno = ("No G6PD deficiency variant found among the sites checked (" + cov_txt +
                 f"); allele(s) of normal function: " + "; ".join(names))
    if unassigned and n_def:
        pheno += f"; additional unassigned CPIC-site variant(s): {', '.join(out['unassigned'])}"
    parts = []
    for x in called:
        z = "hemizygous" if sex == "MALE" else {1: "het", 2: "hom"}.get(x["copies"], str(x["copies"]))
        parts.append(f"{x['name']}({z})")
    out["diplotype"] = "+".join(parts) if parts else "unassigned: " + ", ".join(out["unassigned"])
    out["diplotype"] += " [X-linked, unphased]"
    if qc:
        if n_def:   # a confident deficiency stays the leading phrase (actionable); QC is a note
            pheno += "; QC: " + "; ".join(qc)
        else:       # no confident deficiency: a negative or normal-only call cannot stand over a failed site
            pheno = "Indeterminate - G6PD call(s) failed QC: " + "; ".join(qc) + "; " + pheno
    out["phenotype"] = pheno
    return out


def parse_g6pd_vcf(output_dir: str, gene: str, sample: str) -> CallerResult:
    """Type G6PD from <gene>.vcf.gz against the full CPIC allele table (see above)."""
    result = CallerResult(
        tool="VCF-Check",
        phasing_method="bcftools mpileup/call + CPIC G6PD allele table (171 sites, unphased, X-linked)",
    )
    table = _VCF_CHECK_VARIANTS.get("G6PD")
    if not table:
        result.status = "failed"
        result.diplotype = "G6PD allele table missing"
        return result
    vcf = os.path.join(output_dir, f"{gene}.vcf.gz")
    if not os.path.exists(vcf):
        result.status = "failed"
        result.diplotype = "VCF not found"
        return result
    lo = min(s["pos"] for s in table["sites"])
    hi = max(s["pos"] for s in table["sites"])
    try:
        proc = subprocess.run(["bcftools", "view", "-H", "-r", f"{table['chrom']}:{lo}-{hi}", vcf],
                              capture_output=True, text=True, timeout=60)
    except FileNotFoundError:
        result.status, result.diplotype = "failed", "bcftools not in PATH"
        return result
    except subprocess.TimeoutExpired:
        result.status, result.diplotype = "failed", "bcftools timed out"
        return result
    if proc.returncode != 0:
        result.status, result.diplotype = "failed", "bcftools view failed"
        return result

    records, fmt_by_pos = [], {}
    for line in proc.stdout.splitlines():
        f = line.split("\t")
        if len(f) < 10:
            continue
        fmt = dict(zip(f[8].split(":"), f[9].split(":")))
        gts = [int(g) for g in re.split(r"[/|]", fmt.get("GT", "."))
               if g.isdigit()]
        try:
            qual = float(f[5])
        except ValueError:
            qual = None
        dp = int(fmt["DP"]) if fmt.get("DP", "").isdigit() else None
        records.append((int(f[1]), f[3], f[4].split(","), gts, qual, dp))
        fmt_by_pos[int(f[1])] = (fmt, [_trim_common_suffix(f[3], a)[1] for a in f[4].split(",")])

    sex = (_RUN_CTX.get("sex") or "").upper() or None
    if sex not in ("MALE", "FEMALE"):
        sex = _infer_sex_from_bam(_RUN_CTX["bam"]) if _RUN_CTX.get("bam") else None
    depth = _site_depths(_RUN_CTX.get("bam"), _RUN_CTX.get("ref"), table["chrom"],
                         [s["pos"] for s in table["sites"]])
    call = g6pd_call(records, table, depth, sex)

    result.status = "ok"
    result.diplotype = call["diplotype"]
    result.phenotype = call["phenotype"]
    result.haplotype1 = call["alleles"][0]["name"] if call["alleles"] else "-"
    site_by_pos = {s["pos"]: s for s in table["sites"]}
    for x in call["alleles"]:
        for p, alt in x["keys"]:
            fm, alts = fmt_by_pos.get(p, ({}, []))
            try:
                ad = [int(v) for v in fm["AD"].split(",")]
                af = f"{ad[alts.index(alt) + 1] / sum(ad):.3f}"
            except (KeyError, ValueError, IndexError, ZeroDivisionError):
                af = "-"
            result.supporting_variants.append({
                "allele": x["name"], "chrom": "chrX", "pos": str(p),
                "ref": site_by_pos[p]["ref"],
                "alt": alt, "af": af, "depth": fm.get("DP", "-"),
                "effect": x["function"], "rsid": site_by_pos[p].get("rsid") or "-"})
    result.allele_score = (f"sex={call['sex']}; sites typed "
                           f"{call['n_typed'] if call['n_typed'] is not None else '?'}/{call['n_sites']}")
    return result


# ── Parser: RYR1 VCF-Check (CPIC malignant-hyperthermia alleles) ─────────────
# RYR1 has NO star-allele system; CPIC defines MH risk by a curated list of
# pathogenic / likely-pathogenic missense variants. The star-allele callers
# (PyPGx/Stargazer/Aldy) are the wrong instrument here and routinely miss real
# MH variants the array reports (validation: P_2000 c.1021G>A missed by all 3).
# This VCF-Check is the authoritative RYR1 caller; the 93-variant table lives in
# docker/vcf_check_variants.json (CPIC api.cpicpgx.org + Ensembl-resolved coords),
# NOT hand-embedded — so it is auditable and updatable. GRCh38 / chr19 / +strand.
def parse_ryr1_vcf(output_dir: str, gene: str, sample: str) -> CallerResult:
    """Query RYR1.vcf.gz for CPIC MH-associated variants (data-file driven)."""
    result = CallerResult(
        tool="VCF-Check",
        phasing_method="bcftools mpileup — direct CPIC MH-variant query (unphased)",
    )
    entry = _VCF_CHECK_VARIANTS.get("RYR1", {})
    table = entry.get("variants", [])
    if not table:
        result.status = "failed"
        result.diplotype = "RYR1 variant table missing"
        return result

    vcf = os.path.join(output_dir, f"{gene}.vcf.gz")
    if not os.path.exists(vcf):
        result.status = "failed"
        result.diplotype = "VCF not found"
        return result

    by_pos: dict[int, list] = {}
    for v in table:
        by_pos.setdefault(v["pos"], []).append(v)
    chrom = entry.get("chrom", "chr19")
    region = f"{chrom}:{min(by_pos)}-{max(by_pos)}"
    try:
        proc = subprocess.run(
            ["bcftools", "view", "-r", region, vcf],
            capture_output=True, text=True, timeout=60,
        )
    except FileNotFoundError:
        result.status = "failed"
        result.diplotype = "bcftools not in PATH"
        return result
    except subprocess.TimeoutExpired:
        result.status = "failed"
        result.diplotype = "bcftools timed out"
        return result

    found: list[tuple] = []
    variants_list: list[dict] = []
    for line in proc.stdout.splitlines():
        if line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 10:
            continue
        pos = int(fields[1])
        if pos not in by_pos:
            continue
        ref, alt_str = fields[3], fields[4]
        alts = alt_str.split(",")
        fmt = dict(zip(fields[8].split(":"), fields[9].split(":")))
        gt = fmt.get("GT", "./.").replace("|", "/")
        dp, ad = fmt.get("DP", "-"), fmt.get("AD", "-")
        copies = sum(1 for g in gt.split("/") if g not in ("0", "."))
        if copies == 0:
            continue
        for adef in by_pos[pos]:
            if adef["alt"] not in alts:
                continue
            zyg = "hom" if copies >= 2 else "het"
            found.append((adef["hgvs_c"], zyg, adef.get("status", "")))
            variants_list.append({
                "allele": adef["hgvs_c"], "chrom": chrom, "pos": str(pos),
                "ref": ref, "alt": adef["alt"], "af": ad, "depth": dp,
                "effect": adef.get("note", ""), "rsid": adef.get("rsid") or "-",
            })

    result.status = "ok"
    result.supporting_variants = variants_list
    if not found:
        result.diplotype = "*1/*1"
        result.phenotype = (f"No CPIC MH-associated RYR1 variant found among the {len(table)} substitution "
                            "sites checked (non-substitution alleles are not checked)")
        return result

    parts = [f"{h}({z})" for h, z, _ in found]
    result.diplotype = "+".join(parts) + " [unphased]"
    worst = "Pathogenic" if any(s == "Pathogenic" for _, _, s in found) else "Likely Pathogenic"
    result.phenotype = (f"RYR1 MH-associated variant detected ({worst}) — malignant hyperthermia "
                        "susceptibility; avoid volatile anaesthetics and succinylcholine")
    return result


# ── Parser: UGT1A1 VCF-Check (promoter + coding SNP supplement) ───────────────
# Star-allele callers inconsistently report UGT1A1 *60 (rs45530432) and *80
# (rs887829) — East Asian promoter variants that compound the activity of *28.
# *6 (rs4148323) is missed by some tools in East Asian samples.
# This parser reads the gene VCF produced by bcftools and checks each position
# directly, providing an independent fifth call for UGT1A1.
_UGT1A1_SNP_ALLELES: list[dict] = [
    # (GRCh38 pos, REF, ALT, star-allele, rsID, functional note)
    {"pos": 233757013, "ref": "T", "alt": "G", "allele": "*60", "rsid": "rs45530432",
     "note": "Promoter -3279T>G; reduces UGT1A1 expression; common East Asian allele"},
    {"pos": 233759924, "ref": "C", "alt": "T", "allele": "*80", "rsid": "rs887829",
     "note": "Promoter -401C>T; moderate expression reduction"},
    {"pos": 233760498, "ref": "G", "alt": "A", "allele": "*6",  "rsid": "rs4148323",
     "note": "c.211G>A (p.Gly71Arg); decreased function; common in East Asians"},
]
_UGT1A1_INDEL_POS = 233760233  # *28 TA-repeat locus (rs3064744); detected as indel


def parse_ugt1a1_vcf(output_dir: str, gene: str, sample: str) -> CallerResult:
    """Query UGT1A1.vcf.gz for *28, *60, *80, *6 alleles independently of star-allele callers."""
    result = CallerResult(
        tool="VCF-Check",
        phasing_method="bcftools mpileup — direct variant query (unphased)",
    )

    vcf = os.path.join(output_dir, f"{gene}.vcf.gz")
    if not os.path.exists(vcf):
        result.status = "failed"
        result.diplotype = "VCF not found"
        return result

    # Query the promoter + coding region spanning all target alleles
    region = "chr2:233757013-233760498"
    try:
        proc = subprocess.run(
            ["bcftools", "view", "-r", region, vcf],
            capture_output=True, text=True, timeout=30,
        )
    except FileNotFoundError:
        result.status = "failed"
        result.diplotype = "bcftools not in PATH"
        return result
    except subprocess.TimeoutExpired:
        result.status = "failed"
        result.diplotype = "bcftools timed out"
        return result

    # Parse VCF lines for target positions
    # alleles_found: star-allele -> (gt string, depth, alt_copies, rsid, note)
    alleles_found: dict[str, tuple] = {}
    variants_list: list[dict] = []

    for line in proc.stdout.splitlines():
        if line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 10:
            continue
        pos = int(fields[1])
        ref, alt_str = fields[3], fields[4]
        fmt_keys = fields[8].split(":")
        fmt_vals = fields[9].split(":")
        fmt = dict(zip(fmt_keys, fmt_vals))
        gt_raw = fmt.get("GT", "./.")
        gt = gt_raw.replace("|", "/")
        dp  = fmt.get("DP", "-")
        ad  = fmt.get("AD", "-")
        alt_copies = sum(1 for g in gt.split("/") if g not in ("0", "."))
        if alt_copies == 0:
            continue

        # Check *28 TA-repeat indel (any indel at this position is *28)
        if pos == _UGT1A1_INDEL_POS and len(alt_str) != len(ref):
            alleles_found["*28"] = (gt, dp, alt_copies, "rs3064744",
                                    "TA-repeat insertion; *28 = 7 TA (Gilbert's); "
                                    f"REF={ref[:16]}… ALT={alt_str[:16]}…")
            variants_list.append({
                "allele": "*28", "chrom": "chr2", "pos": str(pos),
                "ref": ref[:20], "alt": alt_str[:20],
                "af": ad, "depth": dp, "effect": "promoter TA-repeat expansion",
                "rsid": "rs3064744",
            })
            continue

        # Check SNP alleles
        alts = alt_str.split(",")
        for adef in _UGT1A1_SNP_ALLELES:
            if adef["pos"] != pos:
                continue
            if adef["alt"] not in alts:
                continue
            alleles_found[adef["allele"]] = (gt, dp, alt_copies, adef["rsid"], adef["note"])
            variants_list.append({
                "allele": adef["allele"], "chrom": "chr2", "pos": str(pos),
                "ref": ref, "alt": adef["alt"],
                "af": ad, "depth": dp, "effect": adef["note"],
                "rsid": adef["rsid"],
            })

    result.status = "ok"
    result.supporting_variants = variants_list

    if not alleles_found:
        result.diplotype = "*1/*1"
        result.phenotype = ("No UGT1A1 variant found among the sites checked "
                            "(*28 TA-repeat, *60, *80, *6); not a Normal Function call")
        return result

    # Report as compound allele list with genotype (unphased; phasing requires LRS).
    # Format: *28(het)+*60(hom)+*80(het)
    # Sorted by allele number so output is deterministic.
    def _sort_key(name: str) -> int:
        m = re.search(r"(\d+)", name)
        return int(m.group(1)) if m else 999

    parts = []
    for aname in sorted(alleles_found, key=_sort_key):
        gt, dp, copies, rsid, note = alleles_found[aname]
        zygosity = "hom" if copies >= 2 else "het"
        parts.append(f"{aname}({zygosity})")
    result.diplotype = "+".join(parts) + " [unphased]"

    # Coarse phenotype
    risk = set(alleles_found)
    is_28_hom = alleles_found.get("*28", ("",))[0] in ("1/1", "1|1")
    if "*28" in risk and is_28_hom and len(risk) == 1:
        result.phenotype = "Poor Function (Gilbert's — *28/*28)"
    elif risk:
        result.phenotype = "Reduced Function (variant allele(s) detected)"
    else:  # unreachable: alleles_found is non-empty here
        result.phenotype = "No UGT1A1 variant found among the sites checked"

    return result


# ── Parser: single-SNP VCF-Check (VKORC1, IFNL3) ─────────────────────────────
# These genes' clinical call IS one SNP (VKORC1 rs9923231 / IFNL3 rs12979860);
# the star-allele callers wrap it in incompatible nomenclature (*H/*S/rsID). Read
# the SNP genotype straight from the VCF. bcftools call -v emits only variant
# sites, so absence of the record at a covered locus = reference homozygous.
# Reported in the clinical sense (report_ref/report_alt) to match the array.
def parse_single_snp_vcf(output_dir: str, gene: str, sample: str) -> CallerResult:
    result = CallerResult(
        tool="VCF-Check",
        phasing_method="bcftools mpileup — single-SNP genotype",
    )
    entry = _VCF_CHECK_VARIANTS.get(gene, {})
    snp = entry.get("snp")
    if not snp:
        result.status = "failed"
        result.diplotype = f"{gene} SNP def missing"
        return result
    vcf = os.path.join(output_dir, f"{gene}.vcf.gz")
    if not os.path.exists(vcf):
        result.status = "failed"
        result.diplotype = "VCF not found"
        return result
    chrom, pos = entry.get("chrom", "chr1"), snp["pos"]
    try:
        proc = subprocess.run(
            ["bcftools", "view", "-r", f"{chrom}:{pos}-{pos}", vcf],
            capture_output=True, text=True, timeout=30,
        )
    except FileNotFoundError:
        result.status = "failed"
        result.diplotype = "bcftools not in PATH"
        return result
    except subprocess.TimeoutExpired:
        result.status = "failed"
        result.diplotype = "bcftools timed out"
        return result

    rref, ralt = snp["report_ref"], snp["report_alt"]
    copies, dp = 0, "-"
    for line in proc.stdout.splitlines():
        if line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 10 or int(f[1]) != pos or snp["alt"] not in f[4].split(","):
            continue
        fmt = dict(zip(f[8].split(":"), f[9].split(":")))
        gt = fmt.get("GT", "./.").replace("|", "/")
        dp = fmt.get("DP", "-")
        copies = sum(1 for g in gt.split("/") if g not in ("0", "."))
        break

    zyg = "hom" if copies >= 2 else ("het" if copies == 1 else "ref")
    result.diplotype = {"hom": f"{ralt}/{ralt}", "het": f"{rref}/{ralt}",
                        "ref": f"{rref}/{rref}"}[zyg]
    result.haplotype1, result.haplotype2 = result.diplotype.split("/")
    result.phenotype = snp.get("phenotypes", {}).get(zyg, "-")
    if copies:
        result.supporting_variants = [{
            "allele": snp.get("rsid"), "chrom": chrom, "pos": str(pos),
            "ref": snp["ref"], "alt": snp["alt"], "af": "-", "depth": dp,
            "effect": snp.get("note", ""), "rsid": snp.get("rsid", "-"),
        }]
    result.status = "ok"
    return result


# ── Parser: ABO blood type (BAM pileup at the defining sites) ─────────────────
# abo_pileup.py (run by the bcftools rule for ABO) writes genes/ABO/abo_pileup.json
# with read counts at the c.261delG O site and four A/B SNPs. GRCh38 carries the O
# allele, so hom-ref at c.261 is O/O. No pileup file or no coverage -> "Not typed".
# Also writes results/<S>/abo/abo_type.tsv (+ .json, with per-site evidence).
def parse_abo_vcf(output_dir: str, gene: str, sample: str) -> CallerResult:
    result = CallerResult(
        tool="ABO-Typer",
        phasing_method="BAM pileup at c.261delG + 4 A/B SNPs (unphased; GRCh38 ref = O)",
    )
    try:
        defs = _oa.load_abo_defs()
    except Exception as exc:  # noqa: BLE001
        result.status = "failed"
        result.diplotype = f"ABO defs missing: {exc}"
        return result

    pileup = None
    pj = os.path.join(output_dir, "abo_pileup.json")
    if os.path.exists(pj):
        try:
            with open(pj) as fh:
                pileup = json.load(fh)
        except Exception:  # noqa: BLE001
            pileup = None
    call = _oa.abo_assign_type(pileup, defs)
    if pileup is None:
        call["flags"].append("abo_pileup.json missing/unreadable (ABO read-level step did not run)")

    sample_root = os.path.dirname(os.path.dirname(output_dir))  # .../genes/ABO -> sample root
    abo_dir = os.path.join(sample_root, "abo")
    try:
        _oa.write_abo_tsv(os.path.join(abo_dir, "abo_type.tsv"), sample, call)
        with open(os.path.join(abo_dir, "abo_type.json"), "w") as fh:
            json.dump({"sample": sample, **call}, fh, indent=2)
    except Exception:  # noqa: BLE001  (contract write is best-effort)
        pass

    typed = call["ABO_type"] != "Not typed"
    result.diplotype = call["alleles"] if typed else "not typed"
    result.phenotype = f"ABO blood group {call['ABO_type']}" if typed else "ABO not typed"
    result.allele_score = f"confidence={call['confidence']}; sites_observed={call['n_sites_observed']}" + (
        "; " + "; ".join(call["flags"]) if call["flags"] else "")
    result.supporting_variants = [
        {"allele": s["rsid"], "chrom": "chr9", "pos": s.get("pos", "-"), "ref": "-", "alt": "-",
         "af": s.get("af") if s.get("af") is not None else "-", "depth": s.get("depth"),
         "effect": f"{s.get('cdna', '')} ({s.get('role', '')}) alt={s.get('n_alt')}/ref={s.get('n_ref')}",
         "rsid": s["rsid"]}
        for s in call["per_site"]
    ]
    result.status = "ok" if typed else "failed"
    return result


# ── Parser: arcasHLA (HLA class II — DQA1/DQB1/DRB1) ──────────────────────────
# arcasHLA types all HLA genes once per sample into hla2/<sample>.genotype.json.
# The per-gene compare reads that shared JSON (sample-level dependency, like
# PharmCAT) and emits the SAME *_comparison.tsv schema OptiType uses for class I.
def parse_arcashla(output_dir: str, gene: str, sample: str) -> CallerResult:
    result = CallerResult(
        tool="arcasHLA",
        phasing_method="arcasHLA — read-partitioned HLA genotyping (class II)",
    )
    # hla2/<sample>.genotype.json is sample-level; search the gene dir, sample root.
    sample_root = os.path.dirname(os.path.dirname(output_dir))
    candidates = [
        os.path.join(sample_root, "hla2", f"{sample}.genotype.json"),
        os.path.join(output_dir, "hla2", f"{sample}.genotype.json"),
    ]
    candidates += glob.glob(os.path.join(sample_root, "hla2", "*.genotype.json"))
    path = next((p for p in candidates if os.path.isfile(p)), None)
    if not path:
        # Distinguish "never ran" from "ran and failed": the Snakefile `hla2` rule
        # writes hla2/hla2.status (exit code). A non-zero code with no genotype JSON
        # is a FAILURE, never a silent not_run (failure is not a clean result).
        for sf in (os.path.join(sample_root, "hla2", "hla2.status"),
                   os.path.join(output_dir, "hla2", "hla2.status")):
            if os.path.isfile(sf):
                rc = open(sf).read().strip()
                if rc not in ("", "0"):
                    result.status = "failed"
                    result.diplotype = (f"arcasHLA typing failed (hla2.status={rc}; "
                                        "see pgx-hla2.sh log / hla2/arcashla_*.log)")
                    return result
        return result  # not_run
    try:
        geno = json.load(open(path))
    except Exception as exc:  # noqa: BLE001
        result.status = "failed"
        result.diplotype = f"parse error: {exc}"
        return result
    haps = _oa.parse_arcashla_genotype(geno, gene)
    if not haps:
        result.status = "failed"
        result.diplotype = f"{gene} absent from arcasHLA output"
        return result
    result.haplotype1, result.haplotype2 = haps
    result.diplotype = f"{haps[0]}/{haps[1]}"
    result.status = "ok"
    return result


# ── VCF-Check dispatch table ───────────────────────────────────────────────────
# Maps gene name → gene-specific VCF parser function.
# Add new genes here; the main() loop calls the appropriate function via
# support.get("vcf_check") + _VCF_CHECK_PARSERS lookup.
_VCF_CHECK_PARSERS: dict = {
    "UGT1A1":  parse_ugt1a1_vcf,
    "CACNA1S": parse_cacna1s_vcf,
    "G6PD":    parse_g6pd_vcf,
    "RYR1":    parse_ryr1_vcf,
    "VKORC1":  parse_single_snp_vcf,
    "IFNL3":   parse_single_snp_vcf,
    "CYP2C-CLUSTER": parse_single_snp_vcf,  # rs12777823 (CPIC warfarin)
    "ABO":     parse_abo_vcf,               # PROVISIONAL — UNVALIDATED
}

# Genes resolved by a single authoritative caller rather than a star-allele vote:
#   VCF-Check for single-SNP / curated-variant-list genes (star-allele callers use
#   incompatible nomenclature there), and Cyrius for CYP2D6 (purpose-built for the
#   2D6/2D7 paralog + SV). When the authority produces a call it IS the verdict;
#   if it declines (no ok result) compute_verdict falls back to the consensus.
AUTHORITATIVE = {
    "RYR1": "VCF-Check", "CACNA1S": "VCF-Check", "G6PD": "VCF-Check",
    "VKORC1": "VCF-Check", "IFNL3": "VCF-Check", "CYP2D6": "Cyrius",
    "CYP2C-CLUSTER": "VCF-Check",  # rs12777823 single SNP (CPIC warfarin)
    "ABO": "VCF-Check",            # provisional ABO typer (UNVALIDATED)
    # PharmCAT (CPIC reference impl) is the authority where star-allele callers
    # genuinely disagree on a CPIC gene: UGT1A1 *28/*80 LD, CYP2B6, CYP4F2.
    "UGT1A1": "PharmCAT", "CYP2B6": "PharmCAT", "CYP4F2": "PharmCAT",
}


# ── SV mode note ──────────────────────────────────────────────────────────────
def _sv_note(gene: str) -> str:
    notes = []
    if gene in PYPGX_SV_GENES:
        notes.append("PyPGx: depth-of-coverage + VDR control stats")
    if gene in STARGAZER_SV_GENES:
        notes.append("Stargazer: GDF from BAM (paralog CN normalisation)")
    return ("SV mode — " + "; ".join(notes)) if notes else "no SVs expected"


# ── Output ────────────────────────────────────────────────────────────────────
def _normalize_diplotype(d: str) -> str:
    """Order-insensitive diplotype so *2/*4 == *4/*2."""
    if not d or d == "-":
        return "-"
    if "/" not in d:
        return d
    return "/".join(sorted(p.strip() for p in d.split("/")))


# ── Clinical phenotype tier ───────────────────────────────────────────────────
# When star-allele callers print DIFFERENT diplotype strings (tool nomenclature:
# Stargazer *S-codes, Aldy rsID tokens, or differing haplotype phasing) but the
# phenotype-emitting callers agree on the CPIC phenotype, the call is clinically
# concordant — CPIC reports by phenotype / gene-activity-score, not by haplotype
# string. "Indeterminate"/unknown phenotypes count as ABSTAIN, not dissent (Aldy
# routinely emits Indeterminate for DPYD even when it called the variants).
_PHENO_ABSTAIN = {"", "-", "indeterminate", "unknown", "n/a", "na", "no result",
                  "not available", "no call", "none", "no function"}
_PHENO_ALIASES = {"extensive": "normal", "rapid metaboliser": "rapid"}


def _norm_phenotype(p: str) -> str | None:
    """Normalise a tool phenotype to a comparable token, or None if it abstains."""
    s = (p or "").strip().lower().replace("_", " ")
    s = re.sub(r"\s+", " ", s).strip()
    if s in _PHENO_ABSTAIN:
        return None
    s = re.sub(r"\s*metaboli[sz]er$", "", s).strip()
    if s in _PHENO_ABSTAIN or not s:
        return None
    s = re.sub(r"^(likely|possible|probable)\s+", "", s)
    s = {"ultra rapid": "ultrarapid"}.get(s, s)
    return _PHENO_ALIASES.get(s, s)


def _phenotype_consensus(results: list[CallerResult]):
    """If >=2 phenotype-emitting callers agree on a single normalised phenotype
    (abstainers ignored), return (display_phenotype, n_agree, display_diplotype);
    else None. Used to rescue a diplotype-STRING discordance/majority that is in
    fact a single clinical phenotype."""
    voting = [(r, _norm_phenotype(r.phenotype)) for r in results
              if r.status == "ok" and r.diplotype not in ("-", "")]
    voting = [(r, p) for r, p in voting if p is not None]
    if len(voting) < 2 or len({p for _, p in voting}) != 1:
        return None
    pheno_disp = Counter(r.phenotype for r, _ in voting).most_common(1)[0][0]
    dip_disp = Counter(r.diplotype for r, _ in voting).most_common(1)[0][0]
    return pheno_disp, len(voting), dip_disp


# ── NAT2: CPIC phenotype from CPIC's own allele function values ───────────────
# PyPGx has no NAT2 diplotype table (reports "Indeterminate") and Stargazer's
# label is not a CPIC term, so no caller yields a CPIC phenotype. Derived here
# from the CPIC NAT2 allele function table and diplotype->phenotype pairing:
# increased+increased = Rapid, increased+decreased = Intermediate,
# decreased+decreased = Poor. An allele CPIC marks uncertain/unknown, or does
# not list, leaves the phenotype Indeterminate (never a guess).
_NAT2_INCREASED = frozenset({"*1", "*4"})
_NAT2_DECREASED = frozenset({"*5", "*6", "*7", "*14", "*15", "*16", "*19", "*29", "*30", "*31", "*32", "*33", "*34", "*35", "*36", "*37", "*38", "*39", "*40", "*46", "*47", "*48", "*49", "*50", "*52", "*53", "*56", "*58", "*59", "*62", "*63", "*64", "*66", "*67", "*68", "*69", "*70", "*72", "*73", "*75"})


def nat2_cpic_phenotype(diplotype: str) -> str:
    """CPIC NAT2 phenotype for an 'a/b' star diplotype, else 'Indeterminate'."""
    parts = [p.strip() for p in (diplotype or "").split("/")]
    if len(parts) != 2:
        return "Indeterminate"
    fn = []
    for p in parts:
        if p in _NAT2_INCREASED:
            fn.append("inc")
        elif p in _NAT2_DECREASED:
            fn.append("dec")
        else:
            return "Indeterminate"
    if fn == ["inc", "inc"]:
        return "Rapid Metabolizer"
    if fn == ["dec", "dec"]:
        return "Poor Metabolizer"
    return "Intermediate Metabolizer"


def _compute_verdict_core(results: list[CallerResult], no_call: bool,
                    gene: str | None = None) -> dict[str, Any]:
    """Single source of truth for the gene-level call.

    Emitted into detail.json so pgx-report.py and the batch summary READ this
    rather than each recomputing concordance (which previously diverged and
    silently resolved 2-vs-2 ties as a winner).

    Agreement is counted on the SYNONYM-COLLAPSED (canonical) diplotype so that
    e.g. DPYD `*9A` / `*S10` / `c.85T>C (*9A)` count as one call, not three
    discordant ones (design A.1). The displayed `consensus_diplotype` stays the
    most common RAW spelling among the agreeing tools (human-readable). For a
    gene with no synonym entry, canonicalization reduces to the old
    order-insensitive normalize, so its verdict is unchanged. `reconciled` is
    True when synonym collapse genuinely raised the agreement count.

    status: no_call | no_data | discordant | majority | concordant
    """
    if no_call:
        return {"consensus_diplotype": "NO_CALL", "consensus_phenotype": "-",
                "n_agree": 0, "n_called": 0, "status": "no_call", "reconciled": False}

    # Authoritative-caller genes: the designated caller (VCF-Check or Cyrius) is
    # the verdict when it produced a call; otherwise fall through to the vote.
    auth_tool = AUTHORITATIVE.get(gene)
    if auth_tool:
        a = next((r for r in results if r.tool == auth_tool and r.status == "ok"
                  and r.diplotype not in ("-", "")), None)
        if a is not None:
            n_called = sum(1 for r in results if r.status == "ok" and r.diplotype not in ("-", ""))
            # G6PD safety net: a "no deficiency" authority call must not stand when an
            # independent caller (PyPGx/Stargazer/Aldy/PharmCAT) reports deficiency.
            if gene == "G6PD" and not re.match(r"G6PD (Deficient|Variable|not typed)|Indeterminate|not typed",
                                               a.phenotype or ""):
                dissent = [r.tool for r in results if r is not a and r.status == "ok"
                           and re.search(r"\bdeficient\b|decreased|\bvariable\b", r.phenotype or "", re.I)]
                if dissent:
                    # never leave the authority's negative text as the phenotype consumers read
                    return {"consensus_diplotype": "DISCORDANT",
                            "consensus_phenotype": (f"Indeterminate - discordant: {a.tool} found no G6PD "
                                                    f"deficiency allele but {', '.join(dissent)} report(s) "
                                                    "deficiency or variable activity; confirm by enzyme assay"),
                            "n_agree": 0, "n_called": n_called, "status": "discordant",
                            "reconciled": False, "reason": f"{a.tool} found no deficiency allele but "
                            f"{', '.join(dissent)} report(s) deficiency"}
            return {"consensus_diplotype": a.diplotype,
                    "consensus_phenotype": a.phenotype if a.phenotype not in ("-", "") else "-",
                    "n_agree": 1, "n_called": n_called, "status": "concordant",
                    "reconciled": False, "authority": auth_tool}

    pairs, phenos = [], []   # pairs: (canonical, raw_normalized)
    for r in results:
        if r.status == "ok" and r.diplotype not in ("-", ""):
            raw = _normalize_diplotype(r.diplotype)
            canon = reconcile.canonicalize_diplotype(gene, r.diplotype) if gene else raw
            pairs.append((canon, raw))
            if r.phenotype not in ("-", ""):
                phenos.append(r.phenotype)

    n_called = len(pairs)
    if n_called == 0:
        return {"consensus_diplotype": "-", "consensus_phenotype": "-",
                "n_agree": 0, "n_called": 0, "status": "no_data", "reconciled": False}

    canon_counts = Counter(c for c, _ in pairs).most_common()
    top_canon, n_agree = canon_counts[0]
    pheno = Counter(phenos).most_common(1)[0][0] if phenos else "-"

    # Did synonym collapse genuinely merge tools that the raw count split apart?
    raw_top_agree = Counter(raw for _, raw in pairs).most_common(1)[0][1]
    reconciled = bool(gene) and n_agree > raw_top_agree

    pheno_rescue = _phenotype_consensus(results)

    # Tie: a second diplotype called by as many tools — no safe winner on the
    # string. Clinically rescue if the phenotype-emitting callers nonetheless agree.
    if len(canon_counts) > 1 and canon_counts[1][1] == n_agree:
        if pheno_rescue:
            p_disp, p_n, d_disp = pheno_rescue
            return {"consensus_diplotype": d_disp, "consensus_phenotype": p_disp,
                    "n_agree": p_n, "n_called": n_called, "status": "concordant",
                    "reconciled": True, "phenotype_tier": True, "authority": "Phenotype"}
        return {"consensus_diplotype": "DISCORDANT", "consensus_phenotype": pheno,
                "n_agree": n_agree, "n_called": n_called, "status": "discordant",
                "reconciled": False}

    # Display the most common RAW spelling among tools whose canonical == winner.
    top_dip = Counter(raw for c, raw in pairs if c == top_canon).most_common(1)[0][0]
    status = "concordant" if n_agree == n_called else "majority"
    # A diplotype-string majority (not unanimous) is clinically concordant when
    # all phenotype-emitting callers agree on the phenotype.
    if status == "majority" and pheno_rescue:
        p_disp, p_n, _ = pheno_rescue
        return {"consensus_diplotype": top_dip, "consensus_phenotype": p_disp,
                "n_agree": p_n, "n_called": n_called, "status": "concordant",
                "reconciled": True, "phenotype_tier": True, "authority": "Phenotype"}
    return {"consensus_diplotype": top_dip, "consensus_phenotype": pheno,
            "n_agree": n_agree, "n_called": n_called, "status": status,
            "reconciled": reconciled}


def compute_verdict(results: list[CallerResult], no_call: bool,
                    gene: str | None = None) -> dict[str, Any]:
    """Gene-level verdict; NAT2's phenotype is overridden by the CPIC lookup."""
    v = _compute_verdict_core(results, no_call, gene)
    if gene == "NAT2" and v.get("status") in ("concordant", "majority"):
        v["consensus_phenotype"] = nat2_cpic_phenotype(v["consensus_diplotype"])
    return v


def pharmcat_cross_check(output_dir: str, gene: str, sample: str,
                         verdict: dict[str, Any]) -> dict[str, Any] | None:
    """Informational PharmCAT call for `gene`, read from the sample-level report.

    PharmCAT runs once per sample and reports ~20 CPIC genes; we surface its call
    as an independent CPIC-reference cross-check even for genes where it is NOT the
    verdict authority. This NEVER feeds the verdict — it only annotates agreement.

    Returns None when: PharmCAT is already this gene's authority (then it's the
    verdict, not a cross-check); the gene is CYP2D6 (PharmCAT does not resolve
    CYP2D6 structurally — Cyrius does); no report exists; or PharmCAT did not call
    the gene. `agrees` compares the CPIC *phenotype* (not the raw diplotype string),
    so tool nomenclature differences don't read as disagreement.
    """
    if AUTHORITATIVE.get(gene) == "PharmCAT" or gene == "CYP2D6":
        return None
    r = parse_pharmcat(output_dir, gene, sample)
    if r.status != "ok" or r.diplotype in ("-", ""):
        return None
    v_ph = _norm_phenotype(verdict.get("consensus_phenotype", ""))
    p_ph = _norm_phenotype(r.phenotype)
    return {
        "diplotype": r.diplotype,
        "phenotype": r.phenotype if r.phenotype not in ("-", "") else "-",
        "agrees": bool(v_ph) and v_ph == p_ph,
    }


def print_table(
    gene: str,
    sample: str,
    results: list[CallerResult],
    output_dir: str,
    region_depth: float | None = None,
    min_depth: float = DEFAULT_MIN_DEPTH,
    alt_block: str = "",
) -> None:
    W = 72
    SEP = "─" * W

    # Gate the call when (a) coverage is too low, or (b) alt_block is set — the
    # BAM is not ALT-aware and this is a paralog/SV gene whose call would be
    # unreliable (e.g. minimap2 -> false CYP2D6 *5 deletion).
    low_cov = region_depth is not None and region_depth < min_depth
    no_call = low_cov or bool(alt_block)
    verdict = compute_verdict(results, no_call, gene)
    if alt_block and verdict.get("status") == "no_call":
        verdict["reason"] = alt_block

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

    for r in results:
        print(f"{r.tool:<14}{r.diplotype:<18}{r.activity_score:<18}{r.phenotype}")

    print(SEP)

    st = verdict["status"]
    if st == "no_call":
        if alt_block:
            print(f"VERDICT: NO CALL — {alt_block}")
        else:
            print(f"VERDICT: NO CALL — insufficient coverage "
                  f"({region_depth:.1f}x < {min_depth:.0f}x minimum)")
        print("Tool calls above are suppressed and retained for audit only.")
    elif st == "discordant":
        print(f"VERDICT: DISCORDANT — tools disagree "
              f"({verdict['n_agree']}/{verdict['n_called']} on the top diplotype); no consensus")
    elif st == "no_data":
        print("Concordance: no calls available")
    else:  # concordant | majority
        print(f"Concordance: {verdict['n_agree']}/{verdict['n_called']} tools agree on "
              f"{verdict['consensus_diplotype']}")

    print("=" * W)
    print()

    # ── TSV ───────────────────────────────────────────────────────────────────
    tsv_path = os.path.join(output_dir, f"{gene}_{sample}_comparison.tsv")
    with open(tsv_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(
            ["Gene", "Sample", "Build", "Tool", "Diplotype",
             "ActivityScore", "Phenotype", "Status", "SVMode"]
        )
        sv_mode = _sv_note(gene)
        for r in results:
            if no_call:
                # Suppress every tool's diplotype so neither this TSV nor the
                # batch-summary headline can present an unreliable call.
                _reason = (alt_block if alt_block
                           else f"Insufficient coverage ({region_depth:.1f}x<{min_depth:.0f}x)")
                dip, score, pheno, status = ("NO_CALL", "-", _reason, "no_call")
            else:
                dip, score, pheno, status = (
                    r.diplotype, r.activity_score, r.phenotype, r.status)
            writer.writerow(
                [gene, sample, "GRCh38", r.tool, dip, score, pheno, status, sv_mode]
            )
    print(f"Full results saved to: {tsv_path}")

    # ── Detail JSON ───────────────────────────────────────────────────────────
    json_path = os.path.join(output_dir, f"{gene}_{sample}_detail.json")
    detail = {
        "gene":    gene,
        "sample":  sample,
        "build":   "GRCh38",
        "sv_mode": _sv_note(gene),
        "coverage": {
            "mean_depth": region_depth,
            "min_depth":  min_depth,
            "verdict":    "NO_CALL" if no_call else "OK",
        },
        # Gene-level verdict — the single authority. pgx-report.py and the batch
        # summary MUST read this instead of recomputing concordance.
        "verdict": verdict,
        # Per-tool dict keeps each caller's raw emitted values for audit even
        # when the gene is NO_CALL (the TSV/headline above are suppressed).
        "tools":   {r.tool: r.to_dict(sample, gene) for r in results},
    }
    # Informational PharmCAT cross-check (never feeds the verdict above).
    _xc = pharmcat_cross_check(output_dir, gene, sample, verdict)
    if _xc:
        detail["cross_check"] = {"PharmCAT": _xc}
    with open(json_path, "w") as fh:
        json.dump(detail, fh, indent=2)
    print(f"Detail JSON saved to:  {json_path}")


# ── Main ──────────────────────────────────────────────────────────────────────
def main() -> None:
    parser = argparse.ArgumentParser(
        description="Parse PGx caller outputs and produce comparison table + detail JSON."
    )
    parser.add_argument("--gene",       required=True, help="Gene name (e.g. CYP2D6)")
    parser.add_argument("--sample",     required=True, help="Sample name")
    parser.add_argument("--output-dir", required=True, help="Results directory")
    parser.add_argument("--region-depth", type=float, default=None,
                        help="Mean read depth over the gene region (drives the coverage gate)")
    parser.add_argument("--min-depth", type=float, default=DEFAULT_MIN_DEPTH,
                        help=f"Min mean depth to attempt a call (default {DEFAULT_MIN_DEPTH:g}x)")
    parser.add_argument("--failed-tools", default="",
                        help="Comma-separated tools that ran but failed (marked 'failed', not 'not_run')")
    parser.add_argument("--alt-aware", default="yes", choices=["yes", "no", "unknown"],
                        help="Whether the BAM was aligned ALT-aware. 'no'/'unknown' forces "
                             "NO_CALL for paralog/SV genes (CYP2D6, CYP2A6/2B6, GSTM1/T1, …).")
    parser.add_argument("--bam", default=None,
                        help="Input BAM/CRAM (lets VCF-Check genes verify per-site depth and infer sex)")
    parser.add_argument("--ref", default=None, help="Reference FASTA (needed to read depth from CRAM)")
    parser.add_argument("--sex", default=None,
                        help="Sample sex M/MALE or F/FEMALE (any case) overrides the chrX/chrY "
                             "inference for X-linked G6PD; any other value (unknown, NA) = infer")
    args = parser.parse_args()
    _RUN_CTX.update(bam=args.bam, ref=args.ref, sex={"M": "MALE", "MALE": "MALE", "F": "FEMALE",
                                                 "FEMALE": "FEMALE"}.get((args.sex or "").strip().upper()))

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
    if support.get("optitype"):
        results.append(parse_optitype(args.output_dir, gene, args.sample))
    if support.get("arcashla"):
        results.append(parse_arcashla(args.output_dir, gene, args.sample))
    if support.get("mutserve"):
        results.append(parse_mutserve(args.output_dir, gene, args.sample))
    if support.get("vcf_check"):
        _vcf_fn = _VCF_CHECK_PARSERS.get(gene)
        if _vcf_fn:
            results.append(_vcf_fn(args.output_dir, gene, args.sample))
    if support.get("cyrius"):
        results.append(parse_cyrius(args.output_dir, gene, args.sample))
    if support.get("pharmcat"):
        results.append(parse_pharmcat(args.output_dir, gene, args.sample))

    # Distinguish "tool ran and crashed" from "tool never ran". A crashed caller
    # leaves no output file, so it would otherwise look identical to not_run.
    failed = {t.strip() for t in args.failed_tools.split(",") if t.strip()}
    for r in results:
        if r.tool in failed and r.status == "not_run":
            r.status = "failed"
            r.diplotype = "tool failed (see log)"

    # ALT-awareness gate: paralog/SV genes are unreliable from a non-ALT-aware
    # alignment (reads mis-distributed across 2D6/2D7 + ALT contigs).
    alt_sensitive = (bool(support.get("cyrius"))
                     or gene in PYPGX_SV_GENES or gene in STARGAZER_SV_GENES)
    alt_block = ""
    if args.alt_aware in ("no", "unknown") and alt_sensitive:
        alt_block = ("aligner not ALT-aware — paralog/SV call unreliable; "
                     "re-align with bwa-mem + GRCh38 .alt")

    print_table(gene, args.sample, results, args.output_dir,
                region_depth=args.region_depth, min_depth=args.min_depth,
                alt_block=alt_block)


if __name__ == "__main__":
    main()
