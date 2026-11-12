#!/usr/bin/env python3
"""Variant-tier reconciliation (design A.1).

Star-allele callers disagree on *nomenclature* far more than on *biology*: PyPGx
calls DPYD `c.85T>C (*9A)`, Stargazer `*S10`, Aldy `*9A` — the same variant, three
names. The plain ensemble (pgx-compare.py) compares diplotype strings and marks
that "discordant". This module collapses known synonyms to a canonical token so
the verdict reflects genuine agreement.

Scoped deliberately: only well-established synonyms in `allele_synonyms.json`.
VKORC1 H-haplotypes and NUDT15 *5/*6 are NOT collapsed (they are not synonyms);
those go to array adjudication instead.

Backward-compatible: a gene with no synonym entry canonicalizes to the same
order-insensitive string the ensemble already used, so its verdict is unchanged.

Usage:
    reconcile.py --selftest
    reconcile.py --prove <SAMPLE> [--results results] [--gene DPYD]
"""
import argparse
import glob
import json
import os
import re
import sys

_HERE = os.path.dirname(os.path.realpath(__file__))
_SYN_PATH = os.path.join(_HERE, "allele_synonyms.json")

with open(_SYN_PATH, encoding="utf-8") as _fh:
    _RAW = json.load(_fh)

# Reference-allele synonyms (collapse to "*1") and per-gene token -> canonical.
_REFERENCE = {t.lower() for t in _RAW.get("_reference", [])}
SYNONYMS: dict[str, dict[str, str]] = {}
for _gene, _entry in _RAW.items():
    if _gene.startswith("_"):
        continue
    token_to_canon: dict[str, str] = {}
    for _canon, _aliases in _entry.items():
        token_to_canon[_canon.lower()] = _canon
        for _a in _aliases:
            token_to_canon[_a.lower()] = _canon
    SYNONYMS[_gene] = token_to_canon


def _strip_allele(tok: str) -> str:
    """Reduce one allele token to a comparable base form (lowercased)."""
    t = (tok or "").strip().lower()
    if not t:
        return ""
    t = re.sub(r"\s*\(.*?\)", "", t)   # drop "(*9a)" parentheticals
    t = re.sub(r"x\d+", "", t)          # copy number  *2xn / *106x2 -> *106
    t = t.split("+")[0]                  # *28+*80 -> *28 (base allele)
    t = t.split(";")[0]                  # pypgx "*17;*1;" -> *17
    return t.strip()


def canonicalize_allele(gene: str, tok: str) -> str:
    """Map one allele token to its canonical form for this gene.

    Copy-number / sub-allele stripping is applied ONLY for genes with a curated
    synonym map. For every other gene we keep the token as-is (lowercased) and
    apply only the universal Reference->*1 collapse — so CNV-meaningful calls
    like CYP2D6 `*1x2` or `*36+*10` are NOT flattened into `*1`/`*36`.
    """
    t = (tok or "").strip().lower()
    if not t or t in ("-", "."):
        return t
    if gene in SYNONYMS:
        base = _strip_allele(tok)
        if base in _REFERENCE:
            return "*1"
        return SYNONYMS[gene].get(base, base)
    return "*1" if t in _REFERENCE else t


# ── Notation normalisation (comparison only) ─────────────────────────────────
# The single place where caller diplotype strings are made comparable. Used by the
# verdict (consensus + n_agree) and by the PharmCAT outside-call builder. The
# per-caller table keeps the ORIGINAL strings; nothing here is written back.
#
# Genes whose trailing letter is a legacy sub-allele of the same star number
# (PharmVar CYP2D6 *4C, *6B, *10B ...). For other genes a letter can be a distinct
# allele (TPMT *3A != *3C, DPYD *2A), so letters are never stripped there.
LETTER_SUBALLELE_GENES = frozenset({"CYP2D6"})
UNRESOLVED = None   # norm_diplotype() result for a call that cannot be compared

_VCF_CHECK_RE = re.compile(r"^(\*[0-9A-Za-z]+)\((het|hom)\)(?:\[unphased\])?$", re.I)


def _clean_allele(gene, tok: str) -> str:
    """One allele token -> comparable spelling (case preserved).

    Rules: (1) Aldy's tandem marker is dropped (`*36X` -> `*36`; a copy number such as
    `*2x2` has digits and is kept); (2) a PharmVar sub-allele suffix is dropped
    (`*4.021` -> `*4`); (3) for genes in LETTER_SUBALLELE_GENES a trailing letter on the
    same star number is dropped (`*6B` -> `*6`, `*4C` -> `*4`); (4) `rsNNN` novel-variant
    tags are NOT stripped: `*2+rs61737946` stays different from `*2`.
    """
    t = re.sub(r"\s+", "", tok or "")
    if re.match(r"^\*[0-9]+x$", t, re.I):          # (1) Aldy tandem marker
        t = t[:-1]
    t = re.sub(r"^(\*[0-9]+)\.[0-9]+", r"\1", t)    # (2) .021 sub-allele
    t = re.sub(r"X(?=[0-9]+$)", "x", t)             # copy number is written x2 (CPIC, PharmVar), never X2
    if gene in LETTER_SUBALLELE_GENES:              # (3) legacy letter sub-allele
        t = re.sub(r"^(\*[0-9]+)[A-Za-z](?=$|x[0-9]+$)", r"\1", t)
    if t.lower() in _REFERENCE:
        t = "*1"
    return t


def clean_diplotype(gene, dip: str) -> str:
    """Display-safe cleaned spelling (order and case kept): used for the PharmCAT call."""
    halves = [("+".join(_clean_allele(gene, a) for a in h.split("+")))
              for h in re.sub(r"\s+", "", dip or "").split("/")]
    return "/".join(halves)


def norm_diplotype(gene, dip: str):
    """Comparison key for a caller's diplotype string, or UNRESOLVED (None).

    Rules, in order:
      1. whitespace is removed and case is ignored;
      2. VCF-Check's `*N(het) [unphased]` is `*1/*N` and `*N(hom)` is `*N/*N`, ONLY for a
         single listed allele; two or more listed alleles (`*6(het)+*60(het) [unphased]`)
         cannot be placed on haplotypes, so the call is UNRESOLVED and the caller must
         be left out of the diplotype comparison (neither agreeing nor disagreeing);
      3. per-allele cleaning as in _clean_allele (Aldy `X` marker, `.NNN`, CYP2D6 letter
         sub-allele; `rsNNN` tags kept as real differences);
      4. alleles inside a tandem (`*36+*10` == `*10+*36`) and the two haplotypes
         (`*2/*4` == `*4/*2`) are sorted.
    `-`/empty returns "-".
    """
    d = re.sub(r"\s+", "", dip or "")
    if not d or d == "-":
        return "-"
    if re.search(r"\((?:het|hom)\)", d, re.I):
        m = _VCF_CHECK_RE.match(d)
        if not m:
            return UNRESOLVED
        allele = _clean_allele(gene, m.group(1)).lower()
        return "/".join(sorted(["*1", allele] if m.group(2).lower() == "het" else [allele, allele]))
    if "/" not in d:
        return "+".join(sorted(_clean_allele(gene, a).lower() for a in d.split("+")))
    halves = ["+".join(sorted(_clean_allele(gene, a).lower() for a in h.split("+")))
              for h in d.split("/")]
    return "/".join(sorted(halves))


def canonicalize_diplotype(gene: str, dip: str):
    """Order-insensitive, notation-normalised, synonym-collapsed diplotype key.

    Genes with a synonym entry (DPYD, CYP2C19, TPMT) additionally collapse synonyms.
    Returns UNRESOLVED (None) for a call that cannot be compared (see norm_diplotype).
    """
    if not dip or dip == "-":
        return "-"
    key = norm_diplotype(gene, dip)
    if key is UNRESOLVED or gene not in SYNONYMS:
        return key
    d = re.sub(r"\s+", "", dip)
    m = _VCF_CHECK_RE.match(d)
    if m:
        d = ("*1/" if m.group(2).lower() == "het" else m.group(1) + "/") + m.group(1)
    if "/" not in d:
        return canonicalize_allele(gene, d)
    halves = [canonicalize_allele(gene, p) for p in d.split("/")]
    return "/".join(sorted(halves))


# ── CYP2D6 activity score (CPIC data, no value in code) ─────────────────────────
_CPIC_PATH = os.path.join(_HERE, "cpic_cyp2d6.json")
PHENO_CLASSES = ("Poor", "Intermediate", "Normal", "Ultrarapid")
_cpic_cache = None


def _cpic():
    """docker/cpic_cyp2d6.json (built by cpic_cyp2d6_build.py from the CPIC API)."""
    global _cpic_cache
    if _cpic_cache is None:
        with open(_CPIC_PATH, encoding="utf-8") as fh:
            _cpic_cache = json.load(fh)
    return _cpic_cache


def _haplotype_activity(hap: str):
    """CPIC activity value of one cleaned haplotype (`*10`, `*2x2`, `*36+*10`), or None.

    A tandem is looked up under every order of its parts (CPIC writes `*36+*10`,
    `*68+*4`). None for: an `rsNNN` novel tag, an unknown copy number (`*2x>=3`), an allele
    absent from the CPIC table, or one whose CPIC activity value is not numeric.
    """
    import itertools
    if re.search(r"rs\d", hap, re.I) or not re.fullmatch(r"\*[0-9]+(x[0-9]+)?(\+\*[0-9]+(x[0-9]+)?)*", hap):
        return None
    table = _cpic()["alleles"]
    for perm in itertools.permutations(hap.split("+")):
        entry = table.get("+".join(perm))
        if entry is not None:
            return entry["activity"]
    return None


def activity_score(dip: str, gene: str = "CYP2D6"):
    """Total CYP2D6 activity score (sum over the two haplotypes) or None when unknown."""
    d = clean_diplotype(gene, dip)
    if "(" in d or d.count("/") != 1:
        return None
    vals = [_haplotype_activity(h) for h in d.split("/")]
    return None if None in vals else round(sum(vals), 2)


def score_phenotype(score):
    """CPIC phenotype for a total activity score (gene_result table), or None."""
    if score is None:
        return None
    sp = _cpic()["score_phenotype"]
    hit = sp["exact"].get(f"{score:.2f}")
    if hit:
        return hit
    ge = [float(k) for k in sp["at_least"] if float(k) <= score]
    return sp["at_least"][f"{max(ge):.2f}"] if ge else None


def phenotype_class(text):
    """Poor/Intermediate/Normal/Ultrarapid for a phenotype string; "Indeterminate" for an
    explicit indeterminate; None when there is no usable phenotype."""
    t = (text or "").strip().lower().replace("_", " ")
    if t in ("", "-", "na", "n/a", "none"):
        return None
    if "indeterminate" in t or "unknown" in t or "no call" in t:
        return "Indeterminate"
    # the phenotype word must lead the string: "Not Normal" or "Normal to Ultrarapid" is not a class
    m = re.match(r"(poor|intermediate|normal|ultrarapid)(\s+metaboli[sz]er\b.*)?$", t)
    return m.group(1).capitalize() if m else None


# ── proof / self-check ────────────────────────────────────────────────────────
def _prove(sample: str, results_dir: str, only_gene: str | None):
    """Show that callers reconcile to the same canonical diplotype on real data."""
    genes = [only_gene] if only_gene else ["DPYD"]
    any_seen = False
    for gene in genes:
        det = glob.glob(os.path.join(results_dir, sample, "genes", gene, "*detail.json"))
        if not det:
            continue
        any_seen = True
        d = json.load(open(det[0]))
        raw, canon = {}, set()
        for t, info in d.get("tools", {}).items():
            dip = info.get("diplotype")
            if dip in (None, "-", ""):
                continue
            c = canonicalize_diplotype(gene, dip)
            raw[t] = (dip, c)
            canon.add(c)
        print(f"[{sample}] {gene}: pipeline verdict = "
              f"{d.get('verdict', {}).get('status', '?')} "
              f"({d.get('verdict', {}).get('consensus_diplotype', '?')})")
        for t, (dip, c) in raw.items():
            print(f"    {t:11} {dip:32} -> {c}")
        verdict = "CONCORDANT after reconciliation" if len(canon) == 1 else \
                  f"still differs: {sorted(canon)}"
        print(f"    => {verdict}\n")
    if not any_seen:
        print(f"(no detail.json found for {sample})")


def _selftest():
    # DPYD: three caller spellings of the same variant collapse to one canonical
    assert canonicalize_diplotype("DPYD", "Reference/c.85T>C (*9A)") == \
           canonicalize_diplotype("DPYD", "*1/*S10") == \
           canonicalize_diplotype("DPYD", "*1/*9A"), "DPYD *9A synonyms must reconcile"
    assert canonicalize_diplotype("DPYD", "*1/*5") == \
           canonicalize_diplotype("DPYD", "*1/*S6") == \
           canonicalize_diplotype("DPYD", "*1/c.1627A>G"), "DPYD *5 synonyms must reconcile"
    # CYP2C19 sub-allele collapse
    assert canonicalize_diplotype("CYP2C19", "*2/*38") == \
           canonicalize_diplotype("CYP2C19", "*2/*2"), "CYP2C19 *38 is a *2 sub-allele"
    # TPMT *1S inflation
    assert canonicalize_diplotype("TPMT", "*1/*1S") == "*1/*1"
    # order-insensitive
    assert canonicalize_diplotype("DPYD", "*9A/*1") == canonicalize_diplotype("DPYD", "*1/*9A")
    # backward-compat: a gene with no entry behaves like plain normalize
    assert canonicalize_diplotype("CYP2D6", "*4/*2") == "*2/*4"
    # we DO NOT over-collapse: VKORC1/NUDT15 are not in the map, stay distinct
    assert canonicalize_diplotype("NUDT15", "*1/*5") != canonicalize_diplotype("NUDT15", "*1/*6")
    assert "VKORC1" not in SYNONYMS, "VKORC1 must be left to adjudication, not collapsed"
    print("reconcile self-test passed")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--prove", metavar="SAMPLE")
    ap.add_argument("--results", default="results")
    ap.add_argument("--gene")
    args = ap.parse_args()
    if args.selftest:
        _selftest()
    elif args.prove:
        _prove(args.prove, args.results, args.gene)
    else:
        ap.error("pass --selftest or --prove <SAMPLE>")


if __name__ == "__main__":
    main()
