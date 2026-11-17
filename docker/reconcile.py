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
import collections
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


# A haplotype is cut into components only at a `+` that precedes `rsNNN`, `*` or an HGVS `c.` name, so an HGVS name keeps its
# own `+` (`c.1905+1G>A`).
_COMPONENT_SPLIT = re.compile(r"\+(?=rs\d|\*|c\.)", re.I)


def _split_tags(tok: str):
    """(comparable base form, dropped tag components) of one allele token. The base is lowercased; the tags are the
    components after the first (`*28+*80` -> `*28` and [`*80`], `*9a+rs1+rs2` -> `*9a` and [`rs1`, `rs2`],
    `*1+c.1905+1g>a` -> `*1` and [`c.1905+1g>a`]); an HGVS name alone keeps its `+` (`c.1905+1g>a`, no tags)."""
    t = (tok or "").strip().lower()
    if not t:
        return "", []
    t = re.sub(r"\s*\(.*?\)", "", t)   # drop "(*9a)" parentheticals
    t = re.sub(r"x\d+", "", t)          # copy number  *2xn / *106x2 -> *106
    tags = []
    if t.startswith("*") or re.match(r"rs\d+\+", t):   # *28+*80 -> *28, rs1+rs2 -> rs1
        parts = _COMPONENT_SPLIT.split(t)
        t, tags = parts[0], [x.split(";")[0].strip() for x in parts[1:]]
    t = t.split(";")[0]                  # pypgx "*17;*1;" -> *17
    return t.strip(), [x for x in tags if x]


def _strip_allele(tok: str) -> str:
    """Reduce one allele token to a comparable base form (lowercased)."""
    return _split_tags(tok)[0]


def _canonical_allele(gene: str, tok: str):
    """(canonical key, blocking tags) of one allele token (see canonicalize_allele)."""
    t = (tok or "").strip().lower()
    if not t or t in ("-", "."):
        return t, []
    if gene in SYNONYMS:
        base, tags = _split_tags(tok)
        canon = "*1" if base in _REFERENCE else SYNONYMS[gene].get(base, base)
        blocked = blocking_tags(gene, tags)
        if blocked:
            # a tag that is not Normal function (or whose function CPIC's table does not give) is evidence the
            # stripped diplotype lacks: the key keeps those tags (not the Normal ones), spelled by the gene's
            # notation / synonym map so one variant is one key, and the call is its own reading (owner 2026-11-16)
            return canon + "+" + "+".join(sorted(k for k, _, _ in blocked)), blocked
        return canon, []
    return ("*1" if t in _REFERENCE else t), []


def canonicalize_allele(gene: str, tok: str) -> str:
    """Map one allele token to its canonical form for this gene.

    Copy-number / sub-allele stripping is applied ONLY for genes with a curated
    synonym map. For every other gene we keep the token as-is (lowercased) and
    apply only the universal Reference->*1 collapse — so CNV-meaningful calls
    like CYP2D6 `*1x2` or `*36+*10` are NOT flattened into `*1`/`*36`.
    A tag the stripping would drop stays in the key when it is not Normal function in CPIC's allele table.
    """
    return _canonical_allele(gene, tok)[0]


# ── Notation normalisation (comparison only) ─────────────────────────────────
# The single place where caller diplotype strings are made comparable. Used by the
# verdict (consensus + n_agree) and by the PharmCAT outside-call builder. The
# per-caller table keeps the ORIGINAL strings; nothing here is written back.
#
# A legacy trailing letter (CYP2D6 `*14A`, `*14B`, `*13C`, `*4N`) is never stripped or mapped.
# PharmVar re-designated alleles that had been defined from partial sequence (Nofziger et al.
# 2020, PharmVar GeneFocus: CYP2D6, PMID 31544239), so a letter-free star number can now be a
# haplotype with another function (CPIC table in cpic_cyp2d6.json: `*14` decreased 0.5,
# `*114` no function). The repo carries no cited legacy->current table, so a lettered name is
# kept verbatim: it compares only equal to itself and has no activity score.
# A legacy lettered allele name (`*13C`, `*4N`, `*4NX`, `*14A`, any case). `x<digit>` is a copy number and a
# bare trailing `x`/`X` (Aldy's tandem marker, `*36X`) is notation that _clean_allele normalises: neither matches.
LEGACY_LETTERED = re.compile(r"\*\d+(?:(?![xX](?:\d|$|[+/]))[A-Za-z])")
UNRESOLVED = None   # norm_diplotype() result for a call that cannot be compared
UNTYPED = {"", "-", "NA", "N/A", "NONE", "UNKNOWN"}


def is_typed_call(status, dip) -> bool:
    """A caller row that typed the gene: status ok and a phased `a/b` diplotype (no no-call marker)."""
    d = (dip or "").strip()
    return (status or "").strip().lower() == "ok" and not (d.upper() in UNTYPED or "NO_CALL" in d.upper()
                                                            or "/" not in d)

_VCF_CHECK_RE = re.compile(r"^(\*[0-9A-Za-z]+)\((het|hom)\)(?:\[unphased\])?$", re.I)


def _clean_allele(gene, tok: str) -> str:
    """One allele token -> comparable spelling (case preserved).

    Rules: (1) Aldy's tandem marker is dropped (`*36X` -> `*36`; a copy number such as
    `*2x2` has digits and is kept); (2) a PharmVar sub-allele suffix is dropped
    (`*4.021` -> `*4`); (3) a legacy trailing letter is NOT stripped (`*14A` stays `*14A`, never `*14`;
    see the note above); (4) `rsNNN` novel-variant
    tags are NOT stripped: `*2+rs61737946` stays different from `*2`.
    """
    t = re.sub(r"\s+", "", tok or "")
    if re.match(r"^\*[0-9]+x$", t, re.I):          # (1) Aldy tandem marker
        t = t[:-1]
    t = re.sub(r"^(\*[0-9]+)\.[0-9]+", r"\1", t)    # (2) .021 sub-allele
    t = re.sub(r"X(?=[0-9]+$)", "x", t)             # copy number is written x2 (CPIC, PharmVar), never X2
    if t.lower() in _REFERENCE:
        t = "*1"
    return t


def _fold_tandem(hap: str) -> str:
    """Fold repeated parts of one CYP2D6 haplotype into CPIC copy-number spelling: `*36+*36+*10` ->
    `*36x2+*10`, `*2+*2` -> `*2x2`, `*36x2+*36` -> `*36x3`. Parts keep first-appearance order. A haplotype
    with any part that is not `*N` or `*NxM` (an `rsNNN` tag) is returned unchanged."""
    counts: dict[str, int] = {}
    for part in hap.split("+"):
        m = re.fullmatch(r"(\*[0-9]+)(?:x([0-9]+))?", part)
        if not m:
            return hap
        counts[m.group(1)] = counts.get(m.group(1), 0) + int(m.group(2) or 1)
    return "+".join(a if n == 1 else f"{a}x{n}" for a, n in counts.items())


def clean_diplotype(gene, dip: str) -> str:
    """Display-safe cleaned spelling (order and case kept): used for the PharmCAT call. For CYP2D6 repeated
    tandem parts are folded into CPIC copy-number spelling (_fold_tandem); comparison keys (norm_diplotype)
    are not folded."""
    halves = [("+".join(_clean_allele(gene, a) for a in h.split("+")))
              for h in re.sub(r"\s+", "", dip or "").split("/")]
    if gene == "CYP2D6":
        halves = [_fold_tandem(h) for h in halves]
    return "/".join(halves)


def norm_diplotype(gene, dip: str):
    """Comparison key for a caller's diplotype string, or UNRESOLVED (None).

    Rules, in order:
      1. whitespace is removed and case is ignored;
      2. VCF-Check's `*N(het) [unphased]` is `*1/*N` and `*N(hom)` is `*N/*N`, ONLY for a
         single listed allele; two or more listed alleles (`*6(het)+*60(het) [unphased]`)
         cannot be placed on haplotypes, so the call is UNRESOLVED and the caller must
         be left out of the diplotype comparison (neither agreeing nor disagreeing);
      3. per-allele cleaning as in _clean_allele (Aldy `X` marker, `.NNN`; a legacy
         letter is kept as written; `rsNNN` tags kept as real differences);
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


def _canonical_diplotype(gene: str, dip: str):
    """(key, blocking tags) behind canonicalize_diplotype; the tags are [(key, tag as written, function)] of the
    components kept in the key because they are not Normal function (see blocking_tags)."""
    if not dip or dip == "-":
        return "-", []
    key = norm_diplotype(gene, dip)
    if key is UNRESOLVED or gene not in SYNONYMS:
        return key, []
    d = re.sub(r"\s+", "", dip)
    m = _VCF_CHECK_RE.match(d)
    if m:
        d = ("*1/" if m.group(2).lower() == "het" else m.group(1) + "/") + m.group(1)
    if "/" not in d:
        return _canonical_allele(gene, d)
    parts = [_canonical_allele(gene, p) for p in d.split("/")]
    return "/".join(sorted(k for k, _ in parts)), [t for _, ts in parts for t in ts]


def canonicalize_diplotype(gene: str, dip: str):
    """Order-insensitive, notation-normalised, synonym-collapsed diplotype key.

    Genes with a synonym entry (DPYD, CYP2C19, TPMT) additionally collapse synonyms.
    Returns UNRESOLVED (None) for a call that cannot be compared (see norm_diplotype).
    """
    return _canonical_diplotype(gene, dip)[0]


def canonical_tags(gene: str, dip: str):
    """[(key, tag, function)] of the tags kept in canonicalize_diplotype's key for `dip` (empty when none)."""
    return _canonical_diplotype(gene, dip)[1]


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


def _haplotype_value(hap: str):
    """CPIC activity of one cleaned haplotype (`*10`, `*2x2`, `*36x2+*10`) as (value, is_lower_bound), or None.

    A tandem is looked up under every order of its parts (CPIC writes `*36+*10`, `*68+*4`). A copy number N>=3
    with no CPIC row of its own (`*1x3`) uses CPIC's `x>=3` row (`*1x≥3`, activity `≥3.0`: a lower bound, or
    `*4x≥3`, activity 0.0: exact). None for: an `rsNNN` novel tag, a copy number CPIC does not list, an allele
    absent from the CPIC table, or one whose CPIC activity value is neither a number nor `≥number`.
    """
    import itertools
    if re.search(r"rs\d", hap, re.I) or not re.fullmatch(r"\*[0-9]+(x[0-9]+)?(\+\*[0-9]+(x[0-9]+)?)*", hap):
        return None
    table = _cpic()["alleles"]

    def lookup(parts):
        for perm in itertools.permutations(parts):
            entry = table.get("+".join(perm))
            if entry is not None:
                return entry
        return None
    parts = hap.split("+")
    entry = lookup(parts)
    if entry is None:
        ge3 = [re.sub(r"x([0-9]+)$", lambda m: "x≥3" if int(m.group(1)) >= 3 else m.group(0), p) for p in parts]
        entry = lookup(ge3) if ge3 != parts else None
    if entry is None:
        return None
    if entry["activity"] is not None:
        return entry["activity"], False
    if entry.get("activity_min") is not None:
        return entry["activity_min"], True
    return None


def _haplotype_activity(hap: str):
    """Exact CPIC activity value of one haplotype, or None (a lower-bound row is not an exact value)."""
    v = _haplotype_value(hap)
    return v[0] if v and not v[1] else None


def _diplotype_values(dip: str, gene: str = "CYP2D6"):
    d = clean_diplotype(gene, dip)
    if "(" in d or d.count("/") != 1:
        return None
    vals = [_haplotype_value(h) for h in d.split("/")]
    return None if None in vals else vals


def activity_score(dip: str, gene: str = "CYP2D6"):
    """Total CYP2D6 activity score (sum over the two haplotypes) or None when unknown or only a lower bound."""
    vals = _diplotype_values(dip, gene)
    if vals is None or any(b for _, b in vals):
        return None
    return round(sum(v for v, _ in vals), 2)


def activity_lower_bound(dip: str, gene: str = "CYP2D6"):
    """Lower bound of the total activity score when a haplotype is a CPIC `x>=3` row with a `≥` activity
    (`*1x3/*10` -> 3.25), else None. Never an exact score."""
    vals = _diplotype_values(dip, gene)
    if vals is None or not any(b for _, b in vals):
        return None
    return round(sum(v for v, _ in vals), 2)


def bound_phenotype(lower_bound):
    """CPIC phenotype for 'activity score >= lower_bound' (the `≥` rows of gene_result), or None."""
    if lower_bound is None:
        return None
    at_least = _cpic()["score_phenotype"]["at_least"]
    ge = [float(k) for k in at_least if float(k) <= lower_bound]
    return at_least[f"{max(ge):.2f}"] if ge else None


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


# ── TPMT / NUDT15: CPIC diplotype -> phenotype (CPIC data, no value in code) ────
_THIO_PATH = os.path.join(_HERE, "cpic_thiopurine.json")
_thio_cache = None


def thiopurine_phenotype(gene: str, dip: str):
    """CPIC phenotype of a TPMT or NUDT15 diplotype from CPIC's own diplotype table
    (docker/cpic_thiopurine.json, built by cpic_thiopurine_build.py), or None when the table has no row.
    The lookup is order-insensitive and TPMT `*1S` is read as `*1` (allele_synonyms.json)."""
    global _thio_cache
    if _thio_cache is None:
        raw = _load_lookup(_THIO_PATH, "CPIC thiopurine table")
        _thio_cache = {g: {canonicalize_diplotype(g, k): v for k, v in raw[g].items()} for g in ("TPMT", "NUDT15")}
    key = canonicalize_diplotype(gene, dip)
    if gene not in _thio_cache or key in (None, "-"):
        return None
    return _thio_cache[gene].get(key)


# ── CPIC diplotype -> phenotype for the other star-allele genes (CPIC data, no value in code) ──
_DIP_PATH = os.path.join(_HERE, "cpic_diplotype.json")
_dip_cache = None
# Legacy lettered names -> the numbered core allele, from an explicit PharmVar-verified table only
# (pharmvar_legacy_suballeles.json, derived by pharmvar_legacy_build.py from the rendered PharmVar gene pages in
# pharmvar_legacy_snapshot.json, PharmVar 6.2.29). A lettered allele that is not listed under `fold` for its gene is
# never folded (it stays unscorable and counts toward the split rule). `redesignated` lists letters whose number
# differs from the core allele they now sit under (CYP2C19 `*1A` is now `*38.001`, PMID 32602114: the old meaning
# was the reference allele); they are never folded, since the caller may mean either (owner decision 2026-11-16:
# CYP2C19 `*1A` stays unscorable rather than being read as `*38`). Such a reading has no CPIC allele function either
# (reading_function_profile reports `*1A` unknown), so it is never set aside by the unscorable-reading rule: it counts
# toward the split.
_LEGACY_PATH = os.path.join(_HERE, "pharmvar_legacy_suballeles.json")
_legacy_fold = None


def _load_lookup(path: str, what: str):
    """Read one lookup JSON on first use; a missing or unreadable file fails loud, naming the file (never an empty table)."""
    try:
        with open(path, encoding="utf-8") as fh:
            return json.load(fh)
    except (OSError, ValueError) as e:
        raise RuntimeError(f"{what}: required lookup file {os.path.basename(path)} ({path}) cannot be read: {e}") from e


def _legacy_table():
    global _legacy_fold
    if _legacy_fold is None:
        raw = _load_lookup(_LEGACY_PATH, "PharmVar legacy sub-allele table")
        _legacy_fold = {g: {k.casefold(): v["core_allele"] for k, v in d.items()} for g, d in raw["fold"].items()}
    return _legacy_fold


def fold_legacy_letter(gene: str, allele: str) -> str:
    """`*2B` -> `*2` only when PharmVar lists `*2B` as a legacy label of a sub-allele of CYP2C19*2; else unchanged.
    Case-insensitive, like every other allele comparison."""
    a = allele.strip()
    return _legacy_table().get(gene, {}).get(a.casefold(), a)


# Caller-specific allele names -> CPIC allele names (docker/dpyd_notation_map.json, built by dpyd_notation_map_build.py
# from Stargazer's shipped allele table and the CPIC allele definitions). Used only to look up a CPIC phenotype
# (cpic_diplotype_phenotype); the verdict's comparison keys (norm_diplotype) keep the caller's own spelling.
_NOTATION_PATH = os.path.join(_HERE, "dpyd_notation_map.json")
_notation_map = None


def _notation():
    global _notation_map
    if _notation_map is None:
        raw = _load_lookup(_NOTATION_PATH, "caller notation map")
        _notation_map = {
            g: {k: v["cpic"] for sec in ("stargazer", "aldy", "cpic_star_names") for k, v in e[sec].items()}
            for g, e in raw.items() if g not in ("retrieved", "source", "cpic_alleles")}
    return _notation_map


def _cpic_allele_key(gene: str, tok: str) -> str:
    """One CPIC-table or caller allele name -> comparable key for the diplotype table lookup. Unlike
    canonicalize_allele it keeps a `+` inside a name (`c.1905+1G>A`), so no two CPIC alleles merge. The synonym
    map of the gene (allele_synonyms.json) is tried on the whole name, then on its parenthetical, then on the
    name without it (`c.85T>C (*9A)`, `*9A`, `c.85T>C` all give `rs1801265`)."""
    t = re.sub(r"\s+", " ", (tok or "")).strip().lower()
    nm = _notation().get(gene, {})
    if t in nm:                                     # Stargazer *S10 / Aldy *rs2297595 -> the CPIC allele name
        t = re.sub(r"\s+", "", nm[t]).lower()   # _clean_allele has already dropped spaces from table keys
    syn = SYNONYMS.get(gene, {})
    m = re.match(r"^(.*?)\s*\((.*)\)$", t)
    for cand in ((t,) + ((m.group(2).strip(), m.group(1).strip()) if m else ())):
        if cand in _REFERENCE:
            return "*1"
        if cand in syn:
            return syn[cand]
    return t


def _dip_key(gene: str, dip: str, collapse: bool = True):
    """Order-insensitive key of an `a/b` diplotype, or None when it is not a two-haplotype call. With
    collapse=False the alleles are only cleaned and lowercased (CPIC's own spelling)."""
    d = re.sub(r"\s+", " ", (dip or "")).strip()
    if d.count("/") != 1 or re.search(r"\((?:het|hom)\)", d, re.I):
        return None
    f = _cpic_allele_key if collapse else (lambda g, t: t.lower())
    return "/".join(sorted(f(gene, _clean_allele(gene, h.strip())) for h in d.split("/")))


def _dip_tables():
    """{gene: (exact, collapsed)}: exact keeps CPIC's own allele spelling, collapsed folds synonyms (a key shared
    by two CPIC rows of different phenotype is None: never scored)."""
    global _dip_cache
    if _dip_cache is None:
        raw = _load_lookup(_DIP_PATH, "CPIC diplotype table")
        _dip_cache = {}
        for g, rows in raw.items():
            if g in ("retrieved", "source", "allele_function"):
                continue
            tabs = ({}, {})
            for k, v in rows.items():
                for tab, collapse in zip(tabs, (False, True)):
                    dk = _dip_key(g, k, collapse)
                    tab[dk] = v if tab.get(dk, v) == v else None
            _dip_cache[g] = tabs
    return _dip_cache


def has_cpic_diplotype_table(gene: str) -> bool:
    return gene in _dip_tables()


def cpic_diplotype_lookup(gene: str, dip: str):
    """(state, phenotype) of one diplotype reading in CPIC's own table. state is "scored" (phenotype set),
    "ambiguous" (the reading's collapsed key is shared by CPIC rows of different phenotype, so none is stored) or
    "no_row" (CPIC has no row). Lookup order and folding are those of cpic_diplotype_phenotype."""
    tabs = _dip_tables()
    if gene not in tabs or gene == "CYP2D6":
        return "no_row", None
    exact, coll = tabs[gene]
    folded = "/".join(fold_legacy_letter(gene, h) for h in (dip or "").split("/"))
    ambiguous = False
    for d, tab, collapse in ((dip, exact, False), (dip, coll, True), (folded, exact, False), (folded, coll, True)):
        k = _dip_key(gene, d, collapse)
        if k is None or k not in tab:
            continue
        if tab[k] is not None:
            return "scored", tab[k]
        ambiguous = True            # an ambiguous (None) row falls through to the folded lookups
    return ("ambiguous" if ambiguous else "no_row"), None


def cpic_diplotype_phenotype(gene: str, dip: str):
    """CPIC phenotype of one diplotype reading from CPIC's own diplotype table (docker/cpic_diplotype.json, built
    by cpic_diplotype_build.py), or None when CPIC has no row for it or its row is ambiguous (see
    cpic_diplotype_lookup, which tells the two apart).
    Lookup order: the reading as spelled (synonym-collapsed, order-insensitive; DPYD `*9A` == `c.85T>C (*9A)`),
    then, when that finds no row, with a legacy lettered name folded to its numbered core allele (`*3B` -> `*3`,
    `*2B` -> `*2`) but only where PharmVar's Legacy Label column lists it (LEGACY_FOLD); any other lettered allele
    is not folded and gets no row.
    Not used for CYP2D6, which has its own activity-score rule."""
    return cpic_diplotype_lookup(gene, dip)[1]


# ── allele function of the components of a reading (for the unscorable-reading rule) ──────────────────
# CPIC's own allele-function table (cpic_diplotype.json `allele_function`, built by cpic_diplotype_build.py from CPIC
# /allele data: CYP2C19, DPYD, TPMT). A gene without an entry has no known component function, so none of its
# readings can be set aside. A component is split from a reading at `+` only before `rsNNN` or `*` (an HGVS name such
# as `c.1905+1G>A` keeps its `+`); `*rsNNN` / `rsNNN` / `*S` codes / HGVS names go through the same notation and
# synonym maps as the diplotype lookup.
_function_cache = None


def _allele_functions(gene: str):
    """{comparable allele key: CPIC function}; a key that two CPIC alleles give different functions is None."""
    global _function_cache
    if _function_cache is None:
        raw = _load_lookup(_DIP_PATH, "CPIC diplotype table").get("allele_function", {})
        _function_cache = {}
        for g, alleles in raw.items():
            tab = {}
            for name, fn in alleles.items():
                k = _cpic_allele_key(g, _clean_allele(g, name))
                tab[k] = fn if tab.get(k, fn) == fn else None
            _function_cache[g] = tab
    return _function_cache.get(gene)


_KNOWN_FUNCTIONS = ("Normal function", "Decreased function", "No function", "Increased function")


def _component_function(gene: str, comp: str):
    """(cleaned component, comparable key, CPIC function) of one allele component (named allele, `rsNNN` / `*rsNNN`
    tag, `*S` code, HGVS name). The function is "Normal function", ... from CPIC's allele-function table, or None
    when the gene has no table or the table does not give it (not listed, ambiguous, Unknown/Uncertain function).
    The cleaned component is the spelling that was looked up (`rsNNN` becomes `*rsNNN`)."""
    comp = _clean_allele(gene, comp.strip())
    if re.fullmatch(r"rs\d+", comp, re.I):         # Aldy's `+rsNNN` tag is its `*rsNNN` name
        comp = "*" + comp
    key = _cpic_allele_key(gene, comp)
    fns = _allele_functions(gene)
    if fns is None:
        return comp, key, None
    fn = fns.get(key)
    if fn is None:
        fn = fns.get(_cpic_allele_key(gene, fold_legacy_letter(gene, comp)))
    return comp, key, (fn if fn in _KNOWN_FUNCTIONS else None)


def blocking_tags(gene: str, tags):
    """[(key, tag as written, function or "unknown")] of the dropped tag components that are not Normal function. A
    tag blocks the stripping when its CPIC function is decreased / no / increased, or unknown (no table for the
    gene, or the table does not give it)."""
    out = []
    for t in tags:
        _, key, fn = _component_function(gene, t)
        if fn != "Normal function":
            out.append((key, t, fn or "unknown"))
    return out


# A caller's own named allele that CPIC's table does not give (a Stargazer `*S` code not in dpyd_notation_map.json, an
# Aldy `*rsNNN` not mapped) is a variant of unknown function, like the same variant written as a `+rsNNN` tag, so it
# also blocks a Normal DPYD verdict (owner 2026-11-16). On the 66 stored DPYD calls (198 callers' diplotypes) none has one.
UNKNOWN_NAMED_ALLELE_BLOCKS_NORMAL = True


def _walk_components(gene: str, dip: str, allow_vcf_form: bool = False):
    """[(component as written, cleaned component, key, CPIC function or None, is_tag)] of one caller diplotype, or None
    when it cannot be read as two haplotypes (or, with allow_vcf_form, as VCF-Check's `*A(het)+*B(het) [unphased]`
    single list of alleles). The diplotype is cut at `/` first; a haplotype written with VCF-Check's `(het)` / `(hom)`
    markers has them removed and all of its components are main alleles. The one place a diplotype is cut into
    components."""
    haps = [h.strip() for h in (dip or "").strip().split("/")]
    if not all(haps) or len(haps) > 2:
        return None
    marker = re.compile(r"\((?:het|hom)\)", re.I)
    if len(haps) == 1 and not (allow_vcf_form and marker.search(haps[0])):
        return None
    out = []
    for hap in haps:
        main_only = bool(marker.search(hap))
        if main_only:
            hap = re.sub(r"\s*\[unphased\]", "", marker.sub("", hap)).strip()
        for i, comp in enumerate(_COMPONENT_SPLIT.split(hap)):
            c, key, fn = _component_function(gene, comp)
            out.append((comp.strip(), c, key, fn, i > 0 and not main_only))
    return out


def nonnormal_components(gene: str, dip: str):
    """[(component as written, function or "unknown")] of one caller's diplotype that argue against a Normal
    phenotype: a named allele rated Decreased or No function in CPIC's allele table, a named allele the table does
    not give (UNKNOWN_NAMED_ALLELE_BLOCKS_NORMAL), and any `+` tag that is not Normal function (decreased, no,
    increased, unknown). A diplotype that cannot be read (not `a/b`, not VCF-Check's allele list) is one unknown
    component. Empty only for a gene without a CPIC allele-function table or an empty call."""
    if _allele_functions(gene) is None or not (dip or "").strip() or dip.strip() == "-":
        return []
    walked = _walk_components(gene, dip, allow_vcf_form=True)
    if walked is None:
        return [(dip.strip(), "unknown")]
    out = []
    for written, _, _, fn, is_tag in walked:
        if is_tag and fn != "Normal function":
            out.append((written, fn or "unknown"))
        elif not is_tag and fn in ("Decreased function", "No function"):
            out.append((written, fn))
        elif not is_tag and fn is None and UNKNOWN_NAMED_ALLELE_BLOCKS_NORMAL:
            out.append((written, "unknown"))
    return list(dict.fromkeys(out))


def reading_function_profile(gene: str, reading: str):
    """(deviating, unknown) for one diplotype reading: `deviating` counts, per comparable allele key, the components
    whose CPIC function is not "Normal function" (decreased, no, increased); `unknown` lists the components whose
    function CPIC's table does not give (not in it, ambiguous, or Unknown/Uncertain function), or the whole
    reading when the gene has no function table or the reading is not an `a/b` pair."""
    walked = _walk_components(gene, reading) if _allele_functions(gene) is not None else None
    if walked is None:
        return collections.Counter(), [reading]
    deviating, unknown = collections.Counter(), []
    for _, c, key, fn, _ in walked:
        if fn is None:
            unknown.append(c)
        elif fn != "Normal function":
            deviating[key] += 1
    return deviating, unknown


# ── CYP2D6: one agreement rule for the verdict and for the PharmCAT outside call ─────
MAJORITY_MIN = 3   # owner decision 2026-11-15: callers sharing one CPIC activity score for a majority-by-score call
_PREFER = ("Cyrius", "PyPGx", "Stargazer", "StellarPGx", "Aldy")


def _fmt_score(x):
    return "" if x is None else f"{x:.2f}".rstrip("0").rstrip(".") if x != int(x) else f"{x:.1f}"


def _e_why(info, cy, same):
    others = [i for i in info if i is not cy and i not in same]
    note = ("; dissenting: " + "; ".join(
        f"{i['tool']}={i['diplotype']} (score {_fmt_score(i['score']) or 'unknown'}, "
        f"{i['own'] or i['class'] or 'no phenotype'})" for i in others)) if others else ""
    return "agreeing: " + ", ".join(i["tool"] for i in same) + note if cy["score"] is None else \
        "agreeing on the score: " + ", ".join(i["tool"] for i in same) + note


def cyp2d6_agreement(rows):
    """Level of agreement among CYP2D6 callers -> dict(diplotype, phenotype, score, level, source, why, callers).

    `rows`: dicts with tool, diplotype, phenotype, one per caller that typed the gene with a comparable call (the
    caller filters; see pgx_outside_calls.typed). Levels, first match wins: a all normalised diplotypes equal ->
    diplotype + CPIC score; b every score known and equal -> score only; e Cyrius has a scored call and >=1 other
    caller has the same score -> Cyrius's diplotype + score; m (after c fails) at least MAJORITY_MIN callers share one
    exact score -> that score only (not when Cyrius typed); c every caller has the same CPIC phenotype class
    (from its score, or from a `≥` lower bound, never from its own label when it has a score) -> phenotype only;
    otherwise no call, with the reason in `why`. The PharmCAT outside call and the verdict's CYP2D6 phenotype
    (when Cyrius gave none) both come from this function, so the gene row and the drug table cannot differ."""
    res = {"gene": "CYP2D6", "diplotype": None, "phenotype": None, "score": None, "level": "-",
           "source": "", "why": "", "callers": [], "bound": None}
    info = []
    for r in rows:
        sc = activity_score(r["diplotype"])
        lb = activity_lower_bound(r["diplotype"])
        own = phenotype_class(r.get("phenotype"))
        mapped = score_phenotype(sc) if sc is not None else bound_phenotype(lb)
        info.append({"tool": r["tool"], "diplotype": r["diplotype"], "key": norm_diplotype("CYP2D6", clean_diplotype("CYP2D6", r["diplotype"])),  # *36+*36 = *36x2
                     "score": sc, "bound": lb, "own": own,
                     "class": phenotype_class(mapped) if mapped else own})
    res["callers"] = [{"tool": i["tool"], "diplotype": i["diplotype"], "score": i["score"],
                       "phenotype": i["class"]} for i in info]
    n = len(info)
    desc = "; ".join(f"{i['tool']}={i['diplotype']} (score {_fmt_score(i['score']) or 'unknown'}, "
                     f"{i['class'] or 'no phenotype'})" for i in info)
    if len({i["key"] for i in info}) == 1:
        pick = min(info, key=lambda i: _PREFER.index(i["tool"]) if i["tool"] in _PREFER else len(_PREFER))
        res.update(level="a", source="all callers agree", diplotype=clean_diplotype("CYP2D6", pick["diplotype"]),
                   score=pick["score"], bound=pick["bound"])
        return res
    scores = {i["score"] for i in info}
    if None not in scores and len(scores) == 1:
        res.update(level="b", source=f"activity score agreed {n}/{n}", score=scores.pop())
        return res
    # e: Cyrius (the read-depth copy-number caller) has a call with a CPIC activity score and at least one other
    # caller's call has the same score. Precedence: a, then b, then e, then c.
    cy = next((i for i in info if i["tool"] == "Cyrius" and i["score"] is not None), None)
    if cy is not None:
        same = [i for i in info if i is not cy and i["score"] == cy["score"]]
        if same:
            res.update(level="e", source=f"Cyrius + {len(same)} other caller(s) on activity score "
                                         f"{_fmt_score(cy['score'])}",
                       diplotype=clean_diplotype("CYP2D6", cy["diplotype"]), score=cy["score"],
                       why=_e_why(info, cy, same))
            return res
    # e (lower bound): Cyrius's call has only a `≥` bound (`*1x3/*10`, CPIC `≥3.25`), which has no exact value to
    # compare, so another caller must report the same canonical diplotype. The diplotype and its bound phenotype are
    # passed, never a number (mirrors pgx-compare._cyp2d6_authority_phenotype).
    cyb = next((i for i in info if i["tool"] == "Cyrius" and i["score"] is None and i["bound"] is not None), None)
    if cyb is not None:
        same = [i for i in info if i is not cyb and i["key"] == cyb["key"]]
        if same and bound_phenotype(cyb["bound"]):
            res.update(level="e", source=f"Cyrius + {len(same)} other caller(s) on diplotype {cyb['diplotype']} "
                                         f"(activity score >= {_fmt_score(cyb['bound'])})",
                       diplotype=clean_diplotype("CYP2D6", cyb["diplotype"]), bound=cyb["bound"],
                       why=_e_why(info, cyb, same))
            return res
    fail = None
    if any("Indeterminate" in (i["own"], i["class"]) for i in info):
        fail = "a caller reports Indeterminate (blocks phenotype agreement): " + desc
    # Every caller's call needs a CPIC activity score or lower bound: the phenotype is then CPIC's for it, never a
    # caller's own label, and no caller is left out of the comparison.
    elif any(i["score"] is None and i["bound"] is None for i in info):
        fail = "callers disagree and a call has no CPIC activity score: " + desc
    elif len({i["class"] for i in info}) == 1 and info[0]["class"]:
        res.update(level="c", source=f"phenotype agreed {n}/{n}", phenotype=f"{info[0]['class']} Metabolizer")
        return res
    else:
        fail = "callers disagree: " + desc
    # m: majority by activity score (owner decision 2026-11-15): at least MAJORITY_MIN callers with an exact CPIC
    # score (a lettered legacy name has none) share one score; that score's CPIC phenotype is used. The agreeing and
    # dissenting callers are named in `why`.
    best = collections.Counter(i["score"] for i in info if i["score"] is not None).most_common(1)
    # Not when Cyrius (the authority) typed the gene: its call is the verdict then and levels a/b/e decide the drug table.
    # An explicit Indeterminate from a caller that has a score blocks m as it blocks c.
    blocked = any(i["own"] == "Indeterminate" and (i["score"] is not None or i["bound"] is not None) for i in info)
    if best and best[0][1] >= MAJORITY_MIN and not blocked and not any(i["tool"] == "Cyrius" for i in info):
        sc = best[0][0]
        agree = [i for i in info if i["score"] == sc]
        others = [i for i in info if i["score"] != sc]
        res.update(level="m", source=f"majority by activity score {_fmt_score(sc)} ({len(agree)}/{n} callers)",
                   score=sc, why="agreeing on the score: " + ", ".join(i["tool"] for i in agree)
                   + ("; dissenting: " + "; ".join(
                       f"{i['tool']}={i['diplotype']} (score {_fmt_score(i['score']) or 'unknown'}, "
                       f"{i['own'] or i['class'] or 'no phenotype'})" for i in others) if others else ""))
        return res
    res["why"] = fail
    return res


def agreement_phenotype(res):
    """CPIC phenotype of a cyp2d6_agreement result (its phenotype, else that of its score or lower bound), or None."""
    if res.get("phenotype"):
        return res["phenotype"]
    if res.get("score") is not None:
        return score_phenotype(res["score"])
    return bound_phenotype(res.get("bound")) if res.get("level") in ("a", "e") else None


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
    # CYP2C19 *38 is its own (normal function) allele in CPIC, never a *2 synonym
    assert canonicalize_diplotype("CYP2C19", "*2/*38") != canonicalize_diplotype("CYP2C19", "*2/*2")
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
