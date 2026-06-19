# Sketch: PharmCAT cross-check column (informational)

**Status:** design sketch (not implemented). **Effort:** ~1–2 h.
**Goal:** surface PharmCAT's call for *every* gene it covers as an **informational
cross-check** in the report — without making it the verdict authority for those genes.

## Why it's near-free

PharmCAT already runs **once per sample** (the `pharmcat` rule) and writes a single
`results/<sample>/pharmcat/*.report.json` containing **all** its CPIC genes
(~20: CYP2C9, CYP2C19, CYP3A4/5, CYP2B6, CYP4F2, DPYD, NUDT15, SLCO1B1, TPMT, UGT1A1,
VKORC1, G6PD, IFNL3, CACNA1S, RYR1, …). We currently consume only 3 (UGT1A1/CYP2B6/CYP4F2,
where PharmCAT is the authority). The expensive part is already paid for; reading more
genes out of the same JSON costs nothing.

`parse_pharmcat(output_dir, gene, sample)` already reads that shared report and looks up
`genes[gene]`, returning an empty/`failed` result when the gene isn't present. So the
cross-check needs **no new parser and no new DAG edges** — just read it opportunistically.

## Design principles

1. **Informational only — never changes the verdict.** Keeps the single-authority model
   intact. Cross-check lives in its own `detail.json` field, not in `tools`/the vote.
2. **Opportunistic read.** If `pharmcat/*.report.json` exists (it does for any whole-sample
   run), use it; else skip. No new rule dependencies; single-gene debug runs simply omit it.
3. **Compare on phenotype, not raw string.** "Agrees?" reuses the phenotype normaliser
   (`_norm_phenotype`) so PharmCAT's nomenclature doesn't spuriously "disagree."
4. **Exclude where it doesn't apply.** Skip genes already PharmCAT-authoritative (it's the
   verdict there, not a cross-check) and **CYP2D6** (PharmCAT doesn't call CYP2D6
   structurally — that's Cyrius).
5. **A disagreement is a signal, not an error.** It flags where the star-allele consensus
   and the CPIC reference matcher diverge — the genes worth manual review or future
   authority promotion (DPYD is the prime candidate).

## Code shape

### `pgx-compare.py` — attach cross-check to detail.json

```python
def pharmcat_cross_check(output_dir, gene, sample, verdict):
    """Informational PharmCAT call for `gene` from the sample-level report.
    Returns None when PharmCAT is already this gene's authority, doesn't cover
    the gene, or no report exists. Never feeds the verdict."""
    if AUTHORITATIVE.get(gene) == "PharmCAT" or gene == "CYP2D6":
        return None
    r = parse_pharmcat(output_dir, gene, sample)          # existing parser
    if r.status != "ok" or r.diplotype in ("-", "", "no PharmCAT report"):
        return None
    agrees = (_norm_phenotype(r.phenotype) is not None
              and _norm_phenotype(r.phenotype)
                  == _norm_phenotype(verdict.get("consensus_phenotype", "")))
    return {"diplotype": r.diplotype, "phenotype": r.phenotype, "agrees": agrees}

# in the detail.json assembly, after the verdict is computed:
cc = pharmcat_cross_check(output_dir, gene, sample, verdict)
detail["cross_check"] = {"PharmCAT": cc} if cc else {}
```

### `pgx-report.py` — show it (landing chip + detail row)

```python
# gene-data dict
gd["cross_check"] = detail.get("cross_check", {})

# landing card (only when present): a small chip, green if agrees, amber if differs
cc = gd["cross_check"].get("PharmCAT")
if cc:
    cls, mark = ("badge-green", "✓") if cc["agrees"] else ("badge-amber", "⚠ differs — review")
    xcheck_html = f'<div class="gene-xcheck"><span class="badge {cls}">PharmCAT {cc["diplotype"]} {mark}</span></div>'
```
Detail page: one row "PharmCAT (cross-check — informational, not part of the verdict):
`<diplotype>` / `<phenotype>` — agrees / differs".

Legend note: "Cross-check = an independent CPIC-reference call shown for confidence; it does
not change the verdict."

## Tests

- `test_pharmcat.py`: `pharmcat_cross_check` returns `None` for an authoritative gene and for
  CYP2D6; returns `{agrees: True}` when phenotypes match, `{agrees: False}` when they differ.
- `test_verdict.py` / report: a gene with a differing cross-check still keeps its own verdict
  (cross-check never overrides).

## Follow-on (separate decision)

Once the cross-check column exists and we can see where PharmCAT diverges across a cohort,
**promote to authority** only the genes where it demonstrably beats the current method — the
strongest candidate is **DPYD** (rigorous CPIC matcher vs the phenotype-tier heuristic). That
is a per-gene `AUTHORITATIVE` decision, made on evidence, not a blanket flip.
