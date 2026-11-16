# PharmCAT outside calls (HLA-A/B, CYP2D6) and pgx.report.json

Moved from README.md on 2026-11-14 without change of content, except the formatting of the level table. Rule: `pharmcat_outside` in `docker/Snakefile`; code: `docker/pgx_outside_calls.py`, `docker/pgx-pharmcat-outside.sh`.

## Outside calls

PharmCAT 3.2.0 does not type HLA from a VCF, and
`pharmcat_pipeline` has no outside-call option (only the core tool does:
`java -jar pharmcat.jar -po <file>`, [Outside Call Format](https://pharmcat.clinpgx.org/using/Outside-Call-Format/)).
After the normal PharmCAT run, rule `pharmcat_outside` (`pgx-pharmcat-outside.sh`) writes a tab-separated
`pharmcat/outside_calls.pharmcat.tsv` (`gene<TAB>diplotype<TAB>phenotype<TAB>activity score`, empty columns allowed) and re-runs PharmCAT on the already
preprocessed VCF with `-po`:
- **HLA-A / HLA-B**: the OptiType diplotype at two-field resolution, e.g. `HLA-B	*15:02/*55:01`.
- **CYP2D6**, by level of agreement among the callers that typed it (at least two with a comparable
  call; notation is normalised first by `reconcile.norm_diplotype`: whitespace and case; Aldy's tandem
  `X` marker; `.NNN` and CYP2D6 letter sub-alleles such as `*6B`; allele order inside a tandem and
  between haplotypes; `+rsNNN` tags, copy numbers and different tandem structures stay real
  differences). Levels are tried in the order a, b, e, c, d and the first match wins. Level e (added 2026-11-13) is not a majority vote: it needs Cyrius, the read-depth copy-number caller, and one other caller to agree on the CPIC activity score. Other callers' disagreement, including an explicit Indeterminate, does not block it and is recorded:

  | Level | Condition | Passed to PharmCAT | `source` in the record |
  |---|---|---|---|
  | a | all normalised diplotypes equal | diplotype + activity score (score omitted if unknown) | `all callers agree` |
  | b | diplotypes differ, every caller has a known activity score, all equal | activity score only | `activity score agreed n/n` |
  | e | not a or b; Cyrius has a call with a known CPIC activity score and at least one other caller's call has the same score | Cyrius's cleaned diplotype + that score | `Cyrius + N other caller(s) on activity score X` (dissenting callers listed in the reason column and provenance `note`) |
  | c | scores differ, every caller has a known activity score, all in the same CPIC phenotype | phenotype only | `phenotype agreed n/n` |
  | d | otherwise | nothing (omitted, reason logged) | |

  The activity score is the sum of the two haplotype values of the CPIC allele table
  (`docker/cpic_cyp2d6.json`, fetched from the CPIC API on 2026-11-12 by `docker/cpic_cyp2d6_build.py`;
  a tandem is looked up in either order, e.g. `*36+*10` = 0.25, `*68+*4` = 0.0). It is unknown, never
  guessed, for an allele absent from the table (e.g. `*79`), an allele with no numeric value, an
  Aldy `+rsNNN` tag or an unknown copy number. The phenotype class always comes from the score through
  CPIC's activity-score table (`gene_result`); a caller's own phenotype label is never passed on.
  Level e keeps level b's precedence (when every score is equal the score-only call is unchanged) and
  goes before level c (a diplotype with a score is more specific than a phenotype). It never applies when
  Cyrius made no call or its call has no CPIC score (e.g. `*134`, `*1/*36+*36+*10`): b, c and d then run as before.
  Level c needs a known score for every caller: no caller is left out of the comparison. A caller
  that reports Indeterminate or unknown agrees with nothing and blocks level c.
  Level a passes the cleaned spelling of one agreeing caller (Cyrius preferred). The file may
  therefore have empty columns: `CYP2D6<TAB><TAB><TAB>0.5` (score only) and
  `CYP2D6<TAB><TAB>Intermediate Metabolizer` (phenotype only) are accepted by PharmCAT 3.2.0.
  **Caution (level c):** PharmCAT answers a phenotype-only call with the annotations of every
  activity score of that phenotype (for Intermediate: 0.25, 0.5, 0.75, 1.0), so a drug can carry
  several annotations. Downstream must treat differing recommendations among them as not assessed.
  In the CPIC recommendation API (retrieved 2026-11-12) the recommendation text, classification or
  implication differs between activity scores of the Intermediate phenotype for atomoxetine,
  nortriptyline and tamoxifen.
- **Fallback**: a missing input or an uncallable gene is omitted and logged; the stage-1 PharmCAT
  step is unchanged and its report is kept as `pharmcat/pgx.report.no_outside.json`. Nothing is
  fabricated, and the rule never fails the run.
- **Record**: `pharmcat/outside_calls.tsv` (gene, diplotype, activity score, phenotype, level, source,
  passed, reason omitted, per-caller `tool=diplotype|score|phenotype`) and `pharmcat_outside_calls`
  in `provenance.json` (same fields, per-caller list).
Re-run only this step for an existing sample (inside the image): `pgx-pharmcat-outside.sh /pgx/results/<sample>`.

## `pgx.report.json` as a downstream contract

Beyond the three PharmCAT-authority
genes, `pharmcat/pgx.report.json` carries PharmCAT's **full per-drug CPIC recommendations** for
the sample — under `drugs["CPIC Guideline Annotation"][<drug>].guidelines[].annotations[]`, each
with `drugRecommendation`, `implications`, `phenotypes`, `classification` (evidence strength) and
`alternateDrugAvailable`. This is the authoritative source for a **medication-response-by-drug**
view (drugs grouped by therapeutic class, each with a dose/caution action). Downstream consumers
(e.g. the OmniGen reporting layer) read this file directly and **render** it — they do not re-call
PGx. Keep the schema stable; it is a consumed contract, not just an internal artifact.
