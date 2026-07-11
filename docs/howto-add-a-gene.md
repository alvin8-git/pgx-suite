# How to add or change a gene

You want the pipeline to call a new gene, or change which tools run for an existing
one. Everything is driven by one file — `docker/genes.tsv` — so this is a one-row
edit plus a rebuild. No Python or Snakefile changes needed for a standard star-allele
gene.

## Prerequisites

- A working pgx-suite checkout (see the [tutorial](tutorial-getting-started.md)).
- The gene's GRCh38 coordinates (`chr:start-end`) and which callers support it.
- Python 3.10+ on the host (to run the config regression test).

## The single source of truth

`docker/genes.tsv` is a tab-separated table; both the Snakemake DAG and
`pgx-compare.py` read it. Columns:

| Column | Meaning |
|--------|---------|
| `gene` | Gene symbol (uppercase; hyphens allowed, e.g. `HLA-A`) |
| `region` | `chr:start-end` on GRCh38, or `ALT_CONTIG` for alt-contig genes (GSTT1) |
| `pypgx` `stargazer` `aldy` `stellarpgx` `optitype` `mutserve` `arcashla` | `1` if that tool calls this gene, else `0` (`arcashla` = HLA class-II typing for DQA1/DQB1/DRB1, runs once per sample like PharmCAT) |
| `vcf_check` | `1` if `pgx-compare.py` runs VCF-Check (direct genotype from the CPIC variant table — RYR1/VKORC1/etc.) |
| `pypgx_sv` | `1` if PyPGx needs SV preprocessing (depth-of-coverage + VDR control) |
| `stargazer_sv` | `1` if Stargazer needs a GDF depth profile (paralog CN) |
| `stargazer_control` | Stargazer control gene (`vdr`), or `-` if none |
| `cyrius` | `1` to run the Cyrius CYP2D6 SV caller (only CYP2D6 today) |
| `pharmcat` | `1` if PharmCAT is authoritative for this gene (UGT1A1/CYP2B6/CYP4F2) |

Column order in the file is: `gene region pypgx stargazer aldy stellarpgx optitype
mutserve vcf_check pypgx_sv stargazer_sv stargazer_control cyrius pharmcat arcashla`.

**Authority wiring (only if a CPIC-reference method should *override* the star-allele
vote).** Setting `cyrius`/`pharmcat`/`vcf_check` to `1` runs the method; to make it the
**verdict** for that gene, also add the gene to the `AUTHORITATIVE` dict in
`docker/pgx-compare.py` (`"GENE": "Cyrius" | "PharmCAT" | "VCF-Check"`). For VCF-Check
genes, add the variant(s) to `docker/vcf_check_variants.json` (guarded by
`test_vcf_check.py`). See [`PGxDocumentation.md`](PGxDocumentation.md#reconciliation--verdict-authority).

## Steps

1. **Add (or edit) the row.** Append a tab-separated line to `docker/genes.tsv`.
   Example — a hypothetical SNP-only gene called by all four star-allele callers:

   ```
   MYGENE	chr1:1000000-1010000	1	1	1	1	0	0	0	0	0	vdr	0	0
   ```

   (The two trailing `0`s are the `cyrius` and `pharmcat` columns.) For a gene with a
   structural-variant component (like CYP2D6), set `pypgx_sv`/`stargazer_sv` to `1`. For
   HLA, set only `optitype=1` and use the MHC region.

2. **Run the config regression test.** This proves the table still parses and that the
   matrix the loaders produce matches what you intended:

   ```bash
   python3 docker/test_genes_config.py
   ```

   It freezes the known-good matrix; if you only *added* a gene, update the frozen
   `ORIG_SUPPORT` / `ORIG_COORDS` dicts in that test to include your new row, then re-run
   until it prints `all genes.tsv regression checks passed`.

3. **Rebuild the image** so the new `genes.tsv` ships:

   ```bash
   docker build -t pgx-suite:latest .
   ```

4. **Run just your gene** to confirm it calls:

   ```bash
   docker run --privileged --rm \
     -v "$(pwd)/StellarPGx:/pgx/stellarpgx" -v "$(pwd)/StellarPGx/containers:/pgx/containers" \
     -v "$(pwd)/pypgx/pypgx-bundle:/pgx/bundle" -v "/path/to/GRCh38:/pgx/ref" \
     -v "/path/to/data:/pgx/data" -v "$(pwd)/results/mygene:/pgx/results" \
     pgx-suite:latest \
     snakemake -s /opt/pgx/Snakefile --cores 4 \
       --config bam=/pgx/data/sample.bam outdir=/pgx/results genes="MYGENE"
   ```

## Verification

```bash
# the verdict block lands in the per-gene detail JSON
python3 -c "import json,glob; d=json.load(open(glob.glob('results/mygene/genes/MYGENE/*_detail.json')[0])); print(d['verdict'])"
```

You should see a `verdict` object with a `status` of `concordant` / `majority` /
`discordant` / `no_call`. The gene also appears in `results/mygene/all_genes_summary.tsv`
and as a card in the HTML report.

## Troubleshooting

| Symptom | Cause / fix |
|---------|-------------|
| `KeyError: 'MYGENE'` from a loader | The row is space-separated, not tab-separated. `genes.tsv` is strict TSV. |
| `test_genes_config.py` fails after adding a gene | Expected until you add the new gene to the frozen `ORIG_*` dicts in the test. |
| Gene always `NO_CALL` | `region` has no/low coverage in your sample, or the coordinates are wrong. Verify with `samtools coverage -r chr:start-end sample.bam`. |
| New gene needs clinical interpretation in the report | Drug/phenotype text comes from the CPIC database in `pgx-report.py`; add an entry there for full clinical reporting. The diplotype call itself works from `genes.tsv` alone. |

## Related

- [Tutorial: from a BAM to your first report](tutorial-getting-started.md)
- [Per-tool reference](ToolsDocumentation.md) — which tools support which genes today
- [README → Architecture](../README.md) — why config lives in one place
