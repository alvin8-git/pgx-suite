# PGx Suite Docker Container

Single container for nine pharmacogenomics callers, orchestrated by Snakemake.
All tools configured for **GRCh38**.

## License Notice

Stargazer (UW) and Aldy (IURTC) are **non-commercial academic use only**.
This image must NOT be pushed to any public registry.

## Prerequisites

- Docker with `--privileged` support (local workstation)
- `StellarPGx/containers/stellarpgx-dev.sif` present (already in repo)
- ~3 GB disk for the built image
- GRCh38 reference FASTA (user-supplied, for StellarPGx full runs)

## Build

```bash
docker build -t pgx-suite:latest .   # run from the repo root
```

First build takes ~15-20 minutes. Subsequent builds use layer cache.

## Smoke Tests

```bash
# PyPGx, Aldy, Stargazer (fast, no volumes needed):
docker run --rm pgx-suite:latest bash /opt/pgx/test.sh

# Full suite including StellarPGx (requires volumes):
./docker/docker-run.sh -- bash /opt/pgx/test.sh
```

## Tools Installed

| Tool | Version | Invocation |
|------|---------|-----------|
| PyPGx | 0.26.0 | `pypgx` |
| Stargazer | 2.0.3 | `stargazer` |
| Aldy | 4.8.3 | `aldy` |
| StellarPGx | 1.2.7 | `nextflow run /pgx/stellarpgx/main.nf` |
| Cyrius | vendored | `cyrius` rule → `/opt/cyrius/star_caller.py` (CYP2D6; authority) |
| PharmCAT | 3.2.0 | `pharmcat` rule → `pharmcat_pipeline` in `pharmcat.sif` (UGT1A1/CYP2B6/CYP4F2; authority) |
| OptiType | 1.3.5 | `optitype` rule (Apptainer SIF, HLA-A/B) |
| mutserve | 2.0.0 | `mutserve` rule (baked-in JAR, MT-RNR1) |
| VCF-Check | in-house | `pgx-compare.py` direct genotype from `vcf_check_variants.json` (RYR1/VKORC1/IFNL3/G6PD/CACNA1S; authority) |
| Snakemake | orchestrator | `snakemake -s /opt/pgx/Snakefile` |

## Volume Mounts (at runtime)

| Host Path | Container Path | Purpose |
|-----------|---------------|---------|
| `./StellarPGx` | `/pgx/stellarpgx` | StellarPGx pipeline scripts, DB, resources |
| `./StellarPGx/containers` | `/pgx/containers` | Apptainer SIFs: `stellarpgx-dev.sif`, `optitype.sif`, `pharmcat.sif` |
| `./pypgx/pypgx-bundle` | `/pgx/bundle` | 1KGP VCFs + CNV data for PyPGx phasing |
| `/path/to/ref` | `/pgx/ref` | GRCh38 reference FASTA + .fai index |
| `/path/to/data` | `/pgx/data` | Input BAM/CRAM files |
| `/path/to/results` | `/pgx/results` | Output directory |

## Why --privileged?

StellarPGx uses Nextflow + Apptainer (Singularity) to run its pipeline steps
inside `stellarpgx-dev.sif`. Apptainer requires kernel namespace support
(`SYS_ADMIN` capability) to unpack SIF overlay filesystems.

## Running the pipeline

The orchestrator is Snakemake (`/opt/pgx/Snakefile`), driven by `genes.tsv`:

```bash
snakemake -s /opt/pgx/Snakefile --cores 4 \
  --config bam=/pgx/data/sample.bam outdir=/pgx/results
# single gene (debug): add  genes="CYP2D6"
```

Each caller runs as a non-fatal rule; `pgx-compare.py` reconciles them into one
verdict per gene (with a NO_CALL coverage gate) and `pgx-report.py` writes the
self-contained HTML report. See the top-level [README](../README.md) for full usage.

## Troubleshooting

| Symptom | Cause | Fix |
|---------|-------|-----|
| `FATAL: could not open image` | SIF not mounted | Add `-v $(pwd)/StellarPGx/containers:/pgx/containers` |
| `FATAL: kernel too old` | Need privilege | Run with `--privileged` |
| Nextflow hangs at startup | JAR download | Ensure internet access during first run |
| `beagle.jar` error | Java not in PATH | Verify `java -version` in container |
