// pgx-suite processes: one per Snakefile rule, same commands, same container paths.
// Every process merges its upstream result trees (inputs staged as up*/res) into ./res,
// points /pgx/results at it (bin/nf_stage.sh) and emits ./res. OUT is therefore
// /pgx/results everywhere, exactly as in docker/Snakefile.
// Caller processes always exit 0 and write the tool's exit code to a .status file
// (a failed tool is reported as failed, never as a negative); a caller killed by the infrastructure
// (134/137/143 killed, 140 time) exits non-zero (bin/nf_kill_gate.sh): retried with more memory, then the
// run stops. The hard-dependency steps stop the run.

// Every process takes the sample map `s` [id, bam, ref, bam_glob, depth_bam, sex, min_depth] as a
// declared input so it is part of the task hash. (Values from helper functions are not hashed:
// an earlier version reused one sample's cached result for another.)
OUT = '/pgx/results'

process PREFLIGHT {
    // the Snakefile's _startup_checks(): truncated input or a contig missing from ref.fai stops the run
    label 'process_low'
    input:
    val s
    val genes
    output: val s
    script:
    """
    samtools quickcheck ${s.bam}
    for r in ${genes.findAll { it.region != 'ALT_CONTIG' }.collect { it.region.tokenize(':')[0] }.unique().join(' ')}; do
        cut -f1 ${s.ref}.fai | grep -qx "\$r" || { echo "contig '\$r' not found in ${s.ref}.fai" >&2; exit 1; }
    done
    """
}

process BCFTOOLS {
    tag "${g.gene}"
    input:  tuple val(g), val(s)
    output: tuple val(g), path('res')
    script:
    def G = "${OUT}/genes/${g.gene}"
    """
    nf_stage.sh
    mkdir -p ${G}/logs; rc=0
    ( set -o pipefail; bcftools mpileup -r ${g.region} -f ${s.ref} -a AD,DP --max-depth 500 -o - ${s.bam} \
      | bcftools call -m -v --output-type z -o ${G}/${g.gene}.vcf.gz && tabix -p vcf ${G}/${g.gene}.vcf.gz ) || rc=\$?
    echo \$rc > ${G}/logs/bcftools.status
    nf_kill_gate.sh \$rc
    """
}

process ALDY {
    tag "${g.gene}"
    cpus 2
    memory { 8.GB * task.attempt }
    input:  tuple val(g), val(s)
    output: tuple val(g), path('res')
    script:
    def G = "${OUT}/genes/${g.gene}"
    """
    nf_stage.sh
    mkdir -p ${G}/logs ${G}/aldy; rc=0
    aldy genotype -g ${g.gene} -p illumina ${s.bam} -o ${G}/aldy/${g.gene}.aldy || rc=\$?
    echo \$rc > ${G}/logs/aldy.status
    nf_kill_gate.sh \$rc
    """
}

process STELLARPGX {
    tag "${g.gene}"
    label 'nested_apptainer'
    cpus 2
    memory { 4.GB * task.attempt }
    input:  tuple val(g), val(s)
    output: tuple val(g), path('res')
    script:
    def G = "${OUT}/genes/${g.gene}"
    def sgene = g.gene.toLowerCase() == 'por' ? 'cypor' : g.gene.toLowerCase()
    """
    nf_stage.sh
    mkdir -p ${G}/logs ${G}/stellarpgx; rc=0
    cd /pgx   # StellarPGx's nextflow.config finds its SIF at \$PWD/containers; /pgx is the Snakemake cwd
    nextflow run /pgx/stellarpgx/main.nf --gene ${sgene} \
        --in_bam '${s.bam_glob}' --ref_file ${s.ref} --out_dir ${G}/stellarpgx \
        --res_init /pgx/stellarpgx/resources --db_init /pgx/stellarpgx/database \
        --caller_init /pgx/stellarpgx/scripts -work-dir ${G}/stellarpgx/.work || rc=\$?
    # the nested nextflow returns 1 when one of its tasks was killed: keep the kill code of that task
    if [ \$rc -ne 0 ]; then k=\$(find ${G}/stellarpgx/.work -name .exitcode 2>/dev/null | xargs -r awk 'FNR==1{print}' | grep -m1 -x -E '134|137|140|143' || true); [ -n "\$k" ] && rc=\$k; fi
    rm -rf ${G}/stellarpgx/.work 2>/dev/null || true
    echo \$rc > ${G}/logs/stellarpgx.status
    nf_kill_gate.sh \$rc
    """
}

process PYPGX_DEPTH {   // hard dependency: a failure stops the run
    tag "${g.gene}"
    input:  tuple val(g), val(s)
    output: tuple val(g), path('res')
    script:
    def G = "${OUT}/genes/${g.gene}"
    """
    nf_stage.sh
    mkdir -p ${G}; rm -f ${G}/depth-of-coverage.zip
    pypgx prepare-depth-of-coverage ${G}/depth-of-coverage.zip ${s.bam} --assembly GRCh38
    """
}

process PYPGX_CONTROL {   // hard dependency
    tag "${g.gene}"
    input:  tuple val(g), val(s)
    output: tuple val(g), path('res')
    script:
    def G = "${OUT}/genes/${g.gene}"
    """
    nf_stage.sh
    mkdir -p ${G}; rm -f ${G}/control-stats-VDR.zip
    pypgx compute-control-statistics VDR ${G}/control-stats-VDR.zip ${s.bam} --assembly GRCh38
    """
}

process STARGAZER_GDF {   // hard dependency
    tag "${g.gene}"
    input:  tuple val(g), val(s)
    output: tuple val(g), path('res')
    script:
    def G = "${OUT}/genes/${g.gene}"
    def gl = g.gene.toLowerCase()
    """
    nf_stage.sh
    mkdir -p ${G}
    stargazer -G ${gl}.gdf -t ${gl} -c vdr -B ${s.bam} -o ${G}/stargazer_gdf -a grc38 -d wgs
    """
}

process PYPGX {
    tag "${g.gene}"
    cpus 2
    memory { 6.GB * task.attempt }
    input:  tuple val(g), val(s), path(ups, stageAs: 'up*/*')
    output: tuple val(g), path('res')
    script:
    def G = "${OUT}/genes/${g.gene}"
    def variants = g.needs_vcf ? "--variants ${G}/${g.gene}.vcf.gz" : ''
    def sv = g.pypgx_sv ? "--depth-of-coverage ${G}/depth-of-coverage.zip --control-statistics ${G}/control-stats-VDR.zip" : ''
    """
    nf_stage.sh
    mkdir -p ${G}/logs; rc=0
    pypgx run-ngs-pipeline ${g.gene} ${G}/pypgx --assembly GRCh38 --force ${variants} ${sv} || rc=\$?
    echo \$rc > ${G}/logs/pypgx.status
    nf_kill_gate.sh \$rc
    """
}

process STARGAZER {
    tag "${g.gene}"
    memory { 4.GB * task.attempt }
    input:  tuple val(g), val(s), path(ups, stageAs: 'up*/*')
    output: tuple val(g), path('res')
    script:
    def G = "${OUT}/genes/${g.gene}"
    def gl = g.gene.toLowerCase()
    def mode = g.stargazer_sv ? "-c vdr -g ${G}/stargazer_gdf/${gl}.gdf"
             : (g.stargazer_control != '-' ? "-c ${g.stargazer_control}" : '')
    """
    nf_stage.sh
    mkdir -p ${G}/logs; rc=0
    stargazer -t ${gl} -d wgs -a grc38 -i ${G}/${g.gene}.vcf.gz -o ${G}/stargazer ${mode} || rc=\$?
    echo \$rc > ${G}/logs/stargazer.status
    nf_kill_gate.sh \$rc
    """
}

process CYRIUS {
    tag "${g.gene}"
    cpus 2
    memory { 4.GB * task.attempt }
    input:  tuple val(g), val(s)
    output: tuple val(g), path('res')
    script:
    def G = "${OUT}/genes/${g.gene}"
    """
    nf_stage.sh
    mkdir -p ${G}/cyrius ${G}/logs; rc=0
    echo ${s.bam} > ${G}/cyrius/manifest.txt
    python3 /opt/cyrius/star_caller.py --manifest ${G}/cyrius/manifest.txt \
        --genome 38 --reference ${s.ref} --prefix ${s.id} --outDir ${G}/cyrius \
        --threads ${task.cpus} || rc=\$?
    echo \$rc > ${G}/logs/cyrius.status
    nf_kill_gate.sh \$rc
    """
}

process OPTITYPE {   // sample level: one run types HLA-A/B/C
    label 'nested_apptainer'
    cpus 2
    memory { 4.GB * task.attempt }
    input:  val s
    output: path 'res'
    script:
    """
    nf_stage.sh
    mkdir -p ${OUT}/optitype; rc=0
    pgx-hla.sh ${s.bam} ${s.id} ${OUT} || rc=\$?
    echo \$rc > ${OUT}/optitype/optitype.status
    nf_kill_gate.sh \$rc
    """
}

process ABO_PILEUP {   // ABO from read-level evidence; the declared input of the ABO compare
    tag "${g.gene}"
    input:  tuple val(g), val(s)
    output: tuple val(g), path('res')
    script:
    def G = "${OUT}/genes/${g.gene}"
    """
    nf_stage.sh
    mkdir -p ${G}; python3 /opt/pgx/abo_pileup.py ${s.bam} ${s.ref} ${G}/abo_pileup.json
    """
}

process MUTSERVE {
    tag "${g.gene}"
    input:  tuple val(g), val(s)
    output: tuple val(g), path('res')
    script:
    def G = "${OUT}/genes/${g.gene}"
    """
    nf_stage.sh
    mkdir -p ${G}/logs; rc=0
    pgx-mt.sh ${s.bam} ${s.id} ${G} --ref ${s.ref} || rc=\$?
    echo \$rc > ${G}/logs/mutserve.status
    nf_kill_gate.sh \$rc
    """
}

process HLA2 {   // sample level: one arcasHLA run types DQA1/DQB1/DRB1
    label 'nested_apptainer'
    cpus 4
    memory { 8.GB * task.attempt }
    input:  val s
    output: path 'res'
    script:
    def h = "${OUT}/hla2"
    """
    nf_stage.sh
    mkdir -p ${h}; rc=0
    pgx-hla2.sh ${s.bam} ${s.id} ${h} --ref ${s.ref} --sif /pgx/containers/arcashla.sif || rc=\$?
    echo \$rc > ${h}/hla2.status
    nf_kill_gate.sh \$rc   # killed (134/137/143/140): the task is retried or the run stops, so the next line only sees ordinary failures
    if [ \$rc -ne 0 ]; then echo "ERROR: HLA class II (arcasHLA) typing FAILED rc=\$rc -- DQA1/DQB1/DRB1 will be reported as failed, not clean. See hla2/ logs." >&2; fi
    """
}

process PHARMCAT {   // sample level
    label 'nested_apptainer'
    cpus 4
    memory { 6.GB * task.attempt }
    input:  val s
    output: path 'res'
    script:
    def pc = "${OUT}/pharmcat"
    """
    nf_stage.sh
    mkdir -p ${pc}; rc=0
    ( set -o pipefail; bcftools mpileup -R /opt/pgx/pharmcat_positions.vcf.bgz -f ${s.ref} --max-depth 500 -Ou ${s.bam} 2>/dev/null \
      | bcftools call -m --output-type z -o ${pc}/pgx.vcf.gz \
      && tabix -f -p vcf ${pc}/pgx.vcf.gz \
      && apptainer exec --cleanenv --bind ${pc}:${pc} /pgx/containers/pharmcat.sif \
         pharmcat_pipeline ${pc}/pgx.vcf.gz --missing-to-ref -reporterJson -o ${pc} ) || rc=\$?
    echo \$rc > ${pc}/pharmcat.status
    nf_kill_gate.sh \$rc
    """
}

process ALTCHECK {   // hard dependency
    input:  val s
    output: path 'res'
    script:
    """
    nf_stage.sh
    python3 /opt/pgx/pgx_altcheck.py ${s.bam} > ${OUT}/alt_aware.json
    """
}

process COMPARE {
    tag "${g.gene}"
    input:  tuple val(g), val(s), path(ups, stageAs: 'up*/*')
    output: tuple val(g), path('res')
    script:
    """
    nf_stage.sh
    mkdir -p ${OUT}/genes/${g.gene}   # Snakemake creates a rule's output folder before it runs
    pgx_compare_wrap.py ${g.gene} ${g.region} ${g.is_alt ? 1 : 0} ${s.id} ${s.min_depth} \
        ${s.bam} ${s.ref} ${OUT} ${s.sex ?: ''}
    """
}

process PROVENANCE {
    input:
    val s
    val n_genes
    output: path 'res'
    script:
    """
    nf_stage.sh
    write_provenance.py ${OUT} ${s.id} ${s.ref} ${n_genes} ${workflow.nextflow.version}
    """
}

process PHARMCAT_OUTSIDE {   // non-fatal; edits provenance.json and pharmcat/ on its own copy of the tree
    label 'nested_apptainer'
    input:
    val s
    path ups, stageAs: 'up*/*'
    output: path 'res'
    script:
    """
    nf_stage.sh
    rc=0; PHARMCAT_SIF=/pgx/containers/pharmcat.sif pgx-pharmcat-outside.sh ${OUT} \
        > ${OUT}/pharmcat/outside.log 2>&1 || rc=\$?
    echo \$rc > ${OUT}/pharmcat/outside.status   # not gated: non-fatal by design; a kill is recorded by the script (pharmcat/outside.killed, provenance.json) and shown in the report
    """
}

process ALL_GENES_SUMMARY {
    input:
    val s
    path ups, stageAs: 'up*/*'
    val genes
    output: path 'res'
    script:
    """
    nf_stage.sh
    all_genes_summary.py ${OUT} ${s.id} ${genes.join(' ')}
    """
}

process BAMSTATS {
    cpus 4
    memory { 4.GB * task.attempt }
    input:  val s
    output: path 'res'
    script:
    """
    nf_stage.sh
    pgx-bamstats.sh ${s.bam} ${s.id} ${OUT} ${s.ref} ${s.depth_bam}
    """
}

process REPORT {   // last step: its tree is the complete result and is the one published
    publishDir "${params.outdir}/${s.id}", mode: params.publish_mode,
               saveAs: { it.replaceFirst(/^res\//, '') }
    input:
    val s
    path ups, stageAs: 'up*/*'
    output: path 'res/*'
    script:
    """
    nf_stage.sh
    pgx-report.py --sample ${s.id} --output ${OUT} \
        --genes-dir ${OUT}/genes --bam ${s.bam} --bam-stats ${OUT}/bam_stats.json \
        --provenance ${OUT}/provenance.json
    """
}

process SANITY_GATE {   // runs after REPORT is published, so a failing gate still leaves the report
    publishDir "${params.outdir}/${s.id}", mode: params.publish_mode
    input:
    val s
    path summary_tree, stageAs: 'up1/*'
    val report_done
    output: path '.sanity_ok', hidden: true
    script:
    """
    sanity_gate.py up1/res/all_genes_summary.tsv .sanity_ok
    """
}
