#!/usr/bin/env nextflow
// pgx-suite on Nextflow DSL2 (port of docker/Snakefile; see docs/nextflow-port-plan.md).
// Pilot scope: every rule except hla2, mutserve and abo_pileup (genes needing those are refused).
//
//   nextflow run main.nf -profile docker --bam /path/HG002.bam --sample HG002 \
//       --ref /path/hg38.fa --genes "CYP2D6 G6PD HLA-A" --outdir results

include { PREFLIGHT; BCFTOOLS; ALDY; STELLARPGX; PYPGX_DEPTH; PYPGX_CONTROL; STARGAZER_GDF;
          PYPGX; STARGAZER; CYRIUS; OPTITYPE; PHARMCAT; ALTCHECK; COMPARE; PROVENANCE;
          PHARMCAT_OUTSIDE; ALL_GENES_SUMMARY; BAMSTATS; REPORT; SANITY_GATE } from './modules/local/pgx'

// docker/genes.tsv is the single source of truth, as for the Snakefile.
def readGenes(String path) {
    def lines = file(path).readLines().findAll { it.trim() }
    def hdr = lines[0].split('\t')
    lines.drop(1).collect { ln ->
        def v = ln.split('\t', -1)
        def r = [:]; hdr.eachWithIndex { h, i -> r[h] = v[i] }
        def on = { String c -> r[c] == '1' }
        def tools = ['pypgx', 'stargazer', 'aldy', 'stellarpgx', 'optitype', 'mutserve',
                     'vcf_check', 'cyrius', 'pharmcat', 'arcashla'].findAll { on(it) }
        def is_alt = r.region == 'ALT_CONTIG'
        [gene: r.gene, region: r.region, tools: tools, pypgx_sv: on('pypgx_sv'),
         stargazer_sv: on('stargazer_sv'), stargazer_control: r.stargazer_control,
         is_alt: is_alt, needs_vcf: (on('pypgx') || on('stargazer')) && !is_alt]
    }
}

// Number of result trees COMPARE waits for, from the gene's flags (never waits on a tool that will not run).
def nCompareInputs(Map g) {
    1 /* altcheck */ +
    ['pypgx', 'stargazer', 'aldy', 'stellarpgx', 'cyrius', 'optitype', 'pharmcat'].count { g.tools.contains(it) } +
    ((g.needs_vcf || g.tools.contains('vcf_check')) ? 1 : 0)
}

// Channel helpers (top-level functions: closures assigned in a workflow block are taken as process names)
def withTool(ch, String t) { ch.filter { g, o -> g.tools.contains(t) } }
def keyed(ch) { ch.map { g, r -> tuple(g.gene, g, r) } }
def withTrees(start, List chans) {   // left-join result trees on the gene name; keep only existing ones
    chans.inject(start.map { g, s -> tuple(g.gene, g, s) }) { acc, c -> acc.join(keyed(c), remainder: true) }
         .filter { it[1] != null }
         .map { row -> tuple(row[1], row[2], row.drop(3).findAll { it instanceof Path }) }
}
def sampleTrees(genes, tree, String t) {
    genes.filter { g, o -> t == null || g.tools.contains(t) }.combine(tree).map { g, o, r -> tuple(g, r) }
}
def inOrder(List chans) {   // concatenate channels' paths in list order (later trees win in nf_stage.sh)
    chans.withIndex().collect { ch, i -> ch.flatten().map { r -> tuple(i, r) } }
         .inject { a, b -> a.mix(b) }
         .toSortedList { a, b -> a[0] <=> b[0] }.map { it*.get(1) }
}

workflow {
    if (!params.bam)    error "--bam is required"
    if (!params.ref)    error "--ref is required"
    if (!params.sample) error "--sample is required"

    def all = readGenes(params.genes_tsv).findAll { it.gene != 'CYPOR' }   // CYPOR is an alias of POR
    def want = params.genes ? params.genes.toString().tokenize(' ,') : all*.gene
    def unknown = want - all*.gene
    if (unknown) error "unknown genes: ${unknown.join(' ')}"
    def target = all.findAll { want.contains(it.gene) }.sort { it.gene }
    def notYet = target.findAll { g -> ['mutserve', 'arcashla'].any { g.tools.contains(it) } || g.gene == 'ABO' }
    if (notYet) error "not yet ported (pilot): ${notYet*.gene.join(' ')}"

    // The sample map is a declared input of every process, so it is part of every task hash.
    def bamName = file(params.bam).name
    def ext = bamName.tokenize('.').last()
    def sample = [id: params.sample.toString(), bam: "/pgx/data/${bamName}", ref: "/pgx/ref/${file(params.ref).name}",
                  bam_glob: "/pgx/data/${bamName.substring(0, bamName.lastIndexOf('.'))}*{${ext},${ext == 'bam' ? 'bai' : 'crai'}}",
                  depth_bam: params.depth_bam ? "/pgx/data/${file(params.depth_bam).name}" : "/pgx/data/${bamName}",
                  sex: params.sex ?: '', min_depth: params.min_depth]
    ok = PREFLIGHT(sample, target)
    genes = ok.combine(channel.fromList(target)).map { s, g -> tuple(g, s) }

    // per-gene callers
    bcf    = BCFTOOLS(genes.filter { g, o -> g.needs_vcf || g.tools.contains('vcf_check') })
    aldy   = ALDY(withTool(genes, 'aldy'))
    stellar = STELLARPGX(withTool(genes, 'stellarpgx'))
    cyrius = CYRIUS(withTool(genes, 'cyrius'))

    // PyPGx / Stargazer: join the VCF and the hard SV-prep steps on the gene name
    pdep  = PYPGX_DEPTH(withTool(genes, 'pypgx').filter { g, o -> g.pypgx_sv })
    pctl  = PYPGX_CONTROL(withTool(genes, 'pypgx').filter { g, o -> g.pypgx_sv })
    gdf   = STARGAZER_GDF(withTool(genes, 'stargazer').filter { g, o -> g.stargazer_sv })
    pyIn  = withTrees(withTool(genes, 'pypgx'), [bcf, pdep, pctl])
    sgIn  = withTrees(withTool(genes, 'stargazer'), [bcf, gdf])
    pypgx = PYPGX(pyIn)
    starg = STARGAZER(sgIn)

    // sample-level steps
    alt  = ALTCHECK(ok)
    opti = OPTITYPE(ok)
    pc   = PHARMCAT(ok)

    // COMPARE: all trees routed to a gene, grouped with an exact expected count.
    // PharmCAT goes only to pharmcat-flag genes, the dependency the Snakefile declares. pgx-compare.py
    // also adds a PharmCAT cross-check to any other gene if pharmcat/pgx.report.json already exists,
    // which under Snakemake is a race (on a full genome PharmCAT usually finishes after the compares,
    // so the cross-check is usually absent). Here it is always absent for those genes: deterministic.
    perGene = bcf.mix(aldy, stellar, cyrius, pypgx, starg)
                 .mix(sampleTrees(genes, alt, null), sampleTrees(genes, opti, 'optitype'), sampleTrees(genes, pc, 'pharmcat'))
                 .map { g, r -> tuple(groupKey(g.gene, nCompareInputs(g)), g, r) }
                 .groupTuple()
                 .map { k, gs, rs -> tuple(gs[0], rs) }
                 .combine(ok).map { g, rs, s -> tuple(g, s, rs) }   // ponytail: one sample per run; key on sample too for multi-sample
    cmp = COMPARE(perGene)
    cmpTrees = cmp.map { g, r -> r }.collect()

    // sample-level finalisation (plan 2.5: PHARMCAT_OUTSIDE edits its own copy; only its result goes on)
    // Trees are merged in list order and later ones win, so the order must be fixed:
    // PHARMCAT_OUTSIDE's edited provenance.json and pharmcat/ must come after any stage-1 copy.
    prov    = PROVENANCE(ok, target.size())
    outside = PHARMCAT_OUTSIDE(ok, inOrder([cmpTrees, pc, prov]))
    summary = ALL_GENES_SUMMARY(ok, cmpTrees, target*.gene)
    bstats  = BAMSTATS(ok)
    report  = REPORT(ok, inOrder([cmpTrees, outside, bstats, summary]))
    SANITY_GATE(ok, summary, report.collect().map { true })
}
