# ----------------------------------------------------------------------------
# Author: Seung-been "Steven" Lee
# The license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

class Locus:
    """Store various information for the given gene.

    Parameters
    ----------
    name : str
        Name of the gene.
    chr : str
        Chromosome where the gene is located.
    start : int
        Start position of the gene.
    end : int
        End position of the gene.
    upstream : int
        Upstream.
    downstream : int
        Downstream.
    masked_starts : list
        Start positions of the masked regions, if present.
    masked_ends : list
        End positions of the masked regions, if present.
    exon_starts : list
        Start positions of the exons.
    exon_ends : list
        End positions of the exons.

    Attributes
    ----------
    name : str
        Name of the gene.
    chr : str
        Chromosome where the gene is located.
    gene_start : int
        Start position of the gene.
    gene_end : int
        End position of the gene.
    upstream : int
        Upstream.
    downstream : int
        Downstream.
    region_start : int
        Start position of the bigger region.
    region_end : int
        End position of the bigger region.
    masked_starts : list
        Start positions of the masked regions, if present.
    masked_ends : list
        End positions of the masked regions, if present.
    region : str
        Genomic region in the `chr:start-end` format.
    region_size : int
        Size of the bigger region.
    gene_size : int
        Size of the gene.
    exon_starts : list
        Start positions of the exons.
    exon_ends : list
        End positions of the exons.

    """
    def __init__(self,
                 name,
                 chr,
                 start,
                 end,
                 upstream,
                 downstream,
                 masked_starts,
                 masked_ends,
                 exon_starts,
                 exon_ends):
        self.name = name
        self.chr = chr
        self.gene_start = start
        self.gene_end = end
        self.upstream = upstream
        self.downstream = downstream
        self.region_start = start - upstream
        self.region_end = end + downstream
        self.masked_starts = masked_starts
        self.masked_ends = masked_ends
        self.exon_starts = exon_starts
        self.exon_ends = exon_ends
        self.region = f'{chr}:{start-upstream}-{end+downstream}'
        self.region_size = (end + downstream) - (start - upstream)
        self.gene_size = end - start

def get_loci(gene_df,
             target_gene,
             control_gene,
             control_genes,
             genome_build):
    """Extract various data from the gene dataframe.

    The 'chr' string will be removed from the chromosome names. That is,
    'chr12' will be converted to '12'.

    Parameters
    ----------
    gene_df : pandas.DataFrame
        Dataframe containing gene information.
    target_gene : str
        Name of the target gene.
    control_gene : str
        Name of the control gene or a genomic region provided by the user.
    control_genes : list
        List of the control genes.
    genome_build : str
        Build of the reference genome assembly used to create the input data.

    Returns
    -------
    Locus
        Locus instance for the target gene.
    Locus
        Locus instance for the control gene.
    Locus
        Locus instance for the target gene's paralog if there is one.
        Otherwise, returns None.
    """
    # Get the target information.
    target_locus = Locus(*_parse_row(
        gene_df[gene_df['name'] == target_gene], genome_build))

    # Get the control information.
    if control_gene in control_genes:
        control_locus = Locus(*_parse_row(
            gene_df[gene_df['name'] == control_gene], genome_build))
    elif control_gene == None:
        control_locus = None
    else:
        name = 'CUSTOM'
        chr = control_gene.split(':')[0].replace('chr', '')
        start = int(control_gene.split(':')[1].split('-')[0])
        end = int(control_gene.split(':')[1].split('-')[1])
        upstream = 0
        downstream = 0
        masked_starts = []
        masked_ends = []
        exon_starts = []
        exon_ends = []
        control_locus = Locus(name, chr, start, end, upstream, downstream,
            masked_starts, masked_ends, exon_starts, exon_ends)

    # Get the paralog information, if present.
    paralog_gene = gene_df[
        gene_df['name'] == target_gene]['paralog'].values[0]

    if paralog_gene == '.':
        paralog_locus = None
    else:
        paralog_locus = Locus(*_parse_row(
            gene_df[gene_df['name'] == paralog_gene], genome_build))

    return target_locus, control_locus, paralog_locus

def _parse_row(row, genome_build):
    """Parse each row in the gene table."""
    name = row['name'].values[0]
    chr = row['chr'].values[0].replace('chr', '')
    start = row[f'{genome_build}_start'].values[0]
    end = row[f'{genome_build}_end'].values[0]
    upstream = row['upstream'].values[0]
    downstream = row['downstream'].values[0]
    masked_starts = row[f'{genome_build}_masked_starts'].values[0]
    masked_ends = row[f'{genome_build}_masked_ends'].values[0]
    exon_starts = row[f'{genome_build}_exon_starts'].values[0]
    exon_ends = row[f'{genome_build}_exon_ends'].values[0]

    def f(x):
        return [int(y) for y in x.strip(',').split(',')]

    _masked_starts = [] if masked_starts == '.' else f(masked_starts)
    _masked_ends = [] if masked_ends == '.' else f(masked_ends)
    _exon_starts = [] if exon_starts == '.' else f(exon_starts)
    _exon_ends = [] if exon_ends == '.' else f(exon_ends)

    return (name, chr, start, end, upstream, downstream,
        _masked_starts, _masked_ends, _exon_starts, _exon_ends)
