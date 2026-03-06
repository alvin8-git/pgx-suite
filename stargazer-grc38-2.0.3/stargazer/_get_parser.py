# ----------------------------------------------------------------------------
# Author: Seung-been "Steven" Lee
# The license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import argparse
from version import __version__

def get_parser(target_genes,
               control_genes,
               data_types,
               program_dir):
    """Create a parser with user-provided arguments.

    Parameters
    ----------
    target_genes : list
        Target gene names.
    control_genes : list
        Control gene names.
    data_types : list
        Input data types.

    Returns
    -------
    argparse.ArgumentParser
        A parser.

    """
    builds = ['hg19', 'grc38']

    parser = argparse.ArgumentParser(
        description=("Stargazer is a bioinformatics tool for calling star "
                     "alleles in pharmacogenes using genomic data from "
                     "next-generation sequencing or single nucleotide "
                     "polymorphism microarray. "
                     "Required arguments are indicated with '[required]'. "
                     "In order to perform structural variation detection, "
                     "you must provide both '-c/--control-gene' and "
                     "'-g/--gdf-file'. Otherwise, Stargazer will run in "
                     "the VCF-only mode. "
                     "Stargazer switches to gdf-creation mode with"
                     "'-G/--create-gdf-file' argument. You must provide "
                     "both -c/--control-gene' and -B/--bam-file"
                     "For more information, visit the official website: "
                     "https://stargazer.gs.washington.edu/stargazerweb.")
    )

    parser.add_argument(
        '-v',
        '--version',
        action='version',
        version=f'Stargazer v{__version__}',
        help="Print the Stargazer version number and exit."
    )

    parser.add_argument(
        '-i',
        '--vcf-file',
        #required=True,  # no longer required with addition of -create-gdf mode
        metavar='PATH',
        help=("Path to the input VCF file (.vcf or .vcf.gz). An empty VCF "
              "file containing only the header and metadata is still "
              "accepted by Stargazer as it means no variants were detected "
              "from the sample(s). [required for standard running]")
    )

    parser.add_argument(
        '-o',
        '--output-dir',
        required=True,
        metavar='PATH',
        help=("Path to the output directory. If the directory with the "
              "same name already exists, Stargazer will overwrite it. "
              "[required for all processes]")
    )

    parser.add_argument(
        '-t',
        '--target-gene',
        required=True,
        metavar=f'TEXT',
        choices=target_genes,
        help=(f"Name of the target gene. Choices: {_func(target_genes)}."
              " [required for all processes]")
    )

    parser.add_argument(
        '-a',
        '--genome-build',
        metavar='TEXT',
        choices=builds,
        default='hg19',
        help=("Build of the reference genome assembly used to create the "
              f"input data. Choices: {_func(builds)}. [default: 'hg19']")
    )

    parser.add_argument(
        '-d',
        '--data-type',
        metavar='TEXT',
        choices=data_types,
        default='wgs',
        help=("Type of the input data: wgs, whole genome sequencing; "
              "ts, targeted sequencing including whole exome sequencing; "
              "or chip, single nucleotide polymorphism microarray. Choices: "
              f"{_func(data_types)}. [default: 'wgs']")
    )

    parser.add_argument(
        '-c',
        '--control-gene',
        metavar='TEXT',
        help=("Name of a preselected control gene. Used for intrasample "
              "normalization during copy number analysis. Choices: "
              f"{_func(control_genes)}. Alternatively, you can provide a "
              "custom genomic region with the 'chr:start-end' format (e.g. "
              "chr12:48232319-48301814). [required with --create-gdf-file]")
    )

    parser.add_argument(
        '-g',
        '--gdf-file',
        metavar='PATH',
        help=("Path to a GDF file containing read depth for both the target "
              "and control genes.")
    )

    parser.add_argument(
        '-r',
        '--hap-panel',
        metavar='PATH',
        help=("Path to a reference haplotype panel file (.vcf or .vcf.gz). "
              "By default, Stargazer uses the haplotype panel from the 1000 "
              "Genomes Project (phase 3).")
    )

    parser.add_argument(
        '-s',
        '--sample-list',
        metavar='TEXT',
        nargs='+',
        help=("List of space-delimited sample IDs to be used as control "
              "for intersample normalization during copy number analysis. "
              "This argument has an effect only when the input data type "
              "is targeted sequencing.")
    )

    parser.add_argument(
        '-p',
        '--include-profiles',
        action='store_true',
        help=("Use this tag to include profiles for read depth, copy number, "
              "and allele fraction. This tag has no effect if Stargazer is "
              "running in the VCF-only mode.")
    )

    parser.add_argument(
        '-b',
        '--impute-ungenotyped-markers',
        action='store_true',
        help=("Use this tag to impute ungenotyped markers during "
              "statistical phasing. Potentially useful for low-density DNA "
              "microarray data.")
    )

    parser.add_argument(
        '-st',
        '--star_table',
        metavar='PATH',
        default=f"{program_dir}/star_table.tsv",
        #action='store_true',
        help=("Use this tag to specify the path for star_table to be used")
    )

    parser.add_argument(
        '-G',
        '--create-gdf-file',
        metavar='NAME',
        help=("This tag causes Stargazer to run as a GDF-file generator. Argument requires a "
              "file name (which will be placed in the output directory).")
    )
    
    parser.add_argument(
        '-B',
        '--bam-file',
        nargs="*", # because this needs to be a list
        metavar='LIST',
        help=("List of paths to BAM or CRAM files. [required for --create-gdf-file runs]")
    )
           
    args = parser.parse_args()

    if not args.create_gdf_file:
        if args.vcf_file is None:
            parser.error('without --create-gdf-file, -v/--vcf-file is required')
    else:
        if args.control_gene is None:
            parser.error('with --create-gdf-file, -c/--control-gene is required')

    return parser

def _func(l):
    """Convert a list to a pretty text."""
    return '{' + ', '.join([f"'{x}'" for x in l]) + '}'
