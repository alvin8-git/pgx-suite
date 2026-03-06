import os
from typing import List

import pysam

from common import logging, sm_tag, get_gene_table
from sglib import sort_regions

logger = logging.getLogger(__name__)

def bam2sdf(
        genome_build: str,
        target_gene: str,
        target_reg: str,
        control_gene: str,
        control_reg: str,
        bam_file: List[str],
        **kwargs
    ) -> str:
    """
    Create SDF file from BAM file(s).

    Returns:
        str: SDF file.

    Args:
        genome_build (str): Genome build (hg19, grc38).
        target_gene (str): Target gene.
        control_gene (str): Control gene or region.
        bam_file (list[str]): BAM file(s).
    """

    gene_table = get_gene_table()

    targets = [k for k, v in gene_table.items() if v["type"] == "target"]

    if target_gene not in targets:
        raise ValueError(f"'{target_gene}' is not among target genes: {targets}")

    # this could be done with all at once with _get_loci.py
    #  AND WE NEED TO CHECK IF STEVEN REALLY MEANT UPSTREAM AND DOWNSTREAM COS HE DOESN'T SEEM TO TAKE STRAND INTO ACCOUNT
    t_chr = gene_table[target_gene]["chr"].replace("chr", "")
    t_start = gene_table[target_gene][f"{genome_build}_start"]
    t_end = gene_table[target_gene][f"{genome_build}_end"]
    t_up = gene_table[target_gene]["upstream"]
    t_dn = gene_table[target_gene]["downstream"]
    t_start_minus_up = int(t_start) - int(t_up)
    t_end_plus_dn = int(t_end) + int(t_dn)
    #tr = gene_table[target_gene][f"{genome_build}_region"].replace("chr", "")
    tr = f"{t_chr}:{t_start_minus_up}-{t_end_plus_dn}"

    if "chr" in control_gene or ":" in control_gene:
        cr = control_gene.replace("chr", "")

    else:
        controls = [k for k, v in gene_table.items() if v["control"] == "yes"]

        if control_gene not in controls:
            raise ValueError(f"'{control_gene}' is not among control genes: {controls}")

        c_chr = gene_table[control_gene]["chr"].replace("chr", "")
        c_start = gene_table[control_gene][f"{genome_build}_start"].replace("chr", "")
        c_end = gene_table[control_gene][f"{genome_build}_end"].replace("chr", "")
        c_up = gene_table[control_gene]["upstream"]
        c_dn = gene_table[control_gene]["downstream"]
        c_start_minus_up = int(c_start) - int(c_up)
        c_end_plus_dn = int(c_end) + int(c_dn)
        cr = f"{c_chr}:{c_start_minus_up}-{c_end_plus_dn}"
        #cr = gene_table[control_gene][f"{genome_build}_region"].replace("chr", "")

    #regions = sort_regions([tr, cr])
    regions = sort_regions([target_reg, control_reg])

    # Get sample and sequence names from BAM headers.
    sm = []
    sn = []

    for x in bam_file: 
        sm.append(sm_tag(x))

        result = pysam.view("-H", x).strip().split("\n")
        
        for line in result:
            fields = line.split("\t")
            if "@SQ" == fields[0]:
                for field in fields:
                    if "SN:" in field:
                        y = field.replace("SN:", "")
                        if y not in sn:
                            sn.append(y)

    logger.info(f"Sample IDs: {sm}")
    # MIGHT WANNA TURN THIS BACK ON FOR RECORD-KEEPING
    #logger.info(f"Contigs: {sn}")

    # Determine whether the "chr" string should be used.
    if any(["chr" in x for x in sn]):
        chr_str = "chr"
    else:
        chr_str = ""

    result = ""

    for region in regions:
        temp = pysam.depth("-a", "-Q", "1", "-r", f"{chr_str}{region}", *bam_file)
        result += temp

    return result
