# ----------------------------------------------------------------------------
# Author: Seung-been "Steven" Lee
# The license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

Y_MAX = 6
Y_MIN = -1.5
STEPS = 500
XLABEL_FONTSIZE = 15
YLABEL_FONTSIZE = 20
XTICKLABELS_FONTSIZE = 15
YTICKLABELS_FONTSIZE = 15
NAME_FONTSIZE = 15
NAME_POSITION = -1.3

def plot(target_df,
         control_df,
         profile_df,
         unmasked_df,
         target_locus,
         control_locus,
         output_dir,
         persons,
         logger,
         paralog_locus=None):
    """Create read depth and copy number plots.

    This method creates a PNG file for each sample, which contains:

    1. Read depth profile for the target gene:
        Gray dots indicate the sample's per-base read depth. Green dots
        indicate average per-base read depth from the entire study samples.

    2. Read depth profile for the control gene:
        Same as above.

    3. Copy number profile for the target gene:
        Gray dots indicate the sample's per-base copy number estimates
        computed from read depth. The navy solid line and the cyan dashed
        line represent copy number profiles for each haplotype. The red line
        represents copy number profile for both haplotypes combined.

    4. Allele fraction profile for the target gene:
        Navy dots and cyan dots indicate allele fraction estimates computed
        from allelic read depth for each haplotype. The navy solid line and
        the cyan dashed line represent copy number profiles for each
        haplotype. The red line represents copy number profile for both
        haplotypes combined.

    Parameters
    ----------
    target_df : pandas.DataFrame
        Dataframe containing read depth for the target gene.
    control_df : pandas.DataFrame
        Dataframe containing read depth for the control gene.
    profile_df : pandas.DataFrame
        Dataframe containing theoretical copy number profiles.
    unmasked_df : pandas.DataFrame
        Dataframe containing copy number for unmasked portion of the
        target gene.
    target_locus : Locus
        Locus instance for the target gene.
    control_locus : Locus
        Locus instance for the control gene.
    output_dir : str
        Path to the output directory.
    persons : list
        List of Person instances.
    logger : logging.RootLogger
        Stargazer logger.
    paralog_locus : Locus, optional
        Locus instance for the target gene's paralog if there is one.

    """
    # Make the plots directory.
    os.mkdir(f"{output_dir}/plots")




    for person in persons:
        filtered = [x for x in person.hap[0].obs if x.td > 10 and x.het]
        if not filtered:
            continue

        pos_list = []
        hap1_al_list = []
        hap2_al_list = []
        hap1_af_list = []
        hap2_af_list = []

        for snp in filtered:
            hap1_snp = [x for x in person.hap[0].obs if x.pos == snp.pos][0]
            hap2_snp = [x for x in person.hap[1].obs if x.pos == snp.pos][0]

            pos_list.append(int(snp.pos))
            hap1_al_list.append(hap1_snp.var)
            hap2_al_list.append(hap2_snp.var)
            hap1_af_list.append(hap1_snp.af)
            hap2_af_list.append(hap2_snp.af)

        df = pd.DataFrame({
            "pos": pos_list,
            "hap1_al": hap1_al_list,
            "hap2_al": hap2_al_list,
            "hap1_af": hap1_af_list,
            "hap2_af": hap2_af_list
        })

        person.af_df = df









    # Get the sequencing coverage.
    tcv = target_df.shape[0] / (target_locus.region_size)
    ccv = control_df.shape[0] / (control_locus.region_size)

    # Set the increment size for x-axis. Increasing the number of steps will
    # increase the density of data points plotted. Increasing the sequencing
    # coverage has the opposite effect.
    tic = round((target_locus.region_size) / STEPS * tcv)
    cic = round((control_locus.region_size) / STEPS * ccv)

    logger.info(f"Step size: {STEPS}")
    logger.info(f"Step size adjusted for target gene: {tic}")
    logger.info(f"Step size adjusted for control gene: {cic}")

    # Create the plots for each sample.
    for person in persons:
        name = person.name

        i = list(target_df).index("Depth_for_" + name)
        hap1_sv = person.get_hap1_sv()
        hap2_sv = person.get_hap2_sv()

        if hap1_sv == '.':
            dip_sv = person.get_dip_sv().split(',')
            hap1_sv = dip_sv[0]
            hap2_sv = dip_sv[1]

        # Determine the y-axis scale for read depth.
        if target_df["Depth_for_" + name].mean() < 60:
            scale1 = 15
        else:
            scale1 = 30

        # Determine the y-axis scale for copy number.
        if (profile_df[hap1_sv] + profile_df[hap2_sv]).max() < 5:
            scale2 = 1
        else:
            scale2 = 2

        fig, [[ax1, ax2], [ax3, ax4]] = plt.subplots(2, 2, figsize=(14, 10))

        # Plot the read depth profile for the target region.
        _add_base_plot(ax1, target_locus.chr, target_locus.region_start, target_locus.region_end, scale1)
        ax1.set_ylabel("Read depth", fontsize=YLABEL_FONTSIZE)
        ax1.set_yticks(np.linspace(0, Y_MAX * scale1, 4))
        _add_gene_model(ax1, target_locus.exon_starts, target_locus.exon_ends, scale1)
        _add_gene_name(ax1, target_locus.gene_start, target_locus.gene_end, scale1, target_locus.name)
        if paralog_locus is not None:
            _add_gene_model(ax1, paralog_locus.exon_starts, paralog_locus.exon_ends, scale1)
            _add_gene_name(ax1, paralog_locus.gene_start, paralog_locus.gene_end, scale1, paralog_locus.name)
        ax1.scatter(
            target_df.iloc[::tic, 0],
            target_df.iloc[::tic, i],
            color="lightgray",
            edgecolors="black"
        )
        ax1.scatter(
            target_df.iloc[::tic, 0],
            target_df.drop(["pos"], axis=1).mean(axis=1).iloc[::tic],
            color="lightgreen",
            edgecolors="black"
        )

        # Plot the read depth profile for the control region.
        _add_base_plot(ax2, control_locus.chr, control_locus.region_start, control_locus.region_end, scale1)
        ax2.set_ylabel("Read depth", fontsize=YLABEL_FONTSIZE)
        ax2.set_yticks(np.linspace(0, Y_MAX * scale1, 4))
        if control_locus.exon_starts:
            _add_gene_model(ax2, control_locus.exon_starts, control_locus.exon_ends, scale1)
        _add_gene_name(ax2, control_locus.region_start, control_locus.region_end, scale1, control_locus.name)
        ax2.scatter(
            control_df.iloc[::cic, 0],
            control_df.iloc[::cic, i],
            color="lightgray",
            edgecolors="black"
        )
        ax2.scatter(
            control_df.iloc[::tic, 0],
            control_df.drop(["pos"], axis=1).mean(axis=1).iloc[::tic],
            color="lightgreen",
            edgecolors="black"
        )



        # Plot the copy number profile for the target region.
        _add_base_plot(ax3, target_locus.chr, target_locus.region_start, target_locus.region_end, scale2)
        ax3.set_ylabel("Copy number", fontsize=YLABEL_FONTSIZE)
        ax3.set_yticks(np.linspace(0, Y_MAX * scale2, 4))

        _add_gene_model(ax3, target_locus.exon_starts, target_locus.exon_ends, scale2)
        _add_gene_name(ax3, target_locus.gene_start, target_locus.gene_end, scale2, target_locus.name)

        if paralog_locus is not None:
            _add_gene_model(ax3, paralog_locus.exon_starts, paralog_locus.exon_ends, scale2)
            _add_gene_name(ax3, paralog_locus.gene_start, paralog_locus.gene_end, scale2, paralog_locus.name)

        ax3.scatter(
            unmasked_df.iloc[::tic, 0],
            unmasked_df.iloc[::tic, i],
            color="lightgray",
            edgecolors="black"
        )

        _add_svcomb_plot(ax3, profile_df, hap1_sv, hap2_sv)


        # Plot the allele fraction profile for the target region.
        _add_base_plot(ax4, target_locus.chr, target_locus.region_start, target_locus.region_end, scale2)
        ax4.set_ylabel("Allele fraction", fontsize=YLABEL_FONTSIZE)
        ax4.set_yticks([0, Y_MAX * scale2 / 2, Y_MAX * scale2])
        ax4.set_yticklabels([0, 0.5, 1])


        _add_gene_model(ax4, target_locus.exon_starts, target_locus.exon_ends, scale2)
        _add_gene_name(ax4, target_locus.gene_start, target_locus.gene_end, scale2, target_locus.name)

        if paralog_locus is not None:
            _add_gene_model(ax4, paralog_locus.exon_starts, paralog_locus.exon_ends, scale2)
            _add_gene_name(ax4, paralog_locus.gene_start, paralog_locus.gene_end, scale2, paralog_locus.name)

        if person.af_df is not None:
            ax4.scatter(
                person.af_df["pos"],
                person.af_df["hap1_af"] * Y_MAX * scale2,
                color="navy",
                edgecolors="black"
            )

            ax4.scatter(
                person.af_df["pos"],
                person.af_df["hap2_af"] * Y_MAX * scale2,
                color="cyan",
                edgecolors="black"
            )

        _add_svcomb_plot(ax4, profile_df, hap1_sv, hap2_sv)


        plt.tight_layout()
        plt.savefig(f"{output_dir}/plots/{name}.png")
        plt.close(fig)

def _add_base_plot(ax, chr, start, end, scale):
    ax.set_xlim([start, end])
    ax.set_ylim([Y_MIN * scale, Y_MAX * scale])
    ax.set_xlabel("Chromosome {}".format(chr.replace("chr", "")), fontsize=XLABEL_FONTSIZE)
    ax.set_xticks(np.linspace(start, end, 5))
    ax.tick_params(axis='x', which='major', labelsize=XTICKLABELS_FONTSIZE)
    ax.tick_params(axis='y', which='major', labelsize=YTICKLABELS_FONTSIZE)

def _add_gene_name(ax, start, end, scale, name):
    ax.text(
        x=(start + end) / 2,
        y=NAME_POSITION * scale,
        s=f"${name.upper()}$",
        horizontalalignment="center",
        fontsize=NAME_FONTSIZE
    )

def _add_gene_model(ax, starts, ends, scale):
    ax.hlines(
        y=Y_MIN * 0.5 * scale,
        xmin=starts[0],
        xmax=ends[-1],
        color="black"
    )

    for i in range(len(starts)):
        ax.add_patch(Rectangle(
            xy=(starts[i], Y_MIN * 0.5 * scale - 0.15 * scale),
            width=ends[i] - starts[i],
            height=0.30 * scale,
            color="black"
        ))

def _add_svcomb_plot(ax, p, s1, s2):
    ax.plot(p["pos"], p[s1] + p[s2], color="red", linewidth=6)
    ax.plot(p["pos"], p[s1], color="navy", linewidth=2)
    ax.plot(p["pos"], p[s2], color="cyan", linewidth=2, linestyle="dashed")
