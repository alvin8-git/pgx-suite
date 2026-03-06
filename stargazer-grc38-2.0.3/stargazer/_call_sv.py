# ----------------------------------------------------------------------------
# Author: Seung-been "Steven" Lee
# The license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np
import pandas as pd
import itertools
from scipy.ndimage import median_filter

def call_sv(persons,
            gdf_file,
            target_locus,
            control_locus,
            data_type,
            sv_dict,
            logger,
            genome_build,
            sample_list=None):
    """Detect structural variants from copy number data.

    Parameters
    ----------
    persons : list
        List of Person instances.
    gdf_file : str
        Path to the GDF file.
    target_locus : Locus
        Locus instance for the target gene.
    control_locus : Locus
        Locus instance for the control gene.
    data_type : str
        Type of the input data.
    sv_dict : dict
        Dictionary of structural variants.
    logger : logging.RootLogger
        Stargazer logger.
    sample_list : list, optional
        List of sample IDs.

    Returns
    -------
    pandas.DataFrame
        Dataframe containing read depth for the target gene.
    pandas.DataFrame
        Dataframe containing read depth for the control gene.
    pandas.DataFrame
        Dataframe containing theoretical copy number profiles.
    pandas.DataFrame
        Dataframe containing copy number for unmasked portion of the
        target gene.

    """
    # Read the input GDF file.
    headers = ['Locus', 'Total_Depth', 'Average_Depth_sample']
    headers += ['Depth_for_' + x.name for x in persons]
    dtypes = ['object', 'int64', 'float64'] + ['int64'] * len(persons)
    input_df = pd.read_table(gdf_file, dtype=dict(zip(headers, dtypes)))

    # Filter out loci with low read depth.
    logger.info(f"Use all samples for normalization: {sample_list is None}")
    if sample_list is None:
        input_df = input_df[input_df.Average_Depth_sample >= 5]
    else:
        i = input_df[['Depth_for' + x for x in sample_list]].mean(axis=1) >= 5
        input_df = input_df[i]

    # Expand the Locus column.
    input_df[['chr', 'pos']] = input_df['Locus'].str.split(':', expand=True)
    input_df = input_df.astype({'pos': 'int64'})

    # Remove the unnecessary columns.
    del input_df['Locus']
    del input_df['Total_Depth']
    del input_df['Average_Depth_sample']

    # Split the dataframe into two.
    def func(row):
        result = (row['chr'].replace("chr","") == target_locus.chr and 
             target_locus.region_start <= row['pos'] <= target_locus.region_end)
        return result
    i = input_df.apply(func, axis=1)
    target_df = input_df[i]
    del target_df["chr"]
    target_df = target_df[["pos"] + [x for x in target_df if x != "pos"]] # Change order.

    # Extract control gene data.
    def func(row):
        result = (row['chr'].replace("chr","") == control_locus.chr and
            control_locus.region_start <= row['pos'] <= control_locus.region_end)
        return result
    i = input_df.apply(func, axis=1)
    control_df = input_df[i]
    del control_df["chr"]
    control_df = control_df[["pos"] + [x for x in control_df if x != "pos"]] # Change order.

    # Check for data missingness.
    if target_df.shape[0] == 0:
        raise ValueError("The input GDF file has no data for the target gene")

    if control_df.shape[0] == 0:
        raise ValueError("The input GDF file has no data for the control gene")

    # Apply the intra-sample normalization.
    norm_df = pd.concat([
        target_df.iloc[:, 0],
        target_df.iloc[:, 1:] / control_df.iloc[:, 1:].mean() * 2], axis=1)

    # Apply the inter-sample normalization.
    if data_type == 'ts':
        if sample_list:
            x = norm_df[['Depth_for_' + x for x in sl]].mean(axis=1)
        else:
            x = norm_df.iloc[:, 1:].mean(axis=1)

        norm_df.iloc[:, 1:] = norm_df.iloc[:, 1:].div(x, axis=0) * 2

    # Divide the dataframe into two again.
    if not target_locus.masked_starts:
        masked_df = pd.DataFrame()
        unmasked_df = norm_df
    else:
        def f(x):
            should_be_masked = False
            for i in range(len(target_locus.masked_starts)):
                if target_locus.masked_starts[i] <= x["pos"] <= target_locus.masked_ends[i]:
                    should_be_masked = True
                    break
            return should_be_masked

        i = norm_df.apply(lambda x: f(x), axis=1)

        masked_df = norm_df[i]
        unmasked_df = norm_df[i == False]

    median_df = pd.concat([
        unmasked_df.iloc[:, 0],
        unmasked_df.iloc[:, 1:].apply(lambda x: median_filter(x, size=1001))], axis=1)

    profile_dict = {}

    for name, dat in sv_dict.items():
        if dat[f"{genome_build}_cp"] == '.':
            changepoints = np.array([target_locus.region_start, target_locus.region_end + 1])
        else:
            changepoints = np.array([target_locus.region_start] + [int(x) for x in dat[f"{genome_build}_cp"].strip(",").split(",")] + [target_locus.region_end + 1])

        copy_numbers = [int(x) for x in dat["cn"].strip(",").split(",")]

        differences = np.diff(changepoints)

        x = np.repeat(copy_numbers, differences, axis=0)
        profile_dict[name] = x

    profile_df = pd.DataFrame.from_dict(profile_dict)
    profile_df.insert(0, "pos", range(target_locus.region_start, target_locus.region_end + 1))

    a = profile_df.pos.tolist()
    b = list(set(profile_df.pos) & set(median_df.pos)) # intersection
    i = [x in b for x in a]

    subset_df = profile_df[i]

    pairs = list(itertools.combinations([*sv_dict], 2)) + [(x, x) for x in [*sv_dict]]

    samples_with_sv = 0

    for person in persons:
        seq1_list = []
        seq2_list = []
        prob1_list = []
        prob2_list = []
        like_list = []
        lambda1_list = []
        lambda2_list = []
        ssr_raw_list = []
        ssr_adj_list = []

        for pair in pairs:
            seq1 = pair[0]
            seq2 = pair[1]
            prob1 = float(sv_dict[pair[0]]["prob"])
            prob2 = float(sv_dict[pair[1]]["prob"])

            if pair == ("no_sv", "no_sv"):
                lambda1 = 1
            else:
                lambda1 = 0.95

            if seq1 == seq2:
                lambda2 = 1
            else:
                lambda2 = 0.99

            like = prob1 * prob2
            ssr_raw = np.square(subset_df[seq1].to_numpy() + subset_df[seq2].to_numpy() - median_df["Depth_for_" + person.name].to_numpy()).sum()

            if ssr_raw == 0:
                ssr_raw = 1

            ssr_adj = ssr_raw / like / lambda1 / lambda2

            seq1_list.append(pair[0])
            seq2_list.append(pair[1])
            prob1_list.append(prob1)
            prob2_list.append(prob2)
            like_list.append(like)
            lambda1_list.append(lambda1)
            lambda2_list.append(lambda2)
            ssr_raw_list.append(ssr_raw)
            ssr_adj_list.append(ssr_adj)

        df = pd.DataFrame({
            "seq1": seq1_list,
            "seq2": seq2_list,
            "prob1": prob1_list,
            "prob2": prob2_list,
            "like": like_list,
            "lambda1": lambda1_list,
            "lambda2": lambda2_list,
            "ssr_raw": ssr_raw_list,
            "ssr_adj": ssr_adj_list,
        })

        df.sort_values(by=["ssr_adj"], inplace=True)
        person.sv = [df["seq1"].iloc[0], df["seq2"].iloc[0]]
        person.ssr = df["ssr_adj"].iloc[0]
        person.ssr_df = df

        if person.sv != ["no_sv", "no_sv"]:
            samples_with_sv += 1

    logger.info(f"Number of samples with SV: {samples_with_sv}")

    return target_df, control_df, profile_df, unmasked_df
