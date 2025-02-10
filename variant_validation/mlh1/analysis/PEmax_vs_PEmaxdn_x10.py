"""
    Project :  Prime editing pilot screens. Revision analysis. Script to investigate difference between variant scores generated in HAP1:PEmax vs HAP1:PEmax+MLH1dn cell lines (x10 experiment).
    Date : 240916
    Python version 3.10

"""

# IMPORT

import pandas as pd
import pathlib as pl
import math
import numpy as np

# VARIABLES

INPUT_DIR = pl.Path("/Users/kajbac/code/230904_S002_final_filter")

OUTPUT_DIR = pl.Path("/Users/kajbac/code/240913_cell_line_differences")
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)

pegRNA_D20_D34_path = INPUT_DIR / "MLH1x10_pegRNA_data_merged_best_pegRNA.csv"

variables_per_sample = [
    {"pegRNA_D20": "CK003", "pegRNA_D34": "CK007"},
    {"pegRNA_D20": "CK005", "pegRNA_D34": "CK009"},
    {"pegRNA_D20": "CK004", "pegRNA_D34": "CK008"},
    {"pegRNA_D20": "CK006", "pegRNA_D34": "CK010"},
]

PEmax_variables_per_sample = [
    {"pegRNA_D20": "CK003", "pegRNA_D34": "CK007"},
    {"pegRNA_D20": "CK005", "pegRNA_D34": "CK009"},
]
PEmaxdn_variables_per_sample = [
    {"pegRNA_D20": "CK004", "pegRNA_D34": "CK008"},
    {"pegRNA_D20": "CK006", "pegRNA_D34": "CK010"},
]

D20_list = ["CK003", "CK005", "CK004", "CK006"]
PEmax_D20_list = ["CK003", "CK005"]
PEmaxdn_D20_list = ["CK004", "CK006"]

PEmax_list = ["CK003", "CK005", "CK007", "CK009"]
PEmaxdn_list = ["CK004", "CK006", "CK008", "CK010"]
pseudo_pre_freq_filter = 1.4e-4
ST_editing_filter_PEmax = 5
ST_editing_filter_PEmaxdn = 25

# Dictionaries to correct for redundancy in pegRNA design between non-sense SNVs and stopGain MNVS: will be used to rename mutant index,variant_keys, GRCh38 start and end locations of stopGain MNV pegRNAs that are similar to non-sense SNV pegRNAs
mutant_index_correction_dict = {
    2184: 1775,
    2185: 1785,
    2190: 1829,
    2193: 1856,
    2196: 1887,
    2199: 1914,
    2202: 1937,
    2209: 2004,
}
GRCh38_correction_dict = {
    "37017517_TAA": 37017517,
    "37017520_TAA": 37017521,
    "37017535_TAA": 37017535,
    "37017544_TAA": 37017544,
    "37017553_TAA": 37017555,
    "37017562_TAA": 37017564,
    "37017571_TAA": 37017571,
    "37017592_TAA": 37017594,
}
var_key_correction_dict = {
    "37017517_TAA": "37017517_T",
    "37017520_TAA": "37017521_A",
    "37017535_TAA": "37017535_T",
    "37017544_TAA": "37017544_T",
    "37017553_TAA": "37017555_A",
    "37017562_TAA": "37017564_A",
    "37017571_TAA": "37017571_T",
    "37017592_TAA": "37017594_A",
}

# high ST threshold: mean ST-editing rate for pegRNAs encoding the same variant (pegRNAs retained after first ST-filter step)
high_ST_threshold = 36

# ----------------------------------------------------------------------------------------------------------------------
#                                                     Functions
# ----------------------------------------------------------------------------------------------------------------------
# Filter out variants with high-editing rates in negative control


def filter_edited_controls_first(df):
    """
    Filter out variants with high correct ST editing in negative control.

    Parameters
    ----------
    df (df): Dataframe to be processed.
    """
    df = df.loc[df["percentage_editing_control_high"] == False]
    return df


def translate_GRCh37locations_to_GRCh38locations(df):
    """
    Lift GRCh37 locations to GRCh38 locations.

    Parameters
    ----------
    df (df): Dataframe to be processed.
    """
    df["start"] = df["start"] - 41491
    df.rename(columns={"start": "GRCh38startlocation"}, inplace=True)
    df["end"] = df["end"] - 41491
    df.rename(columns={"end": "GRCh38endlocation"}, inplace=True)
    return df


def add_var_key_for_dataframe_merging_pegRNA(df):
    """
    Generate a string identifier for each variant used for dataframe merging.

    Parameters
    ----------
    df (df): Dataframe to be processed.
    """
    df["var_key"] = (
        df["GRCh38startlocation"].astype("string")
        + "_"
        + df["mut_allele"].astype("string")
    )
    return df


def correct_stopGain_pegRNA_design_redundancy(df):
    """
    Correct "mutant indices", "GRCh38startlocation", "var_key"s of redundant MNV stopGain pegRNAs.

    Parameters
    ----------
    df (df): Dataframe to be processed.
    """
    df["mutant index"] = df["mutant index"].replace(mutant_index_correction_dict)
    df["GRCh38startlocation"] = (
        df["var_key"].map(GRCh38_correction_dict).fillna(df["GRCh38startlocation"])
    )
    df["GRCh38endlocation"] = (
        df["var_key"].map(GRCh38_correction_dict).fillna(df["GRCh38endlocation"])
    )
    df["var_key"] = df["var_key"].replace(var_key_correction_dict)

    return df


def add_pseudo_count(
    df,
    pegRNA_D20,
    pegRNA_D34,
):
    """
    Add pseudocount of one to pegRNA count.

    Parameters
    ----------
    df (df): Dataframe to be processed.
    pegRNA_D20 (str) : D20 sample identifier
    pegRNA_D34 (str) : D34 sample identifier
    """
    df[f"{pegRNA_D20}_pegRNA_pseudo_count"] = df[f"{pegRNA_D20}_pegRNA_count"] + 1
    df[f"{pegRNA_D34}_pegRNA_pseudo_count"] = df[f"{pegRNA_D34}_pegRNA_count"] + 1
    return df


def calculate_frequencies_ratio_log2_pegRNA(df, pegRNA_D20, pegRNA_D34):
    """
    Calculate pegRNA frequencies, post- over pre-frequency ratios, and log2-ratios from pegRNA pseudocounts.

    Parameters
    ----------
    df (df): Dataframe to be processed.
    pegRNA_D20 (str) : D20 sample identifier
    pegRNA_D34 (str) : D34 sample identifier
    """
    df[f"{pegRNA_D20}_pre_pseudo_freq"] = df[f"{pegRNA_D20}_pegRNA_pseudo_count"] / df[
        f"{pegRNA_D20}_pegRNA_pseudo_count"
    ].sum(axis=0)
    df[f"{pegRNA_D34}_post_pseudo_freq"] = df[f"{pegRNA_D34}_pegRNA_pseudo_count"] / df[
        f"{pegRNA_D34}_pegRNA_pseudo_count"
    ].sum(axis=0)
    df[f"{pegRNA_D34}_{pegRNA_D20}_ratio"] = df[
        f"{pegRNA_D34}_post_pseudo_freq"
    ].astype("float64") / df[f"{pegRNA_D20}_pre_pseudo_freq"].astype("float64")

    df[f"{pegRNA_D20}_{pegRNA_D34}_log2"] = df[
        f"{pegRNA_D34}_{pegRNA_D20}_ratio"
    ].apply(math.log2)
    return df


def process_pegRNA_data_frame(
    df,
):
    """
    Lift GRCh37 locations to GRCh38 locations, filter out variants with high correct ST editing in negative control, generate a string identifier for each variant used for dataframe merging, and correct "mutant index" and "var_key" string .

    Parameters
    ----------
    df (df): Dataframe to be processed.
    """
    df = (
        df.pipe(translate_GRCh37locations_to_GRCh38locations)
        .pipe(filter_edited_controls_first)
        .pipe(add_var_key_for_dataframe_merging_pegRNA)
        .pipe(correct_stopGain_pegRNA_design_redundancy)
    )

    return df


def pegRNA_frequency_and_log2(df, pegRNA_D20, pegRNA_D34):
    """
    Add pseudocount of one to pegRNA count, and calculate pegRNA frequencies, post- over pre-frequency ratios, and log2-ratios from pegRNA pseudocounts.

    Parameters
    ----------
    df (df): Dataframe to be processed.
    """

    df = df.pipe(add_pseudo_count, pegRNA_D20, pegRNA_D34).pipe(
        calculate_frequencies_ratio_log2_pegRNA, pegRNA_D20, pegRNA_D34
    )

    return df


# ----------------------------------------------------------------------------------------------------------------------
# Apply pre-frequency and percentage ST-editing filter


def filter_pegRNAs_on_pre_freq(df, var_dict):
    """
    Set pegRNA log2-ratios and pegRNA percentage editing to nan, if pre-frequency is below specified filter value.

    Parameters
    ----------
    df (df): Dataframe to be processed.
    var_dict (dict): Dictionary used to link sample collection date to sample identifier.
    """
    df.loc[
        df[f"{var_dict['pegRNA_D20']}_pre_pseudo_freq"] < pseudo_pre_freq_filter,
        [
            f"{var_dict['pegRNA_D20']}_{var_dict['pegRNA_D34']}_log2",
            f"{var_dict['pegRNA_D20']}_percentage_editing",
        ],
    ] = np.nan

    return df


def filter_pegRNAs_on_ST_editing_rates(df, var_dict):
    """
    Set pegRNA log2-ratios and pegRNA percentage editing to nan, if surrogate target editing rate is below specified filter value.

    Parameters
    ----------
    df (df): Dataframe to be processed.
    var_dict (dict): Dictionary used to translate sample collection date to sample identifier.
    """
    if var_dict["pegRNA_D20"] in PEmax_list:
        df.loc[
            df[f"{var_dict['pegRNA_D20']}_percentage_editing"]
            < ST_editing_filter_PEmax,
            [
                f"{var_dict['pegRNA_D20']}_{var_dict['pegRNA_D34']}_log2",
                f"{var_dict['pegRNA_D20']}_percentage_editing",
            ],
        ] = np.nan
        df.loc[
            df[f"{var_dict['pegRNA_D20']}_percentage_editing"] == np.nan,
            [
                f"{var_dict['pegRNA_D20']}_{var_dict['pegRNA_D34']}_log2",
                f"{var_dict['pegRNA_D20']}_percentage_editing",
            ],
        ] = np.nan
    elif var_dict["pegRNA_D20"] in PEmaxdn_list:
        df.loc[
            df[f"{var_dict['pegRNA_D20']}_percentage_editing"]
            < ST_editing_filter_PEmaxdn,
            [
                f"{var_dict['pegRNA_D20']}_{var_dict['pegRNA_D34']}_log2",
                f"{var_dict['pegRNA_D20']}_percentage_editing",
            ],
        ] = np.nan
        df.loc[
            df[f"{var_dict['pegRNA_D20']}_percentage_editing"] == np.nan,
            [
                f"{var_dict['pegRNA_D20']}_{var_dict['pegRNA_D34']}_log2",
                f"{var_dict['pegRNA_D20']}_percentage_editing",
            ],
        ] = np.nan
    return df


# ----------------------------------------------------------------------------------------------------------------------
# Calculate variant-level function scores


def collapse_variants(df):
    """
    Average pegRNA log2-ratios of pegRNAs programming the same variant.

    Parameters
    ----------
    df (df): Dataframe to be processed.
    """
    log2_list = []
    for var_dict in variables_per_sample:
        log2_list.append(f"{var_dict['pegRNA_D20']}_{var_dict['pegRNA_D34']}_log2")
    df_log2_mean = df.groupby(["var_key"])[log2_list].mean()
    df_pegRNA_score = pd.merge(
        df,
        df_log2_mean,
        on=["var_key"],
        how="left",
        validate="many_to_one",
        suffixes=("", "_mean_variant"),
    )

    return df_pegRNA_score


def normalise_log2_scores_per_replicate(df):
    """
    Normalise mean pegRNA log2-ratios (function scores) to the median mean pegRNA log2-ratio (function score) of synonymous variants for each condition.

    Parameters
    ----------
    df (df): Dataframe to be processed.
    """
    for var_dict in variables_per_sample:
        df[
            f"{var_dict['pegRNA_D20']}_{var_dict['pegRNA_D34']}_log2_mean_variant_normalised"
        ] = (
            df[f"{var_dict['pegRNA_D20']}_{var_dict['pegRNA_D34']}_log2_mean_variant"]
            - df.loc[df["Consequence"] == "SYNONYMOUS"][
                f"{var_dict['pegRNA_D20']}_{var_dict['pegRNA_D34']}_log2_mean_variant"
            ].median()
        )
    return df


def average_variant_log2_per_cell_line(df):
    """
    Average normalised mean pegRNA log2-ratios (function scores) per cell line.

    Parameters
    ----------
    df (df): Dataframe to be processed.
    """
    PEmax_log2_mean_variant_list = []
    for var_dict in PEmax_variables_per_sample:
        PEmax_log2_mean_variant_list.append(
            f"{var_dict['pegRNA_D20']}_{var_dict['pegRNA_D34']}_log2_mean_variant_normalised"
        )

    df["log2_mean_variant_PEmax"] = df[PEmax_log2_mean_variant_list].mean(axis=1)

    PEmaxdn_log2_mean_variant_list = []
    for var_dict in PEmaxdn_variables_per_sample:
        PEmaxdn_log2_mean_variant_list.append(
            f"{var_dict['pegRNA_D20']}_{var_dict['pegRNA_D34']}_log2_mean_variant_normalised"
        )

    df["log2_mean_variant_PEmaxdn"] = df[PEmaxdn_log2_mean_variant_list].mean(axis=1)

    return df


def normalise_variant_log2_scores_per_cell_line(df):
    """
    Normalise function scores that have been averaged per cell line to median mean function score of synonymous variants per cell line.
    Parameters
    ----------
    df (df): Dataframe to be processed.
    """
    df["log2_mean_variant_PEmax_normalised"] = (
        df["log2_mean_variant_PEmax"]
        - df.loc[df["Consequence"] == "SYNONYMOUS"]["log2_mean_variant_PEmax"].median()
    )

    df["log2_mean_variant_PEmaxdn_normalised"] = (
        df["log2_mean_variant_PEmaxdn"]
        - df.loc[df["Consequence"] == "SYNONYMOUS"][
            "log2_mean_variant_PEmaxdn"
        ].median()
    )

    return df


# ----------------------------------------------------------------------------------------------------------------------
# Calculate variant-level ST-editing rates


def collapse_ST_editing(df):
    """
    Average surrogate target editing rates of pegRNAs programming the same variant.

    Parameters
    ----------
    df (df): Dataframe to be processed.
    """
    D20_editing_list = []
    for sample in D20_list:
        D20_editing_list.append(f"{sample}_percentage_editing")
    df_ST_editing_mean = df.groupby(["var_key"])[D20_editing_list].mean()
    df_pegRNA_score = pd.merge(
        df,
        df_ST_editing_mean,
        on=["var_key"],
        how="left",
        validate="many_to_one",
        suffixes=("", "_mean_variant"),
    )

    return df_pegRNA_score


def average_variant_ST_editing_per_cell_line(df):
    """
    Average mean surrogate target editing rates per variant per cell line.

    Parameters
    ----------
    df (df): Dataframe to be processed.
    """
    PEmax_D20_editing_mean_variant_list = []
    for sample in PEmax_D20_list:
        PEmax_D20_editing_mean_variant_list.append(
            f"{sample}_percentage_editing_mean_variant"
        )

    df["PEmax_percentage_editing_mean_variant"] = df[
        PEmax_D20_editing_mean_variant_list
    ].mean(axis=1)

    PEmaxdn_D20_editing_mean_variant_list = []
    for sample in PEmaxdn_D20_list:
        PEmaxdn_D20_editing_mean_variant_list.append(
            f"{sample}_percentage_editing_mean_variant"
        )

    df["PEmaxdn_percentage_editing_mean_variant"] = df[
        PEmaxdn_D20_editing_mean_variant_list
    ].mean(axis=1)

    return df


# ----------------------------------------------------------------------------------------------------------------------
# Generate column with max SpliceAI score (out of four scores provided by CADD)


def max_SpliceAI_score(df):
    """
    Generate column with max SpliceAI score (out of four scores provided by CADD).

    Parameters
    ----------
    df (df): Dataframe to be processed.
    """
    df["SpliceAI_max"] = df[
        [
            "SpliceAI-acc-gain",
            "SpliceAI-acc-loss",
            "SpliceAI-don-gain",
            "SpliceAI-don-loss",
        ]
    ].max(axis=1)
    return df


# ----------------------------------------------------------------------------------------------------------------------
#                                               MAIN FUNCTION
# ----------------------------------------------------------------------------------------------------------------------
def MLH1x10_pegRNA_analysis(pegRNA_D20_D34_path, variables_per_sample):
    """
    Main pre-processing function to generate a dataframe containing function scores for the MLH1 exon 10 screen.

    Parameters
    ----------
    pegRNA_D20_D34_path (path): Path to input dataframe.
    variables_per_sample (list): List of dictionaries used to link sample collection date to sample identifier.
    """
    ST_filter = f"{ST_editing_filter_PEmax}_PEmax_{ST_editing_filter_PEmaxdn}_PEmaxdn"
    save_dir = OUTPUT_DIR / f"pseudo_pre_freq_{pseudo_pre_freq_filter}" / ST_filter
    save_dir.mkdir(exist_ok=True, parents=True)
    df_pegRNA = pd.read_csv(pegRNA_D20_D34_path, delimiter=",")

    # Data processing
    df_pegRNA_processed = df_pegRNA.copy().pipe(process_pegRNA_data_frame)
    for var_dict in variables_per_sample:
        df_pegRNA_processed = (
            df_pegRNA_processed.pipe(
                pegRNA_frequency_and_log2,
                var_dict["pegRNA_D20"],
                var_dict["pegRNA_D34"],
            )
            .pipe(filter_pegRNAs_on_pre_freq, var_dict)
            .pipe(filter_pegRNAs_on_ST_editing_rates, var_dict)
        )

    df_variant_score = (
        df_pegRNA_processed.pipe(collapse_variants)
        .pipe(normalise_log2_scores_per_replicate)
        .pipe(average_variant_log2_per_cell_line)
        .pipe(normalise_variant_log2_scores_per_cell_line)
        .pipe(collapse_ST_editing)
        .pipe(average_variant_ST_editing_per_cell_line)
        .pipe(max_SpliceAI_score)
    )

    # Data output
    df_variant_score.to_csv(
        (save_dir / f"data_MLH1x10_variant_score_per_cell_line.csv").as_posix()
    )


# ----------------------------------------------------------------------------------------------------------------------

if __name__ == "__main__":
    MLH1x10_pegRNA_analysis(pegRNA_D20_D34_path, variables_per_sample)
