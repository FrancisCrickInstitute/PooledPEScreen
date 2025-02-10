"""
    Project :  Prime editing pilot screens. Revision experiments. Mini-pool screen. PegRNA data analysis.
    Date : 240830
    Python version 3.11.8

"""

# IMPORT

import pandas as pd
import pathlib as pl
import math
import numpy as np

# VARIABLES

INPUT_DIR = pl.Path("/Users/kajbac/code/240830_revision_experiment/")
INPUT_DIR_CADD_CLINVAR = pl.Path("/Users/kajbac/code/230904_S002_final_filter")
INPUT_DIR_VARIANT_SCORE = pl.Path(
    "/Users/kajbac/code/240116_final_scripts_for_PE_paper/supplementary_tables/"
)

OUTPUT_DIR = pl.Path("/Users/kajbac/code/240830_revision_experiment/pegRNA_score")
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)

data_path = INPUT_DIR / "MLH1_pegRNA_freqs_edits_merged.csv"
cadd_clinvar_x10_path = (
    INPUT_DIR_CADD_CLINVAR / "MLH1x10_pegRNA_data_merged_best_pegRNA.csv"
)
cadd_clinvar_non_coding_path = (
    INPUT_DIR_CADD_CLINVAR / "MLH1intronic_pegRNA_data_merged_best_pegRNA.csv"
)
x10_variant_score_path = (
    INPUT_DIR_VARIANT_SCORE / "Supp_table_MLH1_x10_variant_scores.csv"
)
non_coding_variant_score_path = (
    INPUT_DIR_VARIANT_SCORE / "Supp_table_MLH1_non_coding_variant_scores.csv"
)


variables_per_sample = [
    {"pegRNA_D19": "CK001", "pegRNA_D33": "CK006"},
    {"pegRNA_D19": "CK002", "pegRNA_D33": "CK007"},
    {"pegRNA_D19": "CK003", "pegRNA_D33": "CK008"},
    {"pegRNA_D19": "CK004", "pegRNA_D33": "CK009"},
]
HAP1_list = ["CK001", "CK002", "CK006", "CK007"]
K562_list = ["CK003", "CK004", "CK008", "CK009"]
D19_list = ["CK001", "CK002", "CK003", "CK004"]

negative_control_indeces = [1978, 2385]
positive_control_indeces = [1856, 2361]

pseudo_pre_freq_filter = 1.4e-4
ST_editing_filter = 5

# ----------------------------------------------------------------------------------------------------------------------
#                                                     FUNCTIONS
# ----------------------------------------------------------------------------------------------------------------------


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


def add_pseudo_count(
    df,
    pegRNA_D19,
    pegRNA_D33,
):
    """
    Add pseudocount of one to pegRNA count.

    Parameters
    ----------
    df (df): Dataframe to be processed.
    pegRNA_D19 (str) : D19 sample identifier
    pegRNA_D33 (str) : D33 sample identifier
    """
    df[f"{pegRNA_D19}_pegRNA_pseudo_count"] = df[f"{pegRNA_D19}_pegRNA_count"] + 1
    df[f"{pegRNA_D33}_pegRNA_pseudo_count"] = df[f"{pegRNA_D33}_pegRNA_count"] + 1
    return df


def calculate_frequencies_ratio_log2_pegRNA(df, pegRNA_D19, pegRNA_D33):
    """
    Calculate pegRNA frequencies, post- over pre-frequency ratios, and log2-ratios from pegRNA pseudocounts.

    Parameters
    ----------
    df (df): Dataframe to be processed.
    pegRNA_D19 (str) : D19 sample identifier
    pegRNA_D33 (str) : D33 sample identifier
    """
    df[f"{pegRNA_D19}_pre_pseudo_freq"] = df[f"{pegRNA_D19}_pegRNA_pseudo_count"] / df[
        f"{pegRNA_D19}_pegRNA_pseudo_count"
    ].sum(axis=0)
    df[f"{pegRNA_D33}_post_pseudo_freq"] = df[f"{pegRNA_D33}_pegRNA_pseudo_count"] / df[
        f"{pegRNA_D33}_pegRNA_pseudo_count"
    ].sum(axis=0)
    df[f"{pegRNA_D33}_{pegRNA_D19}_ratio"] = df[
        f"{pegRNA_D33}_post_pseudo_freq"
    ].astype("float64") / df[f"{pegRNA_D19}_pre_pseudo_freq"].astype("float64")

    df[f"{pegRNA_D19}_{pegRNA_D33}_log2"] = df[
        f"{pegRNA_D33}_{pegRNA_D19}_ratio"
    ].apply(math.log2)
    return df


def process_pegRNA_data_frame(
    df,
):
    """
    Lift GRCh37 locations to GRCh38 locations, filter out variants with high correct ST editing in negative control, and generate a string identifier for each variant used for dataframe merging.

    Parameters
    ----------
    df (df): Dataframe to be processed.
    """
    df = df.pipe(translate_GRCh37locations_to_GRCh38locations).pipe(
        add_var_key_for_dataframe_merging_pegRNA
    )

    return df


def pegRNA_frequency_and_log2(df, pegRNA_D19, pegRNA_D33):
    """
    Add pseudocount of one to pegRNA count, and calculate pegRNA frequencies, post- over pre-frequency ratios, and log2-ratios from pegRNA pseudocounts.

    Parameters
    ----------
    df (df): Dataframe to be processed.
    """
    df = df.pipe(add_pseudo_count, pegRNA_D19, pegRNA_D33).pipe(
        calculate_frequencies_ratio_log2_pegRNA, pegRNA_D19, pegRNA_D33
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
        df[f"{var_dict['pegRNA_D19']}_pre_pseudo_freq"] < pseudo_pre_freq_filter,
        [
            f"{var_dict['pegRNA_D19']}_{var_dict['pegRNA_D33']}_log2",
            f"{var_dict['pegRNA_D19']}_percentage_editing",
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

    df.loc[
        df[f"{var_dict['pegRNA_D19']}_percentage_editing"] < ST_editing_filter,
        [
            f"{var_dict['pegRNA_D19']}_{var_dict['pegRNA_D33']}_log2",
            f"{var_dict['pegRNA_D19']}_percentage_editing",
        ],
    ] = np.nan

    return df


# ----------------------------------------------------------------------------------------------------------------------
# Calculate pegRNA-level scores


def normalise_log2_scores_per_condition(df):
    """
    Normalise pegRNA log2-ratios to the median pegRNA log2-ratio of control variants for each condition.

    Parameters
    ----------
    df (df): Dataframe to be processed.
    """
    for var_dict in variables_per_sample:
        df[f"{var_dict['pegRNA_D19']}_{var_dict['pegRNA_D33']}_log2_normalised"] = (
            df[f"{var_dict['pegRNA_D19']}_{var_dict['pegRNA_D33']}_log2"]
            - df.loc[df["control"] == "negative"][
                f"{var_dict['pegRNA_D19']}_{var_dict['pegRNA_D33']}_log2"
            ].median()
        )
    return df


# -----------------------------------------------------------------------------------------------------------------------------------------------------------------
# add annotations and data from previous MLH1 x10 and non-coding screens


def add_CADD_and_ClinVar_consequence(
    df, df_cadd_clinvar_x10, df_cadd_clinvar_non_coding
):
    """
    Extracting Cadd consequence and ClinVar pathogenicity information from input dataframes of the big MLH1 screens and adding this information to the output dataframe.

     Parameters
     ----------
     df (df): Dataframe to be processed.
     df_cadd_clinvar_x10 (df): Dataframe containing Cadd and ClinVar information from x10 screen.
     df_cadd_clinvar_non_coding (df): Dataframe containing Cadd and ClinVar information from non coding screen.
    """

    # pre-process x10 and non-coding dataframes and generate single dataframe
    df_cadd_clinvar_x10_unique = df_cadd_clinvar_x10.drop_duplicates(
        subset="mutant index"
    )
    df_cadd_clinvar_non_coding_unique = df_cadd_clinvar_non_coding.drop_duplicates(
        subset="mutant index"
    )
    df_cadd_clinvar = pd.concat(
        [df_cadd_clinvar_x10_unique, df_cadd_clinvar_non_coding_unique],
        axis=0,
        ignore_index=True,
    )

    # merge cadd and clinvar dataframe with pegRNA dataframe of the revision experiment

    df_merged = pd.merge(df, df_cadd_clinvar, on="mutant index", how="left")

    print(df_merged)

    return df_merged


def add_function_scores(df, df_x10_variant_score, df_non_coding_variant_score):
    """
    Extracting PEmax function scores and overall function scores from final variant score dataframes of the big MLH1 screens and adding this information to the output dataframe.

    Parameters
    ----------
    df (df): Dataframe to be processed.
    df_x10_variant_score (df): Final variant score dataframe from MLH1 x10 screen.
    df_non_coding_variant_score (df): Final variant score dataframe from MLH1 non-coding screen.

    """

    # concatenate function score dataframes
    df_function_scores = pd.concat(
        [df_x10_variant_score, df_non_coding_variant_score], axis=0, ignore_index=True
    )

    # merge function score dataframe with pegRNA dataframe of the revision experiment

    df_merged = pd.merge(df, df_function_scores, on="mutant index", how="left")

    print(df_merged)

    return df_merged


def change_start_loss_consequence_annotation(df):
    """
    Changing Cadd "Consequence" annotation for variant 8 from "5' UTR" to "start_loss".

    Parameters
    ----------
    df (df): Dataframe to be processed.

    """
    df.loc[df["mutant index"] == 2353, "Consequence"] = "START_LOSS"

    return df


# ----------------------------------------------------------------------------------------------------------------------
# MAIN FUNCTION
# ----------------------------------------------------------------------------------------------------------------------


def pegRNA_analysis(data_path, variables_per_sample):
    """
    Main pre-processing function to generate a dataframe containing pegRNA scores for the MLH1 mini-pool screen.

    Parameters
    ----------
    data_path (path): Path to input dataframe.
    variables_per_sample (list): List of dictionaries used to link sample collection date to sample identifier.
    """

    save_dir = OUTPUT_DIR
    save_dir.mkdir(exist_ok=True, parents=True)
    df_pegRNA = pd.read_csv(data_path, delimiter=",")
    df_cadd_clinvar_x10 = pd.read_csv(
        cadd_clinvar_x10_path,
        usecols=[
            "mutant index",
            "Consequence",
            "PHRED",
            "ClinVar_Significance_Simple",
            "SpliceAI-acc-gain",
            "SpliceAI-acc-loss",
            "SpliceAI-don-gain",
            "SpliceAI-don-loss",
        ],
        delimiter=",",
    )
    df_cadd_clinvar_non_coding = pd.read_csv(
        cadd_clinvar_non_coding_path,
        usecols=[
            "mutant index",
            "Consequence",
            "PHRED",
            "ClinVar_Significance_Simple",
            "SpliceAI-acc-gain",
            "SpliceAI-acc-loss",
            "SpliceAI-don-gain",
            "SpliceAI-don-loss",
        ],
        delimiter=",",
    )
    df_x10_variant_score = pd.read_csv(
        x10_variant_score_path,
        usecols=["mutant index", "PEmax_Function_Score", "Function_Score"],
        delimiter=",",
    )
    df_non_coding_variant_score = pd.read_csv(
        non_coding_variant_score_path,
        usecols=["mutant index", "PEmax_Function_Score", "Function_Score"],
        delimiter=",",
    )

    # Data processing
    df_pegRNA_processed = df_pegRNA.copy().pipe(process_pegRNA_data_frame)
    for var_dict in variables_per_sample:
        df_pegRNA_processed = (
            df_pegRNA_processed.pipe(
                pegRNA_frequency_and_log2,
                var_dict["pegRNA_D19"],
                var_dict["pegRNA_D33"],
            ).pipe(filter_pegRNAs_on_pre_freq, var_dict)
            #   .pipe(filter_pegRNAs_on_ST_editing_rates,var_dict)
        )

    df_processed_cadd_clinvar = df_pegRNA_processed.pipe(
        add_CADD_and_ClinVar_consequence,
        df_cadd_clinvar_x10,
        df_cadd_clinvar_non_coding,
    ).pipe(change_start_loss_consequence_annotation)

    df_processed_function_scores = df_processed_cadd_clinvar.pipe(
        add_function_scores, df_x10_variant_score, df_non_coding_variant_score
    )

    df_processed_function_scores["control"] = np.where(
        df_processed_function_scores["mutant index"].isin(positive_control_indeces),
        "positive",
        np.where(
            df_processed_function_scores["mutant index"].isin(negative_control_indeces),
            "negative",
            np.nan,
        ),
    )

    df_pegRNA_score_normalised = df_processed_function_scores.pipe(
        normalise_log2_scores_per_condition
    )
    df_pegRNA_score_normalised["HAP1:PE7_score"] = df_pegRNA_score_normalised[
        ["CK001_CK006_log2_normalised", "CK002_CK007_log2_normalised"]
    ].mean(axis=1)
    df_pegRNA_score_normalised["K562:PE7_score"] = df_pegRNA_score_normalised[
        ["CK003_CK008_log2_normalised", "CK004_CK009_log2_normalised"]
    ].mean(axis=1)

    # Data output
    df_pegRNA_score_normalised.to_csv(
        (save_dir / f"data_MLH1_revision_pegRNA_score.csv").as_posix()
    )
    # df_pegRNA_score.to_csv((save_dir / f"data_MLH1_revision_pegRNA_score.csv").as_posix())


# ----------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    df_pegRNA_score = pegRNA_analysis(data_path, variables_per_sample)
