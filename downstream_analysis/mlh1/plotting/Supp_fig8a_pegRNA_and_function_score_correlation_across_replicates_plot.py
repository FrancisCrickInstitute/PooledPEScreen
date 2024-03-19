"""
    Project :  Prime editing pilot screens. Plotting script to generate Supplementary Figure 8a. Correlation of pegRNA scores, before and after ST-editing filter. Correlation of variant scores. MLH1 saturation screen.
    Date : 240316
    Python version 3.10

"""

# IMPORT

import pandas as pd
import pathlib as pl
import math
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.offsetbox import AnchoredText
from scipy.stats import pearsonr

# VARIABLES

pseudo_pre_freq_filter = 0.00014
ST_editing_filter_PEmax = 5
ST_editing_filter_PEmaxdn = 25

# PATHS

INPUT_DIR_PEGRNA_NO_FILT = pl.Path(
    f"/Users/kajbac/code/240116_final_scripts_for_PE_paper/MLH1x10/pegRNA/pegRNA_score/pseudo_pre_freq_{pseudo_pre_freq_filter}/0_PEmax_0_PEmaxdn"
)
INPUT_DIR_PEGRNA = pl.Path(
    f"/Users/kajbac/code/240116_final_scripts_for_PE_paper/MLH1x10/pegRNA/pegRNA_score/pseudo_pre_freq_{pseudo_pre_freq_filter}/{ST_editing_filter_PEmax}_PEmax_{ST_editing_filter_PEmaxdn}_PEmaxdn"
)
INPUT_DIR_VARIANT = pl.Path(
    f"/Users/kajbac/code/240116_final_scripts_for_PE_paper/MLH1x10/pegRNA/variant_score/pseudo_pre_freq_{pseudo_pre_freq_filter}/{ST_editing_filter_PEmax}_PEmax_{ST_editing_filter_PEmaxdn}_PEmaxdn"
)

df_pegRNA_score_no_filt_path = (
    INPUT_DIR_PEGRNA_NO_FILT / "data_MLH1x10_pegRNA_score_replicates.csv"
)
df_pegRNA_score_path = INPUT_DIR_PEGRNA / "data_MLH1x10_pegRNA_score_replicates.csv"
df_variant_score_path = INPUT_DIR_VARIANT / "data_MLH1x10_variant_score_replicates.csv"

OUTPUT_DIR = pl.Path(
    "/Users/kajbac/code/240116_final_scripts_for_PE_paper/plotting/_plots"
)
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)

save_dir = OUTPUT_DIR
save_dir.mkdir(exist_ok=True, parents=True)

variables_per_sample = [
    {"pegRNA_D20": "CK003", "pegRNA_D34": "CK007"},
    {"pegRNA_D20": "CK005", "pegRNA_D34": "CK009"},
    {"pegRNA_D20": "CK004", "pegRNA_D34": "CK008"},
    {"pegRNA_D20": "CK006", "pegRNA_D34": "CK010"},
]

# Renaming dictionaries
consequence_naming_dict = {
    "SYNONYMOUS": "synonymous",
    "STOP_GAINED": "nonsense",
    "CANONICAL_SPLICE": "canonical splice",
    "NON_SYNONYMOUS": "missense",
    "SPLICE_SITE": "splice site",
    "INTRONIC": "intronic",
}

pegRNA_score_naming_dict = {
    "CK003_CK007_log2_normalised": "PEmax pegRNA score",
    "CK005_CK009_log2_normalised": "PEmax-obn pegRNA score",
    "CK004_CK008_log2_normalised": "PEmax+MLH1dn pegRNA score",
    "CK006_CK010_log2_normalised": "PEmax+MLH1dn-obn pegRNA score",
}
function_score_naming_dict = {
    "CK003_CK007_log2_mean_variant_normalised": "PEmax function score",
    "CK005_CK009_log2_mean_variant_normalised": "PEmax-obn function score",
    "CK004_CK008_log2_mean_variant_normalised": "PEmax+MLH1dn function score",
    "CK006_CK010_log2_mean_variant_normalised": "PEmax+MLH1dn-obn function score",
}

# Color palette
consequence_colors = ["#4e77bb", "#e83677", "#f7a83e", "#64bb97", "#243672", "#8acdef"]
clinvar_colors = ["#0433FF", "#882255", "#DCA237"]

# Hue order
consequence_order = [
    "synonymous",
    "nonsense",
    "canonical splice",
    "missense",
    "splice site",
    "intronic",
]
# --------------------------------------------------------------------------------------------------------------------
# FUNCTIONS


def generate_df_pegRNA_correlation_across_replicates(df):
    """
    Generate dataframe containing pegRNA scores of four screen conditions.

    Parameters
    ----------
    df (df): Dataframe to be processed.
    """
    column_list = ["var_key"]
    for var_dict in variables_per_sample:
        column_list.append(
            f"{var_dict['pegRNA_D20']}_{var_dict['pegRNA_D34']}_log2_normalised"
        )
    column_list.append("Consequence")
    df_pegRNA_scores = df[column_list]

    return df_pegRNA_scores


def generate_df_variant_correlation_across_replicates(df):
    """
    Generate dataframe containing function scores of four screen conditions.

    Parameters
    ----------
    df (df): Dataframe to be processed.
    """
    column_list = ["var_key"]
    for var_dict in variables_per_sample:
        column_list.append(
            f"{var_dict['pegRNA_D20']}_{var_dict['pegRNA_D34']}_log2_mean_variant_normalised"
        )
    column_list.append("Consequence")
    df_variant_scores = df[column_list]

    return df_variant_scores


def calculate_correlation_coefficient(df, replicate_1_score, replicate_2_score):
    """
    Calculate Pearson correlation coefficient for pairwise comparisons.

    Parameters
    ----------
    df (df): Dataframe to be processed.
    replicate_1_score (int): Score of replicate 1.
    replicate_2_score (int): Score of replicate 2.
    """
    df[f"{replicate_1_score}_{replicate_2_score}_r"] = df[replicate_1_score].corr(
        df[replicate_2_score]
    )
    print(f"{replicate_1_score}_{replicate_2_score}_r")
    print(df[f"{replicate_1_score}_{replicate_2_score}_r"])

    return df


# --------------------------------------------------------------------------------------------------------------------
# PLOTS


def set_style(font_scale, context="paper", style="ticks"):
    """
    Set plotting style for figures in the manuscript.

    Parameters
    ----------
    context : str, optional
        Context for plotting. Default is "notebook".
        Options are "paper", "notebook", "talk", and "poster".
    style : str, optional
        Style for plotting. Default is "ticks".
        Options are "darkgrid", "whitegrid", "dark", "white", and "ticks".
    font : str, optional
        Font for plotting. Default is "Helvetica".
    font_scale : float, optional
        Font scale for plotting. Default is 1.
    """
    # matplotlib parameters
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42
    plt.rcParams["svg.fonttype"] = "none"
    plt.rcParams["font.family"] = "sans-serif"
    # seaborn settings and palettes
    sns.set(style=style, context=context, font_scale=font_scale)


def correlation_pair_plot(
    df,
    hue,
    hue_order,
    naming_dict,
    color,
    vars=None,
    save_path=None,
):
    """
    Plotting function generating pairplot of all screen conditions.

    Parameters
    ----------
    df (df): Dataframe to be plotted.
    hue: Pairplot parameter.
    hue_order: Pairplot parameter.
    naming_dict (dict): Dictionary for renaming hue annotations.
    color (list): List of hex color codes.
    vars: Pairplot parameter.
    save_path (path): Output path.
    """
    df[hue].replace(to_replace=naming_dict, inplace=True)
    sns.set_palette(sns.color_palette(color))

    sns.pairplot(
        df,
        hue=hue,
        hue_order=hue_order,
        vars=vars,
        corner=True,
        height=1.5,
        diag_kind="kde",
        plot_kws={"s": 10, "linewidth": 0.2, "edgecolor": "black", "alpha": 0.6},
    )
    if save_path is None:
        plt.show()
    else:
        plt.savefig(save_path)
        # plt.savefig(save_path.with_suffix(".pdf"))
        # plt.savefig(save_path.with_suffix(".png"), dpi=300)


# ----------------------------------------------------------------------------------------------------------------------
# MAIN


def MLH1_x10_pegRNA_correlation_analysis_no_filt():
    """
    Plotting function generating pairplot of pegRNA scores before ST-editing filter.
    """
    df_pegRNA_score = pd.read_csv(df_pegRNA_score_no_filt_path)
    df_correlation_across_replicates = generate_df_pegRNA_correlation_across_replicates(
        df_pegRNA_score
    )
    df_correlation_across_replicates = df_correlation_across_replicates.pipe(
        calculate_correlation_coefficient,
        "CK003_CK007_log2_normalised",
        "CK005_CK009_log2_normalised",
    )
    df_correlation_across_replicates = df_correlation_across_replicates.pipe(
        calculate_correlation_coefficient,
        "CK004_CK008_log2_normalised",
        "CK006_CK010_log2_normalised",
    )
    df_correlation_across_replicates = df_correlation_across_replicates.pipe(
        calculate_correlation_coefficient,
        "CK003_CK007_log2_normalised",
        "CK004_CK008_log2_normalised",
    )
    df_correlation_across_replicates = df_correlation_across_replicates.pipe(
        calculate_correlation_coefficient,
        "CK004_CK008_log2_normalised",
        "CK005_CK009_log2_normalised",
    )
    df_correlation_across_replicates = df_correlation_across_replicates.pipe(
        calculate_correlation_coefficient,
        "CK003_CK007_log2_normalised",
        "CK006_CK010_log2_normalised",
    )
    df_correlation_across_replicates = df_correlation_across_replicates.pipe(
        calculate_correlation_coefficient,
        "CK005_CK009_log2_normalised",
        "CK006_CK010_log2_normalised",
    )

    df_correlation_across_replicates.rename(
        columns=pegRNA_score_naming_dict, inplace=True
    )

    correlation_pair_plot(
        df_correlation_across_replicates,
        naming_dict=consequence_naming_dict,
        color=consequence_colors,
        hue="Consequence",
        hue_order=consequence_order,
        vars=[
            "PEmax pegRNA score",
            "PEmax-obn pegRNA score",
            "PEmax+MLH1dn pegRNA score",
            "PEmax+MLH1dn-obn pegRNA score",
        ],
        save_path=save_dir
        / (f"Supp_fig8a_plot_correlation_pegRNA_scores_across_replicates_pre_filt.svg"),
    )


def MLH1_x10_pegRNA_correlation_analysis():
    """
    Plotting function generating pairplot of pegRNA scores after ST-editing filter.
    """
    df_pegRNA_score = pd.read_csv(df_pegRNA_score_path)

    df_correlation_across_replicates = generate_df_pegRNA_correlation_across_replicates(
        df_pegRNA_score
    )
    df_correlation_across_replicates = df_correlation_across_replicates.pipe(
        calculate_correlation_coefficient,
        "CK003_CK007_log2_normalised",
        "CK005_CK009_log2_normalised",
    )
    df_correlation_across_replicates = df_correlation_across_replicates.pipe(
        calculate_correlation_coefficient,
        "CK004_CK008_log2_normalised",
        "CK006_CK010_log2_normalised",
    )
    df_correlation_across_replicates = df_correlation_across_replicates.pipe(
        calculate_correlation_coefficient,
        "CK003_CK007_log2_normalised",
        "CK004_CK008_log2_normalised",
    )
    df_correlation_across_replicates = df_correlation_across_replicates.pipe(
        calculate_correlation_coefficient,
        "CK004_CK008_log2_normalised",
        "CK005_CK009_log2_normalised",
    )
    df_correlation_across_replicates = df_correlation_across_replicates.pipe(
        calculate_correlation_coefficient,
        "CK003_CK007_log2_normalised",
        "CK006_CK010_log2_normalised",
    )
    df_correlation_across_replicates = df_correlation_across_replicates.pipe(
        calculate_correlation_coefficient,
        "CK005_CK009_log2_normalised",
        "CK006_CK010_log2_normalised",
    )

    df_correlation_across_replicates.rename(
        columns=pegRNA_score_naming_dict, inplace=True
    )

    correlation_pair_plot(
        df_correlation_across_replicates,
        naming_dict=consequence_naming_dict,
        color=consequence_colors,
        hue="Consequence",
        hue_order=consequence_order,
        vars=[
            "PEmax pegRNA score",
            "PEmax-obn pegRNA score",
            "PEmax+MLH1dn pegRNA score",
            "PEmax+MLH1dn-obn pegRNA score",
        ],
        save_path=save_dir
        / (
            f"Supp_fig8a_plot_correlation_pegRNA_scores_across_replicates_post_filt.svg"
        ),
    )


def MLH1_x10_variant_correlation_analysis():
    """
    Plotting function generating pairplot of function scores after ST-editing filter.
    """
    df_variant_score = pd.read_csv(df_variant_score_path)
    df_correlation_across_replicates = (
        generate_df_variant_correlation_across_replicates(df_variant_score)
    )
    df_correlation_across_replicates = df_correlation_across_replicates.drop_duplicates(
        subset="var_key"
    )

    df_correlation_across_replicates = df_correlation_across_replicates.pipe(
        calculate_correlation_coefficient,
        "CK003_CK007_log2_mean_variant_normalised",
        "CK005_CK009_log2_mean_variant_normalised",
    )
    df_correlation_across_replicates = df_correlation_across_replicates.pipe(
        calculate_correlation_coefficient,
        "CK004_CK008_log2_mean_variant_normalised",
        "CK006_CK010_log2_mean_variant_normalised",
    )
    df_correlation_across_replicates = df_correlation_across_replicates.pipe(
        calculate_correlation_coefficient,
        "CK003_CK007_log2_mean_variant_normalised",
        "CK004_CK008_log2_mean_variant_normalised",
    )
    df_correlation_across_replicates = df_correlation_across_replicates.pipe(
        calculate_correlation_coefficient,
        "CK004_CK008_log2_mean_variant_normalised",
        "CK005_CK009_log2_mean_variant_normalised",
    )
    df_correlation_across_replicates = df_correlation_across_replicates.pipe(
        calculate_correlation_coefficient,
        "CK003_CK007_log2_mean_variant_normalised",
        "CK006_CK010_log2_mean_variant_normalised",
    )
    df_correlation_across_replicates = df_correlation_across_replicates.pipe(
        calculate_correlation_coefficient,
        "CK005_CK009_log2_mean_variant_normalised",
        "CK006_CK010_log2_mean_variant_normalised",
    )

    df_correlation_across_replicates.rename(
        columns=function_score_naming_dict, inplace=True
    )
    correlation_pair_plot(
        df_correlation_across_replicates,
        naming_dict=consequence_naming_dict,
        color=consequence_colors,
        hue="Consequence",
        hue_order=consequence_order,
        vars=[
            "PEmax function score",
            "PEmax-obn function score",
            "PEmax+MLH1dn function score",
            "PEmax+MLH1dn-obn function score",
        ],
        save_path=save_dir
        / (
            f"Supp_fig8a_plot_correlation_function_scores_across_replicates_post_filt.svg"
        ),
    )


# ----------------------------------------------------------------------------------------------------------------------
# MAIN


def main():
    """
    Main function to plot Supplementary Figure 8 a for PE manuscript. Correlation pairplots of pegRNA scores pre and post St-editing filter, and function scores, color-coded by variant consequence. MLH1 saturation screen.
    """
    set_style(context="paper", font_scale=1, style="ticks")
    MLH1_x10_pegRNA_correlation_analysis_no_filt()
    MLH1_x10_pegRNA_correlation_analysis()
    MLH1_x10_variant_correlation_analysis()


# ----------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    main()
