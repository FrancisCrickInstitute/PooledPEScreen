"""
    Project :  Prime editing pilot screens. Revision experiment. Correlation of variant scores across replicates and cell lines.
    Date : 240925
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
import colorcet as cc

# PATHS

INPUT_DIR = pl.Path("/Users/kajbac/code/240830_revision_experiment/pegRNA_score")

OUTPUT_DIR = pl.Path("/Users/kajbac/code/240830_revision_experiment/_plots")
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)

df_variant_score_path = INPUT_DIR / "data_MLH1_revision_pegRNA_score.csv"

save_dir = OUTPUT_DIR
save_dir.mkdir(exist_ok=True, parents=True)
# VARIABLES

HAP1_sample_list = ["CK001_CK006_log2_normalised", "CK002_CK007_log2_normalised"]
K562_sample_list = ["CK003_CK008_log2_normalised", "CK004_CK009_log2_normalised"]

variant_naming_dict = {
    1966: "V12",
    1980: "V13",
    2029: "V6",
    1856: "V3",
    1978: "V1",
    2289: "V9",
    2297: "V10",
    2353: "V5",
    2358: "V11",
    2818: "V7",
    3128: "V8",
    2361: "V4",
    2385: "V2",
}

function_score_naming_dict = {
    "CK001_CK006_log2_normalised": "HAP1:PE7 R1",
    "CK002_CK007_log2_normalised": "HAP1:PE7 R2",
    "CK003_CK008_log2_normalised": "K562:PE7 R1",
    "CK004_CK009_log2_normalised": "K562:PE7 R2",
}

sort_dict_plot = {
    "V1": 0,
    "V2": 1,
    "V3": 2,
    "V4": 3,
    "V5": 4,
    "V6": 5,
    "V7": 6,
    "V8": 7,
    "V9": 8,
    "V10": 9,
    "V11": 10,
    "V12": 11,
    "V13": 12,
}
# --------------------------------------------------------------------------------------------------------------------
# FUNCTIONS


def generate_df_variant_correlation_across_replicates(df):
    """
    Generate dataframe containing function scores of four screen conditions.

    Parameters
    ----------
    df (df): Dataframe to be processed.
    """
    column_list = [
        "Consequence",
        "CK001_CK006_log2_normalised",
        "CK002_CK007_log2_normalised",
        "CK003_CK008_log2_normalised",
        "CK004_CK009_log2_normalised",
        "variant_name",
    ]
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
    # hue_order,
    # naming_dict,
    # color,
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
    color_index_to_skip = (
        7  # This is the index corresponding to V8 (deep intronic variant)
    )
    # Create a new palette without the color for V8
    original_palette = sns.color_palette(cc.glasbey, n_colors=13)
    new_palette = [
        color for i, color in enumerate(original_palette) if i != color_index_to_skip
    ]
    sns.set_palette(sns.color_palette(new_palette))

    sns.pairplot(
        df,
        hue=hue,
        # hue_order=hue_order,
        vars=vars,
        corner=True,
        height=1.5,
        diag_kind="kde",
        plot_kws={"s": 15, "linewidth": 0.2, "edgecolor": "black", "alpha": 0.6},
    )
    if save_path is None:
        plt.show()
    else:
        plt.savefig(save_path)
        # plt.savefig(save_path.with_suffix(".pdf"))
        # plt.savefig(save_path.with_suffix(".png"), dpi=300)


# ----------------------------------------------------------------------------------------------------------------------
# MAIN


def function_score_correlation_analysis(df_variant_score):
    """
    Plotting function generating pairplot of function scores.
    """
    df_correlation_across_replicates = (
        generate_df_variant_correlation_across_replicates(df_variant_score)
    )

    df_correlation_across_replicates = df_correlation_across_replicates.pipe(
        calculate_correlation_coefficient,
        HAP1_sample_list[0],
        HAP1_sample_list[1],
    )
    df_correlation_across_replicates = df_correlation_across_replicates.pipe(
        calculate_correlation_coefficient,
        HAP1_sample_list[0],
        K562_sample_list[0],
    )
    df_correlation_across_replicates = df_correlation_across_replicates.pipe(
        calculate_correlation_coefficient,
        HAP1_sample_list[0],
        K562_sample_list[1],
    )

    df_correlation_across_replicates = df_correlation_across_replicates.pipe(
        calculate_correlation_coefficient,
        HAP1_sample_list[1],
        K562_sample_list[0],
    )
    df_correlation_across_replicates = df_correlation_across_replicates.pipe(
        calculate_correlation_coefficient,
        HAP1_sample_list[1],
        K562_sample_list[1],
    )
    df_correlation_across_replicates = df_correlation_across_replicates.pipe(
        calculate_correlation_coefficient,
        K562_sample_list[0],
        K562_sample_list[1],
    )

    df_correlation_across_replicates.rename(
        columns=function_score_naming_dict, inplace=True
    )
    correlation_pair_plot(
        df_correlation_across_replicates,
        # naming_dict=consequence_naming_dict,
        # color=consequence_colors,
        hue="variant_name",
        # hue_order=consequence_order,
        vars=[
            "HAP1:PE7 R1",
            "HAP1:PE7 R2",
            "K562:PE7 R1",
            "K562:PE7 R2",
        ],
        save_path=save_dir
        / (
            f"plot_Supp_Figure_14b_HAP1_K562_function_score_correlation.svg"
        ),
    )

# ----------------------------------------------------------------------------------------------------------------------
# MAIN


def main():
    """
    Main function to plot correlation pairplots of function scores obtained in mini-pool experiment in HAP1 and K562 cells color-coded by variant consequence (Supp. Figure 14b).
    """
    df_variant_score = pd.read_csv(df_variant_score_path)
    df_variant_score["variant_name"] = df_variant_score["mutant index"].map(
        variant_naming_dict
    )
    df_variant_score = df_variant_score.sort_values(
        by=["variant_name"], key=lambda x: x.map(sort_dict_plot)
    )
    df_variant_score = df_variant_score[df_variant_score["variant_name"] != "V8"]
    set_style(context="paper", font_scale=1, style="ticks")
    function_score_correlation_analysis(df_variant_score)


# ----------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    main()
