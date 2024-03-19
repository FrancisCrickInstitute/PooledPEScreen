"""
    Project :  Prime editing pilot screens. Plotting script to generate Figure 3 g and Supplementary Figure 9 e. Correlation Scatterplot: functions scores (variant-level scores) against endogenous function scores. MLH1 saturation screen.
    Date : 240316
    Python version 3.10

"""

# IMPORT

import pandas as pd
import pathlib as pl
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy import stats
from scipy.stats import norm
from matplotlib.offsetbox import AnchoredText

# VARIABLES

pseudo_pre_freq_filter = 1.4e-4
ST_editing_filter_PEmax = 5
ST_editing_filter_PEmaxdn = 25

# PATHS

INPUT_DIR_END = pl.Path(
    f"/Users/kajbac/code/240116_final_scripts_for_PE_paper/MLH1x10/endogenous/log2_enriched"
)

INPUT_DIR_PEGRNA = pl.Path(
    f"/Users/kajbac/code/240116_final_scripts_for_PE_paper/MLH1x10/pegRNA/variant_score/pseudo_pre_freq_{pseudo_pre_freq_filter}/{ST_editing_filter_PEmax}_PEmax_{ST_editing_filter_PEmaxdn}_PEmaxdn"
)

df_endogenous_pegRNA_concat_path = (
    INPUT_DIR_END / "data_expected_variants_replicates_log2_filtered.csv"
)
df_variant_score_path = INPUT_DIR_PEGRNA / "data_MLH1x10_variant_score_replicates.csv"

OUTPUT_DIR = pl.Path(
    "/Users/kajbac/code/240116_final_scripts_for_PE_paper/plotting/_plots"
)
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)

# number of conditions in which variant needed to receive score (pegRNA level)
replicate_number = 2

# Renaming dictionary
consequence_naming_dict = {
    "SYNONYMOUS": "synonymous",
    "STOP_GAINED": "nonsense",
    "CANONICAL_SPLICE": "canonical splice",
    "NON_SYNONYMOUS": "missense",
    "SPLICE_SITE": "splice site",
    "INTRONIC": "intronic",
}

# Color palette
consequence_colors = ["#4e77bb", "#e83677", "#f7a83e", "#64bb97", "#243672", "#8acdef"]

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


def determine_fdr(df):
    """
    FDR testing of pegRNA-derived function scores with the following thresholds: 0.01, 0.02, 0.05 using synonymous variants as controls.

    Parameters
    ----------
    df (df): Dataframe to be processed.
    """
    loc, scale = stats.norm.fit(
        df.loc[df["Consequence"] == "SYNONYMOUS"][
            "log2_mean_variant_replicates_normalised"
        ]
    )
    df["p_values"] = 1 - stats.norm.cdf(
        df["log2_mean_variant_replicates_normalised"], loc, scale
    )
    df["p_adjust"] = stats.false_discovery_control(df["p_values"], method="bh")
    df["fdr_0.01"] = df["p_adjust"].map(lambda x: True if x < 0.01 else False)
    df["fdr_0.02"] = df["p_adjust"].map(lambda x: True if x < 0.02 else False)
    df["fdr_0.05"] = df["p_adjust"].map(lambda x: True if x < 0.05 else False)

    return df


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


def pegRNA_variant_score_correlation_scatter(
    df, hue, hue_order, colors, naming_dict, ax, save_path
):
    """
    Plotting function generating correlation scatterplot with function scores against endogenous function scores, marking LoF variants as stars.

    Parameters
    ----------
    df (df): Dataframe to be plotted.
    hue: Scatterplot parameter.
    hue_order: Scatterplot parameter.
    colors (list): List of hex color codes to set color palette.
    naming_dict (dict): Dictionary for renaming hue annotations.
    ax: Scatterplot parameter.
    save_path (path): Output path.
    """
    sns.set_palette(sns.color_palette(colors))

    df.rename(
        columns={
            "log2_mean_variant_replicates_normalised": "Function score",
            "D34_D20_log2_mean_replicates_normalised": "Endogenous function score",
        },
        inplace=True,
    )

    df = df.dropna(subset="Function score")
    df = df.dropna(subset="Endogenous function score")

    end_variants = df["Endogenous function score"].count()
    print(f"number of endogenous variants plotted is {end_variants}")
    peg_variants = df["Function score"].count()
    print(f"number of peg variants plotted is {peg_variants}")

    r = df["Function score"].corr(df["Endogenous function score"], method="pearson")
    r = round(r, 3)
    df[hue].replace(to_replace=naming_dict, inplace=True)

    df_lof = df.copy().loc[df["fdr_0.01"] == True]
    df_null = df.copy().loc[df["fdr_0.01"] == False]

    sns.scatterplot(
        data=df_null,
        x="Endogenous function score",
        y="Function score",
        hue=hue,
        hue_order=hue_order,
        linewidth=0.2,
        edgecolors="black",
        alpha=0.7,
        s=10,
        ax=ax,
    )
    sns.scatterplot(
        data=df_lof,
        x="Endogenous function score",
        y="Function score",
        hue=hue,
        hue_order=hue_order,
        marker=(8, 1, 0),
        s=30,
        linewidth=0.2,
        edgecolors="black",
        alpha=0.6,
        ax=ax,
    )

    ax.spines[["right", "top"]].set_visible(False)
    anc = AnchoredText(
        f"r = {r}",
        loc="upper left",
        frameon=False,
    )
    ax.add_artist(anc)

    ax.get_legend().remove()

    if save_path is None:
        plt.show()
    else:
        plt.savefig(save_path)
        # plt.savefig(save_path.with_suffix(".pdf"))
        # plt.savefig(save_path.with_suffix(".png"), dpi=300)


def pegRNA_variant_score_correlation_scatter_hc(
    df, hue, hue_order, colors, naming_dict, ax, save_path
):
    """
    Plotting function generating correlation scatterplot with function scores against endogenous function scores, marking LoF variants as stars.

    Parameters
    ----------
    df (df): Dataframe to be plotted.
    hue: Scatterplot parameter.
    hue_order: Scatterplot parameter.
    colors (list): List of hex color codes to set color palette.
    naming_dict (dict): Dictionary for renaming hue annotations.
    ax: Scatterplot parameter.
    save_path (path): Output path.
    """
    sns.set_palette(sns.color_palette(colors))

    df.rename(
        columns={
            "log2_mean_variant_replicates_normalised": "Function score",
            "D34_D20_log2_mean_replicates_normalised": "Endogenous function score",
        },
        inplace=True,
    )

    df = df.dropna(subset="Function score")
    df = df.dropna(subset="Endogenous function score")

    r = df["Function score"].corr(df["Endogenous function score"], method="pearson")
    r = round(r, 3)
    df[hue].replace(to_replace=naming_dict, inplace=True)

    end_variants = df["Endogenous function score"].count()
    print(f"number of hc endogenous variants plotted is {end_variants}")
    peg_variants = df["Function score"].count()
    print(f"number of hc peg variants plotted is {peg_variants}")

    sns.scatterplot(
        data=df,
        x="Endogenous function score",
        y="Function score",
        hue=hue,
        hue_order=hue_order,
        linewidth=0.2,
        edgecolors="black",
        alpha=0.7,
        s=10,
        ax=ax,
    )

    ax.spines[["right", "top"]].set_visible(False)
    anc = AnchoredText(
        f"r = {r}",
        loc="upper left",
        frameon=False,
    )
    ax.add_artist(anc)

    ax.get_legend().remove()

    if save_path is None:
        plt.show()
    else:
        plt.savefig(save_path)
        # plt.savefig(save_path.with_suffix(".pdf"))
        # plt.savefig(save_path.with_suffix(".png"), dpi=300)


# --------------------------------------------------------------------------------------------------------------------
# MAIN


def main():
    """
    Main function to plot Figure 3 g and for PE manuscript. Scatterplot: functions scores (variant scores) against endogenous function scores for both relaxed and high-stringency datasets. MLH1 saturation screen.
    """
    save_dir = OUTPUT_DIR
    save_dir.mkdir(exist_ok=True, parents=True)
    save_path = save_dir / (f"Fig_3g_plot_function_endogenous_score_correlation.svg")
    save_path_hc = save_dir / (
        f"Supp_fig9e_plot_function_endogenous_score_correlation_hc.svg"
    )

    set_style(context="paper", font_scale=1, style="ticks")

    df_endogenous = pd.read_csv(
        df_endogenous_pegRNA_concat_path, delimiter=",", low_memory=False
    )

    df_variant_score = pd.read_csv(
        df_variant_score_path, delimiter=",", low_memory=False
    )

    # Replicate analysis to filter pegRNA-derived function scores
    filter_number = replicate_number - 1
    df_variant_score_filtered = df_variant_score.loc[
        df_variant_score["number_replicates_pegRNA"] > filter_number
    ]
    df_variant_score_filtered = df_variant_score_filtered.drop_duplicates(
        subset=["var_key"]
    )  # important! this dataframe now contains variant level information!

    # Determine fdr
    df_variant_score_filtered = df_variant_score_filtered.pipe(determine_fdr)

    # Merge dataframes
    df_endogenous_score_short = df_endogenous[
        ["var_key", "D34_D20_log2_mean_replicates_normalised"]
    ]

    df_merge = pd.merge(
        df_variant_score_filtered,
        df_endogenous_score_short,
        on=["var_key"],
        how="left",
    )

    # Plotting
    fig, ax = plt.subplots(figsize=(3, 1.5))

    pegRNA_variant_score_correlation_scatter(
        df=df_merge,
        hue="Consequence",
        hue_order=consequence_order,
        colors=consequence_colors,
        naming_dict=consequence_naming_dict,
        ax=ax,
        save_path=save_path,
    )
    # High-stringency dataset
    fig, ax = plt.subplots(figsize=(3, 1.5))

    pegRNA_variant_score_correlation_scatter_hc(
        df=df_merge.loc[df_merge["hc"] == True],
        hue="Consequence",
        hue_order=consequence_order,
        colors=consequence_colors,
        naming_dict=consequence_naming_dict,
        ax=ax,
        save_path=save_path_hc,
    )


# ----------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    main()
