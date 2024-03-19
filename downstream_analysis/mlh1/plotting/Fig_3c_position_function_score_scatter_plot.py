"""
    Project :  Prime editing pilot screens. Plotting script to generate Figure 3 c and Supplementary Figure 9 a. Scatterplot: function scores (variant scores) over positions (GRCh38). MLH1 saturation screen.
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

# VARIABLES

pseudo_pre_freq_filter = 1.4e-4
ST_editing_filter_PEmax = 5
ST_editing_filter_PEmaxdn = 25

# Number of replicates in which a variant has to be scored.
replicate_number = 2

# PATHS

INPUT_DIR = pl.Path(
    f"/Users/kajbac/code/240116_final_scripts_for_PE_paper/MLH1x10/pegRNA/variant_score/pseudo_pre_freq_{pseudo_pre_freq_filter}/{ST_editing_filter_PEmax}_PEmax_{ST_editing_filter_PEmaxdn}_PEmaxdn"
)

OUTPUT_DIR = pl.Path(
    "/Users/kajbac/code/240116_final_scripts_for_PE_paper/plotting/_plots"
)
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)

df_variant_score_path = INPUT_DIR / "data_MLH1x10_variant_score_replicates.csv"

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


# Plots


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


def position_scatterplot(df, hue_order, x, y, hue, naming_dict, save_path=None):
    """
    Plotting function generating scatterplot with function scores over GRCh38 locations, marking LoF variants as stars.

    Parameters
    ----------
    df (df): Dataframe to be plotted
    hue_order: Scatterplot parameter.
    x: Scatterplot parameter.
    y: Scatterplot parameter.
    hue: Scatterplot parameter
    naming_dict (dict): Dictionary for renaming hue annotations.
    save_path (path): Output path.
    """
    sns.set_palette(sns.color_palette(consequence_colors))

    variants = df["log2_mean_variant_replicates_normalised"].count()
    print(f"number of variants plotted is {variants}")

    df[hue].replace(to_replace=naming_dict, inplace=True)

    df_lof = df.copy().loc[df["fdr_0.01"] == True]
    df_null = df.copy().loc[df["fdr_0.01"] == False]

    fig, ax = plt.subplots(figsize=(5, 1.6))
    # ax.axhline(y=0, color="#D3D3D3", linestyle='--')
    sns.scatterplot(
        data=df_null,
        x=x,
        y=y,
        hue=hue,
        hue_order=hue_order,
        ax=ax,
        s=10,
        linewidth=0.2,
        edgecolor="black",
        alpha=0.7,
    )
    sns.scatterplot(
        data=df_lof,
        x=x,
        y=y,
        hue=hue,
        hue_order=hue_order,
        ax=ax,
        marker=(8, 1, 0),
        s=30,
        linewidth=0.2,
        edgecolors="black",
        alpha=0.7,
    )
    ax.spines[["right", "top"]].set_visible(False)
    sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1), frameon=False)

    ax.set_xlabel("GRCh38 position")
    ax.set_ylabel("Function score")

    if save_path is None:
        plt.show()
    else:
        plt.savefig(save_path)
        # plt.savefig(save_path.with_suffix(".pdf"))
        # plt.savefig(save_path.with_suffix(".png"), dpi=300)


def position_scatterplot_hc(df, hue_order, x, y, hue, naming_dict, save_path=None):
    """
    Plotting function generating scatterplot with function scores over GRCh38 locations high-stringency subset.

    Parameters
    ----------
    df (df): Dataframe to be plotted
    hue_order: Scatterplot parameter.
    x: Scatterplot parameter.
    y: Scatterplot parameter.
    hue: Scatterplot parameter
    naming_dict (dict): Dictionary for renaming hue annotations.
    save_path (path): Output path.
    """
    sns.set_palette(sns.color_palette(consequence_colors))

    df["Consequence"].replace(to_replace=naming_dict, inplace=True)

    df_not_hc = df.copy().loc[df["hc"] != True]
    df_hc = df.copy().loc[df["hc"] == True]

    variants = df_hc["log2_mean_variant_replicates_normalised"].count()
    print(f"number of hc variants plotted is {variants}")

    fig, ax = plt.subplots(figsize=(5, 1.6))
    sns.scatterplot(
        data=df_not_hc,
        x=x,
        y=y,
        color="#D3D3D3",
        ax=ax,
        s=10,
        linewidth=0.2,
        edgecolor="black",
        alpha=0.7,
    )
    sns.scatterplot(
        data=df_hc,
        x=x,
        y=y,
        hue=hue,
        hue_order=hue_order,
        ax=ax,
        s=10,
        linewidth=0.2,
        edgecolor="black",
        alpha=0.7,
    )
    ax.spines[["right", "top"]].set_visible(False)
    sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1), frameon=False)

    ax.set_xlabel("GRCh38 position")
    ax.set_ylabel("Function score")

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
    Main function to plot Figure 3 c and Supplementary Figure 9 a for PE manuscript. Scatterplot: functions scores (variant scores) over positions (GRCh38) for both relaxed and high-stringency datasets. MLH1 saturation screen.
    """
    save_dir = OUTPUT_DIR
    save_dir.mkdir(exist_ok=True, parents=True)

    set_style(context="paper", font_scale=1, style="ticks")

    df_variant_score = pd.read_csv(df_variant_score_path, delimiter=",")

    # Replicate analysis
    filter_number = replicate_number - 1
    df_variant_score_filtered = df_variant_score.loc[
        df_variant_score["number_replicates_pegRNA"] > filter_number
    ]
    df_variant_score_filtered = df_variant_score_filtered.drop_duplicates(
        subset=["var_key"]
    )  # important! this dataframe now contains variant level information!

    # Determine fdr
    df_variant_score_filtered = df_variant_score_filtered.pipe(determine_fdr)

    # Plotting
    position_scatterplot(
        df_variant_score_filtered.dropna(subset="Consequence"),
        x="GRCh38startlocation",
        y="log2_mean_variant_replicates_normalised",
        hue="Consequence",
        hue_order=consequence_order,
        naming_dict=consequence_naming_dict,
        save_path=save_dir / (f"Fig_3c_plot_position_scatter_function_score.svg"),
    )

    position_scatterplot_hc(
        df_variant_score_filtered.dropna(subset="Consequence"),
        x="GRCh38startlocation",
        y="log2_mean_variant_replicates_normalised",
        hue="Consequence",
        hue_order=consequence_order,
        naming_dict=consequence_naming_dict,
        save_path=save_dir
        / (f"Supp_fig9a_plot_position_scatter_function_score_hc.svg"),
    )


# ----------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    main()
