"""
    Project :  Prime editing pilot screens. Plotting script to generate Supplementary Figure 12 c and d. Endogenous function scores by genomic positions and correlation of endogenous and function scores for HAP1_PEmax condition.
    Date : 240316
    Python version 3.10

"""
# IMPORT

import pandas as pd
import pathlib as pl
import math
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.offsetbox import AnchoredText
from scipy import stats
from scipy.stats import norm

# PATHS

INPUT_DIR = pl.Path(
    f"/Users/kajbac/code/240116_final_scripts_for_PE_paper/MLH1_non_coding/endogenous/log2_enriched"
)

OUTPUT_DIR = pl.Path(
    "/Users/kajbac/code/240116_final_scripts_for_PE_paper/plotting/_plots"
)
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)

u1_path = INPUT_DIR / "u1_data_full.csv"
u2_path = INPUT_DIR / "u2_data_full.csv"
x11_path = INPUT_DIR / "x11_data_full.csv"
i15_path = INPUT_DIR / "i15_data_full.csv"

csv_paths_u = [u1_path, u2_path]
csv_paths_x11_i15 = [x11_path, i15_path]

pegRNA_score = "CK015_CK019_log2_mean_variant"
df_label_list = ["x11", "i15"]

# Color palette
consequence_colors = [
    "#f7a83e",
    "#243672",
    "#8acdef",
    "#B5A695",
    "#d091bf",
    "#784421",
    "#964594",
]
all_SNV_colors = ["#243672"]

# Hue order
consequence_order = [
    "canonical splice",
    "splice site",
    "intronic",
    "5' UTR",
    "3' UTR",
    "upstream",
    "downstream",
]

# Renaming dictionary
consequence_naming_dict = {
    "CANONICAL_SPLICE": "canonical splice",
    "SPLICE_SITE": "splice site",
    "INTRONIC": "intronic",
    "5PRIME_UTR": "5' UTR",
    "3PRIME_UTR": "3' UTR",
    "UPSTREAM": "upstream",
    "DOWNSTREAM": "downstream",
}

# ----------------------------------------------------------------------------------------------------------------------
# Functions


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
    sns.set(style=style, context=context, font_scale=1.2)


def position_scatterplot_non_coding(
    df, hue_order, x, y, xlim_downstr, xlim_upstr, hue, naming_dict, save_path=None
):
    """
    Plotting function generating scatterplot with function scores over GRCh38 locations, marking exons as grey bars.

    Parameters
    ----------
    df (df): Dataframe to be plotted
    hue_order: Scatterplot parameter.
    x: Scatterplot parameter.
    y: Scatterplot parameter.
    xlim_downstr (int): x-axis limit to indicate downtream boundary of genomic region.
    xlim_upstream (int): x-axis limit to indicate uptream boundary of genomic region.
    hue: Scatterplot parameter
    naming_dict (dict): Dictionary for renaming hue annotations.
    save_path (path): Output path.
    """
    sns.set_palette(sns.color_palette(consequence_colors))

    df[hue].replace(to_replace=naming_dict, inplace=True)

    df_null = df.copy().loc[df["fdr_0.01"] == False]
    df_lof = df.copy().loc[df["fdr_0.01"] == True]

    exon_windows = [
        (36993548, 36993663),
        (36996619, 36996709),
        (37000955, 37001053),
        (37004401, 37004474),
        (37006991, 37007063),
        (37008814, 37008905),
        (37011820, 37011862),
        (37012011, 37012099),
        (37014432, 37014544),
        (37017506, 37017599),
        (37020310, 37020463),
        (37025637, 37026007),
        (37028784, 37028932),
        (37040186, 37040294),
        (37042268, 37042331),
        (37047519, 37047683),
        (37048517, 37048609),
        (37048904, 37049017),
        (37050486, 37050653),
    ]

    fig, ax = plt.subplots(figsize=(2, 1))
    for a, b in exon_windows:
        plt.axvspan(a, b, facecolor="#D3D3D3", edgecolor=None, linewidth=0, alpha=1)
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

    ax.set_xlim(xlim_upstr, xlim_downstr)
    ax.set_ylim(-6, 12)
    ax.set_yticks([-5, 0, 5, 10])
    ax.set_xlabel("GRCh38 position")
    ax.set_ylabel("Endogenous function score")

    fig.tight_layout()

    if save_path is None:
        plt.show()
    else:
        plt.savefig(save_path)
        # plt.savefig(save_path.with_suffix(".pdf"))
        # plt.savefig(save_path.with_suffix(".png"), dpi=300)


def pegRNA_variant_score_correlation_scatter(
    df, hue, hue_order, naming_dict, save_path
):
    """
    Plotting function generating scatterplot, correlating endogenous function scores against function scores.

    Parameters
    ----------
    df (df): Dataframe to be plotted.
    hue: Scatterplot parameter.
    hue_order: Scatterplot parameter.
    naming_dict (dict): Dictionary for renaming hue annotations.
    save_path (path): Output path.
    """
    sns.set_palette(sns.color_palette(consequence_colors))

    df.rename(
        columns={
            "CK015_CK019_log2_mean_variant_normalised": "Function score",
            f"D34_D20_log2_normalised": "Endogenous function score",
        },
        inplace=True,
    )

    df = df.copy().dropna(subset="Endogenous function score")
    df = df.copy().dropna(subset="Function score")

    print(df["Endogenous function score"].count())
    print(df["Function score"].count())

    r = df["Function score"].corr(df["Endogenous function score"], method="pearson")
    r = round(r, 3)
    df[hue].replace(to_replace=naming_dict, inplace=True)

    df_no_lof = df.loc[df["fdr_0.01"] == False]
    df_lof = df.loc[df["fdr_0.01"] == True]

    fig, ax = plt.subplots(figsize=(2, 2))

    sns.scatterplot(
        data=df_lof,
        x="Endogenous function score",
        y="Function score",
        hue=hue,
        hue_order=hue_order,
        marker=(8, 1, 0),
        s=50,
        linewidth=0.2,
        edgecolors="black",
        alpha=0.7,
        ax=ax,
    )
    sns.scatterplot(
        data=df_no_lof,
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
        # \nno. variants = {variant_count} <----add this line to displau variant count
        loc="upper left",
        frameon=False,
    )
    ax.add_artist(anc)
    sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1), frameon=False)
    ax.set_ylim(-2, 15)

    fig.tight_layout()

    if save_path is None:
        plt.show()
    else:
        plt.savefig(save_path)
        # plt.savefig(save_path.with_suffix(".pdf"))
        # plt.savefig(save_path.with_suffix(".png"), dpi=300)


# ----------------------------------------------------------------------------------------------------------------------
# MAIN FUNCTION


def main():
    """
    Main function to plot Supplementary Figures 12 c and d for PE manuscript. Endogenous function scores by genomic positions and correlation of endogenous and function scores for HAP1_PEmax condition.
    """
    save_dir = OUTPUT_DIR
    save_dir.mkdir(exist_ok=True, parents=True)
    save_path_gene_map_all_variants_up = (
        save_dir / "Supp_fig12c_endogenous_function_score_position_scatter_up.svg"
    )
    save_path_gene_map_all_variants_x11 = (
        save_dir / "Supp_fig12c_endogenous_function_score_position_scatter_x11.svg"
    )
    save_path_gene_map_all_variants_i15 = (
        save_dir / "Supp_fig12c_endogenous_function_score_position_scatter_i15.svg"
    )
    save_path_score_correlation = (
        save_dir / "Supp_fig12d_endogenous_function_score_correlation_scatter.svg"
    )

    set_style(context="paper", font_scale=1, style="ticks")

    df_u1 = pd.read_csv(csv_paths_u[0])
    df_u2 = pd.read_csv(csv_paths_u[1])
    df_u1[["u1"]] = True
    df_u2[["u2"]] = True

    df_u = pd.concat([df_u1, df_u2], axis=0, join="outer")

    df_u["variant_scored"] = df_u[pegRNA_score].notna()
    df_u[["u"]] = True
    df_u = df_u.loc[df_u["variant_scored"] == True]
    df_u = df_u.loc[(df_u["expected_variant"] == True)]
    df_u = df_u[~(df_u["expected_SNV"] & (df_u["D20_neg_log2"] <= 1))]
    df_all_expected_variants_list = [df_u]

    df_x11 = pd.read_csv(x11_path)
    df_x11["variant_scored"] = df_x11[pegRNA_score].notna()
    df_x11["x11"] = True
    df_x11 = df_x11.loc[(df_x11["expected_variant"] == True)]
    df_x11 = df_x11[~(df_x11["expected_SNV"] & (df_x11["D20_neg_log2"] <= 1))]
    df_x11 = df_x11.loc[df_x11["variant_scored"] == True]
    df_all_expected_variants_list.append(df_x11)

    df_i15 = pd.read_csv(i15_path)
    df_i15["variant_scored"] = df_i15[pegRNA_score].notna()
    df_i15["i15"] = True
    df_i15 = df_i15.loc[(df_i15["expected_variant"] == True)]
    df_i15 = df_i15[~(df_i15["expected_SNV"] & (df_i15["D20_neg_log2"] <= 1))]
    df_i15 = df_i15.loc[df_i15["variant_scored"] == True]
    df_all_expected_variants_list.append(df_i15)

    df_concat = pd.concat([df_u, df_x11, df_i15], axis=0, join="outer")

    # normalising endogenous function scores to median score of splice neutral intronic variants
    df_concat["D34_D20_log2_normalised"] = (
        df_concat["D34_D20_log2"]
        - df_concat.loc[
            (df_concat["Consequence"] == "INTRONIC")
            & (df_concat["mutant index"] != 3128)
        ]["D34_D20_log2"].median()
    )

    df_u_normalised = df_concat.copy().loc[df_concat["u"] == True]
    df_x11_normalised = df_concat.loc[df_concat["x11"] == True]
    df_i15_normalised = df_concat.loc[df_concat["i15"] == True]

    position_scatterplot_non_coding(
        df_u_normalised,
        x="GRCh38startlocation",
        y="D34_D20_log2_normalised",
        hue="Consequence",
        hue_order=consequence_order,
        xlim_upstr=36993150,
        xlim_downstr=36993700,
        naming_dict=consequence_naming_dict,
        save_path=save_path_gene_map_all_variants_up,
    )

    position_scatterplot_non_coding(
        df_x11_normalised,
        x="GRCh38startlocation",
        y="D34_D20_log2_normalised",
        hue="Consequence",
        hue_order=consequence_order,
        xlim_upstr=37020100,
        xlim_downstr=37020650,
        naming_dict=consequence_naming_dict,
        save_path=save_path_gene_map_all_variants_x11,
    )

    position_scatterplot_non_coding(
        df_i15_normalised,
        x="GRCh38startlocation",
        y="D34_D20_log2_normalised",
        hue="Consequence",
        hue_order=consequence_order,
        xlim_upstr=37047095,
        xlim_downstr=37047645,
        naming_dict=consequence_naming_dict,
        save_path=save_path_gene_map_all_variants_i15,
    )

    pegRNA_variant_score_correlation_scatter(
        df_concat,
        hue="Consequence",
        hue_order=consequence_order,
        naming_dict=consequence_naming_dict,
        save_path=save_path_score_correlation,
    )


# ----------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    main()
