"""
    Project :  Prime editing pilot screens. Plotting script to generate Figure 4 e. Scatterplot: function scores (variant scores) over positions (GRCh38) for intron 15 region. MLH1 non-coding screen.
    Date : 240316
    Python version 3.10

"""


# IMPORT

import pandas as pd
import pathlib as pl
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.offsetbox import AnchoredText
from scipy import stats
from scipy.stats import norm

# VARIABLES

pseudo_pre_freq_filter = 1.0e-4
ST_editing_filter_PEmax = 5
ST_editing_filter_PEmaxdn = 25

# Number of replicates in which a variant has to be observed
replicate_number = 2

# PATHS

INPUT_DIR = pl.Path(
    f"/Users/kajbac/code//240116_final_scripts_for_PE_paper/MLH1_non_coding/variant_score/pseudo_pre_freq_{pseudo_pre_freq_filter}/{ST_editing_filter_PEmax}_PEmax_{ST_editing_filter_PEmaxdn}_PEmaxdn"
)

df_non_coding_variant_score_path = (
    INPUT_DIR / f"data_MLH1intronic_variant_score_replicates.csv"
)

OUTPUT_DIR = pl.Path(
    "/Users/kajbac/code/240116_final_scripts_for_PE_paper/plotting/_plots"
)
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)

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

clinvar_naming_dict = {
    "Uncertain": "Uncertain / Conflicting",
    "Conflicting": "Uncertain / Conflicting",
    "Benign": "Benign / Likely benign",
    "Likely benign": "Benign / Likely benign",
    "Pathogenic": "Pathogenic / Likely pathogenic",
    "Likely pathogenic": "Pathogenic / Likely pathogenic",
}

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
clinvar_colors = ["#0433FF", "#882255", "#DCA237"]

# hue order
consequence_order = [
    "canonical splice",
    "splice site",
    "intronic",
    "5' UTR",
    "3' UTR",
    "upstream",
    "downstream",
]

hue_clin_sig_order = [
    "Benign / Likely benign",
    "Pathogenic / Likely pathogenic",
    "Uncertain / Conflicting",
]

# --------------------------------------------------------------------------------------------------------------------
# FUNCTIONS


def determine_fdr(df):
    """
    FDR testing of pegRNA-derived function scores with the following thresholds: 0.01, 0.02, 0.05 using intronic variants as controls.

    Parameters
    ----------
    df (df): Dataframe to be processed.
    """
    loc, scale = stats.norm.fit(
        df.loc[df["Consequence"] == "INTRONIC"][
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


# plotting


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


def position_scatterplot_non_coding(
    df,
    hue_order,
    x,
    y,
    xlim_downstr,
    xlim_upstr,
    hue,
    naming_dict,
    colors,
    save_path=None,
):
    """
    Plotting function generating scatterplot with function scores over GRCh38 locations, marking LoF variants as stars and exons as grey bars.

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
    colors (list): List of hex color codes.
    save_path (path): Output path.
    """
    sns.set_palette(sns.color_palette(colors))

    df[hue].replace(to_replace=naming_dict, inplace=True)

    df_null = df.copy().loc[df["fdr_0.01"] != True]
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
    ax.set_xlabel("GRCh38 position")
    ax.set_ylabel("Function score")

    if save_path is None:
        plt.show()
    else:
        plt.savefig(save_path)
        # plt.savefig(save_path.with_suffix(".pdf"))
        # plt.savefig(save_path.with_suffix(".png"), dpi=300)


# ----------------------------------------------------------------------------------------------------------------------
# MAIN


def main():
    """
    Main function to plot Figure 4 e. Scatterplot: function scores (variant-level scores) over positions (GRCh38) for intron 15 region, with both "consequence" and "ClinVar" color-coding. MLH1 non-coding screen.
    """
    save_dir = OUTPUT_DIR
    save_dir.mkdir(exist_ok=True, parents=True)
    save_path_consequence = save_dir / (
        f"Fig_4e_plot_position_scatter_function_score_i15_consequence.svg"
    )
    save_path_clinvar = save_dir / (
        f"Fig_4e_plot_position_scatter_function_score_i15_clinvar.svg"
    )
    set_style(context="paper", font_scale=1, style="ticks")

    df_variant_score = pd.read_csv(
        df_non_coding_variant_score_path, delimiter=",", low_memory=False
    )

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

    position_scatterplot_non_coding(
        df_variant_score_filtered.dropna(subset="Consequence"),
        x="GRCh38startlocation",
        y="log2_mean_variant_replicates_normalised",
        hue="Consequence",
        hue_order=consequence_order,
        xlim_upstr=37047095,
        xlim_downstr=37047645,
        naming_dict=consequence_naming_dict,
        colors=consequence_colors,
        save_path=save_path_consequence,
    )
    position_scatterplot_non_coding(
        df_variant_score_filtered.dropna(subset="ClinVar_Significance_Simple"),
        x="GRCh38startlocation",
        y="log2_mean_variant_replicates_normalised",
        hue="ClinVar_Significance_Simple",
        hue_order=hue_clin_sig_order,
        xlim_upstr=37047095,
        xlim_downstr=37047645,
        naming_dict=clinvar_naming_dict,
        colors=clinvar_colors,
        save_path=save_path_clinvar,
    )


# ----------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    main()
