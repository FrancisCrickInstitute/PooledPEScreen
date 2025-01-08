"""
    Project :  Prime editing pilot screens. Plotting script to generate Supplementary Figure 10a using only data generated in the HAP1:PEmax cell line. Correlation Scatterplot: functions scores (variant scores) against SpliceAI scores. MLH1 non-coding screen.
    Date : 241016
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

# PATHS

pseudo_pre_freq_filter = 1.0e-4
ST_editing_filter_PEmax = 5
ST_editing_filter_PEmaxdn = 25

INPUT_DIR = pl.Path(
    f"/Users/kajbac/code/240913_cell_line_differences/pseudo_pre_freq_{pseudo_pre_freq_filter}/{ST_editing_filter_PEmax}_PEmax_{ST_editing_filter_PEmaxdn}_PEmaxdn"
)

OUTPUT_DIR = pl.Path(f"/Users/kajbac/code/240913_cell_line_differences/_plots")
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)

df_cell_line_path = INPUT_DIR / f"data_MLH1_noncoding_variant_score_per_cell_line.csv"

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

# --------------------------------------------------------------------------------------------------------------------
# FUNCTIONS

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


def SpliceAI_correlation_scatter(
    df, x, y, hue, hue_order, ax, naming_dict, save_path=None
):
    """
    Plotting function generating correlation scatterplot with function scores against SpliceAI scores.

    Parameters
    ----------
    df (df): Dataframe to be plotted.
    x: Scatterplot parameter.
    y: Scatterplot parameter.
    hue: Scatterplot parameter.
    hue_order: Scatterplot parameter.
    ax: Scatterplot parameter.
    naming_dict (dict): Dictionary for renaming hue annotations.
    save_path (path): Output path.
    """

    df[hue].replace(to_replace=naming_dict, inplace=True)
    sns.set_palette(sns.color_palette(consequence_colors))

    df = df.dropna(subset="SpliceAI_max")
    df = df.dropna(subset="log2_mean_variant_PEmax_normalised")
    variants = df["SpliceAI_max"].count()
    print(f"number of variants plotted is {variants}")

    r = df[x].corr(df[y])
    r = round(r, 3)


    sns.scatterplot(
        data=df,
        x=x,
        y=y,
        hue=hue,
        hue_order=hue_order,
        ax=ax,
        s=10,
        linewidth=0.2,
        edgecolor="black",
        alpha=0.6,
    )

    anc = AnchoredText(
        f"r = {r}",
        loc="upper left",
        frameon=False,
    )
    ax.add_artist(anc)
    ax.spines[["right", "top"]].set_visible(False)
    ax.set_xlabel("Splice AI")

    ax.set_ylabel(f"Function score")
    ax.get_legend().remove()

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
    Main function to plot Supplementary Figure 10a using Function Scores generated in PEmax. Scatterplot: function scores (variant-level scores) against SpliceAI scores. MLH1 non-coding screen.
    """
    save_dir = OUTPUT_DIR
    save_dir.mkdir(exist_ok=True, parents=True)
    save_path = save_dir / (
        f"Supp_Fig10a_plot_SpliceAI_score_correlation_PEmax_only_no_stars.svg"
    )

    set_style(context="paper", font_scale=1, style="ticks")

    df_cell_line = pd.read_csv(df_cell_line_path, delimiter=",", low_memory=False)

    df_PEmax = df_cell_line.copy()[
        [
            "mutant index",
            "var_key",
            "GRCh38startlocation",
            "Consequence",
            "log2_mean_variant_PEmax_normalised",
            "PEmax_percentage_editing_mean_variant",
            "PHRED",
            "SpliceAI_max",
        ]
    ]

    df_PEmax = df_PEmax.drop_duplicates(subset=["var_key"])
    # important! this dataframe now contains variant level information!

    # Plots
    fig, ax = plt.subplots(figsize=(2, 1.5))

    SpliceAI_correlation_scatter(
        df=df_PEmax.dropna(subset="SpliceAI_max"),
        x="SpliceAI_max",
        y="log2_mean_variant_PEmax_normalised",
        hue="Consequence",
        ax=ax,
        naming_dict=consequence_naming_dict,
        hue_order=consequence_order,
        save_path=save_path,
    )


# ----------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    main()
