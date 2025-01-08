"""
    Project :  Prime editing pilot screens. Revisions. Plotting script to generate Supplementary Figure 10b using only data generated in the HAP1:PEmax+MLH1dn cell line. Correlation Scatterplot: functions scores (variant scores) against Cadd (Phred) scores. MLH1x10 screen.
    Date : 240918
    Python version 3.10

"""

# IMPORT

import pandas as pd
import pathlib as pl
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.offsetbox import AnchoredText
from scipy import stats
from scipy.stats import norm

# VARIABLES

pseudo_pre_freq_filter = 1.4e-4
ST_editing_filter_PEmax = 5
ST_editing_filter_PEmaxdn = 25

# PATHS

INPUT_DIR = pl.Path(
    f"/Users/kajbac/code/240913_cell_line_differences/pseudo_pre_freq_{pseudo_pre_freq_filter}/{ST_editing_filter_PEmax}_PEmax_{ST_editing_filter_PEmaxdn}_PEmaxdn"
)

OUTPUT_DIR = pl.Path(f"/Users/kajbac/code/240913_cell_line_differences/_plots")
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)

df_cell_line_path = INPUT_DIR / f"data_MLH1x10_variant_score_per_cell_line.csv"

# Order for plotting
sort_dict_consequence = {
    "SYNONYMOUS": 0,
    "STOP_GAINED": 1,
    "CANONICAL_SPLICE": 2,
    "NON_SYNONYMOUS": 3,
    "SPLICE_SITE": 4,
    "INTRONIC": 5,
}

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


def Correlation_scatter(df, x, y, hue, hue_order, ax, save_path=None):
    """
    Plotting function generating correlation scatterplot with function scores against CADD (Phred) scores.

    Parameters
    ----------
    df (df): Dataframe to be plotted.
    x: Scatterplot parameter.
    y: Scatterplot parameter.
    hue: Scatterplot parameter.
    hue_order: Scatterplot parameter.
    ax: Scatterplot parameter.
    save_path (path): Output path.
    """
    sns.set_palette(sns.color_palette(consequence_colors))

    df["Consequence"].replace(to_replace=consequence_naming_dict, inplace=True)

    variants = df["log2_mean_variant_PEmaxdn_normalised"].count()
    print(f"number of variants plotted is {variants}")

    r = df[x].corr(df[y], method="pearson")
    r = round(r, 3)

    sns.scatterplot(
        data=df,
        x=x,
        y=y,
        hue=hue,
        hue_order=hue_order,
        ax=ax,
        linewidth=0.2,
        edgecolors="black",
        alpha=0.7,
        s=10,
    )

    ax.spines[["right", "top"]].set_visible(False)
    anc = AnchoredText(
        f"r = {r}",
        loc="upper left",
        frameon=False,
    )
    ax.add_artist(anc)

    ax.set_ylabel(f"Function score")
    ax.set_xlabel(f"CADD score")
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
    Main function to plot Supplementary Figure 10b using only Function Score generated in PEmaxdn. Scatterplot: functions scores (variant scores) against CADD scores for both relaxed and high-stringency datasets. MLH1x10 screen.
    """
    save_dir = OUTPUT_DIR
    save_dir.mkdir(exist_ok=True, parents=True)

    save_path = save_dir / (
        f"Supp_Fig_10b_plot_replicates_cadd_score_correlation_PEmaxdn_only_no_stars.svg"
    )


    set_style(context="paper", font_scale=1, style="ticks")

    df_cell_line = pd.read_csv(df_cell_line_path, delimiter=",", low_memory=False)

    df_PEmax = df_cell_line.copy()[
        [
            "mutant index",
            "var_key",
            "GRCh38startlocation",
            "Consequence",
            "log2_mean_variant_PEmaxdn_normalised",
            "PEmaxdn_percentage_editing_mean_variant",
            "PHRED",
        ]
    ]

    df_PEmax = df_PEmax.drop_duplicates(subset=["var_key"])
    # important! this dataframe now contains variant level information!

    fig, ax = plt.subplots(figsize=(2, 1.5))

    Correlation_scatter(
        df=df_PEmax.dropna(subset="PHRED").loc[
            df_PEmax["Consequence"].isin(["NON_SYNONYMOUS"])
        ],
        x="PHRED",
        y="log2_mean_variant_PEmaxdn_normalised",
        hue="Consequence",
        ax=ax,
        hue_order=consequence_order,
        save_path=save_path,
    )


# ----------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    main()
