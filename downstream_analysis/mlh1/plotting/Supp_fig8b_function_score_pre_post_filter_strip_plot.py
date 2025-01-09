"""
    Project :  Project :  Prime editing pilot screens. Plotting script to generate Supplementary Figure 8 b function score strip plots color-coded by variant consequence. MLH1 non-coding screen.
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
variables_per_sample = [
    {"pegRNA_D20": "CK015", "pegRNA_D34": "CK019"},
    {"pegRNA_D20": "CK017", "pegRNA_D34": "CK021"},
    {"pegRNA_D20": "CK016", "pegRNA_D34": "CK020"},
    {"pegRNA_D20": "CK018", "pegRNA_D34": "CK022"},
]

PEmax_list = ["CK015", "CK019", "CK017", "CK021"]
PEmaxdn_list = ["CK016", "CK020", "CK018", "CK022"]
D20_list = ["CK015", "CK017", "CK016", "CK018"]

pseudo_pre_freq_filter = 1e-4

# PATHS
OUTPUT_DIR = pl.Path(
    "/Users/kajbac/code/240116_final_scripts_for_PE_paper/plotting/_plots"
)
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)

# order for plotting
sort_dict_consequence = {
    "intronic": 0,
    "canonical splice": 1,
    "splice site": 2,
    "5' UTR": 3,
    "3' UTR": 4,
    "upstream": 5,
    "downstream": 6,
}

# renaming dictionary
consequence_naming_dict = {
    "CANONICAL_SPLICE": "canonical splice",
    "SPLICE_SITE": "splice site",
    "INTRONIC": "intronic",
    "5PRIME_UTR": "5' UTR",
    "3PRIME_UTR": "3' UTR",
    "UPSTREAM": "upstream",
    "DOWNSTREAM": "downstream",
}

# color palette
consequence_colors = [
    "#f7a83e",
    "#243672",
    "#8acdef",
    "#B5A695",
    "#d091bf",
    "#784421",
    "#964594",
]

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
# --------------------------------------------------------------------------------------------------------------------
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
    sns.set(style=style, context=context, font_scale=font_scale)


def stripplot_pegRNA(df, df_pegRNA, x, y, hue, hue_order, ax, label, legend=False):
    """
    Plotting function generating stripplots of function scores before and after surrogate target editing filter.

    Parameters
    ----------
    df (df): Dataframe to be plotted.
    df_pegRNA (df): Dataframe filtered on pegRNA scores.
    x: Stripplot parameter.
    y: Stripplot parameter.
    hue: Stripplot parameter.
    hue_order: Stripplot parameter.
    ax: Stripplot parameter.
    label (str): Set as y-label.
    legend (boolean): Decides whether legend is added.
    """
    sns.set_palette(sns.color_palette(consequence_colors))

    df["Consequence"].replace(to_replace=consequence_naming_dict, inplace=True)
    df = df.sort_values(by=["Consequence"], key=lambda x: x.map(sort_dict_consequence))

    variant_count = df[y].count()
    pegRNA_count = df_pegRNA[y].count()

    boxprops = {"edgecolor": "k", "linewidth": 1, "facecolor": "None"}
    lineprops = {"color": "k", "linewidth": 1}
    medianprops = {"linewidth": 1, "solid_capstyle": "butt", "color": "k"}

    boxplot_kwargs = dict(
        {
            "boxprops": boxprops,
            "medianprops": medianprops,
            "whiskerprops": lineprops,
            "capprops": lineprops,
            "width": 0.6,
        },
    )

    sns.stripplot(
        data=df,
        x=x,
        y=y,
        hue=hue,
        hue_order=hue_order,
        legend=legend,
        ax=ax,
        zorder=0,
        alpha=0.7,
        linewidth=0.2,
        s=3,
        edgecolor="black",
    )
    sns.boxplot(
        data=df, x=x, y=y, ax=ax, dodge=False, showfliers=False, **boxplot_kwargs
    )
    anc = AnchoredText(
        f"variants = {variant_count}\npegRNAs = {pegRNA_count}",
        loc="upper left",
        frameon=False,
    )
    ax.add_artist(anc)
    ax.spines[["right", "top"]].set_visible(False)
    ax.set_ylabel(label)
    if legend == True:
        sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1), frameon=False)

    ax.tick_params(axis="x", rotation=25)


# --------------------------------------------------------------------------------------------------------------------
# MAIN


def main():
    """
    Main function to plot Supplementary Figure 8 b for PE manuscript. Stripplot: function scores (variant scores) color-coded by variant consequence. MLH1 non-coding screen.
    """
    save_dir = OUTPUT_DIR
    save_dir.mkdir(exist_ok=True, parents=True)
    set_style(context="paper", font_scale=1, style="ticks")
    save_path = save_dir / "Supp_fig8b_plot_function_score_pre_post_filter_strip.svg"

    fig, axes = plt.subplots(4, 2, sharex=True, sharey=True, figsize=(5, 8))

    df_variant_score_path = pl.Path(
        f"/Users/kajbac/code/240116_final_scripts_for_PE_paper/MLH1_non_coding/variant_score/pseudo_pre_freq_0.0001/0_PEmax_0_PEmaxdn/data_MLH1intronic_variant_score_replicates.csv"
    )

    df_variant_score = pd.read_csv(
        df_variant_score_path, delimiter=",", low_memory=False
    )

    df_variant_score_single = df_variant_score.copy()

    for i, var_dict in enumerate(variables_per_sample):
        # plotting
        stripplot_pegRNA(
            df=df_variant_score_single.dropna(
                subset=[f"{var_dict['pegRNA_D20']}_{var_dict['pegRNA_D34']}_log2"]
            ).drop_duplicates(subset=["var_key"]),
            df_pegRNA=df_variant_score_single.dropna(
                subset=[f"{var_dict['pegRNA_D20']}_{var_dict['pegRNA_D34']}_log2"]
            ),
            x="Consequence",
            y=f"{var_dict['pegRNA_D20']}_{var_dict['pegRNA_D34']}_log2_mean_variant_normalised",
            hue="Consequence",
            ax=axes[i, 0],
            hue_order=consequence_order,
            label="Function score",
            legend=True,
        )

    df_variant_score_path = pl.Path(
        f"/Users/kajbac/code/240116_final_scripts_for_PE_paper/MLH1_non_coding/variant_score/pseudo_pre_freq_0.0001/5_PEmax_25_PEmaxdn/data_MLH1intronic_variant_score_replicates.csv"
    )

    df_variant_score = pd.read_csv(
        df_variant_score_path, delimiter=",", low_memory=False
    )

    df_variant_score_single = df_variant_score.copy()

    for i, var_dict in enumerate(variables_per_sample):
        # plotting
        stripplot_pegRNA(
            df=df_variant_score_single.dropna(
                subset=[f"{var_dict['pegRNA_D20']}_{var_dict['pegRNA_D34']}_log2"]
            ).drop_duplicates(subset=["var_key"]),
            df_pegRNA=df_variant_score_single.dropna(
                subset=[f"{var_dict['pegRNA_D20']}_{var_dict['pegRNA_D34']}_log2"]
            ),
            x="Consequence",
            y=f"{var_dict['pegRNA_D20']}_{var_dict['pegRNA_D34']}_log2_mean_variant_normalised",
            hue="Consequence",
            ax=axes[i, 1],
            hue_order=consequence_order,
            label="Function score",
            legend=True,
        )

        if save_path is None:
            plt.show()
        else:
            plt.savefig(save_path)
            # plt.savefig(save_path.with_suffix(".pdf"))
            # plt.savefig(save_path.with_suffix(".png"), dpi=300)


# ----------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    main()
