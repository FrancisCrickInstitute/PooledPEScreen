"""
    Project :  Prime editing pilot screens. Plotting script to generate Supplementary Figure 8 a function score strip plots color-coded by variant consequence. MLH1 saturation screen.
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
    {"pegRNA_D20": "CK003", "pegRNA_D34": "CK007"},
    {"pegRNA_D20": "CK005", "pegRNA_D34": "CK009"},
    {"pegRNA_D20": "CK004", "pegRNA_D34": "CK008"},
    {"pegRNA_D20": "CK006", "pegRNA_D34": "CK010"},
]

PEmax_list = ["CK003", "CK005", "CK007", "CK009"]
PEmaxdn_list = ["CK004", "CK006", "CK008", "CK010"]
D20_list = ["CK003", "CK005", "CK004", "CK006"]

pseudo_pre_freq_filter = 1.4e-4

# PATHS
OUTPUT_DIR = pl.Path(
    "/Users/kajbac/code/240116_final_scripts_for_PE_paper/plotting/_plots"
)
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)

# Order for plotting
sort_dict_consequence = {
    "synonymous": 0,
    "nonsense": 1,
    "canonical splice": 2,
    "missense": 3,
    "splice site": 4,
    "intronic": 5,
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


def stripplot(
    df,
    df_pegRNA,
    x,
    y,
    hue,
    hue_order,
    ax,
    label,
    legend=False,
):
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
    Main function to plot Supplementary Figure 8 a for PE manuscript. Stripplot: functions scores (variant scores) color-coded by variant consequence. MLH1 saturation screen.
    """
    save_dir = OUTPUT_DIR
    save_dir.mkdir(exist_ok=True, parents=True)
    set_style(context="paper", font_scale=1, style="ticks")
    save_path = save_dir / f"Supp_fig8a_plot_function_score_pre_post_filter_strip.svg"

    fig, axes = plt.subplots(4, 2, sharex=True, sharey=True, figsize=(5, 8))

    ST_editing_filter_PEmax = 0
    ST_editing_filter_PEmaxdn = 0

    df_variant_score_path = pl.Path(
        f"/Users/kajbac/code/240116_final_scripts_for_PE_paper/MLH1x10/pegRNA/variant_score/pseudo_pre_freq_{pseudo_pre_freq_filter}/{ST_editing_filter_PEmax}_PEmax_{ST_editing_filter_PEmaxdn}_PEmaxdn/data_MLH1x10_variant_score_replicates.csv"
    )

    df_variant_score = pd.read_csv(
        df_variant_score_path, delimiter=",", low_memory=False
    )

    df_variant_score_single = df_variant_score.copy()

    for i, var_dict in enumerate(variables_per_sample):
        # plotting
        stripplot(
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

    ST_editing_filter_PEmax = 5
    ST_editing_filter_PEmaxdn = 25

    df_variant_score_path = pl.Path(
        f"/Users/kajbac/code/240116_final_scripts_for_PE_paper/MLH1x10/pegRNA/variant_score/pseudo_pre_freq_{pseudo_pre_freq_filter}/{ST_editing_filter_PEmax}_PEmax_{ST_editing_filter_PEmaxdn}_PEmaxdn/data_MLH1x10_variant_score_replicates.csv"
    )

    df_variant_score = pd.read_csv(
        df_variant_score_path, delimiter=",", low_memory=False
    )

    df_variant_score_single = df_variant_score.copy()

    for i, var_dict in enumerate(variables_per_sample):
        # plotting
        stripplot(
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
