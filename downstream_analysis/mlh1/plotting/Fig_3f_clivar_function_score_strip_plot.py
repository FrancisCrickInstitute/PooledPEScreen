"""
    Project :  Prime editing pilot screens. Plotting script to generate Figure 3 f and Supplementary Figure 9 d. Stripplot of function scores (variant-level scores) by ClinVar pathogenicity annotation. MLH1 saturation screen.
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

pseudo_pre_freq_filter = 1.4e-4
ST_editing_filter_PEmax = 5
ST_editing_filter_PEmaxdn = 25

# Number of replicates in which a variant has to be observed
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

# Order for plotting
sort_dict_clinvar = {
    "Benign / Likely benign": 0,
    "Pathogenic / Likely pathogenic": 1,
    "Uncertain / Conflicting": 2,
}

# Renaming dictionary
clinvar_naming_dict = {
    "Uncertain": "Uncertain / Conflicting",
    "Conflicting": "Uncertain / Conflicting",
    "Benign": "Benign / Likely benign",
    "Likely benign": "Benign / Likely benign",
    "Pathogenic": "Pathogenic / Likely pathogenic",
    "Likely pathogenic": "Pathogenic / Likely pathogenic",
}

# Color palette
clinvar_colors = ["#0433FF", "#882255", "#DCA237"]

# Hue order
hue_clin_sig_order = [
    "Benign / Likely benign",
    "Pathogenic / Likely pathogenic",
    "Uncertain / Conflicting",
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
    x,
    y,
    hue,
    hue_order,
    ax,
    save_path=None,
    legend=False,
):
    """
    Plotting function generating stripplot of function scores grouped by ClinVar pathogenicity annotation.

    Parameters
    ----------
    df (df): Dataframe to be plotted.
    x: Stripplot parameter.
    y: Stripplot parameter.
    hue: Stripplot parameter.
    hue_order: Stripplot parameter.
    ax: Stripplot parameter.
    save_path (path): Output path.
    """
    sns.set_palette(sns.color_palette(clinvar_colors))

    df["ClinVar_Significance_Simple"].replace(
        to_replace=clinvar_naming_dict, inplace=True
    )
    df = df.sort_values(
        by=["ClinVar_Significance_Simple"], key=lambda x: x.map(sort_dict_clinvar)
    )

    variants = df["log2_mean_variant_replicates_normalised"].count()
    print(f"number of variants plotted is {variants}")

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
        alpha=0.8,
        linewidth=0.2,
        s=4,
        edgecolor="black",
    )
    sns.boxplot(
        data=df, x=x, y=y, ax=ax, dodge=False, showfliers=False, **boxplot_kwargs
    )
    ax.spines[["right", "top"]].set_visible(False)
    ax.set_ylabel(f"Function score")
    if legend == True:
        sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1), frameon=False)

    ax.tick_params(axis="x", rotation=25)

    if save_path is None:
        plt.show()
    else:
        plt.savefig(save_path)
        # plt.savefig(save_path.with_suffix(".pdf"))
        # plt.savefig(save_path.with_suffix(".png"), dpi=300)

    plt.tight_layout()


# --------------------------------------------------------------------------------------------------------------------
# MAIN


def main():
    """
    Main function to plot Figure 3 f and Supplemantary Figure 9 d for PE manuscript. Stripplot: functions scores (variant scores) by ClinVar pathogenicity annotation. MLH1 saturation screen.
    """
    save_dir = OUTPUT_DIR
    save_dir.mkdir(exist_ok=True, parents=True)
    save_path = save_dir / (f"Fig_3f_clinvar_strip_plot.svg")
    save_path_hc = save_dir / (f"Supp_fig9d_clinvar_strip_plot_hc.svg")

    set_style(context="paper", font_scale=1, style="ticks")

    df_variant_score = pd.read_csv(
        df_variant_score_path, delimiter=",", low_memory=False
    )

    # Replicate analysis
    filter_number = replicate_number - 1
    df_variant_score_filtered = df_variant_score.loc[
        df_variant_score["number_replicates_pegRNA"] > filter_number
    ]
    df_variant_score_filtered = df_variant_score_filtered.drop_duplicates(
        subset=["var_key"]
    )
    # important! this dataframe now contains variant level information!
    df_variant_score_filtered = df_variant_score_filtered.dropna(
        subset=["ClinVar_Significance_Simple"]
    )

    # plotting
    fig, ax = plt.subplots(figsize=(2, 2))

    stripplot(
        df_variant_score_filtered.dropna(subset="ClinVar_Significance_Simple"),
        x="ClinVar_Significance_Simple",
        y="log2_mean_variant_replicates_normalised",
        hue="ClinVar_Significance_Simple",
        ax=ax,
        hue_order=hue_clin_sig_order,
        legend=True,
        save_path=save_path,
    )

    fig, ax = plt.subplots(figsize=(2, 2))

    stripplot(
        df_variant_score_filtered.dropna(subset="ClinVar_Significance_Simple").loc[
            df_variant_score_filtered["hc"] == True
        ],
        x="ClinVar_Significance_Simple",
        y="log2_mean_variant_replicates_normalised",
        hue="ClinVar_Significance_Simple",
        ax=ax,
        hue_order=hue_clin_sig_order,
        legend=True,
        save_path=save_path_hc,
    )


# ----------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    main()
