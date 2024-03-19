"""
    Project :  Prime editing pilot screens. Plotting script to generate Supplementary Figure 6 b. Endogenous variant frequencies across the 4 conditions. MLH1 saturation screen.
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
import matplotlib.ticker as ticker
from scipy import stats
from scipy.stats import norm

# PATHS

INPUT_DIR_ENDOGENOUS = pl.Path(
    f"/Users/kajbac/code/240116_final_scripts_for_PE_paper/MLH1x10/endogenous/log2_enriched"
)

df_endogenous_score_path = (
    INPUT_DIR_ENDOGENOUS / "data_expected_variants_replicates_log2_filtered.csv"
)

OUTPUT_DIR = pl.Path(
    "/Users/kajbac/code/240116_final_scripts_for_PE_paper/plotting/_plots"
)
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)

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


def variant_frequencies_across_replicates(df, save_path):
    """
    Plotting function generating violinplots of endogenous variant frequencies for all conditions.

    Parameters
    ----------
    df (df): Dataframe to be plotted.
    save_path (path): Output path.
    """
    figure, (ax, ax2) = plt.subplots(2, 1, sharex=True, figsize=(2, 2))

    df["variant_frequency"] = df["variant_frequency"] * 1000
    df["selection"] = df["condition"].map(
        lambda x: True if x in ["PEmax + obn", "PEmax-MLH1dn + obn"] else False
    )

    color_palette = ["#1e76b4", "#ff7f0e"]
    sns.violinplot(
        data=df,
        x="condition",
        y="variant_frequency",
        inner="quart",
        hue="selection",
        palette=color_palette,
        ax=ax,
    )

    sns.violinplot(
        data=df,
        x="condition",
        y="variant_frequency",
        inner="quart",
        hue="selection",
        palette=color_palette,
        ax=ax2,
    )

    for l in ax.lines:
        l.set_linestyle("-")
        l.set_linewidth(0.6)
        l.set_color("black")
        l.set_alpha(1)
    for l in ax.lines[1::3]:
        l.set_linestyle("-")
        l.set_linewidth(1.2)
        l.set_color("black")
        l.set_alpha(1)

    for l in ax2.lines:
        l.set_linestyle("-")
        l.set_linewidth(0.6)
        l.set_color("black")
        l.set_alpha(1)
    for l in ax2.lines[1::3]:
        l.set_linestyle("-")
        l.set_linewidth(1.2)
        l.set_color("black")
        l.set_alpha(1)

    ax.set_ylim(3.5, 10)  # outliers only
    ax2.set_ylim(0, 2)

    ax.spines[["bottom", "top", "right"]].set_visible(False)
    ax2.spines[["top", "right"]].set_visible(False)
    ax.tick_params(
        axis="x",  # changes apply to the x-axis
        which="both",  # both major and minor ticks are affected
        bottom=False,  # ticks along the bottom edge are off
        top=False,  # ticks along the top edge are off
        labelbottom=False,
    )
    ax2.xaxis.tick_bottom()
    ax2.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    ax.set_ylabel([])
    ax2.set_ylabel(f"Variant Frequency ($10^-3$)")

    d = 0.015

    kwargs = dict(transform=ax.transAxes, color="k", clip_on=False)
    ax.plot((-d, +d), (-d, +d), **kwargs)
    ax.plot((1 - d, 1 + d), (-d, +d), **kwargs)

    kwargs.update(transform=ax2.transAxes)
    ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)
    ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)

    ax.get_legend().remove()
    ax2.get_legend().remove()

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
    Main function to plot Supplementary Figure 6 b and for PE manuscript. Violinplot: Endogenous variant frequencies across 4 conditions. MLH1 saturation screen.
    """
    save_dir = OUTPUT_DIR
    save_dir.mkdir(exist_ok=True, parents=True)
    save_path_ET = (
        save_dir
        / f"Supp_fig6b_plot_endogenous_variant_frequencies_log2_filtered_violin.svg"
    )
    set_style(context="paper", font_scale=1, style="ticks")

    df_endogenous = pd.read_csv(df_endogenous_score_path, delimiter=",")

    column_list = ["var_key"]
    condition_list = ["S1", "S1O", "S2", "S2O"]
    for condition in condition_list:
        column_list.append(f"D20_freq_{condition}")
    df_var_freq = df_endogenous[column_list]
    df_var_freq.rename(
        columns={
            "D20_freq_S1": "PEmax",
            "D20_freq_S1O": "PEmax + obn",
            "D20_freq_S2": "PEmax-MLH1dn",
            "D20_freq_S2O": "PEmax-MLH1dn + obn",
        },
        inplace=True,
    )

    df_var_freq = df_var_freq.melt(
        id_vars="var_key", var_name="condition", value_name="variant_frequency"
    )

    variant_frequencies_across_replicates(
        df_var_freq,
        save_path=save_path_ET,
    )


# ----------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    main()
