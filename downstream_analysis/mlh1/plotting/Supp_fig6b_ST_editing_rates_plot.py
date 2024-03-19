"""
    Project :  Prime editing pilot screens. Plotting script to generate Supplementary Figure 6 b. ST editing rates across the 4 conditions.
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

# VARIABLES
variables_per_sample = [
    {"pegRNA_D20": "CK003", "pegRNA_D34": "CK007"},
    {"pegRNA_D20": "CK004", "pegRNA_D34": "CK008"},
    {"pegRNA_D20": "CK005", "pegRNA_D34": "CK009"},
    {"pegRNA_D20": "CK006", "pegRNA_D34": "CK010"},
]

PEmax_list = ["CK003", "CK005", "CK007", "CK009"]
PEmaxdn_list = ["CK004", "CK006", "CK008", "CK010"]
D20_list = ["CK003", "CK005", "CK004", "CK006"]

pseudo_pre_freq_filter = 1.4e-4
ST_editing_filter_PEmax = 0
ST_editing_filter_PEmaxdn = 0

INPUT_DIR = pl.Path(
    f"/Users/kajbac/code/240116_final_scripts_for_PE_paper/MLH1x10/pegRNA/pegRNA_score/pseudo_pre_freq_{pseudo_pre_freq_filter}/{ST_editing_filter_PEmax}_PEmax_{ST_editing_filter_PEmaxdn}_PEmaxdn"
)

df_pegRNA_path = INPUT_DIR / "data_MLH1x10_pegRNA_score_replicates.csv"

OUTPUT_DIR = pl.Path(
    "/Users/kajbac/code/240116_final_scripts_for_PE_paper/plotting/_plots"
)
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)
# --------------------------------------------------------------------------------------------------------------------
# FUNCTIONS
# --------------------------------------------------------------------------------------------------------------------
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


def ST_editing_violin_plot(df, x, y, hue, save_path=None):
    """
    Plotting function generating violinplots of ST-editing rates for all conditions.

    Parameters
    ----------
    df (df): Dataframe to be plotted.
    save_path (path): Output path.
    """
    fig, ax = plt.subplots(figsize=(2, 2))
    color_palette = ["#1e76b4", "#ff7f0e"]
    sns.violinplot(data=df, x=x, y=y, inner="quart", hue=hue, palette=color_palette)

    for l in ax.lines:
        l.set_linestyle("-")
        l.set_linewidth(0.6)
        l.set_color("black")
        l.set_alpha(0.8)
    for l in ax.lines[1::3]:
        l.set_linestyle("-")
        l.set_linewidth(1.2)
        l.set_color("black")
        l.set_alpha(0.8)

    ax.set_ylim(0, 100)
    # ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    ax.set_ylabel(f"% Correct ST Edits")

    ax.spines[["right", "top"]].set_visible(False)
    ax.set(xlabel=None)
    # ax.tick_params(axis="x", rotation=25)
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
    Main function to plot Supplementary Figure 6 b for PE manuscript. Violinplot: Surrogate target editing rates across 4 conditions. MLH1 saturation screen.
    """
    save_dir = OUTPUT_DIR
    save_dir.mkdir(exist_ok=True, parents=True)
    save_path = save_dir / f"Supp_fig6b_plot_ST_editing_rates_violin.svg"
    df = pd.read_csv(df_pegRNA_path, delimiter=",")

    column_list = ["var_key"]
    for condition in D20_list:
        column_list.append(f"{condition}_percentage_editing")
    df_ST_editing = df[column_list]
    df_ST_editing.rename(
        columns={
            "CK003_percentage_editing": "PEmax",
            "CK005_percentage_editing": "PEmax + obn",
            "CK004_percentage_editing": "PEmax-MLH1dn",
            "CK006_percentage_editing": "PEmax-MLH1dn + obn",
        },
        inplace=True,
    )

    df_ST_editing = df_ST_editing.melt(
        id_vars="var_key", var_name="condition", value_name="ST_editing_rate"
    )
    df_ST_editing["selection"] = df_ST_editing["condition"].map(
        lambda x: True if x in ["PEmax + obn", "PEmax-MLH1dn + obn"] else False
    )

    set_style(context="paper", font_scale=1, style="ticks")

    ST_editing_violin_plot(
        df=df_ST_editing,
        x="condition",
        y="ST_editing_rate",
        hue="selection",
        save_path=save_path,
    )


# ----------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    main()
