"""
    Project :  Prime editing pilot screens. Plotting script to generate Supplementary Figure 6 a. 6TG dilution series in HAP1:PEmax, HAP1:PEmax+MLH1dn, and HAP1:PE2+MLH1KO.
    Date : 240316
    Python version 3.10

"""

# IMPORT

import pandas as pd
import pathlib as pl
import matplotlib.pyplot as plt
import seaborn as sns

# PATHS

INPUT_DIR = pl.Path(f"/Users/kajbac/code/240116_final_scripts_for_PE_paper")
dilution_series_path = INPUT_DIR / "MLH1_6TG_dilution_series_results.csv"

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


def dilution_plot(
    data,
    x,
    y1,
    y2,
    y3,
    ax,
    save_path,
):
    """
    Plotting function generating lineplot: Viable cell counts by 6TG dose.

    Parameters
    ----------
    data (df): Dataframe to be plotted
    x: Lineplot parameter.
    y1: Lineplot parameter.
    y2: Lineplot parameter.
    y3: Lineplot parameter.
    hue: Lineplot parameter
    ax: Lineplot parameter.
    save_path (path): Output path.
    """
    sns.lineplot(data=data, x=x, y=y1, ax=ax, color="#1e76b4", label="HAP1:PEmax")
    sns.lineplot(
        data=data, x=x, y=y2, ax=ax, color="#2ba02b", label="HAP1:PEmax+MLH1dn"
    )
    sns.lineplot(data=data, x=x, y=y3, ax=ax, color="#ff7f0e", label="HAP1:PE2+MLH1KO")
    ax.set_xlabel(f"6TG [µg / ml]")
    ax.set_ylabel(f"Viable cells ($10^6$)")
    ax.spines[["right", "top"]].set_visible(False)
    plt.legend(frameon=False)

    plt.savefig(save_path)
    # plt.savefig(save_path.with_suffix(".pdf"))
    # plt.savefig(save_path.with_suffix(".png"), dpi=300)


# --------------------------------------------------------------------------------------------------------------------
# MAIN


def main():
    """
    Main function to plot Supplementary Figure 6 a for PE manuscript. Lineplot: Viable cell counts by 6TG dosage.
    """
    set_style(context="paper", font_scale=1, style="ticks")
    df = pd.read_csv(dilution_series_path, delimiter=",")

    save_dir = OUTPUT_DIR
    save_dir.mkdir(exist_ok=True, parents=True)
    save_path = save_dir / (f"Supp_fig6a_plot_6TG_dilution_series.svg")

    fig, ax = plt.subplots(figsize=(4, 3))
    dilution_plot(
        data=df,
        x="Concentration",
        y1="HAP1-PEmax-B4",
        y2="HAP1-PEmax-MLH1dn-C4",
        y3="HAP1-PE2-MLH1KO-C1",
        ax=ax,
        save_path=save_path,
    )


# ----------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    main()
