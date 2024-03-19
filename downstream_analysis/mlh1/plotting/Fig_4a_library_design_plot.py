"""
    Project :  Prime editing pilot screens. Plotting script to generate Figure 4 a. Intronic library design: psoitions of variants.
    Date : 240316
    Python version 3.10

"""

# IMPORT

import pandas as pd
import pathlib as pl
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import matplotlib.offsetbox
from matplotlib.offsetbox import AnchoredText
from matplotlib.lines import Line2D
from scipy import stats
from scipy.stats import norm

# VARIABLES

pseudo_pre_freq_filter = 1e-4
ST_editing_filter_PEmax = 5
ST_editing_filter_PEmaxdn = 25

# PATHS

INPUT_DIR = pl.Path("/Users/kajbac/code/230904_S002_final_filter")

OUTPUT_DIR = pl.Path(
    "/Users/kajbac/code/240116_final_scripts_for_PE_paper/plotting/_plots"
)
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)

df_non_coding_library_path = (
    INPUT_DIR / "MLH1intronic_pegRNA_data_merged_best_pegRNA.csv"
)

# --------------------------------------------------------------------------------------------------------------------
# FUNCTIONS


def translate_GRCh37locations_to_GRCh38locations(df):
    """
    Lift GRCh37 locations to GRCh38 locations.

    Parameters
    ----------
    df (df): Dataframe to be processed.
    """
    df["start"] = df["start"] - 41491
    df.rename(columns={"start": "GRCh38startlocation"}, inplace=True)
    df["end"] = df["end"] - 41491
    df.rename(columns={"end": "GRCh38endlocation"}, inplace=True)


# Intronic library variant positions
def intronic_position_plot(df, save_path):
    """
    Plotting function to visualise variants present in the library as ticks across the MLH1 sequence space.

    Parameters
    ----------
    df (df): Dataframe to be plotted
    save_path (path): Output path.
    """

    greys = sns.color_palette("Greys_r", 4)
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
    fig, ax = plt.subplots(figsize=(6, 1.5))
    ax.axhline(
        y=0,
        color=greys[0],
        linestyle="-",
    )
    ax.plot([36990548, 36990548 + 2000], [-0.9, -0.9], "-", color="k")
    for x_pos in df["GRCh38startlocation"].values:
        ax.plot([x_pos, x_pos], [-0.25, 0.25], lw=1, color="#4878cf", alpha=0.3)
    ax.set_ylim(-1, 1)
    for a, b in exon_windows:
        plt.axvspan(a, b, facecolor=greys[0], edgecolor=None)
    fig.subplots_adjust(left=0.025, right=0.975, bottom=0.2, top=0.8)
    ax.spines[["left", "right", "top", "bottom"]].set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_yticklabels([])

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
    Main function to plot Figure 4 a for PE manuscript. Variants across MLH1 sequence space. MLH1 non-coding screen.
    """
    save_dir = OUTPUT_DIR
    save_dir.mkdir(exist_ok=True, parents=True)
    save_path = save_dir / (f"Fig_4a_library_design_variant_positions_plot.svg")
    df_intronic_library = pd.read_csv(
        df_non_coding_library_path, delimiter=",", low_memory=False
    )
    translate_GRCh37locations_to_GRCh38locations(df_intronic_library)
    intronic_position_plot(df_intronic_library, save_path=save_path)


# ----------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    main()
