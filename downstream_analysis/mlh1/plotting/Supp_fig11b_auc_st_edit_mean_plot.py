"""
    Project :  Prime editing pilot screens. Plotting script to generate Supplementary Figure 11 b. AUC by mean ST editing filters of increasing stringency. Filter applied on the variant level, AUC calculated on the variant level. MLH1 saturation screen.
    Date : 240316
    Python version 3.10

"""

# IMPORT

import pandas as pd
import pathlib as pl
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from sklearn.metrics import roc_auc_score

# VARIABLES

pseudo_pre_freq_filter = 1.4e-4
ST_editing_filter_PEmax = 5
ST_editing_filter_PEmaxdn = 25

replicate_number = 2

# PATHS

INPUT_DIR = pl.Path(
    f"/Users/kajbac/code/240116_final_scripts_for_PE_paper/MLH1x10/pegRNA/variant_score/pseudo_pre_freq_{pseudo_pre_freq_filter}/{ST_editing_filter_PEmax}_PEmax_{ST_editing_filter_PEmaxdn}_PEmaxdn"
)

df_variant_score_path = INPUT_DIR / "data_MLH1x10_variant_score_replicates.csv"

OUTPUT_DIR = pl.Path(
    "/Users/kajbac/code/240116_final_scripts_for_PE_paper/plotting/_plots"
)
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)
# --------------------------------------------------------------------------------------------------------------------
# FUNCTIONS

# derive scores on the variant level to calculate AUC


# generate AUC dataframe
def generate_auc_dataframe(df):
    """
    Calculate AUC values as a function of continues mean ST editing thresholds with synonymous variants defined as pNeutral and non-sense and canonical splice variants as pLoF. Filter applied on the variant level, AUC calculated on the variant level. Generate data frame containing AUC values, corresponding ST-editing thresholds, and variant numbers.

    Parameters
    ----------
    df (df): Dataframe to be processed.
    """
    ST_editing_percentages = range(0, 101)
    aucs = []
    variant_numbers = []

    for percentage in ST_editing_percentages:
        df_consequence = df.copy().loc[
            df["percentage_editing_mean_variant_replicates"] >= percentage
        ]
        df_variant = df_consequence.drop_duplicates(subset=["var_key"])
        df_variant_filtered = (
            df.copy()
            .loc[df["percentage_editing_mean_variant_replicates"] >= percentage]
            .drop_duplicates(subset=["var_key"])
        )

        # map and dataframe to be used in auc_consequence calculation
        dict_consequence_map = {
            "CANONICAL_SPLICE": 1,
            "STOP_GAINED": 1,
            "SYNONYMOUS": 0,
        }
        df_variant = df_variant.loc[
            df_variant["Consequence"].isin(
                ["CANONICAL_SPLICE", "STOP_GAINED", "SYNONYMOUS"]
            )
        ]
        df_variant = df_variant.assign(
            expected_lof=df_variant["Consequence"].map(dict_consequence_map)
        )

        # calculate fpr, tpr, tnr, lof threshold
        consequence_unique = df_consequence["Consequence"].unique()
        if all([key in consequence_unique for key in dict_consequence_map.keys()]):
            auc_value_consequence = roc_auc_score(
                y_true=df_variant["expected_lof"],
                y_score=df_variant["log2_mean_variant_replicates_normalised"],
            )

            aucs.append(auc_value_consequence)
            variant_number = df_variant_filtered[
                "log2_mean_variant_replicates_normalised"
            ].count()
            variant_numbers.append(variant_number)
            print(
                f"editing_{percentage}_auc_{auc_value_consequence}_variant_{variant_number}"
            )

        else:
            aucs.append(np.nan)
            variant_numbers.append(np.nan)

    df_auc = pd.DataFrame()
    df_auc.loc[:, "AUC"] = aucs
    df_auc.loc[:, "percentage_editing"] = ST_editing_percentages
    df_auc.loc[:, "number_variants"] = variant_numbers

    return df_auc


# PLOTS


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


def auc_plot(
    data,
    x,
    y1,
    y2,
    ax,
    ST_threshold,
):
    """
    Plotting function generating lineplots with AUC values and variant numbers over a continuous range of mean ST-editing thresholds.

    Parameters
    ----------
    data (df): Dataframe to be plotted
    x: Scatterplot parameter.
    y1: Scatterplot parameter.
    y2: Scatterplot parameter.
    ax: Scatterplot parameter
    ST_threshold (int): ST-editing filters set as high-stringency filter.
    """
    sns.lineplot(data=data, x=x, y=y2, ax=ax, color="#d48640")
    if y2 == "number_variants":
        ax.set_ylim(0, 450)
        ax.set_ylabel(f"No. variants")
    elif y2 == "number_pegRNAs":
        ax.set_ylim(0, 2000)
        ax.set_ylabel(f"No. pegRNAs ($10^3$)")
        new_labels = [x / 1000 for x in ax.get_yticks()]
        ax.set_yticklabels(new_labels)
    ax.set_xlabel(f"% Correct Mean ST Edits")

    ax2 = ax.twinx()
    sns.lineplot(data=data, x=x, y=y1, ax=ax2, color="#44729d")
    ax2.set_ylim(0, 1.05)
    ax.spines[["right", "top"]].set_visible(False)
    ax2.spines[["left", "top", "bottom"]].set_visible(False)

    ax.axvline(x=ST_threshold, color="#D3D3D3", linestyle="--")

    plt.legend([], [], frameon=False)


# --------------------------------------------------------------------------------------------------------------------
# MAIN


def main():
    """
    Main function to plot Supplementary Figure 11 b for PE manuscript. AUC values as a function of continues mean ST editing thresholds with synonymous variants defined as pNeutral and non-sense and canonical splice variants as pLoF (MLH1 saturation screen).
    """
    save_dir = OUTPUT_DIR
    save_dir.mkdir(exist_ok=True, parents=True)
    save_path = save_dir / (f"Supp_fig11b_plot_auc_mean_ST_editing.svg")

    set_style(context="paper", font_scale=1, style="ticks")
    filter_number = replicate_number - 1

    df_variant_score = pd.read_csv(df_variant_score_path, delimiter=",")

    df_variant_score_filtered = df_variant_score.loc[
        df_variant_score["number_replicates_pegRNA"] > filter_number
    ]

    df_auc = generate_auc_dataframe(df_variant_score_filtered)

    # Plots
    fig, ax = plt.subplots(figsize=(2, 2))
    auc_plot(
        data=df_auc,
        # replicate = None,
        x=f"percentage_editing",
        y1="AUC",
        y2="number_variants",
        ax=ax,
        ST_threshold=36,
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
