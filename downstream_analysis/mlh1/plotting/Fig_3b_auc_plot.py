"""
    Project :  Prime editing pilot screens. Plotting script to generate Figure 3 b. AUC values as a function of continues ST editing thresholds with synonymous variants defined as pNeutral and non-sense and canonical splice variants as pLoF (MLH1 saturation screen). Filter applied on the pegRNA level, AUC calculated on the variant level.
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

threshold_list = [5, 25, 5, 25]

# PATHS
INPUT_DIR = pl.Path(
    f"/Users/kajbac/code/240116_final_scripts_for_PE_paper/MLH1x10/pegRNA/pegRNA_score/pseudo_pre_freq_{pseudo_pre_freq_filter}/{ST_editing_filter_PEmax}_PEmax_{ST_editing_filter_PEmaxdn}_PEmaxdn"
)

OUTPUT_DIR = pl.Path(
    "/Users/kajbac/code/240116_final_scripts_for_PE_paper/plotting/_plots"
)
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)

df_pegRNA_path = INPUT_DIR / "data_MLH1x10_pegRNA_score_replicates.csv"
# --------------------------------------------------------------------------------------------------------------------
# FUNCTIONS


# Derive variant-level scores (function scores) to calculate AUC
def collapse_variants(df):
    """
    Average pegRNA log2-ratios of pegRNAs programming the same variant.

    Parameters
    ----------
    df (df): Dataframe to be processed.
    """
    log2_list = []
    for var_dict in variables_per_sample:
        log2_list.append(f"{var_dict['pegRNA_D20']}_{var_dict['pegRNA_D34']}_log2")
    df_log2_mean = df.groupby(["var_key"])[log2_list].mean()
    df_pegRNA_score = pd.merge(
        df,
        df_log2_mean,
        on=["var_key"],
        how="left",
        validate="many_to_one",
        suffixes=("", "_mean_variant"),
    )

    return df_pegRNA_score


def normalise_log2_scores_per_replicate_variant(df):
    """
    Normalise mean pegRNA log2-ratios (function scores) to the median mean pegRNA log2-ratio (function score) of synonymous variants for each condition.

    Parameters
    ----------
    df (df): Dataframe to be processed.
    """
    for var_dict in variables_per_sample:
        df[
            f"{var_dict['pegRNA_D20']}_{var_dict['pegRNA_D34']}_log2_mean_variant_normalised"
        ] = (
            df[f"{var_dict['pegRNA_D20']}_{var_dict['pegRNA_D34']}_log2_mean_variant"]
            - df.loc[df["Consequence"] == "SYNONYMOUS"][
                f"{var_dict['pegRNA_D20']}_{var_dict['pegRNA_D34']}_log2_mean_variant"
            ].median()
        )
    return df


# generate AUC dataframe
def generate_auc_dataframe(df, var_dict):
    """
    Calculate AUC values as a function of continues ST editing thresholds with synonymous variants defined as pNeutral and non-sense and canonical splice variants as pLoF. Filter applied on the pegRNA level, AUC calculated on the variant level. Generate data frame containing AUC values, corresponding ST-editing thresholds, and pegRNA and variant numbers.

    Parameters
    ----------
    df (df): Dataframe to be processed.
    var_dict (dict): Dictionary used to translate sample collection date to sample identifier.
    """
    ST_editing_percentages = range(0, 101)
    aucs = []
    variant_numbers = []
    pegRNA_numbers = []

    for percentage in ST_editing_percentages:
        df_consequence = df.copy().loc[
            df[f"{var_dict['pegRNA_D20']}_percentage_editing"] >= percentage
        ]
        df_variant = df_consequence.pipe(collapse_variants)
        df_variant = df_variant.pipe(normalise_log2_scores_per_replicate_variant)
        df_variant = df_variant.drop_duplicates(subset=["var_key"])
        df_filtered = df.copy().loc[
            df[f"{var_dict['pegRNA_D20']}_percentage_editing"] >= percentage
        ]

        # Map and dataframe to be used in auc_consequence and cm_consequence calculation
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

        # Calculate fpr, tpr, tnr, lof threshold
        consequence_unique = df_consequence["Consequence"].unique()
        if all([key in consequence_unique for key in dict_consequence_map.keys()]):
            auc_value_consequence = roc_auc_score(
                y_true=df_variant["expected_lof"],
                y_score=df_variant[
                    f"{var_dict['pegRNA_D20']}_{var_dict['pegRNA_D34']}_log2_mean_variant"
                ],
            )

            aucs.append(auc_value_consequence)
            variant_numbers.append(len(df_filtered.groupby("var_key")))
            pegRNA_number = df_filtered[
                f"{var_dict['pegRNA_D20']}_{var_dict['pegRNA_D34']}_log2"
            ].count()
            pegRNA_numbers.append(pegRNA_number)
            print(
                f"editing_{percentage}_auc_{auc_value_consequence}_pegRNA_{pegRNA_number}"
            )

        else:
            aucs.append(np.nan)
            variant_numbers.append(np.nan)
            pegRNA_numbers.append(np.nan)

    df_auc = pd.DataFrame()
    df_auc.loc[:, "AUC"] = aucs
    df_auc.loc[:, "percentage_editing"] = ST_editing_percentages
    df_auc.loc[:, "number_pegRNAs"] = pegRNA_numbers
    df_auc.loc[:, "number_variants"] = variant_numbers

    return df_auc


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


def auc_plot(
    data,
    x,
    y1,
    y2,
    ax,
    ST_threshold,
):
    """
    Plotting function generating lineplots with AUC values, and pegRNA numbers over a continuous range of ST-editing thresholds.

    Parameters
    ----------
    data (df): Dataframe to be plotted
    x: Scatterplot parameter.
    y1: Scatterplot parameter.
    y2: Scatterplot parameter.
    ax: Scatterplot parameter
    ST_threshold (int): ST-editing filters set for different cell lines.
    """
    sns.lineplot(data=data, x=x, y=y2, ax=ax, color="#d48640")
    if y2 == "number_variants":
        ax.set_ylim(0, 450)
        ax.set_ylabel(f"No. variants")
    elif y2 == "number_pegRNAs":
        ax.set_ylim(0, 2100)
        ax.set_ylabel(f"No. pegRNAs ($10^3$)")
        new_labels = [x / 1000 for x in ax.get_yticks()]
        ax.set_yticklabels(new_labels)
    ax.set_xlabel(f"Minimum % Correct ST Edits")

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
    Main function to plot Figure 3b for PE manuscript. AUC values as a function of continues ST editing thresholds with synonymous variants defined as pNeutral and non-sense and canonical splice variants as pLoF (MLH1x10 screen).
    """
    save_dir = OUTPUT_DIR
    save_dir.mkdir(exist_ok=True, parents=True)
    save_path = save_dir / (f"Fig_3b_plot_auc_variants.svg")

    set_style(context="paper", font_scale=1, style="ticks")

    df_pegRNA = pd.read_csv(df_pegRNA_path, delimiter=",")

    # Plots (individual replicates)
    fig, ax_list = plt.subplots(2, 2, sharex=True, figsize=(3, 2))

    ax_list = ax_list.flatten()

    for i, var_dict in enumerate(variables_per_sample):
        df_auc = generate_auc_dataframe(df_pegRNA, var_dict)
        print(df_auc)

        auc_plot(
            data=df_auc,
            # replicate = None,
            x=f"percentage_editing",
            y1="AUC",
            y2="number_pegRNAs",
            ax=ax_list[i],
            ST_threshold=threshold_list[i],
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
