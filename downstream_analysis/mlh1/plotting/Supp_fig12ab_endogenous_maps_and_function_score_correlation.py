"""
    Project :  Prime editing pilot screens. Plotting script to generate Supplementary Figure 12 a and b. Endogenous function scores by genomic positions for each screen condition and correlation of endogenous and function scores for each condition.
    Date : 240316
    Python version 3.10

"""
# IMPORT

import pandas as pd
import pathlib as pl
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.offsetbox import AnchoredText
from scipy import stats
from scipy.stats import norm

# PATHS
INPUT_DIR = pl.Path(
    f"/Users/kajbac/code/240116_final_scripts_for_PE_paper/MLH1x10/endogenous/log2_enriched"
)

OUTPUT_DIR = pl.Path(
    "/Users/kajbac/code/240116_final_scripts_for_PE_paper/plotting/_plots"
)
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)

PE_max_path = INPUT_DIR / "S1_data_full.csv"
PEmax_O_path = INPUT_DIR / "S1O_data_full.csv"
PEmax_MLH1dn_path = INPUT_DIR / "S2_data_full.csv"
PEmax_MLH1dn_O_path = INPUT_DIR / "S2O_data_full.csv"

csv_paths = [PE_max_path, PEmax_O_path, PEmax_MLH1dn_path, PEmax_MLH1dn_O_path]

score_list = [
    "CK003_CK007_log2_mean_variant",
    "CK005_CK009_log2_mean_variant",
    "CK004_CK008_log2_mean_variant",
    "CK006_CK010_log2_mean_variant",
]
sample_list = ["CK003_CK007", "CK005_CK009", "CK004_CK008", "CK006_CK010"]

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
# ----------------------------------------------------------------------------------------------------------------------
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
    sns.set(style=style, context=context, font_scale=1.2)


def position_scatterplot(df, hue_order, x, y, hue, naming_dict, ax):
    """
    Plotting function generating scatterplot with endogenous function scores over GRCh38 locations.

    Parameters
    ----------
    df (df): Dataframe to be plotted
    hue_order: Scatterplot parameter.
    x: Scatterplot parameter.
    y: Scatterplot parameter.
    hue: Scatterplot parameter.
    naming_dict (dict): Dictionary for renaming hue annotations.
    ax = Scatterplot parameter.
    """
    sns.set_palette(sns.color_palette(consequence_colors))

    variant_count = df[y].count()
    print(f"variant_count here: {variant_count}")

    df[hue].replace(to_replace=naming_dict, inplace=True)

    df = df.dropna(subset="fdr_0.01")
    df_no_lof = df.loc[df["fdr_0.01"] == False]
    df_lof = df.loc[df["fdr_0.01"] == True]

    sns.scatterplot(
        data=df_no_lof,
        x=x,
        y=y,
        hue=hue,
        hue_order=hue_order,
        ax=ax,
        s=10,
        linewidth=0.2,
        edgecolor="black",
        alpha=0.7,
    )
    sns.scatterplot(
        data=df_lof,
        x=x,
        y=y,
        hue=hue,
        hue_order=hue_order,
        ax=ax,
        marker=(8, 1, 0),
        s=30,
        linewidth=0.2,
        edgecolors="black",
        alpha=0.7,
    )
    ax.spines[["right", "top"]].set_visible(False)

    ax.set_xlabel("GRCh38 position")
    ax.set_ylabel("Endogenous function score")


def pegRNA_variant_score_correlation_scatter(
    df, hue, hue_order, naming_dict, pegRNA, ax
):
    """
    Plotting function generating scatterplot, correlating endogenous function scores against function scores.

    Parameters
    ----------
    df (df): Dataframe to be plotted.
    hue: Scatterplot parameter.
    hue_order: Scatterplot parameter.
    naming_dict (dict): Dictionary for renaming hue annotations.
    pegRNA (str): Sample identifier.
    ax = Scatterplot parameter.
    """
    sns.set_palette(sns.color_palette(consequence_colors))

    df.rename(
        columns={
            f"{pegRNA}_log2_mean_variant_normalised": "Function score",
            f"D34_D20_log2_normalised": "Endogenous function score",
        },
        inplace=True,
    )
    df = df.dropna(subset="Endogenous function score")
    df = df.dropna(subset="Function score")

    r = df["Function score"].corr(df["Endogenous function score"], method="pearson")
    r = round(r, 3)
    df[hue].replace(to_replace=naming_dict, inplace=True)

    df = df.dropna(subset="fdr_0.01")
    df_no_lof = df.loc[df["fdr_0.01"] == False]
    df_lof = df.loc[df["fdr_0.01"] == True]

    sns.scatterplot(
        data=df_no_lof,
        x="Endogenous function score",
        y="Function score",
        hue=hue,
        hue_order=hue_order,
        linewidth=0.2,
        edgecolors="black",
        alpha=0.7,
        s=10,
        ax=ax,
    )
    sns.scatterplot(
        data=df_lof,
        x="Endogenous function score",
        y="Function score",
        hue=hue,
        hue_order=hue_order,
        marker=(8, 1, 0),
        s=30,
        linewidth=0.2,
        edgecolors="black",
        alpha=0.7,
        ax=ax,
    )
    ax.spines[["right", "top"]].set_visible(False)
    anc = AnchoredText(
        f"r = {r}\n",
        # \nno. variants = {variant_count} <----add this line to display variant count
        loc="upper left",
        frameon=False,
    )
    ax.add_artist(anc)


def panel_position_scatter_plot(
    df_expected_variants_list, save_path, hue, hue_order, naming_dict
):
    """
    Plotting function generating panel of with endogenous function scores over GRCh38 locations for each screen condition.

    Parameters
    ----------
    df_expected_variants_list (list): List of dataframes to be plotted.
    save_path (Path): Output path.
    hue: Scatterplot parameter.
    hue_order: Scatterplot parameter.
    naming_dict (dict): Dictionary for renaming hue annotations.
    """
    fig, axes = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(6, 4))

    position_scatterplot(
        df_expected_variants_list[0],
        hue_order,
        x=f"GRCh38startlocation",
        y=f"D34_D20_log2_normalised",
        hue=hue,
        naming_dict=naming_dict,
        ax=axes[0, 0],
    )

    position_scatterplot(
        df_expected_variants_list[1],
        hue_order,
        x=f"GRCh38startlocation",
        y=f"D34_D20_log2_normalised",
        hue=hue,
        naming_dict=naming_dict,
        ax=axes[1, 0],
    )

    position_scatterplot(
        df_expected_variants_list[2],
        hue_order,
        x=f"GRCh38startlocation",
        y=f"D34_D20_log2_normalised",
        hue=hue,
        naming_dict=naming_dict,
        ax=axes[0, 1],
    )

    position_scatterplot(
        df_expected_variants_list[3],
        hue_order,
        x=f"GRCh38startlocation",
        y=f"D34_D20_log2_normalised",
        hue=hue,
        naming_dict=naming_dict,
        ax=axes[1, 1],
    )

    axes[0, 0].get_legend().remove()
    axes[1, 0].get_legend().remove()
    axes[0, 1].get_legend().remove()
    axes[1, 1].get_legend().remove()

    fig.tight_layout()
    if save_path is None:
        plt.show()
    else:
        plt.savefig(save_path)
        # plt.savefig(save_path.with_suffix(".pdf"))
        # plt.savefig(save_path.with_suffix(".png"), dpi=300)


def panel_pegRNA_variant_score_correlation_scatter(
    df_list, save_path, hue, hue_order, naming_dict
):
    """
    Plotting function generating panel of scatterplots, correlating endogenous function scores against function scores for each screen condition.

    Parameters
    ----------
    df_list (list): List of dataframe to be plotted.
    save_path (Path): Output path.
    hue: Scatterplot parameter.
    hue_order: Scatterplot parameter.
    naming_dict (dict): Dictionary for renaming hue annotations.
    """
    fig, axes = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(4, 4))

    pegRNA_variant_score_correlation_scatter(
        df=df_list[0],
        hue=hue,
        hue_order=hue_order,
        naming_dict=naming_dict,
        pegRNA=sample_list[0],
        ax=axes[0, 0],
    )

    pegRNA_variant_score_correlation_scatter(
        df=df_list[1],
        hue=hue,
        hue_order=hue_order,
        naming_dict=naming_dict,
        pegRNA=sample_list[1],
        ax=axes[1, 0],
    )

    pegRNA_variant_score_correlation_scatter(
        df=df_list[2],
        hue=hue,
        hue_order=hue_order,
        naming_dict=naming_dict,
        pegRNA=sample_list[2],
        ax=axes[0, 1],
    )

    pegRNA_variant_score_correlation_scatter(
        df=df_list[3],
        hue=hue,
        hue_order=hue_order,
        naming_dict=naming_dict,
        pegRNA=sample_list[3],
        ax=axes[1, 1],
    )

    axes[0, 0].get_legend().remove()
    axes[1, 0].get_legend().remove()
    axes[0, 1].get_legend().remove()
    axes[1, 1].get_legend().remove()

    fig.tight_layout()
    if save_path is None:
        plt.show()
    else:
        plt.savefig(save_path)
        # plt.savefig(save_path.with_suffix(".pdf"))
        # plt.savefig(save_path.with_suffix(".png"), dpi=300)


# ----------------------------------------------------------------------------------------------------------------------
# MAIN FUNCTION
# ----------------------------------------------------------------------------------------------------------------------
def main():
    """
    Main function to plot Supplementary Figures 12 a and b for PE manuscript. Endogenous function scores by genomic positions for each screen condition and correlation of endogenous and function scores for each condition.
    """
    save_dir = OUTPUT_DIR
    save_dir.mkdir(exist_ok=True, parents=True)
    save_path_gene_map_log2_pegRNA_filtered = (
        save_dir / "Supp_fig12a_endogenous_function_score_position_scatter_panel.svg"
    )
    save_path_pegRNA_variant_score_correlation = (
        save_dir / "Supp_fig12b_endogenous_function_score_correlation_scatter_panel.svg"
    )

    set_style(context="paper", font_scale=1, style="ticks")

    df_expected_variants_list_log2_pegRNA_filtered = []

    for i, csv_path in enumerate(csv_paths):
        df = pd.read_csv(csv_path)

        # filter on expected variants, D4 to D20 enrichment, variants that received a pegRNA-derived function score
        df_expected_variants_log2_filtered = df.loc[
            (df["expected_variant"] == True) & (df["D20_D4_log2"] > 1)
        ]
        df_expected_variants_log2_filtered[
            "variant_scored"
        ] = df_expected_variants_log2_filtered[score_list[i]].notna()
        df_expected_variants_log2_pegRNA_filtered = (
            df_expected_variants_log2_filtered.loc[
                df_expected_variants_log2_filtered["variant_scored"] == True
            ]
        )
        df_expected_variants_log2_pegRNA_filtered["D34_D20_log2_normalised"] = (
            df_expected_variants_log2_pegRNA_filtered["D34_D20_log2"]
            - df_expected_variants_log2_pegRNA_filtered.loc[
                (
                    df_expected_variants_log2_pegRNA_filtered["Consequence"]
                    == "SYNONYMOUS"
                )
            ]["D34_D20_log2"].median()
        )
        df_expected_variants_list_log2_pegRNA_filtered.append(
            df_expected_variants_log2_pegRNA_filtered
        )

    panel_position_scatter_plot(
        df_expected_variants_list_log2_pegRNA_filtered,
        save_path_gene_map_log2_pegRNA_filtered,
        hue="Consequence",
        hue_order=consequence_order,
        naming_dict=consequence_naming_dict,
    )
    panel_pegRNA_variant_score_correlation_scatter(
        df_expected_variants_list_log2_pegRNA_filtered,
        save_path_pegRNA_variant_score_correlation,
        hue="Consequence",
        hue_order=consequence_order,
        naming_dict=consequence_naming_dict,
    )


# ----------------------------------------------------------------------------------------------------------------------

if __name__ == "__main__":
    main()
