"""
    Project :  Prime editing pilot screens. Plotting script to generate Figure 3 e and Supplementary Figure 9 c. Correlation Scatterplot: functions scores (variant scores) against Cadd (Phred) scores. MLH1 saturation screen.
    Date : 240316
    Python version 3.10

"""

# IMPORT

import pandas as pd
import pathlib as pl
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.offsetbox import AnchoredText
from scipy import stats
from scipy.stats import norm

# VARIABLES

pseudo_pre_freq_filter = 1.4e-4
ST_editing_filter_PEmax = 5
ST_editing_filter_PEmaxdn = 25

# Number of replicates in which a variant has to be scored.
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
sort_dict_consequence = {
    "SYNONYMOUS": 0,
    "STOP_GAINED": 1,
    "CANONICAL_SPLICE": 2,
    "NON_SYNONYMOUS": 3,
    "SPLICE_SITE": 4,
    "INTRONIC": 5,
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


def determine_fdr(df):
    """
    FDR testing of pegRNA-derived function scores with the following thresholds: 0.01, 0.02, 0.05 using synonymous variants as controls.

    Parameters
    ----------
    df (df): Dataframe to be processed.
    """
    loc, scale = stats.norm.fit(
        df.loc[df["Consequence"] == "SYNONYMOUS"][
            "log2_mean_variant_replicates_normalised"
        ]
    )
    df["p_values"] = 1 - stats.norm.cdf(
        df["log2_mean_variant_replicates_normalised"], loc, scale
    )
    df["p_adjust"] = stats.false_discovery_control(df["p_values"], method="bh")
    df["fdr_0.01"] = df["p_adjust"].map(lambda x: True if x < 0.01 else False)
    df["fdr_0.02"] = df["p_adjust"].map(lambda x: True if x < 0.02 else False)
    df["fdr_0.05"] = df["p_adjust"].map(lambda x: True if x < 0.05 else False)

    return df


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


def Correlation_scatter(df, x, y, hue, hue_order, ax, save_path=None):
    """
    Plotting function generating correlation scatterplot with function scores against CADD (Phred) scores, marking LoF variants as stars.

    Parameters
    ----------
    df (df): Dataframe to be plotted.
    x: Scatterplot parameter.
    y: Scatterplot parameter.
    hue: Scatterplot parameter.
    hue_order: Scatterplot parameter.
    ax: Scatterplot parameter.
    save_path (path): Output path.
    """
    sns.set_palette(sns.color_palette(consequence_colors))

    df["Consequence"].replace(to_replace=consequence_naming_dict, inplace=True)

    variants = df["log2_mean_variant_replicates_normalised"].count()
    print(f"number of variants plotted is {variants}")

    r = df[x].corr(df[y], method="pearson")
    r = round(r, 3)

    df_lof = df.copy().loc[df["fdr_0.01"] == True]
    df_null = df.copy().loc[df["fdr_0.01"] != True]

    sns.scatterplot(
        data=df_null,
        x=x,
        y=y,
        hue=hue,
        hue_order=hue_order,
        ax=ax,
        linewidth=0.2,
        edgecolors="black",
        alpha=0.7,
        s=10,
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
        alpha=0.6,
    )
    ax.spines[["right", "top"]].set_visible(False)
    anc = AnchoredText(
        f"r = {r}",
        loc="upper left",
        frameon=False,
    )
    ax.add_artist(anc)

    ax.set_ylabel(f"Function score")
    ax.set_xlabel(f"CADD score")
    ax.get_legend().remove()

    if save_path is None:
        plt.show()
    else:
        plt.savefig(save_path)
        # plt.savefig(save_path.with_suffix(".pdf"))
        # plt.savefig(save_path.with_suffix(".png"), dpi=300)


def Correlation_scatter_hc(df, x, y, hue, hue_order, ax, save_path=None):
    """
    Plotting function generating correlation scatterplot with function scores against CADD (Phred) scores, marking LoF variants as stars.

    Parameters
    ----------
    df (df): Dataframe to be plotted.
    x: Scatterplot parameter.
    y: Scatterplot parameter.
    hue: Scatterplot parameter.
    hue_order: Scatterplot parameter.
    ax: Scatterplot parameter.
    save_path (path): Output path.
    """
    sns.set_palette(sns.color_palette(consequence_colors))

    df["Consequence"].replace(to_replace=consequence_naming_dict, inplace=True)

    variants = df["log2_mean_variant_replicates_normalised"].count()
    print(f"number of hc variants plotted is {variants}")

    r = df[x].corr(df[y], method="pearson")
    r = round(r, 3)

    sns.scatterplot(
        data=df,
        x=x,
        y=y,
        hue=hue,
        hue_order=hue_order,
        ax=ax,
        linewidth=0.2,
        edgecolors="black",
        alpha=0.7,
        s=10,
    )
    ax.spines[["right", "top"]].set_visible(False)
    anc = AnchoredText(
        f"r = {r}",
        loc="upper left",
        frameon=False,
    )
    ax.add_artist(anc)

    ax.set_ylabel(f"Function score")
    ax.set_xlabel(f"CADD score")
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
    Main function to plot Figure 3 e and Supplemantary Figure 9 c for PE manuscript. Scatterplot: functions scores (variant scores) against CADD scores for both relaxed and high-stringency datasets. MLH1 saturation screen.
    """
    save_dir = OUTPUT_DIR
    save_dir.mkdir(exist_ok=True, parents=True)

    save_path = save_dir / (f"Fig_3e_plot_replicates_cadd_score_correlation.svg")
    save_path_hc = save_dir / (
        f"Supp_fig9c_plot_replicates_cadd_score_correlation_hc.svg"
    )

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
    )  # important! this dataframe now contains variant level information!

    # determine fdr
    df_variant_score_filtered = df_variant_score_filtered.pipe(determine_fdr)

    fig, ax = plt.subplots(figsize=(2, 1.5))

    Correlation_scatter(
        df=df_variant_score_filtered.dropna(subset="PHRED").loc[
            df_variant_score_filtered["Consequence"].isin(["NON_SYNONYMOUS"])
        ],
        x="PHRED",
        y="log2_mean_variant_replicates_normalised",
        hue="Consequence",
        ax=ax,
        hue_order=consequence_order,
        save_path=save_path,
    )

    # high-stringency dataset
    fig, ax = plt.subplots(figsize=(3, 1.5))

    Correlation_scatter_hc(
        df=df_variant_score_filtered.dropna(subset="PHRED")
        .loc[df_variant_score_filtered["Consequence"].isin(["NON_SYNONYMOUS"])]
        .loc[df_variant_score_filtered["hc"] == True],
        x="PHRED",
        y="log2_mean_variant_replicates_normalised",
        hue="Consequence",
        ax=ax,
        hue_order=consequence_order,
        save_path=save_path_hc,
    )


# ----------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    main()
