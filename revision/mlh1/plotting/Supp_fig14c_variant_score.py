"""
    Project :  Prime editing pilot screens. Revision experiment. Plotting script to visualise function scores (variant scores) over variant IDs (Supp. Fig. 14c). 
    Date : 240925
    Python version 3.11.3

"""

# IMPORT

import pandas as pd
import pathlib as pl
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy import stats
from scipy.stats import norm
import colorcet as cc


# PATHS

INPUT_DIR = pl.Path("/Users/kajbac/code/240830_revision_experiment/pegRNA_score")

OUTPUT_DIR = pl.Path("/Users/kajbac/code/240830_revision_experiment/_plots")
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)

df_variant_score_path = INPUT_DIR / "data_MLH1_revision_pegRNA_score.csv"

# VARIABLES
variant_naming_dict = {
    1966: "V12",
    1980: "V13",
    2029: "V6",
    1856: "V3",
    1978: "V1",
    2289: "V9",
    2297: "V10",
    2353: "V5",
    2358: "V11",
    2818: "V7",
    3128: "V8",
    2361: "V4",
    2385: "V2",
}


HAP1_hue_naming_dict = {
    "CK001_CK006_log2_normalised": "R1",
    "CK002_CK007_log2_normalised": "R2",
}

K562_hue_naming_dict = {
    "CK003_CK008_log2_normalised": "R1",
    "CK004_CK009_log2_normalised": "R2",
}

sort_dict_plot = {
    "V1": 0,
    "V2": 1,
    "V3": 2,
    "V4": 3,
    "V5": 4,
    "V6": 5,
    "V7": 6,
    "V8": 7,
    "V9": 8,
    "V10": 9,
    "V11": 10,
    "V12": 11,
    "V13": 12,
}


# --------------------------------------------------------------------------------------------------------------------
# FUNCTIONS


def normalise_log2_scores_per_cell_line(df):
    """
    Normalise pegRNA log2-ratios to the median pegRNA log2-ratio of control variants for each cell line.

    Parameters
    ----------
    df (df): Dataframe to be processed.
    """
    df["HAP1:PE7_score_normalised"] = (
        df["HAP1:PE7_score"]
        - df.loc[df["control"] == "negative"]["HAP1:PE7_score"].median()
    )

    df["K562:PE7_score_normalised"] = (
        df["K562:PE7_score"]
        - df.loc[df["control"] == "negative"]["K562:PE7_score"].median()
    )

    return df


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


def variant_scatterplot(df, x, y, hue, ax, save_path=None):
    """
    Plotting function generating a scatterplot with function scores over variant IDs.

    Parameters
    ----------
    df (df): Dataframe to be plotted
    x: Scatterplot parameter.
    y: Scatterplot parameter.
    hue: Scatterplot parameter
    save_path (path): Output path.
    """
    color_index_to_skip = (
        7  # This is the index corresponding to V8 (deep intronic variant)
    )
    # Create a new palette without the color for V8
    original_palette = sns.color_palette(cc.glasbey, n_colors=13)
    new_palette = [
        color for i, color in enumerate(original_palette) if i != color_index_to_skip
    ]
    sns.set_palette(sns.color_palette(new_palette))

    sns.scatterplot(
        data=df,
        x=x,
        y=y,
        hue=hue,
        ax=ax,
        style="Measure",
        s=10,
        linewidth=0.2,
        edgecolor="black",
        alpha=0.3,
    )

    # Compute the average of the two replicates for each variant
    df_avg = df.groupby(["variant_name"])[y].mean().reset_index()

    # Ensure the averages follow the same x-axis order as the scatterplot
    # Sort df_avg to match the order of `variant_name` in the original df
    df_avg["order"] = df_avg["variant_name"].map(
        lambda name: df[df["variant_name"] == name].index[0]
    )
    df_avg_sorted = df_avg.sort_values(by="order").reset_index(drop=True)

    # Get the color map for the scatter plot's hue from the legend
    handles, labels = ax.get_legend_handles_labels()  # Get legend handles and labels
    color_map = {label: handle.get_color() for handle, label in zip(handles, labels)}

    # Plot horizontal lines for the average values, matching the hue color
    for i, row in df_avg_sorted.iterrows():
        variant = row["variant_name"]  # The variant name

        # Find the corresponding color based on the hue value
        hue_value = df[df["variant_name"] == variant][hue].iloc[
            0
        ]  # Get the correct hue for this variant
        color = color_map.get(
            hue_value, "black"
        )  # Get the color from the palette, default to black if not found

        # Draw the horizontal line for the average score
        ax.hlines(
            y=row[y],  # The average function score (y-value)
            xmin=i - 0.3,  # Start of the line
            xmax=i + 0.3,  # End of the line
            color=color,  # Color of the line matching the hue
            linewidth=2,  # Thickness of the line
            label="Average" if i == 0 else "",  # Add label only for the first line
        )

    ax.spines[["right", "top"]].set_visible(False)
    # sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1), frameon=False)

    ax.set_xlabel("Variants")
    ax.set_ylabel("Function score")
    plt.xticks(rotation=45)


# --------------------------------------------------------------------------------------------------------------------
# MAIN


def main():
    """
    Function to plot function scores over variant IDs (Supp. Figure 14c).
    """
    save_dir = OUTPUT_DIR
    save_dir.mkdir(exist_ok=True, parents=True)
    save_path = save_dir / (
        f"Supp_Fig_14c_function_score_HAP1_K562.svg"
    )

    set_style(context="paper", font_scale=1, style="ticks")

    df_variant_score = pd.read_csv(df_variant_score_path, delimiter=",")
    df_variant_score = normalise_log2_scores_per_cell_line(df_variant_score)

    df_variant_score["variant_name"] = df_variant_score["mutant index"].map(
        variant_naming_dict
    )
    df_variant_score = df_variant_score[df_variant_score["variant_name"] != "V8"]
    df_variant_score = df_variant_score.sort_values(
        by=["variant_name"], key=lambda x: x.map(sort_dict_plot)
    )

    df_HAP1 = df_variant_score.copy()
    df_K562 = df_variant_score.copy()

    df_HAP1_melted = df_HAP1.melt(
        id_vars="variant_name",
        value_vars=["CK001_CK006_log2_normalised", "CK002_CK007_log2_normalised"],
        var_name="Measure",
        value_name="Function Score",
    )
    df_HAP1_melted["Measure"] = df_HAP1_melted["Measure"].map(HAP1_hue_naming_dict)

    print(df_HAP1_melted)

    df_K562_melted = df_K562.melt(
        id_vars="variant_name",
        value_vars=["CK003_CK008_log2_normalised", "CK004_CK009_log2_normalised"],
        var_name="Measure",
        value_name="Function Score",
    )
    df_K562_melted["Measure"] = df_K562_melted["Measure"].map(K562_hue_naming_dict)

    print(df_K562_melted)

    # Plotting
    # HAP1 and K562

    fig, axes = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(8, 2))

    variant_scatterplot(
        df=df_HAP1_melted,
        x="variant_name",
        y="Function Score",
        hue="variant_name",
        ax=axes[0],
    )

    variant_scatterplot(
        df=df_K562_melted,
        x="variant_name",
        y="Function Score",
        hue="variant_name",
        ax=axes[1],
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
