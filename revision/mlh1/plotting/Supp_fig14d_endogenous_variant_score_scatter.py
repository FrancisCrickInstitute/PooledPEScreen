"""
    Project :  Prime editing pilot screens. Revision experiment. Plotting script to visualise endogenous function scores over variant IDs (Supp. Fig. 14d). 
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

INPUT_DIR = pl.Path(
    "/Users/kajbac/code/240830_revision_experiment/preprocessing_output"
)

OUTPUT_DIR = pl.Path("/Users/kajbac/code/240830_revision_experiment/_plots")
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)

df_variant_score_path = INPUT_DIR / "data_MLH1_revision_pegRNA_score.csv"

df_u1_R1_path = INPUT_DIR / f"CK_u1_HAP1_L_R1_data_full.csv"
df_u1_R2_path = INPUT_DIR / f"CK_u1_HAP1_L_R2_data_full.csv"

df_u2_R1_path = INPUT_DIR / f"CK_u2_HAP1_L_R1_data_full.csv"
df_u2_R2_path = INPUT_DIR / f"CK_u2_HAP1_L_R2_data_full.csv"

df_i1_R1_path = INPUT_DIR / f"CK_i1_HAP1_L_R1_data_full.csv"
df_i1_R2_path = INPUT_DIR / f"CK_i1_HAP1_L_R2_data_full.csv"

df_x10_R1_path = INPUT_DIR / f"CK_x10_HAP1_L_R1_data_full.csv"
df_x10_R2_path = INPUT_DIR / f"CK_x10_HAP1_L_R2_data_full.csv"

df_x11_R1_path = INPUT_DIR / f"CK_x11_HAP1_L_R1_data_full.csv"
df_x11_R2_path = INPUT_DIR / f"CK_x11_HAP1_L_R2_data_full.csv"

df_i15_R1_path = INPUT_DIR / f"CK_i15_HAP1_L_R1_data_full.csv"
df_i15_R2_path = INPUT_DIR / f"CK_i15_HAP1_L_R2_data_full.csv"

# VARIABLES
u1_variants = [2289, 2297]
u2_variants = [2353, 2358]
i1_variants = "{36993971: 'C'}"
x10_variants = [1966, 1980, 2029, 1856, 1978]
x11_variants = [2818]
# i15_variants = [3128]
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
neg_controls = [1978, 2385]
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


def average_replicates(df):
    """
    Average function scores per variant across replicates.

    Parameters
    ----------
    df (df): Dataframe to be processed.
    """
    df_log2_mean = df.groupby(["variant_name"])["D33_D19_log2_normalised"].mean()
    df_merge = pd.merge(
        df,
        df_log2_mean,
        on=["variant_name"],
        how="left",
        validate="many_to_one",
        suffixes=("", "_mean_replicates"),
    )

    return df_merge


def normalise_log2_scores(df):
    """
    Normalise endogenous log2-ratios to the median endogenous log2-ratio of control variants for each replicate.

    Parameters
    ----------
    df (df): Dataframe to be processed.
    """
    df["D33_D19_log2_normalised"] = (
        df["D33_D19_log2"]
        - df.loc[df["mutant index"].isin(neg_controls)]["D33_D19_log2"].median()
    )

    return df


def normalise_averaged_log2_scores(df):
    """
    Normalise endogenous log2-ratios to the median endogenous log2-ratio of control variants of scores averaged across replicates.

    Parameters
    ----------
    df (df): Dataframe to be processed.
    """
    df["D33_D19_log2_normalised_mean_replicates_normalised"] = (
        df["D33_D19_log2_normalised_mean_replicates"]
        - df.loc[df["mutant index"].isin(neg_controls)][
            "D33_D19_log2_normalised_mean_replicates"
        ].median()
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
    Plotting function generating scatterplot with function scores over variant names.

    Parameters
    ----------
    df (df): Dataframe to be plotted
    x: Scatterplot parameter.
    y: Scatterplot parameter.
    hue: Scatterplot parameter
    save_path (path): Output path.
    """
    color_index_to_skip = [
        7,
        3,
    ]  # This is the index corresponding to V8 (deep intronic variant) and V4 (canonical splice control)
    # Create a new palette without the color for V8 and V4
    original_palette = sns.color_palette(cc.glasbey, n_colors=13)
    new_palette = [
        color
        for i, color in enumerate(original_palette)
        if i not in color_index_to_skip
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
    Function to plot function scores over variant IDs (Supp. Fig. 14d).
    """
    save_dir = OUTPUT_DIR
    save_dir.mkdir(exist_ok=True, parents=True)
    save_path = save_dir / (
        f"Supp_Fig_14d_HAP1_endogenous_variant_score.svg"
    )

    set_style(context="paper", font_scale=1, style="ticks")

    df_u1_R1 = pd.read_csv(df_u1_R1_path, delimiter=",")
    df_u1_R1 = df_u1_R1[["mutant index", "D19_freq", "D33_D19_log2"]]
    df_u1_R1 = df_u1_R1.loc[df_u1_R1["mutant index"].isin(u1_variants)]

    df_u1_R2 = pd.read_csv(df_u1_R2_path, delimiter=",")
    df_u1_R2 = df_u1_R2[["mutant index", "D19_freq", "D33_D19_log2"]]
    df_u1_R2 = df_u1_R2.loc[df_u1_R2["mutant index"].isin(u1_variants)]

    df_u2_R1 = pd.read_csv(df_u2_R1_path, delimiter=",")
    df_u2_R1 = df_u2_R1[["mutant index", "D19_freq", "D33_D19_log2"]]
    df_u2_R1 = df_u2_R1.loc[df_u2_R1["mutant index"].isin(u2_variants)]

    df_u2_R2 = pd.read_csv(df_u2_R2_path, delimiter=",")
    df_u2_R2 = df_u2_R2[["mutant index", "D19_freq", "D33_D19_log2"]]
    df_u2_R2 = df_u2_R2.loc[df_u2_R2["mutant index"].isin(u2_variants)]

    df_i1_R1 = pd.read_csv(df_i1_R1_path, delimiter=",")
    df_i1_R1 = df_i1_R1[["mutant index", "variant_dict", "D19_freq", "D33_D19_log2"]]
    df_i1_R1 = df_i1_R1.loc[df_i1_R1["variant_dict"] == i1_variants]

    df_i1_R2 = pd.read_csv(df_i1_R2_path, delimiter=",")
    df_i1_R2 = df_i1_R2[["mutant index", "variant_dict", "D19_freq", "D33_D19_log2"]]
    df_i1_R2 = df_i1_R2.loc[df_i1_R2["variant_dict"] == i1_variants]

    df_x10_R1 = pd.read_csv(df_x10_R1_path, delimiter=",")
    df_x10_R1 = df_x10_R1[["mutant index", "D19_freq", "D33_D19_log2"]]
    df_x10_R1 = df_x10_R1.loc[df_x10_R1["mutant index"].isin(x10_variants)]

    df_x10_R2 = pd.read_csv(df_x10_R2_path, delimiter=",")
    df_x10_R2 = df_x10_R2[["mutant index", "D19_freq", "D33_D19_log2"]]
    df_x10_R2 = df_x10_R2.loc[df_x10_R2["mutant index"].isin(x10_variants)]

    df_x11_R1 = pd.read_csv(df_x11_R1_path, delimiter=",")
    df_x11_R1 = df_x11_R1[["mutant index", "D19_freq", "D33_D19_log2"]]
    df_x11_R1 = df_x11_R1.loc[df_x11_R1["mutant index"].isin(x11_variants)]

    df_x11_R2 = pd.read_csv(df_x11_R2_path, delimiter=",")
    df_x11_R2 = df_x11_R2[["mutant index", "D19_freq", "D33_D19_log2"]]
    df_x11_R2 = df_x11_R2.loc[df_x11_R2["mutant index"].isin(x11_variants)]

    dfs_R1 = [df_u1_R1, df_u2_R1, df_i1_R1, df_x10_R1, df_x11_R1]
    dfs_R2 = [df_u1_R2, df_u2_R2, df_i1_R2, df_x10_R2, df_x11_R2]

    df_concat_R1 = pd.concat(dfs_R1, axis=0, ignore_index=True)
    df_concat_R1["mutant index"] = df_concat_R1["mutant index"].astype(int)
    df_concat_R1 = normalise_log2_scores(df_concat_R1)

    df_concat_R2 = pd.concat(dfs_R2, axis=0, ignore_index=True)
    df_concat_R2["mutant index"] = df_concat_R2["mutant index"].astype(int)
    df_concat_R2 = normalise_log2_scores(df_concat_R2)

    df_merge = pd.merge(
        df_concat_R1, df_concat_R2, on="mutant index", suffixes=("_R1", "_R2")
    )

    df_merge["mutant index"] = df_merge["mutant index"].astype(int)
    df_merge["variant_name"] = df_merge["mutant index"].map(variant_naming_dict)
    # df_merge = df_merge.sort_values(by='variant_name')
    print(df_merge)
    # df_concat = average_replicates(df_concat)

    df_merge.to_csv((save_dir / f"data_test_merge.csv").as_posix())
    # df_concat = df_concat.drop_duplicates(subset = "variant_name")
    # df_concat = normalise_averaged_log2_scores(df_concat)
    # df_concat = df_concat.sort_values(by='variant_name')
    df_merge = df_merge.sort_values(
        by=["variant_name"], key=lambda x: x.map(sort_dict_plot)
    )

    df_melted = df_merge.melt(
        id_vars="variant_name",
        value_vars=["D33_D19_log2_normalised_R1", "D33_D19_log2_normalised_R2"],
        var_name="Measure",
        value_name="Endogenous Function Score",
    )

    # df_melted = df_melted.sort_values(by='variant_name')

    print(df_melted)

    # Plotting
    # HAP1 endogenous scores
    fig, ax = plt.subplots(figsize=(3.5, 2))

    variant_scatterplot(
        df=df_melted,
        x="variant_name",
        y="Endogenous Function Score",
        hue="variant_name",
        ax=ax,
    )

    if save_path is None:
        plt.show()
    else:
        plt.savefig(save_path)
        #     # # plt.savefig(save_path.with_suffix(".pdf"))
        plt.savefig(save_path.with_suffix(".png"), dpi=300)


# ----------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    main()
