"""
    Project :  Prime editing pilot screens. Script to generate a heatmap to visualise score differences across cell lines. Non coding experiment.
    Date : 241113
    Python version 3.10
"""

# IMPORT

import pandas as pd
import pathlib as pl
import math
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# VARIABLES

INPUT_DIR = pl.Path(
    "/Users/kajbac/code/240913_cell_line_differences/supplementary_table_alteration"
)

OUTPUT_DIR = pl.Path(
    "/Users/kajbac/code/240913_cell_line_differences/supplementary_table_alteration"
)
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)

supp_path = INPUT_DIR / "Table_supp_8_MLH1_non_coding_variant_scores.xlsx"


# ----------------------------------------------------------------------------------------------------------------------
# MAIN FUNCTION
# ----------------------------------------------------------------------------------------------------------------------
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


def heatmap(df, save_path):
    """
    Produce heatmap to visualize differences in the scores between HAP1:PEmax and HAP1:PEmax+MLH1dn cell lines.

    Parameters
    ----------
    df (df): Dataframe to be processed.
    """
    df.set_index("HGVSc", inplace=True)
    ax = sns.heatmap(
        data=df,
        # annot=True,         # Annotates cells with the score values
        cmap="viridis_r",  # Choose a color palette that suits your data
        cbar=True,  # Show color bar for reference
        linewidths=0.5,  # Add grid lines for separation
        linecolor="grey",  # Color of grid lines
        fmt=".2f",  # Format for scores (optional)
        vmin=-3.55,  # Set the minimum color scale value
        vmax=12.31,  # Set the maximum color scale value
    )

    # Add labels and title
    plt.title("Function Scores of Loss-of-Function Variants Across Conditions")
    plt.xlabel("Experimental Conditions")
    plt.ylabel("Variants")

    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")

    if save_path is None:
        plt.show()
    else:
        plt.savefig(save_path, bbox_inches="tight")
        plt.savefig(save_path.with_suffix(".pdf"), bbox_inches="tight")
        plt.savefig(save_path.with_suffix(".png"), dpi=300, bbox_inches="tight")

    return df


def main():
    """
    Main pre-processing function to generate a dataframe containing variant scores for the MLH1 non coding screen.

    Parameters
    ----------
    :
    :
    """
    save_dir = OUTPUT_DIR
    save_dir.mkdir(exist_ok=True, parents=True)
    # save_path = save_dir / (f"MLH1non_coding_score_differences_heatmap.svg")
    # save_path_neutral = save_dir / (f"MLH1non_coding_score_differences_heatmap_neutral.svg")

    set_style(context="paper", font_scale=1, style="ticks")

    df_supp = pd.read_excel(supp_path, header=1)
    df_supp = df_supp[df_supp["fdr_0.01"] == True]
    print(len(df_supp))
    df_supp = df_supp[
        [
            "HGVSc",
            "PEmax_Function_Score",
            "PEmax_obn_Function_Score",
            "PEmax_MLH1dn_Function_Score",
            "PEmax_MLH1dn_obn_Function_Score",
        ]
    ]
    df_supp["HGVSc"] = df_supp["HGVSc"].str.split(":").str[1]
    df_supp.columns = df_supp.columns.str.replace("_Function_Score", "", regex=False)

    # Split DataFrame into three parts
    df1 = df_supp.iloc[:22]
    df2 = df_supp.iloc[22:44]
    df3 = df_supp.iloc[44:]

    # Generate heatmaps for each subset
    fig, ax = plt.subplots(figsize=(4, 6))
    heatmap(
        df1, save_path=save_dir / "MLH1non_coding_score_differences_heatmap_part1.svg"
    )

    fig, ax = plt.subplots(figsize=(4, 6))
    heatmap(
        df2, save_path=save_dir / "MLH1non_coding_score_differences_heatmap_part2.svg"
    )

    fig, ax = plt.subplots(figsize=(4, 3.81))
    heatmap(
        df3, save_path=save_dir / "MLH1non_coding_score_differences_heatmap_part3.svg"
    )

    # df_neutral= pd.read_excel(
    #     supp_path,header=1)
    # df_neutral= df_neutral[df_neutral['fdr_0.01'] == False]
    # df_neutral=df_neutral[["HGVSc", "PEmax_Function_Score", "PEmax_obn_Function_Score", "PEmax_MLH1dn_Function_Score","PEmax_MLH1dn_obn_Function_Score"]]
    # df_neutral["HGVSc"] = df_neutral["HGVSc"].str.split(':').str[1]

    # fig, ax = plt.subplots(figsize=(4, 40))
    # heatmap(df_neutral, save_path_neutral)


# ----------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    main()
