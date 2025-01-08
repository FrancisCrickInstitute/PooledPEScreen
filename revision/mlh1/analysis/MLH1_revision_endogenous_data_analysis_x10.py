"""
    Project :  Prime editing pilot screens. Revision experiment. Pre-processing script for endogenous MLH1 x10 data. Input are DIMSUM output dataframes containing variant counts. Output dataframes are: "x10_data_full.csv". IMPORTANT: no wildtype control sequenced!
    Date : 240905
    Python version 3.10

"""
# IMPORT

import pandas as pd
import pathlib as pl
import numpy as np
from scipy import stats
from scipy.stats import norm

# VARIABLES

INPUT_DIR = pl.Path(
    "/Users/kajbac/code/240830_revision_experiment/dimsum_output_all_samples"
)

OUTPUT_DIR = pl.Path(
    "/Users/kajbac/code/240830_revision_experiment/preprocessing_output"
)
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)

samples = [
    "CK_x10_HAP1_L_R1",
    "CK_x10_HAP1_L_R2",
]

# Endogenous data paths (Dimsum output)
# D19 and D33 (variants)

sample_1_path = INPUT_DIR / f"{samples[0]}_bq5_fitness_replicates.csv"

sample_2_path = INPUT_DIR / f"{samples[1]}_bq5_fitness_replicates.csv"


# D19 and D33 (indels)

sample_1_path_indel_path = INPUT_DIR / f"{samples[0]}_bq5_indel_variant_data_merge.tsv"
sample_2_path_indel_path = INPUT_DIR / f"{samples[1]}_bq5_indel_variant_data_merge.tsv"


# D19 and D33 (rejected variants)


sample_1_rejected_path = INPUT_DIR / f"{samples[0]}_bq5_rejected_variant_data_merge.tsv"
sample_2_rejected_path = INPUT_DIR / f"{samples[1]}_bq5_rejected_variant_data_merge.tsv"


# Endogenous data csv_path lists

csv_paths = [
    sample_1_path,
    sample_2_path,
]

csv_paths_indels_u = [
    sample_1_path_indel_path,
    sample_2_path_indel_path,
]

csv_paths_rejected_u = [
    sample_1_rejected_path,
    sample_2_rejected_path,
]


csv_path_function_score = pl.Path(
    "/Users/kajbac/code/240830_revision_experiment/pegRNA_score/data_MLH1_revision_pegRNA_score.csv"
)

csv_path_u2_non_SNV = pl.Path(
    "/Users/kajbac/code/240830_revision_experiment/u2_DEL.csv"
)

pegRNA_score = "HAP1:PE7_score"

# offset is used to translate amplicon position to GRCh38Location
offset_x10 = 37017448

amplicon_length = 216


# ----------------------------------------------------------------------------------------------------------------------
#                                                    FUNCTIONS
# ----------------------------------------------------------------------------------------------------------------------
def generate_SNV_white_list_new():
    """
    Generate a "SNV white list", i.e, a list of variant key dictionaries with combinations of GRCh38Location and SNV for all expected SNVs.

    """
    df = pd.read_csv(csv_path_function_score, delimiter=",")
    df = df.loc[df["variant_type"] == "SNP"]

    SNV_white_list = []
    location_list = df.GRCh38startlocation.values.tolist()
    mut_allele_list = df.mut_allele.values.tolist()
    for i, location in enumerate(location_list):
        SNV_white_list += [{location: mut_allele_list[i]}]

    print(f"SNV_whitelist here {SNV_white_list}")
    return SNV_white_list


def create_var_key_df_pegRNA_score(df):
    """
    Generate a string identifier for each variant used for dataframe merging.

    Parameters
    ----------
    df (df): Dataframe to be processed.
    """
    df["var_key"] = df["GRCh38startlocation"].astype("string") + "_" + df["mut_allele"]
    return df


# ----------------------------------------------------------------------------------------------------------------------
# Process dimsum "fitness_replicates.RData" dataframe:


def extract_variant_positions_from_dimsum_sequence(df):
    """
    Extract variant position from nucleotide sequences.

    Parameters
    ----------
    df (df): Dataframe to be processed.
    """
    wild_type_seq = df.loc[df["Nham_nt"] == 0, "nt_seq"].iloc[0]
    sequences = df["nt_seq"].values.tolist()
    position_list = []
    for seq in sequences:
        positions = []
        for i in range(len(wild_type_seq)):
            if seq[i] != wild_type_seq[i]:
                positions.append(i + 1)
        position_list.append(positions)
    df["variant_position"] = position_list
    return df


def convert_variant_positions_to_GRCh38location(positions, offset):
    """
    Convert variant position in nucleotide sequence to GRCh38 location.

    Parameters
    ----------
    position (int): Position in nucleotide sequence to be converted.
    offset (int): Offset used to convert positions.
    """
    locations = []
    for position in positions:
        locations.append(position + offset - 1)
    return locations


def extract_variants_from_dimsum_sequence(df):
    """
    Extract variants from nucleotide sequences.

    Parameters
    ----------
    df (df): Dataframe to be processed.
    """
    wild_type_seq = df.loc[df["Nham_nt"] == 0, "nt_seq"].iloc[0]
    sequences = df["nt_seq"].values.tolist()
    variant_list = []
    for seq in sequences:
        variants = []
        for i in range(len(wild_type_seq)):
            if seq[i] != wild_type_seq[i]:
                variants.append(seq[i])
        variant_list.append(variants)
    df["variants"] = variant_list
    return df


def create_variant_dictionaries_from_dimsum_sequence(df, offset):
    """
    Compare wildtype sequence to nucleotide sequence to generate a list of dictionaries, each linking mutated allele to GRCh38 location.

    Parameters
    ----------
    df (df): Dataframe to be processed.
    offset (int): Offset used to translate variant positions to GRCh38Locations.
    """
    wild_type_seq = df.loc[df["Nham_nt"] == 0, "nt_seq"].iloc[0]
    variants_dict_list = []
    sequences = df["nt_seq"].values.tolist()
    for seq in sequences:
        variants_in_seq = {}
        for i, nt in enumerate(seq):
            if nt != wild_type_seq[i]:
                variants_in_seq[i + offset] = nt.upper()
        variants_dict_list.append(variants_in_seq)
    df["variant_dict"] = variants_dict_list
    return df


def add_white_list_annotation(variant_key, SNV_white_list):
    """
    Identify SNVs that are expected to be present based on pegRNA library design.

    Parameters
    ----------
    variant_key (dict): Dictionary linking mutated allele to GRCh38 location.
    SNV_white_list (list): List of variant key dictionaries with combinations of GRCh38Location and SNVs for all expected SNVs.
    """
    if variant_key in SNV_white_list:
        return True
    else:
        return False


def calculate_pre_and_post_frequencies(df, pre, post):
    """
    Calculate pre- and post-variant frequencies.

    Parameters
    ----------
    df (df): Dataframe to be processed.
    pre (str) : pre-sample identifier
    post (str) : post-sample identifier
    """
    df[f"{pre}_freq"] = df["count_e1_s0"] / df["count_e1_s0"].sum(axis=0)
    df[f"{post}_freq"] = df["count_e1_s1"] / df["count_e1_s1"].sum(axis=0)
    return df


def calculate_ratio_and_log2(df, pre, post):
    """
    Calculate post- over pre-frequency ratios, and log2-ratios from variant frequencies.

    Parameters
    ----------
    df (df): Dataframe to be processed.
    pre (str) : D20 sample identifier
    post (str) : D34 sample identifier
    """
    df[f"{post}_{pre}_ratio"] = df[f"{post}_freq"] / df[f"{pre}_freq"]
    df[f"{post}_{pre}_log2"] = np.where(
        df[f"{post}_{pre}_ratio"] > 0, np.log2(df[f"{post}_{pre}_ratio"]), np.nan
    )

    return df


def create_var_key_df_endogenous(df):
    """
    Generate a string identifier for each variant used for dataframe merging.

    Parameters
    ----------
    df (df): Dataframe to be processed.
    """
    df["var_key"] = df["GRCh38location"] + df["variants"]
    df["var_key"] = df["var_key"].apply(lambda list_: "_".join([str(i) for i in list_]))
    return df


def split_var_key_string(string):
    """
    Split string on _.

    Parameters
    ----------
    string (str): String to be split.
    """
    new_string = string.split("_")[0]
    return new_string


def GRCh38startlocation(df):  # this column is used in position scatterplots
    """
    Extract GRCh38 location from var_key string.

    Parameters
    ----------
    df (df): Dataframe to be processed.
    """
    df["GRCh38startlocation"] = df["var_key"].apply(split_var_key_string)
    return df


def mark_validation_variants(df, validation_list):
    """
    Mark variants selected for validation.

    Parameters
    ----------
    df (df): Dataframe to be processed.
    """
    df["validation_var"] = df["var_key"].map(
        lambda x: True if x in validation_list else False
    )
    return df


def process_variant_score_dataframe(df):
    """
    Processing and filtering dataframe containing pegRNA-derived function scores to contain mutant indices, "GRCh38startlocation", CADD consequence annotations, HAP1:PE7 function scores.

    Parameters
    ----------
    df (df): Dataframe to be processed.
    """
    df_pegRNA_score = df.copy()
    df_pegRNA_score = create_var_key_df_pegRNA_score(df_pegRNA_score)

    df_short = df_pegRNA_score[
        [
            "var_key",
            "mutant index",
            "GRCh38startlocation",
            "Consequence",
            "ClinVar_Significance_Simple",
            "CK001_percentage_editing",
            "CK002_percentage_editing",
            "HAP1:PE7_score",
        ]
    ]

    return df_short


# ----------------------------------------------------------------------------------------------------------
# process sample dataframe
# This function carries out all of the above preprocessing steps:


def process_dimsum_RData_file(
    csv_path_D19_D33,
    df_D19_D33,
    csv_path_indels,
    csv_path_rejected,
    csv_path_non_SNV,
    df_variant_score,
    save_dir,
    correct_amplicon_length,
    offset,
):
    """
    Processing individual endogenous samples with dimsum variant count dataframes as input (D33 vs D19).

    Parameters
    ----------
    csv_path (path): Path to CSV file containing variant counts on D33 and D19.
    df_D19_D33 (df): Dataframe containing variant counts on D33 and 19 to be processed.
    csv_path_D19_indels (path): Path to CSV file containing D33 and D19 indel counts.
    csv_path_D19_rejected (path): Path to CSV file containing D33 and D19 counts of rejected variants.
    csv_path_non_SNV (path): Path to CSV file containinge expected non_SNV variants.
    df_variant_score (df): Dataframe with pegRNA-derived function scores.
    save_dir (path): Path to output directory.
    correct_amplicon_length (int): Correct length of amplicon processed.
    offset (int): Offset used to translate nucleotide position in amplicon to GRCh38location.
    """
    # Import dataframes (D19/D33 indels, D19/D33 rejected)
    df_D19_D33_indels = pd.read_csv(csv_path_indels, delimiter="\t")
    df_D19_D33_indels.rename(
        columns={
            f"input_e1_s0_bNA_count": "count_e1_s0",
            "output_e1_s1_b1_count": "count_e1_s1",
        },
        inplace=True,
    )
    df_D19_D33_indels["nt_seq"] = df_D19_D33_indels["nt_seq"].str.upper()

    df_D19_D33_rejected = pd.read_csv(csv_path_rejected, delimiter="\t")
    df_D19_D33_rejected.rename(
        columns={
            f"input_e1_s0_bNA_count": "count_e1_s0",
            "output_e1_s1_b1_count": "count_e1_s1",
        },
        inplace=True,
    )
    df_D19_D33_rejected["rejected"] = True

    # Concat D19/D33 dataframes
    df_D19_D33_concat = pd.concat(
        [df_D19_D33, df_D19_D33_rejected, df_D19_D33_indels], axis=0, join="outer"
    )

    df_merge_full = df_D19_D33_concat.copy()

    # Calculate D19 and D33 frequencies
    df_merge_full = df_merge_full.pipe(
        calculate_pre_and_post_frequencies, pre="D19", post="D33"
    )

    # Split datframe into indel and no_indel dataframe
    df_merge_indels = df_merge_full.copy().loc[df_merge_full["indel_discarded"] == True]
    df_merge_indels.to_csv(
        (
            save_dir / f"{csv_path_D19_D33.name.split('bq5')[0]}data_indels.csv"
        ).as_posix()
    )
    df_merge_no_indels = df_merge_full.copy().loc[
        df_merge_full["indel_discarded"] != True
    ]
    df_merge_no_indels.to_csv(
        (
            save_dir / f"{csv_path_D19_D33.name.split('bq5')[0]}data_no_indels.csv"
        ).as_posix()
    )

    # Process indel dataframe (add ratios and log2)
    df_merge_indels = (
        df_merge_indels.pipe(calculate_ratio_and_log2, pre="D19", post="D33")
        #   .pipe(calculate_ratio_and_log2, pre = "neg", post = "D19")
    )

    # Process no_indel dataframe (add ratios and log2)
    df_merge_no_indels = (
        df_merge_no_indels.pipe(calculate_ratio_and_log2, pre="D19", post="D33")
        # .pipe(calculate_ratio_and_log2, pre = "neg", post = "D19")
    )
    df_merge_no_indels["nt_seq"] = df_merge_no_indels["nt_seq"].str.upper()

    # Split datframe into correct and incorrect amplicon length
    df_subset_correct_length = df_merge_no_indels.copy()[
        df_merge_no_indels["nt_seq"].str.len() == correct_amplicon_length
    ]
    df_subset_incorrect_length = df_merge_no_indels.copy()[
        df_merge_no_indels["nt_seq"].str.len() != correct_amplicon_length
    ]

    # Process correct length dataframe
    df_subset_correct_length = df_subset_correct_length.pipe(
        extract_variant_positions_from_dimsum_sequence
    )
    df_subset_correct_length["GRCh38location"] = df_subset_correct_length[
        "variant_position"
    ].apply(convert_variant_positions_to_GRCh38location, args=(offset,))
    df_subset_correct_length = df_subset_correct_length.pipe(
        extract_variants_from_dimsum_sequence
    ).pipe(create_variant_dictionaries_from_dimsum_sequence, offset)
    SNV_white_list = generate_SNV_white_list_new()  # Generate SNV white list

    df_subset_correct_length["SNV_in_minipool"] = df_subset_correct_length[
        "variant_dict"
    ].apply(add_white_list_annotation, args=(SNV_white_list,))
    df_subset_correct_length = df_subset_correct_length.pipe(
        create_var_key_df_endogenous
    ).pipe(  # add this point only SNVs are present
        GRCh38startlocation
    )
    # Merge correct length_dataframe with incorrect length dataframe
    df_no_indels_final = pd.concat(
        [df_subset_correct_length, df_subset_incorrect_length], axis=0, join="outer"
    )

    # Merge no_indel_dataframe with variant score dataframe
    df_variant_score = process_variant_score_dataframe(df_variant_score)

    df_no_indels_final = pd.merge(
        df_no_indels_final,
        df_variant_score,
        on=["var_key"],
        how="left",
    )

    # Concat indel and no_indel_dataframe
    df_processed = pd.concat(
        [df_no_indels_final, df_merge_indels], axis=0, join="outer"
    )

    # Add expected non_SNV annotation
    df_expected_non_SNV = pd.read_csv(csv_path_non_SNV, sep=",")
    df_expected_non_SNV = df_expected_non_SNV[["nt_seq", "mutant index"]].copy()
    df_expected_non_SNV_merge = pd.merge(
        df_expected_non_SNV,
        df_variant_score,
        on=["mutant index"],
        how="left",
    )

    df_expected_non_SNV_short = df_expected_non_SNV_merge[
        ["nt_seq", "Consequence", "GRCh38startlocation", "HAP1:PE7_score"]
    ].copy()
    df_expected_non_SNV_short["variant_in_minipool"] = True

    df_processed = pd.merge(
        df_processed,
        df_expected_non_SNV_short,
        on=["nt_seq"],
        how="left",
        suffixes=("_x", "_y"),
    )

    if "variant_in_minipool" in df_processed.columns:
        df_processed["variant_in_minipool"] = df_processed[
            "variant_in_minipool"
        ].fillna(df_processed["SNV_in_minipool"])
    else:
        df_processed["variant_in_minipool"] = df_processed["SNV_in_minipool"]

    df_processed["GRCh38startlocation"] = (
        df_processed[
            [
                "GRCh38startlocation_x",
                "GRCh38startlocation",
            ]
        ]
        .ffill(axis=1)
        .bfill(axis=1)["GRCh38startlocation_x"]
    )
    df_processed["Consequence"] = (
        df_processed[
            [
                "Consequence_x",
                "Consequence_y",
            ]
        ]
        .ffill(axis=1)
        .bfill(axis=1)["Consequence_x"]
    )
    df_processed["HAP1:PE7_score"] = (
        df_processed[
            [
                "HAP1:PE7_score_x",
                "HAP1:PE7_score_y",
            ]
        ]
        .ffill(axis=1)
        .bfill(axis=1)["HAP1:PE7_score_x"]
    )

    return df_processed


# ----------------------------------------------------------------------------------------------------------------------
# MAIN FUNCTION
# ----------------------------------------------------------------------------------------------------------------------
# Endogenous data preprocessing


def MLH1_endogenous_preprocessing(
    csv_path,
    csv_path_D19_indels,
    csv_path_D19_rejected,
    csv_path_non_SNV,
    save_dir,
    correct_amplicon_length,
    offset,
):
    """
    Processing individual endogenous samples with dimsum variant count dataframes as input (D33 vs D19).

    Parameters
    ----------
    csv_path (path): Path to CSV file containing variant counts on D33 and D19.
    csv_path_D19_indels (path): Path to CSV file containing D33 and D19 indel counts.
    csv_path_D19_rejected (path): Path to CSV file containing D33 and D19 counts of rejected variants.
    csv_path_non_SNV (path): Path to CSV file containing expected non_SNV variants.
    save_dir (path): Path to output directory.
    correct_amplicon_length (int): Correct length of amplicon processed.
    offset (int): Offset used to translate nucleotide position in amplicon to GRCh38location.
    """

    # Data processing
    df_D19_D33 = pd.read_csv(csv_path)
    df_variant_score = pd.read_csv(csv_path_function_score)

    (df_processed) = process_dimsum_RData_file(
        csv_path,
        df_D19_D33,
        csv_path_D19_indels,
        csv_path_D19_rejected,
        csv_path_non_SNV,
        df_variant_score,
        save_dir,
        correct_amplicon_length,
        offset,
    )

    return df_processed


def main(csv_paths):
    """
    Processing endogenous data across all amplicons sequenced.

    Parameters
    ----------
    csv_paths (list): List containing paths to CSV files with variant counts on D33 and D19.
    """
    save_dir = OUTPUT_DIR
    save_dir.mkdir(exist_ok=True, parents=True)

    # Data processing

    for i, csv_path in enumerate(csv_paths):
        (df) = MLH1_endogenous_preprocessing(
            csv_path,
            csv_paths_indels_u[i],
            csv_paths_rejected_u[i],
            csv_path_u2_non_SNV,
            save_dir,
            amplicon_length,
            offset_x10,
        )

        df.to_csv(
            (save_dir / f"{csv_path.name.split('bq5')[0]}data_full.csv").as_posix()
        )


# ----------------------------------------------------------------------------------------------------------------------
# Run pre-processing script on all samples

if __name__ == "__main__":
    main(csv_paths)
