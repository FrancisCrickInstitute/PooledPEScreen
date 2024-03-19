"""
    Project :  Prime editing pilot screens. MLH1x10 endogenous analysis script processing DIMSUM output dataframes containing variant counts. Only data processing: Output dataframe is: "data_expected_variants_replicates_log2_filtered.csv".
    Date : 240208
    Python version 3.10

"""
# IMPORT

import pandas as pd
import pathlib as pl
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
from scipy import stats
from scipy.stats import norm

# VARIABLES

OFFSET = 37017448  # used to translate amplicon position to GRCh38Location

INPUT_DIR = pl.Path("/Users/kajbac/code/230904_S002_final_filter")

OUTPUT_DIR = pl.Path("/Users/kajbac/code/240116_final_scripts_for_PE_paper/MLH1x10/endogenous/log2_enriched")
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)

# PegRNA analysis parameter
pseudo_pre_freq_filter = 1.4e-4
ST_editing_filter_PEmax = 5
ST_editing_filter_PEmaxdn = 25

# score list pegRNA annotation
score_list = ["CK003_CK007_log2_mean_variant","CK005_CK009_log2_mean_variant","CK004_CK008_log2_mean_variant","CK006_CK010_log2_mean_variant"]

# CADD and ClinVar path
cadd_1_path = "/Users/kajbac/code/230807_S002_MLH1_endogenous/CADD-GRCh38-v1.6-3-37017453-37017553.tsv"
cadd_2_path = "/Users/kajbac/code/230807_S002_MLH1_endogenous/CADD-GRCh38-v1.6-3-37017554-37017652.tsv"

# CADD and ClinVar dataframes
df_cadd_1 = pd.read_table(cadd_1_path)
df_cadd_2 = pd.read_table(cadd_2_path)

# PegRNA-derived function score 
variant_score_path = pl.Path(f"/Users/kajbac/code/240116_final_scripts_for_PE_paper/MLH1x10/pegRNA/variant_score/pseudo_pre_freq_{pseudo_pre_freq_filter}/{ST_editing_filter_PEmax}_PEmax_{ST_editing_filter_PEmaxdn}_PEmaxdn/data_MLH1x10_variant_score_replicates.csv")

# PegRNA-derived function score dataframe
df_variant_score = pd.read_csv(variant_score_path)

# Endogenous data paths (Dimsum output)
# D20 and D34 (variants)
sample_S1_D20_D34 = (
    INPUT_DIR / "S1_bq5_fitness_replicates.csv"
)
sample_S1_O_D20_D34 = (
    INPUT_DIR / "S1O_bq5_fitness_replicates.csv"
)
sample_S2_D20_D34 = (
    INPUT_DIR / "S2_bq5_fitness_replicates.csv"
)
sample_S2_O_D20_D34 = (
    INPUT_DIR / "S2O_bq5_fitness_replicates.csv"
)

# D4 and D20 (variants)
sample_S1_D4_D20 = (
    INPUT_DIR / "S1_bq5_D4_D20_fitness_replicates.csv"
)
sample_S1_O_D4_D20 = (
    INPUT_DIR / "S1O_bq5_D4_D20_fitness_replicates.csv"
)
sample_S2_D4_D20 = (
    INPUT_DIR / "S2_bq5_D4_D20_fitness_replicates.csv"
)
sample_S2_O_D4_D20 = (
    INPUT_DIR / "S2O_bq5_D4_D20_fitness_replicates.csv"
)

# negative control and D20 (variants)
sample_S1_neg_D20 = (
    INPUT_DIR / "S1_bq5_neg_fitness_replicates.csv"
)
sample_S1_O_neg_D20 = (
    INPUT_DIR / "S1O_bq5_neg_fitness_replicates.csv"
)
sample_S2_neg_D20 = (
    INPUT_DIR / "S2_bq5_neg_fitness_replicates.csv"
)
sample_S2_O_neg_D20 = (
    INPUT_DIR / "S2O_bq5_neg_fitness_replicates.csv"
)

# D20 and D34 (indels)
sample_S1_D20_D34_indel = (
    INPUT_DIR / "S1_D20_D34_bq5_indel_variant_data_merge.tsv"
)

sample_S1_O_D20_D34_indel = (
    INPUT_DIR / "S1O_D20_D34_bq5_indel_variant_data_merge.tsv"
)

sample_S2_D20_D34_indel = (
    INPUT_DIR / "S2_D20_D34_bq5_indel_variant_data_merge.tsv"
)

sample_S2_O_D20_D34_indel = (
    INPUT_DIR / "S2O_D20_D34_bq5_indel_variant_data_merge.tsv"
)

# D4 and D20 (indels)
sample_S1_D4_D20_indel = (
    INPUT_DIR / "S1_D4_D20_bq5_indel_variant_data_merge.tsv"
)

sample_S1_O_D4_D20_indel = (
    INPUT_DIR / "S1O_D4_D20_bq5_indel_variant_data_merge.tsv"
)

sample_S2_D4_D20_indel = (
    INPUT_DIR / "S2_D4_D20_bq5_indel_variant_data_merge.tsv"
)

sample_S2_O_D4_D20_indel = (
    INPUT_DIR / "S2O_D4_D20_bq5_indel_variant_data_merge.tsv"
)

# Negative control and D20 (indels)
sample_S1_neg_D20_indel = (
    INPUT_DIR / "S1_neg_D20_bq5_indel_variant_data_merge.tsv"
)

sample_S1_O_neg_D20_indel = (
    INPUT_DIR / "S1O_neg_D20_bq5_indel_variant_data_merge.tsv"
)

sample_S2_neg_D20_indel = (
    INPUT_DIR / "S2_neg_D20_bq5_indel_variant_data_merge.tsv"
)

sample_S2_O_neg_D20_indel = (
    INPUT_DIR / "S2O_neg_D20_bq5_indel_variant_data_merge.tsv"
)

# D20 and D34 (rejected variants)
sample_S1_D20_D34_rejected = (
    INPUT_DIR / "S1_D20_D34_bq5_rejected_variant_data_merge.tsv"
)

sample_S1_O_D20_D34_rejected = (
    INPUT_DIR / "S1O_D20_D34_bq5_rejected_variant_data_merge.tsv"
)

sample_S2_D20_D34_rejected = (
    INPUT_DIR / "S2_D20_D34_bq5_rejected_variant_data_merge.tsv"
)

sample_S2_O_D20_D34_rejected = (
    INPUT_DIR / "S2O_D20_D34_bq5_rejected_variant_data_merge.tsv"
)

# D4 and D20 (rejected variants)
sample_S1_D4_D20_rejected = (
    INPUT_DIR / "S1_D4_D20_bq5_rejected_variant_data_merge.tsv"
)

sample_S1_O_D4_D20_rejected = (
    INPUT_DIR / "S1O_D4_D20_bq5_rejected_variant_data_merge.tsv"
)

sample_S2_D4_D20_rejected = (
    INPUT_DIR / "S2_D4_D20_bq5_rejected_variant_data_merge.tsv"
)

sample_S2_O_D4_D20_rejected = (
    INPUT_DIR / "S2O_D4_D20_bq5_rejected_variant_data_merge.tsv"
)

# Negative control and D20 (rejected variants)
sample_S1_neg_D20_rejected = (
    INPUT_DIR / "S1_neg_D20_bq5_rejected_variant_data_merge.tsv"
)

sample_S1_O_neg_D20_rejected = (
    INPUT_DIR / "S1O_neg_D20_bq5_rejected_variant_data_merge.tsv"
)

sample_S2_neg_D20_rejected = (
    INPUT_DIR / "S2_neg_D20_bq5_rejected_variant_data_merge.tsv"
)

sample_S2_O_neg_D20_rejected = (
    INPUT_DIR / "S2O_neg_D20_bq5_rejected_variant_data_merge.tsv"
)

# Endogenous data csv_path lists
csv_paths = [
    sample_S1_D20_D34,
    sample_S1_O_D20_D34,
    sample_S2_D20_D34,
    sample_S2_O_D20_D34,
]

csv_paths_D4 = [
    sample_S1_D4_D20,
    sample_S1_O_D4_D20,
    sample_S2_D4_D20,
    sample_S2_O_D4_D20
]

csv_paths_neg = [
    sample_S1_neg_D20,
    sample_S1_O_neg_D20,
    sample_S2_neg_D20,
    sample_S2_O_neg_D20,
]

csv_paths_indels = [
    sample_S1_D20_D34_indel,
    sample_S1_O_D20_D34_indel,
    sample_S2_D20_D34_indel,
    sample_S2_O_D20_D34_indel
]

csv_paths_indels_D4 = [
    sample_S1_D4_D20_indel,
    sample_S1_O_D4_D20_indel,
    sample_S2_D4_D20_indel,
    sample_S2_O_D4_D20_indel
]

csv_paths_indels_neg = [
    sample_S1_neg_D20_indel,
    sample_S1_O_neg_D20_indel,
    sample_S2_neg_D20_indel,
    sample_S2_O_neg_D20_indel
]

csv_paths_rejected = [
    sample_S1_D20_D34_rejected,
    sample_S1_O_D20_D34_rejected,
    sample_S2_D20_D34_rejected,
    sample_S2_O_D20_D34_rejected
]

csv_paths_rejected_D4 = [
    sample_S1_D4_D20_rejected,
    sample_S1_O_D4_D20_rejected,
    sample_S2_D4_D20_rejected,
    sample_S2_O_D4_D20_rejected
]

csv_paths_rejected_neg = [
    sample_S1_neg_D20_rejected,
    sample_S1_O_neg_D20_rejected,
    sample_S2_neg_D20_rejected,
    sample_S2_O_neg_D20_rejected
]

#----------------------------------------------------------------------------------------------------------------------
#                                                           FUNCTIONS
#----------------------------------------------------------------------------------------------------------------------
# Create dataframe containing cadd annotations for SNV window GRCh38 37017453-37017652
# Filter out entries for RPL29P11, which have also been included when downloading cadd data for this window
# Create variant key string for data frame merging

def generate_cadd_reference(df_cadd_1, df_cadd_2):
    """
    Generating dataframe containing CADD annotations for variants in entire mutated region by combining two subset dataframes and filtering on entries for MLH1.

    Parameters
    ----------
    df_cadd_1 (df): First subset CADD dataframe.
    df_cadd_2 (df): Second subset CADD dataframe.
    """
    df_cadd = pd.concat([df_cadd_1, df_cadd_2])
    df_cadd = df_cadd.loc[df_cadd["GeneName"] == "MLH1"]
    df_cadd["var_key"] = df_cadd["Pos"].astype("string") + "_" + df_cadd["Alt"]
    assert len(df_cadd) == 600
    return df_cadd

#---------------------------------------------------------------------------------------------------------------------
# Generate "SNV white list" – variant key dictionary with combinations of GRCh38Location and SNV for all expected SNVs
# Exclude SNVs between GRCh38Location 370717503 and 37017507, and 370717630 and 37017507632 -> these are pegRNA design gaps


def generate_SNV_white_list():
    """
    Generate a "SNV white list", i.e, a liest of variant key dictionaries with combinations of GRCh38Location and SNV for all expected SNVs excluding pegRNA design gaps

    """
    wild_type_sequence = "CCTGTGACCTCACCCCTCAGGACAGTTTTGAACTGGTTGCTTTCTTTTTATTGTTTAGATCGTCTGGTAGAATCAACTTCCTTGAGAAAAGCCATAGAAACAGTGTATGCAGCCTATTTGCCCAAAAACACACACCCATTCCTGTACCTCAGGTAATGTAGCACCAAACTCCTCAACCAAGACTCACAAGGAACAGATGTTCTATCAGGCTCTCCT"
    wild_type_list = []
    for i, nt in enumerate(wild_type_sequence):
        wild_type_list += [{i + OFFSET: nt.upper()}]
    SNV_list = []
    SNV_windows = [(37017453, 37017503), (37017508, 37017630), (37017633, 37017653)]
    for window in SNV_windows:
        for nt in range(window[0], window[1]):
            SNV_list += [{nt: "A"}, {nt: "C"}, {nt: "G"}, {nt: "T"}]

    SNV_white_list = []
    for SNV_dic in SNV_list:
        if SNV_dic not in wild_type_list:
            SNV_white_list.append(SNV_dic)
    assert len(SNV_white_list) == 576

    return SNV_white_list

#---------------------------------------------------------------------------------------------------------------------
# Generate "StopGain" white list – variant key dictionary with combinations of GRCh38Location and variants for all expected StopGain variants, includes 2nt and 3nt changes
# Generate stopGain_var_key_translation_dict - dictionary for translating 2nt change var_key strings to 3nt change variant key strings, additionally 3nt variant key strings of 3nt changes introduced by stopGain pegRNAs are reformatted
# variant key strings are needed for dataframe merging, specifically for mapping pegRNAs onto endogenous variants
# -> more broadly, here I am dealing with the fact that 3nt pegRNAs can introduce 1nt, 2nt or 3nt changes in the genome depending on the target sequence.


def generate_stopGain_white_list_and_stopGain_var_key_translation_dict():
    """
    Generate "StopGain" white list, i.e, variant key dictionary with combinations of GRCh38Location and variants for all expected StopGain variants including 2nt and 3nt changes.

    """
    # stopGain_white_list
    wild_type_seq = "CGTCTGGTAGAATCAACTTCCTTGAGAAAAGCCATAGAAACAGTGTATGCAGCCTATTTGCCCAAAAACACACACCCATTCCTGTACCTC"
    T_in_position1_list = []
    A_in_position2_list = []
    A_in_position3_list = []
    for i, nt in enumerate(wild_type_seq):
        if nt == "T":
            if i % 3 == 0:
                T_in_position1_list.append(i + OFFSET + 60)
        if nt == "A":
            if (i - 1) % 3 == 0:
                A_in_position2_list.append((i - 1) + OFFSET + 60)
        if nt == "A":
            if (i - 2) % 3 == 0:
                A_in_position3_list.append((i - 2) + OFFSET + 60)

    T_in_position1_stopGain_dict = []
    A_in_position2_stopGain_dict = []
    A_in_position3_stopGain_dict = []

    for position in T_in_position1_list:
        T_in_position1_stopGain_dict += [{position + 1: "A", position + 2: "A"}]
    for position in A_in_position2_list:
        A_in_position2_stopGain_dict += [{position - 1: "T", position + 1: "A"}]
    for position in A_in_position3_list:
        A_in_position3_stopGain_dict += [{position - 2: "T", position - 1: "A"}]

    stopGain_white_list = []
    for nt in range(OFFSET + 60, OFFSET + 150, 3):
        stopGain_white_list += [{nt: "T", nt + 1: "A", nt + 2: "A"}]
    stopGain_white_list = (
        stopGain_white_list
        + T_in_position1_stopGain_dict
        + A_in_position2_stopGain_dict
        + A_in_position3_stopGain_dict
    )

    # stopGain_var_key_translation_dict:
    T_in_position1_stopGain_2nt_var_key = []
    A_in_position2_stopGain_2nt_var_key = []
    A_in_position3_stopGain_2nt_var_key = []
    TAA_wrong_format_var_key = []

    for position in T_in_position1_list:
        T_in_position1_stopGain_2nt_var_key += [f"{position+ 1}_{position + 2}_A_A"]

    for position in A_in_position2_list:
        A_in_position2_stopGain_2nt_var_key += [f"{position - 1}_{position + 1}_T_A"]

    for position in A_in_position3_list:
        A_in_position3_stopGain_2nt_var_key += [f"{position - 2}_{position - 1}_T_A"]

    for nt in range(OFFSET + 60, OFFSET + 150, 3):
        TAA_wrong_format_var_key += [f"{nt}_{nt + 1}_{nt + 2}_T_A_A"]

    T_in_position1_stopGain_3nt_var_key = []
    A_in_position2_stopGain_3nt_var_key = []
    A_in_position3_stopGain_3nt_var_key = []
    TAA_right_format_var_key = []

    for position in T_in_position1_list:
        T_in_position1_stopGain_3nt_var_key += [f"{position}_TAA"]
    for position in A_in_position2_list:
        A_in_position2_stopGain_3nt_var_key += [f"{position -1}_TAA"]
    for position in A_in_position3_list:
        A_in_position3_stopGain_3nt_var_key += [f"{position -2}_TAA"]

    for nt in range(OFFSET + 60, OFFSET + 150, 3):
        TAA_right_format_var_key += [f"{nt}_TAA"]

    two_nt_to_three_nt_variant_translation_dict_T1 = dict(
        zip(T_in_position1_stopGain_2nt_var_key, T_in_position1_stopGain_3nt_var_key)
    )

    two_nt_to_three_nt_variant_translation_dict_A2 = dict(
        zip(A_in_position2_stopGain_2nt_var_key, A_in_position2_stopGain_3nt_var_key)
    )

    two_nt_to_three_nt_variant_translation_dict_A3 = dict(
        zip(A_in_position3_stopGain_2nt_var_key, A_in_position3_stopGain_3nt_var_key)
    )

    TAA_change_format_var_key_dict = dict(
        zip(TAA_wrong_format_var_key, TAA_right_format_var_key)
    )

    stopGain_var_key_translation_dict = {
        **two_nt_to_three_nt_variant_translation_dict_T1,
        **two_nt_to_three_nt_variant_translation_dict_A2,
        **two_nt_to_three_nt_variant_translation_dict_A3,
        **TAA_change_format_var_key_dict,
    }

    print(stopGain_white_list)
    print(stopGain_var_key_translation_dict)
    return stopGain_white_list, stopGain_var_key_translation_dict

#---------------------------------------------------------------------------------------------------------------------
# Process dimsum "fitness_replicates.RData" dataframe:
# Extract variant type and location from nucleotide sequences
# Translate variant positions into GRCh38startlocations
# Add variant key dictionaries - used to iedntify expected, unexpected, and stopGain variants
# Identify expected SNVs
# Identify unexpected variants including unexpected SNVs (i.e. SNVs in the non-mutated region, SNVs for which pegRNAs have not been designed)
# Identify expected stopGain variants (3nt or 2nt changes introduced by stopGain pegRNAs)
# Translate 2nt stopGain variant key strings to 3nt stopGain variant key strings, needed for data frame merging and downstream analysis
# Add cadd and clinvar annotations (currently only SNVs)
# Add "STOP_GAINED" annotation to cadd consequence column for variants made by 3nt stopGain pegRNAs
# Calculate pre_ (D20) and post_ (D34) frequencies
# Calculate log2 fold changes
# Calculate normalised log2 fold changes (currently normalised to the median log2 fold change of synonymous variants)


def extract_variant_positions_from_nucleotide_sequence(df):
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
        for i in range(len(seq)):
            if seq[i] != wild_type_seq[i]:
                positions.append(i + 1)
        position_list.append(positions)
    df["variant_position"] = position_list
    return df


def convert_variant_positions_to_GRCh38location(positions):
    """
    Convert variant position in nucleotide sequence to GRCh38 location.

    Parameters
    ----------
    position (int): Position in nucleotide sequence to be converted.
    """
    locations = []
    for position in positions:
        locations.append(position + OFFSET - 1)
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


def create_variant_dictionaries_from_nucleotide_sequence(df):
    """
    Compare wildtype sequence to nucleotide sequence to generate a list of dictionaries, each linking mutated allele to GRCh38 location.

    Parameters
    ----------
    df (df): Dataframe to be processed.
    """
    wild_type_seq = df.loc[df["Nham_nt"] == 0, "nt_seq"].iloc[0]
    variants_dict_list = []
    sequences = df["nt_seq"].values.tolist()
    for seq in sequences:
        variants_in_seq = {}
        for i, nt in enumerate(seq):
            if nt != wild_type_seq[i]:
                variants_in_seq[i + OFFSET] = nt.upper()
        variants_dict_list.append(variants_in_seq)
    df["variant_dict"] = variants_dict_list
    return df


def add_white_list_annotation(variant_key):
    """
    Identify SNVs that are expected to be present based on pegRNA library design.

    Parameters
    ----------
    variant_key (dict): Dictionary linking mutated allele to GRCh38 location.
    """
    SNV_white_list = generate_SNV_white_list()
    (
        stopGain_white_list,
        _,
    ) = generate_stopGain_white_list_and_stopGain_var_key_translation_dict()
    variant_white_list = stopGain_white_list + SNV_white_list
    if variant_key in variant_white_list:
        return True
    else:
        return False


def identify_stopGain_variants(variant_key):
    (
        stopGain_white_list,
        _,
    ) = generate_stopGain_white_list_and_stopGain_var_key_translation_dict()
    """
    Identify SNVs and MNVs introduced by stopGain pegRNAs.

    Parameters
    ----------
    variant_key (dict): Dictionary linking mutated allele to GRCh38 location.
    """
    if variant_key in stopGain_white_list:
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
    df[f"{pre}_freq"] = df["count_e1_s0"] / df["count_e1_s0"].sum(
        axis=0
    )
    df[f"{post}_freq"] = df["count_e1_s1"] / df["count_e1_s1"].sum(
        axis=0
    )
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
    df[f"{post}_{pre}_log2"] = np.where(df[f"{post}_{pre}_ratio"] > 0, np.log2(df[f"{post}_{pre}_ratio"]), np.nan)

    return df

def create_var_key_for_dataframe_merging(df):
    """
    Generate a string identifier for each variant used for dataframe merging.

    Parameters
    ----------
    df (df): Dataframe to be processed.
    """
    df["var_key"] = df["GRCh38location"] + df["variants"]
    df["var_key"] = df["var_key"].apply(
        lambda list_: "_".join([str(i) for i in list_])
    )
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

def add_cadd_annotations(
    df, df_cadd_reference, 
):
    """
    Merge dataframe with CADD dataframe on var_key string.

    Parameters
    ----------
    df (df): Dataframe to be processed.
    """
    df = pd.merge(
        df,
        df_cadd_reference, 
        on=["var_key"], 
        how="left", 
        validate="one_to_one"
    )
    return df

def add_stopGain_peg_to_cadd_consequence(
    df,
):  
    """
    Add "STOP_GAINED" annotation to CADD consequence column for 2nt stopGain or 3nt stopGain variants introduced by stopGain pegRNAs.

    Parameters
    ----------
    df (df): Dataframe to be processed.
    """
    for i in df.index:
        additional_stopGain_variants = df.loc[i, "stopGain_peg"]
        if additional_stopGain_variants == True:
            df.loc[i, "Consequence"] = "STOP_GAINED"
    return df

def translate_stopGain_var_keys(df):
    """
    Add "STOP_GAINED" annotation to CADD consequence column for 2nt stopGain or 3nt stopGain variants introduced by stopGain pegRNAs.

    Parameters
    ----------
    df (df): Dataframe to be processed.
    """
    (
        _,
        stopGain_var_key_translation_dict,
    ) = generate_stopGain_white_list_and_stopGain_var_key_translation_dict()

    df["var_key"] = (
        df["var_key"]
        .map(stopGain_var_key_translation_dict)
        .fillna(df["var_key"])
    )
    return df


def normalise_log2_score(
    df,
):  
    """
    Normalise log2-ratios to median log2-ratio of synonymous variants.
    
    Parameters
    ----------
    df (df): Dataframe to be processed.
    """
    df["D34_D20_log2_normalised"] = (
        df["D34_D20_log2"]
        - df.loc[df["Consequence"] == "SYNONYMOUS"]["D34_D20_log2"].median()
    )
    return df

#----------------------------------------------------------------------------------------------------------
# Process D4 dataframe

def process_D4_dataframe(csv_path_D4_D20, df, df_indels, df_rejected):
    """
    Generate dataframe containing counts of variants, indels, and rejected variants on D4.
    
    Parameters
    ----------
    csv_path_D4_D20: Path to CSV file containing D4 and D20 variant counts.
    df (df): Dataframe to be processed.
    df_indels (df): Dataframe containing D4 and D20 indel counts.
    df_rejected (df): Dataframe containing D4 and D20 counts of rejected variants.
    """
    save_dir = OUTPUT_DIR / csv_path_D4_D20.name.split(".")[0]
    save_dir.mkdir(exist_ok=True, parents=True)

    df_concat = pd.concat([df, df_rejected, df_indels], axis = 0, join = "outer")
    df_concat = (df_concat.pipe(calculate_pre_and_post_frequencies, pre = "D4", post = "D20"))

    df_concat.to_csv(
        (save_dir / f"data_D4_endogenous.csv").as_posix())

    df_D4 = df_concat[["nt_seq", "count_e1_s0", "D4_freq", "rejected"]]
    df_D4.rename(
        columns={"count_e1_s0": "D4_count", "rejected": "D4_rejected"},
        inplace=True,
    )

    return df_D4

#----------------------------------------------------------------------------------------------------------
# Process negative control dataframe

def process_negative_control(csv_path_neg_D20, df, df_indels, df_rejected):
    """
    Generate dataframe containing counts of variants, indels, and rejected variants in the negative control.
    
    Parameters
    ----------
    csv_path_D4_D20: Path to CSV file containing negative control and D20 variant counts.
    df (df): Dataframe to be processed.
    df_indels (df): Dataframe containing negative control and D20 indel counts.
    df_rejected (df): Dataframe containing negative control and D20 counts of rejected variants.
    """
    save_dir = OUTPUT_DIR / csv_path_neg_D20.name.split(".")[0]
    save_dir.mkdir(exist_ok=True, parents=True)

    df_concat = pd.concat([df, df_rejected, df_indels], axis = 0, join = "outer")
    df_concat = (df_concat.pipe(calculate_pre_and_post_frequencies, pre = "neg", post = "D20"))

    df_concat.to_csv(
        (save_dir / f"data_neg_endogenous.csv").as_posix())
    
    df_neg = df_concat[["nt_seq", "count_e1_s0", "neg_freq", "rejected"]]
    df_neg.rename(
        columns={"count_e1_s0": "neg_count", "rejected": "neg_rejected"},
        inplace=True,
    )

    return df_neg

#----------------------------------------------------------------------------------------------------------
# Process variant score dataframe
def determine_fdr(df):
    """
    FDR testing of pegRNA-derived function scores with the following thresholds: 0.01, 0.02, 0.05 using synonymous variants as controls.
    
    Parameters
    ----------
    df (df): Dataframe to be processed.
    """
    loc, scale = stats.norm.fit(df.loc[df["Consequence"]=="SYNONYMOUS"]["log2_mean_variant_replicates_normalised"])
    df["p_values"] = 1- stats.norm.cdf(df["log2_mean_variant_replicates_normalised"],loc,scale)
    df["p_adjust"] = stats.false_discovery_control(df["p_values"], method='bh')
    df["fdr_0.01"] = df["p_adjust"].map(lambda x: True if x < 0.01 else False)    
    df["fdr_0.02"] = df["p_adjust"].map(lambda x: True if x < 0.02 else False)  
    df["fdr_0.05"] = df["p_adjust"].map(lambda x: True if x < 0.05 else False)

    return df  


def process_variant_score_dataframe(df):
    """
    Processing and filtering dataframe containing pegRNA-derived function scores to contain ClinVar annotations, function scores, normalised function scores, high-confidence annotations, and FDR annotations.
    
    Parameters
    ----------
    df (df): Dataframe to be processed.
    """
    # Replicate analysis
    replicate_number = 2
    filter_number = replicate_number - 1
    df_filtered = df.copy().loc[df['number_replicates_pegRNA']>filter_number]
    df_filtered = df_filtered.drop_duplicates(subset=['var_key']) # important! this dataframe now contains variant level information!
    df_fdr = df_filtered.pipe(determine_fdr)
    df_fdr = df_fdr[["var_key", "fdr_0.01"]]

    df_variants = df.copy()
    df_variants = df_variants.drop_duplicates(subset=['var_key'])
    
    df_variants = pd.merge(
        df_variants,
        df_fdr,
        on=["var_key"],
        how="left",
    )

    df_short = df_variants[["var_key", "ClinVar_Significance_Simple", "CK003_percentage_editing_mean_variant", "CK005_percentage_editing_mean_variant", "CK004_percentage_editing_mean_variant", "CK006_percentage_editing_mean_variant", "CK003_CK007_log2_mean_variant","CK005_CK009_log2_mean_variant","CK004_CK008_log2_mean_variant","CK006_CK010_log2_mean_variant", "CK003_CK007_log2_mean_variant_normalised","CK005_CK009_log2_mean_variant_normalised","CK004_CK008_log2_mean_variant_normalised","CK006_CK010_log2_mean_variant_normalised","log2_mean_variant_replicates_normalised","hc", "fdr_0.01"]]
    df_short.to_csv(OUTPUT_DIR / "df_short.csv")
    return df_short
#----------------------------------------------------------------------------------------------------------
# Process sample dataframe
# This function carries out all of the above preprocessing steps:

def process_dimsum_RData_file(

        df_D20_D34,
        csv_path_D20_indels,
        csv_path_D20_rejected,
        csv_path_D4_D20,
        csv_path_D4_D20_indels,
        csv_path_D4_D20_rejected, 
        csv_path_neg_D20,
        csv_path_neg_D20_indels,
        csv_path_neg_D20_rejected, 

        df_cadd_reference, 
        # df_clinvar_reference,
        df_variant_score,
        save_dir 

        ):
    """
    Processing individual endogenous samples with dimsum variant count dataframes as input (D34 vs D20, D20 vs D4, D20 vs negative control).
    
    Parameters
    ----------
    df_D20_D34 (df): Dataframe containing variant counts on D34 and D20 to be processed.
    csv_path_D20_indels (path): Path to CSV file containing D34 and D20 indel counts. 
    csv_path_D20_rejected (path): Path to CSV file containing D34 and D20 counts of rejected variants. 
    csv_path_D4_D20 (path): Path to CSV file containing variant counts on D20 and D4.
    csv_path_D4_D20_indels (path): Path to CSV file containing D20 and D4 indel counts.
    csv_path_D4_D20_rejected (path): Path to CSV file containing D20 and D4 counts of rejected variants.
    csv_path_neg_D20 (path): Path to CSV file containing D20 and negative control variant counts.
    csv_path_neg_D20_indels (path): Path to CSV file containing D20 and negative control indel counts.
    csv_path_neg_D20_rejected (path): Path to CSV file containing D20 and negative control counts of rejected variants.
    df_cadd_reference (df): CADD reference dataframe.
    df_clinvar_reference (df): ClinVar reference dataframe.
    df_variant_score (df): Dataframe with pegRNA-derived function scores.
    save_dir (path): Path to output directory. 
    """
    # Import dataframes (D20/D34 indels, D20/D34 rejected)

    df_D20_D34_indels = pd.read_csv(csv_path_D20_indels,delimiter="\t")
    df_D20_D34_indels.rename(columns={f"input_e1_s0_bNA_count": "count_e1_s0", "output_e1_s1_b1_count" : "count_e1_s1" },
        inplace=True)     

    df_D20_D34_rejected = pd.read_csv(csv_path_D20_rejected,delimiter="\t") 
    df_D20_D34_rejected.rename(columns={f"input_e1_s0_bNA_count": "count_e1_s0", "output_e1_s1_b1_count" : "count_e1_s1"},
        inplace=True)
    df_D20_D34_rejected["rejected"] = True
    
    # Import dataframes (D4/D20, D4/D20 indels, D4/D20 rejected)
    df_D4_D20 = pd.read_csv(csv_path_D4_D20)

    df_D4_D20_indels = pd.read_csv(csv_path_D4_D20_indels,delimiter="\t")
    df_D4_D20_indels.rename(columns={f"input_e1_s0_bNA_count": "count_e1_s0", "output_e1_s1_b1_count" : "count_e1_s1" },
        inplace=True)     

    df_D4_D20_rejected = pd.read_csv(csv_path_D4_D20_rejected,delimiter="\t")
    df_D4_D20_rejected.rename(columns={f"input_e1_s0_bNA_count": "count_e1_s0", "output_e1_s1_b1_count" : "count_e1_s1"},
        inplace=True)
    df_D4_D20_rejected["rejected"] = True

    # Import dataframes (neg/D20, neg/D20 indels, neg/D20 rejected)

    df_neg_D20 = pd.read_csv(csv_path_neg_D20)

    df_neg_D20_indels = pd.read_csv(csv_path_neg_D20_indels,delimiter="\t")
    df_neg_D20_indels.rename(columns={f"input_e1_s0_bNA_count": "count_e1_s0", "output_e1_s1_b1_count" : "count_e1_s1" },
        inplace=True)     

    df_neg_D20_rejected = pd.read_csv(csv_path_neg_D20_rejected,delimiter="\t")
    df_neg_D20_rejected.rename(columns={f"input_e1_s0_bNA_count": "count_e1_s0", "output_e1_s1_b1_count" : "count_e1_s1"},
        inplace=True)
    df_neg_D20_rejected["rejected"] = True

    # Concat D20/D34 dataframes

    df_D20_D34_concat = pd.concat(
        [df_D20_D34, 
         df_D20_D34_rejected, 
         df_D20_D34_indels], 
         axis = 0, 
         join = "outer")

    # Process D4 dataframe
    df_D4 = process_D4_dataframe (csv_path_D4_D20, df = df_D4_D20, df_indels = df_D4_D20_indels, df_rejected = df_D4_D20_rejected)

    # Process negative ctrl dataframe
    df_neg = process_negative_control (csv_path_neg_D20, df = df_neg_D20, df_indels = df_neg_D20_indels, df_rejected = df_neg_D20_rejected)

    # Merge dataframes (D20/D34 with D4)
    df_merge = pd.merge(
        df_D20_D34_concat,
        df_D4,
        on=["nt_seq"],
        how="outer",
    )

    # Merge dataframes (D20/D34/D4 with neg)
    df_merge_full = pd.merge(
        df_merge,
        df_neg,
        on=["nt_seq"],
        how="outer",
    )

    # Calculate D20 and D34 frequencies
    df_merge_full = df_merge_full.pipe(calculate_pre_and_post_frequencies, pre = "D20", post = "D34")
    # df_merge_full = df_merge_full.dropna(subset = "nt_seq")

    # Split datframe into indel and no_indel dataframe
    df_merge_indels = df_merge_full.copy().loc[df_merge_full["indel_discarded"] == True]
    df_merge_indels.to_csv((save_dir / f"{csv_path_D4_D20.name.split('_')[0]}_data_indels.csv").as_posix())
    df_merge_no_indels = df_merge_full.copy().loc[df_merge_full["indel_discarded"] != True]
    df_merge_no_indels.to_csv((save_dir / f"{csv_path_D4_D20.name.split('_')[0]}_data_no_indels.csv").as_posix())

    # Process indel dataframe (add ratios)
    df_merge_indels = (df_merge_indels.pipe(calculate_ratio_and_log2, pre = "D20", post = "D34")
                          .pipe(calculate_ratio_and_log2, pre = "D4", post = "D20")
                          .pipe(calculate_ratio_and_log2, pre = "neg", post = "D20")
                          )
    
    # Process no indel dataframe
    df_merge_no_indels = (
        df_merge_no_indels.pipe(calculate_ratio_and_log2, pre = "D20", post = "D34")
        .pipe(calculate_ratio_and_log2, pre = "D4", post = "D20")
        .pipe(calculate_ratio_and_log2, pre = "neg", post = "D20"))
    df_merge_no_indels["nt_seq"] = df_merge_no_indels["nt_seq"].str.upper()

    # Split no indel dataframe into correct and incorrect length
    correct_amplicon_length = 216
    df_subset_correct_length = df_merge_no_indels.copy()[df_merge_no_indels["nt_seq"].str.len() == correct_amplicon_length]
    df_subset_incorrect_length = df_merge_no_indels.copy()[df_merge_no_indels["nt_seq"].str.len() != correct_amplicon_length]

    # Process correct length dataframe
    df_subset_correct_length = df_subset_correct_length.pipe(extract_variant_positions_from_nucleotide_sequence)
    df_subset_correct_length["GRCh38location"] = df_subset_correct_length["variant_position"].apply(
        convert_variant_positions_to_GRCh38location
    )
    df_subset_correct_length = df_subset_correct_length.pipe(extract_variants_from_dimsum_sequence).pipe(
        create_variant_dictionaries_from_nucleotide_sequence
    )
    df_subset_correct_length["expected_variant"] = df_subset_correct_length["variant_dict"].apply(
        add_white_list_annotation
    )
    df_subset_correct_length["stopGain_peg"] = df_subset_correct_length["variant_dict"].apply(
        identify_stopGain_variants
    )
    df_subset_correct_length = (
        df_subset_correct_length
        .pipe(create_var_key_for_dataframe_merging)
        .pipe(GRCh38startlocation)
        .pipe(add_cadd_annotations, df_cadd_reference)
        .pipe(add_stopGain_peg_to_cadd_consequence)
        .pipe(translate_stopGain_var_keys)
    )
    # Merge correct length_dataframe with incorrect length dataframe
    df_no_indels_final = pd.concat([df_subset_correct_length, 
                                   df_subset_incorrect_length], 
                                   axis = 0, 
                                   join = "outer")

    # Process variant score dataframe
    df_variant_score = process_variant_score_dataframe(df_variant_score)

    # Merge no indel dataframe with variant score dataframe
    df_no_indels_final = pd.merge(
    df_no_indels_final,
    df_variant_score,
    on=["var_key"],
    how="left",
    )

    # Concat indel and no indel dataframe
    df_processed = pd.concat(
        [df_no_indels_final, 
         df_merge_indels], 
         axis = 0, 
         join = "outer")


    return df_processed

# ----------------------------------------------------------------------------------------------------------------------
#  Functions to generate final scores across replicates

def average_replicate_log2 (df):
    """
    Average log2-ratios across conditions.
    
    Parameters
    ----------
    df (df): Dataframe to be processed.
    """
    df["D34_D20_log2_mean_replicates"] = df[["D34_D20_log2_normalised_S1","D34_D20_log2_normalised_S1O","D34_D20_log2_normalised_S2","D34_D20_log2_normalised_S2O"]].mean(axis=1)

    return df

def number_replicates (df):
    """
    Count number of conditions, in which a variant received an edogenous function score.
    
    Parameters
    ----------
    df (df): Dataframe to be processed.
    """
    df["number_replicates"] = df[["D34_D20_log2_S1","D34_D20_log2_S1O","D34_D20_log2_S2","D34_D20_log2_S2O"]].count(1)

    return df

def normalise_replicate_log2_mean(
    df,
):  
    """
    Normalise endogenous function scores that have been averaged across conditions to median mean endogenous function score of synonymous variants.
    Parameters
    ----------
    df (df): Dataframe to be processed.
    """   
    # currently set to median log2 of synonymous variants
    df["D34_D20_log2_mean_replicates_normalised"] = (
        df["D34_D20_log2_mean_replicates"]
        - df.loc[df["Consequence"] == "SYNONYMOUS"]["D34_D20_log2_mean_replicates"].median()
    )
    return df
#---------------------------------------------------------------------------------------------------------------------
#                                                               MAIN FUNCTION
#---------------------------------------------------------------------------------------------------------------------
# Endogenous data pre-processing

def MLH1_x10_endogenous_preprocessing(
        csv_path,
        csv_path_D20_indels,
        csv_path_D20_rejected,
        csv_path_D4_D20,
        csv_path_D4_D20_indels,
        csv_path_D4_D20_rejected, 
        csv_path_neg_D20,
        csv_path_neg_D20_indels,
        csv_path_neg_D20_rejected,
        df_variant_score,
        save_dir
        ):
    """
    Processing individual endogenous samples with dimsum variant count dataframes as input (D34 vs D20, D20 vs D4, D20 vs negative control). Generates and adds CADD references to processed dataframes.
    
    Parameters
    ----------
    csv_path (path): Path to CSV file containing variant counts on D34 and D20.
    csv_path_D20_indels (path): Path to CSV file containing D34 and D20 indel counts. 
    csv_path_D20_rejected (path): Path to CSV file containing D34 and D20 counts of rejected variants. 
    csv_path_D4_D20 (path): Path to CSV file containing variant counts on D20 and D4.
    csv_path_D4_D20_indels (path): Path to CSV file containing D20 and D4 indel counts.
    csv_path_D4_D20_rejected (path): Path to CSV file containing D20 and D4 counts of rejected variants.
    csv_path_neg_D20 (path): Path to CSV file containing D20 and negative control variant counts.
    csv_path_neg_D20_indels (path): Path to CSV file containing D20 and negative control indel counts.
    csv_path_neg_D20_rejected (path): Path to CSV file containing D20 and negative control counts of rejected variants.
    df_variant_score (df): Dataframe with pegRNA-derived function scores.
    save_dir (path): Path to output directory. 
    """
    # Data processing
    df_D20_D34 = pd.read_csv(csv_path)
    df_cadd_reference = generate_cadd_reference(df_cadd_1, df_cadd_2)
    # df_clinvar_reference = generate_clinvar_reference(df_clinvar)

    (df_processed
    ) = process_dimsum_RData_file(
        df_D20_D34,
        csv_path_D20_indels,
        csv_path_D20_rejected,
        csv_path_D4_D20,
        csv_path_D4_D20_indels,
        csv_path_D4_D20_rejected, 
        csv_path_neg_D20,
        csv_path_neg_D20_indels,
        csv_path_neg_D20_rejected, 

        df_cadd_reference, 
        # df_clinvar_reference,
        df_variant_score,
        save_dir 
    )

    return df_processed

def MLH1x10_replicate_analyis(csv_paths, df_variant_score):
        """
        Processing endogenous data across all screen conditions.
    
        Parameters
        ----------
        csv_paths (list): List contoining paths to CSV files with variant counts on D34 and D20.
        df_variant_score (df): Dataframe with pegRNA-derived function scores.
        """
        save_dir = OUTPUT_DIR 
        save_dir.mkdir(exist_ok=True, parents=True)

    # Data processing
        MLH1x10_replicate_list = []
        for i, csv_path in enumerate(csv_paths):
           (df
            ) = MLH1_x10_endogenous_preprocessing(
                csv_path, 
                csv_paths_indels [i],
                csv_paths_rejected [i],
                csv_paths_D4[i],
                csv_paths_indels_D4[i],
                csv_paths_rejected_D4[i], 
                csv_paths_neg[i],
                csv_paths_indels_neg[i],
                csv_paths_rejected_neg[i],
                df_variant_score,
                save_dir 
                )

           df.to_csv((save_dir / f"{csv_path.name.split('_')[0]}_data_full.csv").as_posix())

           # Filter expected variants / apply D4 to D20 enrichment filter
           df_expected_variants = df.copy().loc[(df["expected_variant"] == True) & (df["D20_D4_log2"] > 1)]
           # Normalise on replicate level
           df_expected_variants = normalise_log2_score(df_expected_variants)

           df_expected_variants = df_expected_variants.rename(columns={c: c+f"_{csv_path.name.split('_')[0]}" for c in df_expected_variants.columns if c != 'var_key'})
           MLH1x10_replicate_list.append(df_expected_variants)
            

        df_concat = pd.concat(
          [df_expected_variants.set_index('var_key') for df_expected_variants in MLH1x10_replicate_list], axis=1, join='outer').reset_index()
        
        df_concat['Consequence'] = df_concat[['Consequence_S1','Consequence_S1O','Consequence_S2','Consequence_S2O']].ffill(axis=1).bfill(axis=1)['Consequence_S1']
        df_concat["ClinVar_Significance_Simple"] = df_concat[['ClinVar_Significance_Simple_S1','ClinVar_Significance_Simple_S1O','ClinVar_Significance_Simple_S2','ClinVar_Significance_Simple_S2O']].ffill(axis=1).bfill(axis=1)['ClinVar_Significance_Simple_S1']
        df_concat["GRCh38startlocation"] = df_concat[['GRCh38startlocation_S1','GRCh38startlocation_S1O','GRCh38startlocation_S2','GRCh38startlocation_S2O']].ffill(axis=1).bfill(axis=1)['GRCh38startlocation_S1']
        df_concat["log2_mean_variant_replicates_normalised"] = df_concat[['log2_mean_variant_replicates_normalised_S1','log2_mean_variant_replicates_normalised_S1O','log2_mean_variant_replicates_normalised_S2','log2_mean_variant_replicates_normalised_S2O']].ffill(axis=1).bfill(axis=1)['log2_mean_variant_replicates_normalised_S1']
        df_concat["hc"] = df_concat[['hc_S1','hc_S1O','hc_S2','hc_S2O']].ffill(axis=1).bfill(axis=1)['hc_S1']
        df_concat["fdr_0.01"] = df_concat[['fdr_0.01_S1','fdr_0.01_S1O','fdr_0.01_S2','fdr_0.01_S2O']].ffill(axis=1).bfill(axis=1)['fdr_0.01_S1']

        df_concat_filtered = (df_concat.pipe(number_replicates)
                                        .pipe(average_replicate_log2)
                                        .pipe(normalise_replicate_log2_mean))
        df_concat_filtered.to_csv(
           (save_dir / f"data_expected_variants_replicates_log2_filtered.csv").as_posix())
  
# ----------------------------------------------------------------------------------------------------------------------
# Run pre-processing script on all samples

if __name__ == "__main__":
    MLH1x10_replicate_analyis(csv_paths, df_variant_score)    