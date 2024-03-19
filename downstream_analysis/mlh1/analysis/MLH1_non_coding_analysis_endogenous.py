"""
    Project :  Prime editing pilot screens. Pre-processing script for endogenous MLH1 u1/u2/x11/i15 data. Input are DIMSUM output dataframes containing variant counts. Output dataframes are: "u1_data_full.csv", "u2_data_full.csv", "x11_data_full.csv", "i15_data_full.csv". IMPORTANT: log2-enrichment filter (D20 over negative control) has not been applied yet and is applied later in plotting scripts!
    Date : 240211
    Python version 3.10

"""
# IMPORT

import pandas as pd
import pathlib as pl
import numpy as np
from scipy import stats
from scipy.stats import norm

# VARIABLES

INPUT_DIR = pl.Path("/Users/kajbac/code")

OUTPUT_DIR = pl.Path("/Users/kajbac/code/240116_final_scripts_for_PE_paper/MLH1_non_coding/endogenous/log2_enriched")
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)

u1_save_path = (OUTPUT_DIR / "data_MLH1_u1_expected_SNV.csv")
u2_save_path = (OUTPUT_DIR / "data_MLH1_u2_expected_SNV.csv")
x11_save_path = (OUTPUT_DIR / "data_MLH1_x11_expected_SNV.csv")
i15_save_path = (OUTPUT_DIR / "data_MLH1_i15_expected_SNV.csv")

save_path_list = [u1_save_path, u2_save_path, x11_save_path, i15_save_path]
save_extension_list = ["u1", "u2", "x11", "i15" ]

# Endogenous data paths (Dimsum output)
# D20 and D34 (variants)
u1_D20_D34_path = (
    INPUT_DIR / "240211_MLH1_non_coding_endogenous_pre_processing/dimsum_fitness_data/u1_bq5_fitness_replicates.csv"
)
u2_D20_D34_path = (
    INPUT_DIR / "240211_MLH1_non_coding_endogenous_pre_processing/dimsum_fitness_data/u2_bq5_fitness_replicates.csv"
)
x11_D20_D34_path = (
    INPUT_DIR / "240211_MLH1_non_coding_endogenous_pre_processing/dimsum_fitness_data/x11_bq5_fitness_replicates.csv"
)
i15_D20_D34_path = (
    INPUT_DIR / "240211_MLH1_non_coding_endogenous_pre_processing/dimsum_fitness_data/i15_bq5_fitness_replicates.csv"
)

# Negative control and D20 (variants)
u1_D20_D34_neg_path = (
    INPUT_DIR / "240211_MLH1_non_coding_endogenous_pre_processing/dimsum_fitness_data/u1_neg_bq5_fitness_replicates.csv"
)
u2_D20_D34_neg_path = (
    INPUT_DIR / "240211_MLH1_non_coding_endogenous_pre_processing/dimsum_fitness_data/u2_neg_bq5_fitness_replicates.csv"
)
x11_D20_D34_neg_path = (
    INPUT_DIR / "240211_MLH1_non_coding_endogenous_pre_processing/dimsum_fitness_data/x11_neg_bq5_fitness_replicates.csv"
)
i15_D20_D34_neg_path = (
    INPUT_DIR / "240211_MLH1_non_coding_endogenous_pre_processing/dimsum_fitness_data/i15_neg_bq5_fitness_replicates.csv"
)

# D20 and D34 (indels)
u1_D20_D34_indel_path = (
    INPUT_DIR / "240211_MLH1_non_coding_endogenous_pre_processing/dimsum_fitness_data/u1_bq5_indel_variant_data_merge.tsv"
    )
u2_D20_D34_indel_path = (
    INPUT_DIR / "240211_MLH1_non_coding_endogenous_pre_processing/dimsum_fitness_data/u2_bq5_indel_variant_data_merge.tsv"
    )
x11_D20_D34_indel_path = (
    INPUT_DIR / "240211_MLH1_non_coding_endogenous_pre_processing/dimsum_fitness_data/x11_bq5_indel_variant_data_merge.tsv"
    )
i15_D20_D34_indel_path = (
    INPUT_DIR / "240211_MLH1_non_coding_endogenous_pre_processing/dimsum_fitness_data/i15_bq5_indel_variant_data_merge.tsv"
    )

# Negative control and D20 (indels)
u1_neg_indel_path = (
    INPUT_DIR / "240211_MLH1_non_coding_endogenous_pre_processing/dimsum_fitness_data/u1_neg_bq5_indel_variant_data_merge.tsv"
)
u2_neg_indel_path = (
    INPUT_DIR / "240211_MLH1_non_coding_endogenous_pre_processing/dimsum_fitness_data/u2_neg_bq5_indel_variant_data_merge.tsv"
)
x11_neg_indel_path = (
    INPUT_DIR / "240211_MLH1_non_coding_endogenous_pre_processing/dimsum_fitness_data/x11_neg_bq5_indel_variant_data_merge.tsv"
)
i15_neg_indel_path = (
    INPUT_DIR / "240211_MLH1_non_coding_endogenous_pre_processing/dimsum_fitness_data/i15_neg_bq5_indel_variant_data_merge.tsv"
)

# D20 and D34 (rejected variants)
u1_D20_D34_rejected_path = (
    INPUT_DIR / "240211_MLH1_non_coding_endogenous_pre_processing/dimsum_fitness_data/u1_bq5_rejected_variant_data_merge.tsv"
)
u2_D20_D34_rejected_path = (
    INPUT_DIR / "240211_MLH1_non_coding_endogenous_pre_processing/dimsum_fitness_data/u2_bq5_rejected_variant_data_merge.tsv"
)
x11_D20_D34_rejected_path = (
    INPUT_DIR / "240211_MLH1_non_coding_endogenous_pre_processing/dimsum_fitness_data/x11_bq5_rejected_variant_data_merge.tsv"
)
i15_D20_D34_rejected_path = (
    INPUT_DIR / "240211_MLH1_non_coding_endogenous_pre_processing/dimsum_fitness_data/i15_bq5_rejected_variant_data_merge.tsv"
)

# Negative control and D20 (rejected variants)
u1_neg_rejected_path = (
    INPUT_DIR / "240211_MLH1_non_coding_endogenous_pre_processing/dimsum_fitness_data/u1_neg_bq5_rejected_variant_data_merge.tsv"
)
u2_neg_rejected_path = (
    INPUT_DIR / "240211_MLH1_non_coding_endogenous_pre_processing/dimsum_fitness_data/u2_neg_bq5_rejected_variant_data_merge.tsv"
)
x11_neg_rejected_path = (
    INPUT_DIR / "240211_MLH1_non_coding_endogenous_pre_processing/dimsum_fitness_data/x11_neg_bq5_rejected_variant_data_merge.tsv"
)
i15_neg_rejected_path = (
    INPUT_DIR / "240211_MLH1_non_coding_endogenous_pre_processing/dimsum_fitness_data/i15_neg_bq5_rejected_variant_data_merge.tsv"
)

# Endogenous data csv_path lists
csv_paths = [
    u1_D20_D34_path,
    u2_D20_D34_path,
    x11_D20_D34_path,
    i15_D20_D34_path,
]

csv_paths_neg = [
    u1_D20_D34_neg_path,
    u2_D20_D34_neg_path,
    x11_D20_D34_neg_path,
    i15_D20_D34_neg_path,
]

csv_paths_indels = [
    u1_D20_D34_indel_path,
    u2_D20_D34_indel_path,
    x11_D20_D34_indel_path,
    i15_D20_D34_indel_path,
]

csv_paths_neg_indels = [
    u1_neg_indel_path,
    u2_neg_indel_path,
    x11_neg_indel_path,
    i15_neg_indel_path,
]

csv_paths_rejected = [
    u1_D20_D34_rejected_path,
    u2_D20_D34_rejected_path,
    x11_D20_D34_rejected_path,
    i15_D20_D34_rejected_path,
]

csv_paths_neg_rejected = [
    u1_neg_rejected_path,
    u2_neg_rejected_path,
    x11_neg_rejected_path,
    i15_neg_rejected_path,
]

csv_path_library = (INPUT_DIR / f"230904_S002_final_filter/MLH1intronic_pegRNA_data_merged.csv")
csv_path_pegRNA = (INPUT_DIR / f"240116_final_scripts_for_PE_paper/MLH1_intronic/variant_score/pseudo_pre_freq_0.0001/5_PEmax_25_PEmaxdn/data_MLH1intronic_variant_score_replicates.csv")

csv_path_u1_non_SNV = pl.Path("/Users/kajbac/code/240215_non_coding_endogenous_all_variants/u1_INS_DEL_MNV.csv")
csv_path_u2_non_SNV = pl.Path("/Users/kajbac/code/240215_non_coding_endogenous_all_variants/u2_INS_DEL_MNV.csv")
csv_path_x11_non_SNV = pl.Path("/Users/kajbac/code/240215_non_coding_endogenous_all_variants/x11_INS_DEL_MNV.csv")
csv_path_i15_non_SNV = pl.Path("/Users/kajbac/code/240215_non_coding_endogenous_all_variants/i15_INS_DEL_MNV.csv")

csv_paths_non_SNV =[
    csv_path_u1_non_SNV,
    csv_path_u2_non_SNV,
    csv_path_x11_non_SNV,
    csv_path_i15_non_SNV
]

# variants to validate 
validation_list =["36993329_T","36993335_A", "36993545_T", "37020307_G", "37047255_T"]

# pegRNA annotation
pegRNA_score = "CK015_CK019_log2_mean_variant"

offset_u1 = 36993191  # used to translate amplicon position to GRCh38Location
offset_u2 = 36993418
offset_x11 = 37020286
offset_i15 = 37047216

offset_list = [offset_u1, offset_u2, offset_x11, offset_i15]

u1_region = [36993191, 36993385]
u2_region = [36993418, 36993586]
x11_region = [37020286, 37020482]
i15_region = [37047216, 37047386]

amplicon_length_list = [195,169,197,171]

# Color palette
consequence_colors = ["#f7a83e", "#243672", "#8acdef", "#B5A695", "#d091bf", "#784421", "#964594"]

consequence_order = [
        "canonical splice",
        "splice site",
        "intronic",
        "5' UTR",
        "3' UTR",
        "upstream",
        "downstream", 
    ]

consequence_naming_dict = {"CANONICAL_SPLICE" : "canonical splice", 
               "SPLICE_SITE" : "splice site", 
               "INTRONIC" : "intronic",
               "5PRIME_UTR" : "5' UTR", 
               "3PRIME_UTR" : "3' UTR", 
               "UPSTREAM" : "upstream",
               "DOWNSTREAM": "downstream"
               ,
               }

hue_order_scored = [False, True]

# ----------------------------------------------------------------------------------------------------------------------
#                                                    FUNCTIONS
# ----------------------------------------------------------------------------------------------------------------------
# Generate variant white lists, i.e., variants present in the pegRNA library

def translate_GRCh37locations_to_GRCh38locations(df):
    """
    Lift GRCh37 locations to GRCh38 locations.

    Parameters
    ----------
    df (df): Dataframe to be processed.
    """
    df["start"] = df["start"].astype("Int64") - 41491
    df.rename(columns={"start": "GRCh38startlocation"}, inplace=True)
    df["end"] = df["end"].astype("Int64") - 41491
    df.rename(columns={"end": "GRCh38endlocation"}, inplace=True)
    return df

def create_var_key_for_library_dataframe(df):
    """
    Generate a string identifier for each variant used for dataframe merging.

    Parameters
    ----------
    df (df): Dataframe to be processed.
    """
    df["var_key"] = df["GRCh38startlocation"].astype("string") + "_" + df["ref_allele"] + "_" + df["mut_allele"]
    return df

def generate_SNV_white_list():
    """
    Generate a "SNV white list", i.e, a list of variant key dictionaries with combinations of GRCh38Location and SNV for all expected SNVs.

    """
    df = pd.read_csv(csv_path_library, delimiter=",")
    df = translate_GRCh37locations_to_GRCh38locations(df)
    df = create_var_key_for_library_dataframe(df)

    df = df

    df_list =[]
    df_u1 = df.copy().loc[(df["GRCh38startlocation"].astype("Int64") >= u1_region[0]) & (df["GRCh38startlocation"].astype("Int64")<= u1_region[1])]
    df_u1 = df_u1.loc[df_u1["variant_type"]=="SNP"]
    df_u1 = df_u1.drop_duplicates(subset=['var_key'])
    df_list.append(df_u1)

    df_u2 = df.copy().loc[(df["GRCh38startlocation"].astype("Int64") >= u2_region[0]) & (df["GRCh38startlocation"].astype("Int64") <= u2_region[1])]
    df_u2 = df_u2.loc[df_u2["variant_type"]=="SNP"]
    df_u2 = df_u2.drop_duplicates(subset=['var_key'])
    df_list.append(df_u2)

    df_x11 = df.copy().loc[(df["GRCh38startlocation"].astype("Int64") >= x11_region[0]) & (df["GRCh38startlocation"].astype("Int64") <= x11_region[1])]
    df_x11 = df_x11.loc[df_x11["variant_type"]=="SNP"]
    df_x11 = df_x11.drop_duplicates(subset=['var_key'])
    df_list.append(df_x11)

    df_i15 = df.copy().loc[(df["GRCh38startlocation"].astype("Int64") >= i15_region[0]) & (df["GRCh38startlocation"].astype("Int64") <= i15_region[1])] 
    df_i15 = df_i15.loc[df_i15["variant_type"]=="SNP"]
    df_list.append(df_i15)

    SNV_white_list = []
    for df in df_list:
        location_list = df.GRCh38startlocation.values.tolist()
        mut_allele_list = df.mut_allele.values.tolist()
        for i, location in enumerate(location_list):
            SNV_white_list += [{location: mut_allele_list[i]}]

    return SNV_white_list

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

def mark_validation_variants (df):
    """
    Mark variants selected for validation.

    Parameters
    ----------
    df (df): Dataframe to be processed.
    """
    df["validation_var"]= df["var_key"].map(lambda x: True if x in validation_list else False)
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
# process variant score dataframe
def determine_fdr(df):
    """
    FDR testing of pegRNA-derived function scores with the following thresholds: 0.01, 0.02, 0.05 using synonymous variants as controls.
    
    Parameters
    ----------
    df (df): Dataframe to be processed.
    """
    loc, scale = stats.norm.fit(df.loc[df["Consequence"]=="INTRONIC"]["log2_mean_variant_replicates_normalised"])
    df["p_values"] = 1- stats.norm.cdf(df["log2_mean_variant_replicates_normalised"],loc,scale)
    df["p_adjust"] = stats.false_discovery_control(df["p_values"], method='bh')
    df["fdr_0.01"] = df["p_adjust"].map(lambda x: True if x < 0.01 else False)    
    df["fdr_0.02"] = df["p_adjust"].map(lambda x: True if x < 0.02 else False)  
    df["fdr_0.05"] = df["p_adjust"].map(lambda x: True if x < 0.05 else False)

    return df  


def process_variant_score_dataframe(df):
    """
    Processing and filtering dataframe containing pegRNA-derived function scores to contain mutant indices, "GRCh38startlocation", CADD consequence annotations, function scores fro PEmax cell line, normalised function scores for PEmax cell line, and FDR annotations.
    
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

    df_short = df_variants[["var_key", "mutant index", "GRCh38startlocation", "Consequence", "CK015_percentage_editing_mean_variant", "CK015_CK019_log2_mean_variant", "CK015_CK019_log2_mean_variant_normalised","log2_mean_variant_replicates_normalised","fdr_0.01"]]

    return df_short
#----------------------------------------------------------------------------------------------------------
# process sample dataframe
# This function carries out all of the above preprocessing steps:

def process_dimsum_RData_file(
        csv_path_D20_D34,
        df_D20_D34,
        csv_path_D20_indels,
        csv_path_D20_rejected,
        csv_path_neg_D20,
        csv_path_neg_D20_indels,
        csv_path_neg_D20_rejected,
        csv_path_non_SNV, 
        df_variant_score,
        save_dir,
        correct_amplicon_length,
        offset

        ):
    """
    Processing individual endogenous samples with dimsum variant count dataframes as input (D34 vs D20, D20 vs D4, D20 vs negative control).
    
    Parameters
    ----------
    csv_path_D20_D34 (path): Path to CSV file containing variant counts on D34 and D20.
    df_D20_D34 (df): Dataframe containing variant counts on D34 and D20 to be processed.
    csv_path_D20_indels (path): Path to CSV file containing D34 and D20 indel counts. 
    csv_path_D20_rejected (path): Path to CSV file containing D34 and D20 counts of rejected variants. 
    csv_path_neg_D20 (path): Path to CSV file containing D20 and negative control variant counts.
    csv_path_neg_D20_indels (path): Path to CSV file containing D20 and negative control indel counts.
    csv_path_neg_D20_rejected (path): Path to CSV file containing D20 and negative control counts of rejected variants.
    csv_path_non_SNV (path): Path to CSV file containinge expected non_SNV variants.
    df_variant_score (df): Dataframe with pegRNA-derived function scores.
    save_dir (path): Path to output directory. 
    correct_amplicon_length (int): Correct length of amplicon processed.
    offset (int): Offset used to translate nucleotide position in amplicon to GRCh38location.
    """
    # Import dataframes (D20/D34 indels, D20/D34 rejected)
    df_D20_D34_indels = pd.read_csv(csv_path_D20_indels,delimiter="\t")
    df_D20_D34_indels.rename(columns={f"input_e1_s0_bNA_count": "count_e1_s0", "output_e1_s1_b1_count" : "count_e1_s1" },
        inplace=True)
    df_D20_D34_indels["nt_seq"] = df_D20_D34_indels["nt_seq"].str.upper()     

    df_D20_D34_rejected = pd.read_csv(csv_path_D20_rejected,delimiter="\t") 
    df_D20_D34_rejected.rename(columns={f"input_e1_s0_bNA_count": "count_e1_s0", "output_e1_s1_b1_count" : "count_e1_s1"},
        inplace=True)
    df_D20_D34_rejected["rejected"] = True

    # Import dataframes (neg/D20, neg/D20 indels, neg/D20 rejected)
    df_neg_D20 = pd.read_csv(csv_path_neg_D20)
    df_neg_D20_indels = pd.read_csv(csv_path_neg_D20_indels,delimiter="\t")
    df_neg_D20_indels.rename(columns={f"input_e1_s0_bNA_count": "count_e1_s0", "output_e1_s1_b1_count" : "count_e1_s1" },
        inplace=True) 
    df_neg_D20_indels["nt_seq"] = df_neg_D20_indels["nt_seq"].str.upper()     

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

    # Process negative ctrl dataframe
    df_neg = process_negative_control (csv_path_neg_D20, df = df_neg_D20, df_indels = df_neg_D20_indels, df_rejected = df_neg_D20_rejected)

    # Merge dataframes (D20/D34/D4 with neg)
    df_merge_full = pd.merge(
        df_D20_D34_concat,
        df_neg,
        on=["nt_seq"],
        how="outer",
    )
    # Calculate D20 and D34 frequencies
    df_merge_full = df_merge_full.pipe(calculate_pre_and_post_frequencies, pre = "D20", post = "D34")

    # Split datframe into indel and no_indel dataframe
    df_merge_indels = df_merge_full.copy().loc[df_merge_full["indel_discarded"] == True]
    df_merge_indels.to_csv((save_dir / f"{csv_path_D20_D34.name.split('_')[0]}_data_indels.csv").as_posix())
    df_merge_no_indels = df_merge_full.copy().loc[df_merge_full["indel_discarded"] != True]
    df_merge_no_indels.to_csv((save_dir / f"{csv_path_D20_D34.name.split('_')[0]}_data_no_indels.csv").as_posix())

    # Process indel dataframe (add ratios and log2)
    df_merge_indels = (df_merge_indels.pipe(calculate_ratio_and_log2, pre = "D20", post = "D34")
                          .pipe(calculate_ratio_and_log2, pre = "neg", post = "D20")
                          )
    
    # Process no_indel dataframe (add ratios and log2)
    df_merge_no_indels = (
        df_merge_no_indels.pipe(calculate_ratio_and_log2, pre = "D20", post = "D34")
        .pipe(calculate_ratio_and_log2, pre = "neg", post = "D20"))
    df_merge_no_indels["nt_seq"] = df_merge_no_indels["nt_seq"].str.upper()

    # Split datframe into correct and incorrect amplicon length
    df_subset_correct_length = df_merge_no_indels.copy()[df_merge_no_indels["nt_seq"].str.len() == correct_amplicon_length]
    df_subset_incorrect_length = df_merge_no_indels.copy()[df_merge_no_indels["nt_seq"].str.len() != correct_amplicon_length]

    # Process correct length dataframe
    df_subset_correct_length = df_subset_correct_length.pipe(extract_variant_positions_from_dimsum_sequence)
    df_subset_correct_length["GRCh38location"] = df_subset_correct_length["variant_position"].apply(
        convert_variant_positions_to_GRCh38location, args = (offset,)
    )
    df_subset_correct_length = df_subset_correct_length.pipe(extract_variants_from_dimsum_sequence).pipe(
        create_variant_dictionaries_from_dimsum_sequence, offset
    )
    SNV_white_list = generate_SNV_white_list() # Generate SNV white list

    df_subset_correct_length["expected_SNV"] = df_subset_correct_length["variant_dict"].apply(
        add_white_list_annotation, args = (SNV_white_list,)
    )
    df_subset_correct_length = (
        df_subset_correct_length
        .pipe(create_var_key_for_dataframe_merging)
        .pipe(GRCh38startlocation)
    )
    # Merge correct length_dataframe with incorrect length dataframe
    df_no_indels_final = pd.concat([df_subset_correct_length, 
                                   df_subset_incorrect_length], 
                                   axis = 0, 
                                   join = "outer")

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
        [df_no_indels_final, 
         df_merge_indels], 
         axis = 0, 
         join = "outer")

    # Add expected non_SNV annotation
    df_expected_non_SNV = pd.read_csv(csv_path_non_SNV, sep=',')
    df_expected_non_SNV = df_expected_non_SNV[["nt_seq", "mutant index"]].copy()
    df_expected_non_SNV_merge = pd.merge(
    df_expected_non_SNV,
    df_variant_score,
    on=["mutant index"],
    how="left",
    )

    df_expected_non_SNV_short = df_expected_non_SNV_merge[["nt_seq", "Consequence", "GRCh38startlocation", "CK015_CK019_log2_mean_variant", "CK015_CK019_log2_mean_variant_normalised","fdr_0.01"]].copy()
    df_expected_non_SNV_short["expected_variant"] = True

    df_processed = pd.merge(
        df_processed,
        df_expected_non_SNV_short,
        on=["nt_seq"],
        how="left",
        suffixes=('_x', '_y')
    ) 

    if "expected_variant" in df_processed.columns:
        df_processed["expected_variant"] = df_processed["expected_variant"].fillna(df_processed["expected_SNV"])
    else:
        df_processed["expected_variant"] = df_processed["expected_SNV"]

    df_processed['GRCh38startlocation'] = df_processed[['GRCh38startlocation_x','GRCh38startlocation',]].ffill(axis=1).bfill(axis=1)['GRCh38startlocation_x']
    df_processed['Consequence'] = df_processed[['Consequence_x','Consequence_y',]].ffill(axis=1).bfill(axis=1)['Consequence_x']
    df_processed['CK015_CK019_log2_mean_variant'] = df_processed[['CK015_CK019_log2_mean_variant_x','CK015_CK019_log2_mean_variant_y',]].ffill(axis=1).bfill(axis=1)['CK015_CK019_log2_mean_variant_x']
    df_processed['CK015_CK019_log2_mean_variant_normalised'] = df_processed[['CK015_CK019_log2_mean_variant_normalised_x','CK015_CK019_log2_mean_variant_normalised_y',]].ffill(axis=1).bfill(axis=1)['CK015_CK019_log2_mean_variant_normalised_x']
    df_processed['fdr_0.01'] = df_processed[['fdr_0.01_x','fdr_0.01_y',]].ffill(axis=1).bfill(axis=1)['fdr_0.01_x']

    return df_processed

# ----------------------------------------------------------------------------------------------------------------------
# MAIN FUNCTION
# ----------------------------------------------------------------------------------------------------------------------
# Endogenous data preprocessing

def MLH1_endogenous_preprocessing(
        csv_path,
        csv_path_D20_indels,
        csv_path_D20_rejected, 
        csv_path_neg_D20,
        csv_path_neg_D20_indels,
        csv_path_neg_D20_rejected,
        csv_path_non_SNV,
        save_dir,
        correct_amplicon_length,
        offset
        ):
    """
    Processing individual endogenous samples with dimsum variant count dataframes as input (D34 vs D20, D20 vs D4, D20 vs negative control).
    
    Parameters
    ----------
    csv_path (path): Path to CSV file containing variant counts on D34 and D20.
    csv_path_D20_indels (path): Path to CSV file containing D34 and D20 indel counts. 
    csv_path_D20_rejected (path): Path to CSV file containing D34 and D20 counts of rejected variants. 
    csv_path_neg_D20 (path): Path to CSV file containing D20 and negative control variant counts.
    csv_path_neg_D20_indels (path): Path to CSV file containing D20 and negative control indel counts.
    csv_path_neg_D20_rejected (path): Path to CSV file containing D20 and negative control counts of rejected variants.
    csv_path_non_SNV (path): Path to CSV file containing expected non_SNV variants.
    save_dir (path): Path to output directory.
    correct_amplicon_length (int): Correct length of amplicon processed. 
    offset (int): Offset used to translate nucleotide position in amplicon to GRCh38location. 
    """

    # Data processing
    df_D20_D34 = pd.read_csv(csv_path)
    df_variant_score = pd.read_csv(csv_path_pegRNA)

    (df_processed
    ) = process_dimsum_RData_file(
        csv_path,
        df_D20_D34,
        csv_path_D20_indels,
        csv_path_D20_rejected,
        csv_path_neg_D20,
        csv_path_neg_D20_indels,
        csv_path_neg_D20_rejected, 
        csv_path_non_SNV,

        df_variant_score,
        save_dir,
        correct_amplicon_length,
        offset
    )

    return df_processed

def main(csv_paths):
        """
        Processing endogenous data across all amplicons sequenced.
    
        Parameters
        ----------
        csv_paths (list): List containing paths to CSV files with variant counts on D34 and D20.
        """
        save_dir = OUTPUT_DIR 
        save_dir.mkdir(exist_ok=True, parents=True)

    # Data processing

        for i, csv_path in enumerate(csv_paths):
           (df
            ) = MLH1_endogenous_preprocessing(
                csv_path, 
                csv_paths_indels [i],
                csv_paths_rejected [i], 
                csv_paths_neg[i],
                csv_paths_neg_indels[i],
                csv_paths_neg_rejected[i],
                csv_paths_non_SNV[i],
                save_dir,
                amplicon_length_list[i],
                offset_list[i]
                )

           df.to_csv((save_dir / f"{csv_path.name.split('_')[0]}_data_full.csv").as_posix())

# ----------------------------------------------------------------------------------------------------------------------
# Run pre-processing script on all samples

if __name__ == "__main__":
    main(csv_paths)    