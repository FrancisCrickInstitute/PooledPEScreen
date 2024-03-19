"""
    Project :  Prime editing pilot screens. Script to generate supplementary variant score table for MLH1 non-coding screen.
    Date : 240315
    Python version 3.10

"""
# IMPORT

import pandas as pd
import pathlib as pl
import numpy as np
from scipy import stats
from scipy.stats import norm

# VARIABLES

pseudo_pre_freq_filter = 1.0e-4
ST_editing_filter_PEmax = 5
ST_editing_filter_PEmaxdn = 25

D20_list = ["CK015","CK016","CK017","CK018"]

# PEmax function score
pegRNA_score = "CK015_CK019_log2_mean_variant"

# Number of replicates in which a variant has to be scored.
replicate_number = 2

# Columns for variant-level output dataframe
variant_score_df_columns = ['mutant index', 'gene', 'chr', 'GRCh38startlocation', 'GRCh38endlocation', 'variant_type', 'ref_allele',
                             'mut_allele', 'Consequence', 'Accession', 'ClinVar_Significance_Simple', 'PHRED', 'SpliceAI-acc-gain', 'SpliceAI-acc-loss', 'SpliceAI-don-gain', 'SpliceAI-don-loss','SpliceAI_max', 
                             'CK015_CK019_log2_mean_variant_normalised', 
                             'CK017_CK021_log2_mean_variant_normalised', 
                             'CK016_CK020_log2_mean_variant_normalised', 
                             'CK018_CK022_log2_mean_variant_normalised', 
                             'log2_mean_variant_replicates_normalised', 
                             'p_adjust', 
                             'fdr_0.01',
                             "D34_D20_log2_normalised",
                             ]


# Dictionary to rename variant-level output dataframe
column_renaming_dict_variant_score = { 'CK015_CK019_log2_mean_variant_normalised' : "PEmax_Function_Score",         
                                      'CK017_CK021_log2_mean_variant_normalised' : "PEmax_obn_Function_Score", 'CK016_CK020_log2_mean_variant_normalised' : "PEmax_MLH1dn_Function_Score", 'CK018_CK022_log2_mean_variant_normalised' : "PEmax_MLH1dn_obn_Function_Score", 'log2_mean_variant_replicates_normalised' : "Function_Score", 
                                      "D34_D20_log2_normalised" : "PEmax_endogenous_Function_Score",
                                      "p_adjust" : "q_value"}


# Dictionary to rename variant types
variant_type_renaming_dict = {"SNP" : "SNV", "ONP" : "MNV"}


# PATHS

# variant score
INPUT_DIR_VAR_SCORE = pl.Path(
    f"/Users/kajbac/code//240116_final_scripts_for_PE_paper/MLH1_non_coding/variant_score/pseudo_pre_freq_{pseudo_pre_freq_filter}/{ST_editing_filter_PEmax}_PEmax_{ST_editing_filter_PEmaxdn}_PEmaxdn"
)
df_variant_score_path = INPUT_DIR_VAR_SCORE / "data_MLH1intronic_variant_score_replicates.csv"

# endogenous score
INPUT_DIR_END = pl.Path(
    f"/Users/kajbac/code/240116_final_scripts_for_PE_paper/MLH1_non_coding/endogenous/log2_enriched")

u1_path = INPUT_DIR_END / "u1_data_full.csv"
u2_path = INPUT_DIR_END / "u2_data_full.csv"
x11_path = INPUT_DIR_END / "x11_data_full.csv"
i15_path = INPUT_DIR_END / "i15_data_full.csv"

csv_paths_u = [u1_path, u2_path]
csv_paths_x11_i15 = [x11_path, i15_path]

# VEP annotaion
VEP_path = pl.Path(
    f"/Users/kajbac/code/240116_final_scripts_for_PE_paper/MLH1intronic_vep_output.txt")

# output
OUTPUT_DIR = pl.Path("/Users/kajbac/code/240116_final_scripts_for_PE_paper/supplementary_tables")
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)

# non-SNV mapping dict
non_SNV_var_key_mapping_dict = {"CGCTGAAGGGTGGGGCTGGATGGCGTAAGCTACAGCTGAAGGAAGAACGTGAGCACGAGGCACTGAGGTGATTGGCTGAAGGCACTTCCGTTGAGCATCTAGACGTTTCCTTGGCTCTTCTGGTGTCGTTCGTGGCAGGGGTTATTCGGCGGCTGGACGAG" : "36993540_G", "GCACACAAGCCCGGTTCCGGCATCTCTGCTCCTATTGGCTGGATATTTCGTATTCCCCGAGCTCCTAAAAACGAACCAATAGGAAGAGCGGACAGCGATCTCTAACGCGCAAGCGCATATCCTTCTAGGTAGCGGGCAGTAGCGCTTCAGGGAGGGACGAAGAGACCCAGCAACCCACAGAGTTGAGAAATTTG" : "36993332_G"}

#--------------------------------------------------------------------------------------------------------------------
# FUNCTIONS

def determine_fdr(df):
    """
    FDR testing of pegRNA-derived function scores with the following thresholds: 0.01, 0.02, 0.05 using intronic variants as controls.

    Parameters
    ----------
    df (df): Dataframe to be processed.
    """
    loc, scale = stats.norm.fit(
        df.loc[df["Consequence"] == "INTRONIC"][
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

#--------------------------------------------------------------------------------------------------------------------
# generate endogenous dataframe

def generate_endogenous_function_score_dataframe ():
    """
    Generate dataframe containing endogenous function scores from all expected variants (per library design) in all regions sequenced for validation in HAP1:PEmax.

    Parameters
    ----------
    df (df): Dataframe to be processed.
    """
    df_u1 = pd.read_csv(csv_paths_u[0])
    df_u2 = pd.read_csv(csv_paths_u[1])
    df_u1[["u1"]] = True
    df_u2[["u2"]] = True

    df_u = pd.concat([df_u1, df_u2], axis=0, join="outer")

    df_u["variant_scored"] = df_u[pegRNA_score].notna()
    df_u[["u"]] = True
    df_u = df_u.loc[df_u["variant_scored"] == True]
    df_u = df_u.loc[(df_u["expected_variant"] == True)]
    df_u = df_u[~(df_u["expected_SNV"] & (df_u["D20_neg_log2"] <= 1))]
    df_all_expected_variants_list = [df_u]

    df_x11 = pd.read_csv(x11_path)
    df_x11["variant_scored"] = df_x11[pegRNA_score].notna()
    df_x11["x11"] = True
    df_x11 = df_x11.loc[(df_x11["expected_variant"] == True)]
    df_x11 = df_x11[~(df_x11["expected_SNV"] & (df_x11["D20_neg_log2"] <= 1))]
    df_x11 = df_x11.loc[df_x11["variant_scored"] == True]
    df_all_expected_variants_list.append(df_x11)

    df_i15 = pd.read_csv(i15_path)
    df_i15["variant_scored"] = df_i15[pegRNA_score].notna()
    df_i15["i15"] = True
    df_i15 = df_i15.loc[(df_i15["expected_variant"] == True)]
    df_i15 = df_i15[~(df_i15["expected_SNV"] & (df_i15["D20_neg_log2"] <= 1))]
    df_i15 = df_i15.loc[df_i15["variant_scored"] == True]
    df_all_expected_variants_list.append(df_i15)

    df_concat = pd.concat([df_u, df_x11, df_i15], axis=0, join="outer")

    # normalising endogenous function scores to median score of splice neutral intronic variants
    df_concat["D34_D20_log2_normalised"] = (
        df_concat["D34_D20_log2"]
        - df_concat.loc[
            (df_concat["Consequence"] == "INTRONIC")
            & (df_concat["mutant index"] != 3128)
        ]["D34_D20_log2"].median()
    )
    return df_concat

#--------------------------------------------------------------------------------------------------------------------
# MAIN

def main ():
    """
    Main function to generate supplementary table for MLH1 saturation screen.
    """
    save_dir = OUTPUT_DIR

    # Import df
    df_variant_score = pd.read_csv(
        df_variant_score_path,
        delimiter=",") 

    df_VEP =  pd.read_table(
        VEP_path) 

    # Generate variant score supplementary table
    # Reduce variant-score dataframe to variant level
    filter_number = replicate_number -1
    df_variant_score_filtered = df_variant_score.loc[df_variant_score['number_replicates_pegRNA']>filter_number]
    df_variant_score_filtered = df_variant_score_filtered.drop_duplicates(subset=['var_key']) # important! this dataframe now contains variant level information!

    # Determine fdr
    df_variant_score_filtered = df_variant_score_filtered.pipe(determine_fdr) 

    # Add endogenous function score
    df_endogenous = generate_endogenous_function_score_dataframe()
    df_endogenous["var_key_original"] = df_endogenous["var_key"]
    df_endogenous['var_key'] = df_endogenous['nt_seq'].map(non_SNV_var_key_mapping_dict)
    df_endogenous['var_key'] = df_endogenous['var_key'].fillna(df_endogenous["var_key_original"])

    df_endogenous_short = df_endogenous[['var_key', "D34_D20_log2_normalised"]]
    print(list(df_endogenous.columns.values))
    df_variant_score_filtered = pd.merge(
        df_variant_score_filtered,
        df_endogenous_short,
        on=['var_key'],
        how="left",
    )

    # Re-name columns
    df_variant_supp_table = df_variant_score_filtered[variant_score_df_columns]
    df_variant_supp_table['variant_type'] = df_variant_supp_table['variant_type'].replace(variant_type_renaming_dict)
    df_variant_supp_table.rename(columns=column_renaming_dict_variant_score, inplace=True)

    # Add HGVSc column
    df_VEP_HGVSc = df_VEP[["mutant index", "HGVSc",'HGVSp']]
    df_variant_supp_table = pd.merge(
        df_variant_supp_table,
        df_VEP_HGVSc,
        on=["mutant index"],
        how="left",
    )
    print(list(df_variant_supp_table.columns.values))
    df_variant_supp_table.to_csv((save_dir /  f"Supp_table_MLH1_non_coding_variant_scores.csv").as_posix())

# ----------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
   main ()