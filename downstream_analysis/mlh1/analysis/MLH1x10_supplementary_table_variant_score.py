"""
    Project :  Prime editing pilot screens. Script to generate supplementary variant score table for MLH1 saturation screen.
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

pseudo_pre_freq_filter = 1.4e-4
ST_editing_filter_PEmax = 5
ST_editing_filter_PEmaxdn = 25

D20_list = ["CK003","CK005","CK004","CK006"]

# Number of replicates in which a variant has to be scored.
replicate_number = 2


# Columns for variant-level output dataframe
variant_score_df_columns = ['mutant index', 'gene', 'chr', 'GRCh38startlocation', 'GRCh38endlocation', 'variant_type', 'ref_allele',
                             'mut_allele', 'Consequence', 'Accession', 'ClinVar_Significance_Simple', 'PHRED', 'SpliceAI-acc-gain', 'SpliceAI-acc-loss', 'SpliceAI-don-gain', 'SpliceAI-don-loss','SpliceAI_max', 'percentage_editing_mean_variant_replicates','hc', 
                             'CK003_CK007_log2_mean_variant_normalised', 
                             'CK005_CK009_log2_mean_variant_normalised', 
                             'CK004_CK008_log2_mean_variant_normalised', 
                             'CK006_CK010_log2_mean_variant_normalised', 
                             'log2_mean_variant_replicates_normalised', 
                             'p_adjust', 
                             'fdr_0.01',
                             "D34_D20_log2_mean_replicates_normalised",]


# Dictionary to rename variant-level output dataframe
column_renaming_dict_variant_score = {'percentage_editing_mean_variant_replicates': "mean_correct_ST_editing", 
                                      "hc" :'high-stringency', 
                                      'CK003_CK007_log2_mean_variant_normalised' : "PEmax_Function_Score", 'CK005_CK009_log2_mean_variant_normalised' : "PEmax_obn_Function_Score", 'CK004_CK008_log2_mean_variant_normalised' : "PEmax_MLH1dn_Function_Score", 'CK006_CK010_log2_mean_variant_normalised' : "PEmax_MLH1dn_obn_Function_Score", 'log2_mean_variant_replicates_normalised' : "Function_Score", 
                                      "D34_D20_log2_mean_replicates_normalised" : "Endogenous_Function_Score",
                                      "p_adjust" : "q_value"}


# Dictionary to rename variant types
variant_type_renaming_dict = {"SNP" : "SNV", "ONP" : "MNV"}


# PATHS

# variant score
INPUT_DIR_VAR_SCORE = pl.Path(f"/Users/kajbac/code/240116_final_scripts_for_PE_paper/MLH1x10/pegRNA/variant_score/pseudo_pre_freq_{pseudo_pre_freq_filter}/{ST_editing_filter_PEmax}_PEmax_{ST_editing_filter_PEmaxdn}_PEmaxdn")
df_variant_score_path = INPUT_DIR_VAR_SCORE / "data_MLH1x10_variant_score_replicates.csv"

# endogenous variant score
INPUT_DIR_END_SCORE = pl.Path(
    f"/Users/kajbac/code/240116_final_scripts_for_PE_paper/MLH1x10/endogenous/log2_enriched"
)
df_endogenous_pegRNA_concat_path = (
    INPUT_DIR_END_SCORE / "data_expected_variants_replicates_log2_filtered.csv"
)

VEP_path = pl.Path(
    f"/Users/kajbac/code/240116_final_scripts_for_PE_paper/MLH1x10_vep_output.txt")

# output
OUTPUT_DIR = pl.Path("/Users/kajbac/code/240116_final_scripts_for_PE_paper/supplementary_tables")
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)

#--------------------------------------------------------------------------------------------------------------------
# FUNCTIONS

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
    
    df_endogenous_pegRNA_concat = pd.read_csv(
        df_endogenous_pegRNA_concat_path,
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
    df_endogenous_pegRNA_concat.dropna(subset = "D34_D20_log2_mean_replicates_normalised")
    df_endogenous_pegRNA_short = df_endogenous_pegRNA_concat[["var_key","D34_D20_log2_mean_replicates_normalised"]]

    df_variant_score_filtered = pd.merge(
        df_variant_score_filtered,
        df_endogenous_pegRNA_short,
        on=["var_key"],
        how="left",
    )

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
    df_variant_supp_table.to_csv((save_dir / f"Supp_table_MLH1_x10_variant_scores.csv").as_posix())

# ----------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
   main ()