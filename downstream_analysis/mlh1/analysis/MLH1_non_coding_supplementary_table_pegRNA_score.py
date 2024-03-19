"""
    Project :  Prime editing pilot screens. Script to generate supplementary pegRNA score table for MLH1 non-coding screen.
    Date : 240315
    Python version 3.10
"""

# IMPORT

import pandas as pd
import pathlib as pl
import math
import numpy as np

# VARIABLES

INPUT_DIR = pl.Path("/Users/kajbac/code/230904_S002_final_filter")

OUTPUT_DIR = pl.Path("/Users/kajbac/code/240116_final_scripts_for_PE_paper/supplementary_tables")
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)

pegRNA_D20_D34 = INPUT_DIR / "MLH1intronic_pegRNA_data_merged_best_pegRNA.csv"
VEP_path = pl.Path(
    f"/Users/kajbac/code/240116_final_scripts_for_PE_paper/MLH1intronic_vep_output.txt")

variables_per_sample = [
    {"pegRNA_D20": "CK015", "pegRNA_D34": "CK019"},
    {"pegRNA_D20": "CK017", "pegRNA_D34": "CK021"},
    {"pegRNA_D20": "CK016", "pegRNA_D34": "CK020"},
    {"pegRNA_D20": "CK018", "pegRNA_D34": "CK022"},
]

PEmax_list = ["CK015","CK017","CK019","CK021"]
PEmaxdn_list = ["CK016","CK018","CK020","CK022"]
D20_list = ["CK015","CK017","CK016","CK018"]

pseudo_pre_freq_filter = 1.0e-4
ST_editing_filter_PEmax = 5
ST_editing_filter_PEmaxdn = 25

# only best performing pegRNA
only_best_performing_pegRNA = False


# Columns for pegRNA design dataframe
pegRNA_score_df_columns = ['mutant index', 'gene', 'chr',  'GRCh38startlocation', 'GRCh38endlocation', 'variant_type', 'ref_allele',
                            'mut_allele', 'PAM location', 'PAM strand', 'PAM', 'protospacer', 'PBS',
                            'PBS length', 'RTT', 'RTT length', 'PBS_RTT_5to3', "distance mut to 5' RTT",
                            'distance to nick', 'PBS GC content', 'CRISPick_score',
                            'PRIDICT_editing_Score_deep', 'target sequence', 'target length', 'pegRNA_BC', 
                            'Library_ID', 'Library_BC', 'epegRNA_oligo_final', 'ref_search_term',
                            'mut_search_term', 'edited target sequence', 'Consequence',
                            'Accession', 'ClinVar_Significance_Simple',
                            
                             'CK015_percentage_editing', 'CK016_percentage_editing', 'CK017_percentage_editing',  'CK018_percentage_editing', 'CK012_percentage_editing', 
                            'percentage_editing_control_high',
                            
                            'CK015_pre_pseudo_freq', 'CK019_post_pseudo_freq', 
                            'CK017_pre_pseudo_freq', 'CK021_post_pseudo_freq', 
                            'CK016_pre_pseudo_freq', 'CK020_post_pseudo_freq',
                            'CK018_pre_pseudo_freq', 'CK022_post_pseudo_freq',
                            
                            'CK015_CK019_log2', 
                            'CK017_CK021_log2', 
                            'CK016_CK020_log2',  
                            'CK018_CK022_log2', 
                            'CK015_CK019_log2_normalised', 
                            'CK017_CK021_log2_normalised', 
                            'CK016_CK020_log2_normalised', 
                            'CK018_CK022_log2_normalised', 
                            'number_replicates']



high_in_negative_ctrl_columns = ['mutant index', 'gene', 'chr',  'GRCh38startlocation', 'GRCh38endlocation', 'variant_type', 'ref_allele',
                            'mut_allele', 'PAM location', 'PAM strand', 'PAM', 'protospacer', 'PBS',
                            'PBS length', 'RTT', 'RTT length', 'PBS_RTT_5to3', "distance mut to 5' RTT",
                            'distance to nick', 'PBS GC content', 'CRISPick_score',
                            'PRIDICT_editing_Score_deep', 'target sequence', 'target length', 'pegRNA_BC',
                            'Library_ID', 'Library_BC', 'epegRNA_oligo_final', 'ref_search_term',
                            'mut_search_term', 'edited target sequence', 'Consequence',
                            'Accession', 'ClinVar_Significance_Simple',  'CK015_percentage_editing', 'CK016_percentage_editing', 'CK017_percentage_editing',  'CK018_percentage_editing', 'CK012_percentage_editing', 
                            'percentage_editing_control_high']

# Dictionary to rename pegRNA-level output dataframe
column_renaming_dict_pegRNA_score = {'CK015_percentage_editing' : "PEmax_D20_percentage_editing", 
                            'CK016_percentage_editing' : "PEmax_MLH1dn_D20_percentage_editing", 
                            'CK017_percentage_editing' : "PEmax_obn_D20_percentage_editing", 
                            'CK018_percentage_editing' : "PEmax_MLH1dn_obn_D20_percentage_editing", 
                            'CK012_percentage_editing' : "neg_control_percentage_editing",

                            'CK015_pre_pseudo_freq' : "PEmax_D20_pseudo_freq", 
                            'CK019_post_pseudo_freq' : "PEmax_D34_pseudo_freq", 
                            'CK017_pre_pseudo_freq' : "PEmax_obn_D20_pseudo_freq", 
                            'CK021_post_pseudo_freq' : "PEmax_obn_D34_pseudo_freq", 
                            'CK016_pre_pseudo_freq' : "PEmax_MLH1dn_D20_pseudo_freq", 
                            'CK020_post_pseudo_freq' : "PEmax_MLH1dn_D34_pseudo_freq", 
                            'CK018_pre_pseudo_freq' : "PEmax_MLH1dn_obn_D20_pseudo_freq", 
                            'CK022_post_pseudo_freq' : "PEmax_MLH1dn_obn_D34_pseudo_freq", 

                            'CK015_CK019_log2' : "PEmax_pegRNA_score_unnormalized", 
                            'CK017_CK021_log2' : "PEmax_obn_pegRNA_score_unnormalized", 
                            'CK016_CK020_log2' : "PEmax_MLH1dn_pegRNA_score_unnormalized",  
                            'CK018_CK022_log2' : "PEmax_MLH1dn_obn_pegRNA_score_unnormalized", 

                            'CK015_CK019_log2_normalised' : "PEmax_pegRNA_score", 
                            'CK017_CK021_log2_normalised' : "PEmax_obn_pegRNA_score", 
                            'CK016_CK020_log2_normalised'  : "PEmax_MLH1dn_pegRNA_score", 
                            'CK018_CK022_log2_normalised' : "PEmax_MLH1dn_obn_pegRNA_score",
                            'number_replicates' : 'number_conditions'}

# Dictionary to rename variant types
variant_type_renaming_dict = {"SNP" : "SNV", "ONP" : "MNV"}
# ----------------------------------------------------------------------------------------------------------------------
#                                                     FUNCTIONS
# ----------------------------------------------------------------------------------------------------------------------

def filter_edited_controls(df):
   """
    Filter out variants with high correct ST editing in negative control.

    Parameters
    ----------
    df (df): Dataframe to be processed.
    """
   df_out = df.loc[df["percentage_editing_control_high"] == True]
   df= df.loc[df["percentage_editing_control_high"] == False]

   return df, df_out

def translate_GRCh37locations_to_GRCh38locations(df):
    """
    Lift GRCh37 locations to GRCh38 locations.

    Parameters
    ----------
    df (df): Dataframe to be processed.
    """
    df["start"] = df["start"] - 41491
    df.rename(columns={"start": "GRCh38startlocation"}, inplace=True)
    df["end"] = df["end"] - 41491
    df.rename(columns={"end": "GRCh38endlocation"}, inplace=True)
    return df

def add_var_key_for_dataframe_merging_pegRNA(df):
    """
    Generate a string identifier for each variant used for dataframe merging.

    Parameters
    ----------
    df (df): Dataframe to be processed.
    """
    df["var_key"] = (
        df["GRCh38startlocation"].astype("string")
        + "_"
        + df["mut_allele"].astype("string")
    )
    return df

def add_pseudo_count(
    df,
    pegRNA_D20,
    pegRNA_D34,
):
    """
    Add pseudocount of one to pegRNA count.

    Parameters
    ----------
    df (df): Dataframe to be processed.
    pegRNA_D20 (str) : D20 sample identifier
    pegRNA_D34 (str) : D34 sample identifier
    """
    df[f"{pegRNA_D20}_pegRNA_pseudo_count"] = (
        df[f"{pegRNA_D20}_pegRNA_count"] + 1
    )
    df[f"{pegRNA_D34}_pegRNA_pseudo_count"] = (
        df[f"{pegRNA_D34}_pegRNA_count"] + 1
    )
    return df

def calculate_frequencies_ratio_log2_pegRNA(
    df, pegRNA_D20, pegRNA_D34
):
    """
    Calculate pegRNA frequencies, post- over pre-frequency ratios, and log2-ratios from pegRNA pseudocounts.

    Parameters
    ----------
    df (df): Dataframe to be processed.
    pegRNA_D20 (str) : D20 sample identifier
    pegRNA_D34 (str) : D34 sample identifier
    """
    df[f"{pegRNA_D20}_pre_pseudo_freq"] = df[f"{pegRNA_D20}_pegRNA_pseudo_count"] / df[f"{pegRNA_D20}_pegRNA_pseudo_count"
    ].sum(axis=0)
    df[f"{pegRNA_D34}_post_pseudo_freq"] = df[f"{pegRNA_D34}_pegRNA_pseudo_count"] / df[
        f"{pegRNA_D34}_pegRNA_pseudo_count"
    ].sum(axis=0)
    df[f"{pegRNA_D34}_{pegRNA_D20}_ratio"] = df[f"{pegRNA_D34}_post_pseudo_freq"].astype("float64") / df[
        f"{pegRNA_D20}_pre_pseudo_freq"
    ].astype("float64")

    df[f"{pegRNA_D20}_{pegRNA_D34}_log2"] = df[f"{pegRNA_D34}_{pegRNA_D20}_ratio"].apply(math.log2)
    return df

def process_pegRNA_data_frame(
    df,
):
    """
    Lift GRCh37 locations to GRCh38 locations, filter out variants with high correct ST editing in negative control, and generate a string identifier for each variant used for dataframe merging.

    Parameters
    ----------
    df (df): Dataframe to be processed.
    """
    df = (
        df.pipe(translate_GRCh37locations_to_GRCh38locations)
        .pipe(add_var_key_for_dataframe_merging_pegRNA))

    return df

def pegRNA_frequency_and_log2(df,
    pegRNA_D20,
    pegRNA_D34):
    """
    Add pseudocount of one to pegRNA count, and calculate pegRNA frequencies, post- over pre-frequency ratios, and log2-ratios from pegRNA pseudocounts.

    Parameters
    ----------
    df (df): Dataframe to be processed.
    """
    df= (df
        .pipe(add_pseudo_count,pegRNA_D20, pegRNA_D34)
        .pipe(
            calculate_frequencies_ratio_log2_pegRNA, pegRNA_D20, pegRNA_D34))

    return df

# ----------------------------------------------------------------------------------------------------------------------
# Apply pre-frequency and percentage ST-editing filter

def filter_pegRNAs_on_pre_freq(
    df, var_dict
):
   """
   Set pegRNA log2-ratios and pegRNA percentage editing to nan, if pre-frequency is below specified filter value.

   Parameters
   ----------
   df (df): Dataframe to be processed.
   var_dict (dict): Dictionary used to link sample collection date to sample identifier.
   """
   df.loc[df[f"{var_dict['pegRNA_D20']}_pre_pseudo_freq"] < pseudo_pre_freq_filter, [f"{var_dict['pegRNA_D20']}_{var_dict['pegRNA_D34']}_log2",f"{var_dict['pegRNA_D20']}_percentage_editing"]] = np.nan

   return df

def filter_pegRNAs_on_ST_editing_rates(
    df, var_dict
):
   """
   Set pegRNA log2-ratios and pegRNA percentage editing to nan, if surrogate target editing rate is below specified filter value.

   Parameters
   ----------
   df (df): Dataframe to be processed.
   var_dict (dict): Dictionary used to translate sample collection date to sample identifier.
   """
   if var_dict['pegRNA_D20'] in PEmax_list:
      df.loc[df[f"{var_dict['pegRNA_D20']}_percentage_editing"] < ST_editing_filter_PEmax, [f"{var_dict['pegRNA_D20']}_{var_dict['pegRNA_D34']}_log2",f"{var_dict['pegRNA_D20']}_percentage_editing"]] = np.nan      
   elif var_dict['pegRNA_D20'] in PEmaxdn_list:
      df.loc[df[f"{var_dict['pegRNA_D20']}_percentage_editing"] < ST_editing_filter_PEmaxdn, 
             [f"{var_dict['pegRNA_D20']}_{var_dict['pegRNA_D34']}_log2",f"{var_dict['pegRNA_D20']}_percentage_editing"]] = np.nan

   return df

# ----------------------------------------------------------------------------------------------------------------------
# Calculate pegRNA-level scores

def normalise_log2_scores_per_replicate (df):
    """
    Normalise pegRNA log2-ratios to the median pegRNA log2-ratio of intronic variants for each condition.
    
    Parameters
    ----------
    df (df): Dataframe to be processed.
    """
    for var_dict in variables_per_sample:
        df[f"{var_dict['pegRNA_D20']}_{var_dict['pegRNA_D34']}_log2_normalised"] = (
        df[f"{var_dict['pegRNA_D20']}_{var_dict['pegRNA_D34']}_log2"]
        - df.loc[df["Consequence"] == "INTRONIC"][f"{var_dict['pegRNA_D20']}_{var_dict['pegRNA_D34']}_log2"].median()
    )
    return df

def number_replicates (df):
    """
    Count number of conditions, in which at least one pegRNA passed the specified pre-frequency and ST-filter.
    
    Parameters
    ----------
    df (df): Dataframe to be processed.
    """
    log2_list = []
    for var_dict in variables_per_sample:
        log2_list.append(f"{var_dict['pegRNA_D20']}_{var_dict['pegRNA_D34']}_log2")
    df["number_replicates"] = df[log2_list].count(1)

    return df

# ----------------------------------------------------------------------------------------------------------------------
# MAIN FUNCTION
# ----------------------------------------------------------------------------------------------------------------------
def pegRNA_design_supplementary_table(pegRNA_D20_D34, variables_per_sample):
    """
    Main pre-processing function to generate a dataframe containing pegRNA scores for the MLH1 non-coding variant screen.

    Parameters
    ----------
    pegRNA_D20_D34_path (path): Path to input dataframe.
    variables_per_sample (list): List of dictionaries used to link sample collection date to sample identifier. 
    """
    save_dir = OUTPUT_DIR 
    save_dir.mkdir(exist_ok=True, parents=True)
    
    df_pegRNA = pd.read_csv(
        pegRNA_D20_D34,
        delimiter=",") 
    df_VEP =  pd.read_table(
        VEP_path) 

    # Data processing
    df_pegRNA_processed = df_pegRNA.copy().pipe(process_pegRNA_data_frame)
    df_pegRNA_processed, df_edited_controls = filter_edited_controls (df_pegRNA_processed)
    for var_dict in variables_per_sample:
        df_pegRNA_processed = (df_pegRNA_processed.pipe(pegRNA_frequency_and_log2,var_dict["pegRNA_D20"], var_dict["pegRNA_D34"])
                                                  .pipe(filter_pegRNAs_on_pre_freq, var_dict)
                                                  .pipe(filter_pegRNAs_on_ST_editing_rates,var_dict))   
    df_pegRNA_processed = (df_pegRNA_processed.pipe(normalise_log2_scores_per_replicate)
                     .pipe(number_replicates))

    df_pegRNA_design_supp_table = df_pegRNA_processed[pegRNA_score_df_columns]
    df_edited_controls = df_edited_controls[high_in_negative_ctrl_columns]

    df_pegRNA_design_supp_table_final = pd.concat([df_pegRNA_design_supp_table, 
                                   df_edited_controls], 
                                   axis = 0, 
                                   join = "outer")
    df_pegRNA_design_supp_table_final.rename(columns= column_renaming_dict_pegRNA_score, inplace=True)
    df_pegRNA_design_supp_table_final['variant_type'] = df_pegRNA_design_supp_table_final['variant_type'].replace(variant_type_renaming_dict)


    # Add HGVSc column
    df_VEP_HGVSc = df_VEP[["mutant index", "HGVSc",'HGVSp']]
    df_pegRNA_design_supp_table_final = pd.merge(
        df_pegRNA_design_supp_table_final,
        df_VEP_HGVSc,
        on=["mutant index"],
        how="left",
    )
    
    print(len(df_pegRNA_design_supp_table_final))
    print(df_pegRNA_design_supp_table_final["pegRNA_BC"].nunique())
    print(list(df_pegRNA_design_supp_table_final.columns.values))
    df_pegRNA_design_supp_table_final.to_csv((save_dir / f"Supp_table_MLH1_non_coding_pegRNA_scores.csv").as_posix())

# ----------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
   pegRNA_design_supplementary_table (pegRNA_D20_D34, variables_per_sample)