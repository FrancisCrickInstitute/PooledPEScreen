"""
    Project :  Prime editing pilot screens. MLH1 intronic pegRNA analysis script processing "MLH1intronic_pegRNA_data_merged_best_pegRNA.csv". Output dataframe is: "data_MLH1intronic_pegRNA_score_replicates.csv"
    Date : 230820
    Python version 3.10

"""

# IMPORT

import pandas as pd
import pathlib as pl
import math
import numpy as np

# VARIABLES

INPUT_DIR = pl.Path("/Users/kajbac/code/230904_S002_final_filter")

OUTPUT_DIR = pl.Path("/Users/kajbac/code/240116_final_scripts_for_PE_paper/MLH1_non_coding/pegRNA_score")
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)

pegRNA_D20_D34_path = INPUT_DIR / "MLH1intronic_pegRNA_data_merged_best_pegRNA.csv"

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
   df= df.loc[df["percentage_editing_control_high"] == False]
   return df

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
        .pipe(filter_edited_controls)
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
# Best-performing pegRNA filter

def filter_best_performing_pegRNA(df, var_dict):
    """
    Set pegRNA log2-ratios and pegRNA percentage editing to nan, if not pegRNA with highest activity.

    Parameters
    ----------
    df (df): Dataframe to be processed.
    var_dict (dict): Dictionary used to translate sample collection date to sample identifier.
    """
    df.loc[df[f"{var_dict['pegRNA_D20']}_{var_dict['pegRNA_D34']}_log2_max_pegRNA"] == False, [f"{var_dict['pegRNA_D20']}_{var_dict['pegRNA_D34']}_log2",f"{var_dict['pegRNA_D20']}_percentage_editing"]] = np.nan 

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

def MLH1_intronic_pegRNA_analysis (pegRNA_D20_D34_path,variables_per_sample):
    """
    Main pre-processing function to generate a dataframe containing pegRNA scores for the MLH1 non-coding screen.

    Parameters
    ----------
    pegRNA_D20_D34_path (path): Path to input dataframe.
    variables_per_sample (list): List of dictionaries used to link sample collection date to sample identifier. 
    """
    ST_filter= f"{ST_editing_filter_PEmax}_PEmax_{ST_editing_filter_PEmaxdn}_PEmaxdn"
    save_dir = OUTPUT_DIR / f"pseudo_pre_freq_{pseudo_pre_freq_filter}" / ST_filter
    save_dir.mkdir(exist_ok=True, parents=True)
    df_pegRNA = pd.read_csv(
        pegRNA_D20_D34_path,
        delimiter=",")   

    # Data processing
    df_pegRNA_processed = df_pegRNA.copy().pipe(process_pegRNA_data_frame)
    for var_dict in variables_per_sample:
        df_pegRNA_processed = (df_pegRNA_processed.pipe(pegRNA_frequency_and_log2,var_dict["pegRNA_D20"], var_dict["pegRNA_D34"])
                                                  .pipe(filter_pegRNAs_on_pre_freq, var_dict)
                                                  .pipe(filter_pegRNAs_on_ST_editing_rates,var_dict))   

    if only_best_performing_pegRNA == True:
        for var_dict in variables_per_sample:
            df_pegRNA_processed = filter_best_performing_pegRNA(df_pegRNA_processed, var_dict)
            save_dir = OUTPUT_DIR / f"only_best_performing_pegRNA" / f"pseudo_pre_freq_{pseudo_pre_freq_filter}" / ST_filter 
            save_dir.mkdir(exist_ok=True, parents=True)   

    df_pegRNA_score = (df_pegRNA_processed.pipe(normalise_log2_scores_per_replicate)
                     .pipe(number_replicates))
  
    # Data output
    df_pegRNA_score.to_csv((save_dir / f"data_MLH1intronic_pegRNA_score_replicates.csv").as_posix())
# ----------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
   df_pegRNA_score = MLH1_intronic_pegRNA_analysis (pegRNA_D20_D34_path, variables_per_sample)