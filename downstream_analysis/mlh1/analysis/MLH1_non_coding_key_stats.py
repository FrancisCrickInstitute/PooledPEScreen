"""
    Project :  Prime editing pilot screens. Analysis script to determine key numbers and statistics for MLH1 non-coding screen.
    Date : 240415
    Python version 3.10

"""
# IMPORT

import pandas as pd
import pathlib as pl
import numpy as np
from scipy import stats
from scipy.stats import norm
from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay

# VARIABLES

pseudo_pre_freq_filter = 1e-4
ST_editing_filter_PEmax = 5
ST_editing_filter_PEmaxdn = 25


# Number of replicates in which a variant has to be scored.
replicate_number = 2

# Renaming dictionary
clinvar_naming_dict = {"Uncertain" : "Uncertain / Conflicting", 
               "Conflicting" : "Uncertain / Conflicting", 
               "Benign" : "Benign / Likely benign", 
               "Likely benign" : "Benign / Likely benign", 
               "Pathogenic" : "Pathogenic / Likely pathogenic",
               "Likely pathogenic" : "Pathogenic / Likely pathogenic"}
# PATHS

INPUT_DIR_VAR_SCORE = pl.Path(f"/Users/kajbac/code/240116_final_scripts_for_PE_paper/MLH1_non_coding/variant_score/pseudo_pre_freq_{pseudo_pre_freq_filter}/{ST_editing_filter_PEmax}_PEmax_{ST_editing_filter_PEmaxdn}_PEmaxdn")
df_variant_score_path = INPUT_DIR_VAR_SCORE / "data_MLH1intronic_variant_score_replicates.csv"

INPUT_DIR_PEGRNA_NO_FILT = INPUT_DIR = pl.Path(f"/Users/kajbac/code/240116_final_scripts_for_PE_paper/MLH1_non_coding/pegRNA_score/pseudo_pre_freq_{pseudo_pre_freq_filter}/0_PEmax_0_PEmaxdn")
df_pegRNA_no_filt_path = INPUT_DIR_PEGRNA_NO_FILT / "data_MLH1intronic_pegRNA_score_replicates.csv"

#--------------------------------------------------------------------------------------------------------------------
# FUNCTIONS

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

# --------------------------------------------------------------------------------------------------------------------------------------------
# from stripplot

def clinvar_stats(df):
    """
    Function to print out benchmarking statistics.
    
    Parameters
    ----------
    df (df): Dataframe to be used.
    """
    df["ClinVar_Significance_Simple"].replace(to_replace = clinvar_naming_dict, inplace = True)
    df_path = df.loc[df["ClinVar_Significance_Simple"] == "Pathogenic / Likely pathogenic"]
    df_benign = df.loc[df["ClinVar_Significance_Simple"] == "Benign / Likely benign"]
    df_VUS = df.loc[df["ClinVar_Significance_Simple"] == "Uncertain / Conflicting" ] 
    path_median = df_path["log2_mean_variant_replicates_normalised"].median()
    benign_median = df_benign["log2_mean_variant_replicates_normalised"].median()
    VUS_median = df_VUS["log2_mean_variant_replicates_normalised"].median() 
    path_std = df_path["log2_mean_variant_replicates_normalised"].std()
    benign_std = df_benign["log2_mean_variant_replicates_normalised"].std()
    VUS_std = df_VUS["log2_mean_variant_replicates_normalised"].std()
    print ("look here for statistics, median and SD")
    print(f"pathogenic_median is {path_median}")
    print(f"benign_median is {benign_median}")
    print(f"VUS_median is {VUS_median}")
    print(f"pathogenic_std is {path_std}")
    print(f"benign_std is {benign_std}")
    print(f"VUS_std is {VUS_std}")

# from SpliceAI correlation plot 
    
def calculate_fraction_of_lof_spliceAI_0_variants (df):
    """
    Function to calculate intronic variant statistics.

    Parameters
    ----------
    df (df): Dataframe to be processed.
    """
    df = df.copy().loc[df["Consequence"]=="INTRONIC"]
    
    df_spliceAI_05 = df.copy().loc[df["SpliceAI_max"]< 0.05]
    df_lof = df.copy().loc[df["fdr_0.01"]==True]
    df_lof_spliceAI_05 = df_lof.loc[df_lof["SpliceAI_max"]< 0.05]

    n_spliceAI_05 = len(df_spliceAI_05)
    n_lof = len(df_lof)
    n_lof_splice_AI_05 = len(df_lof_spliceAI_05)
    fraction = n_lof_splice_AI_05 / n_spliceAI_05

    print("look here for variants with spliceAI score up to 0.5 and fraction of lof variants")
    print(f"number of spliceAI_05 is {n_spliceAI_05}")
    print(f"number of lof is {n_lof}")
    print(f"number of lof_spliceAI_05 is {n_lof_splice_AI_05}")
    print(f"fraction of lof_spliceAI_05 is {fraction}")


#--------------------------------------------------------------------------------------------------------------------
# MAIN

def main ():
    """
    Main function to print key statistics for MLH1 non-coding screen.
    """

    df_variant_score = pd.read_csv(
        df_variant_score_path,
        delimiter=",") 

    df_pegRNA_no_filt = pd.read_csv(
        df_pegRNA_no_filt_path,
        delimiter=",") 
    
    # Determine: Number of pegRNAs encoding all variants in library, number of variants (all, SNVs, MNVs), percentage of designed SNVs over all possible SNVs in 200bp region 
    library_pegRNAs = len(df_variant_score["mutant index"]) + 78 # high in negative CTRL
    library_variants = df_variant_score["mutant index"].nunique()
    library_SNVs = df_variant_score.loc[df_variant_score["variant_type"] == "SNP"]["mutant index"].nunique()
    library_MNVs = df_variant_score.loc[df_variant_score["variant_type"] == "ONP"]["mutant index"].nunique()
    library_INS = df_variant_score.loc[df_variant_score["variant_type"] == "INS"]["mutant index"].nunique()
    library_DEL = df_variant_score.loc[df_variant_score["variant_type"] == "DEL"]["mutant index"].nunique()
    print(f"Number of pegRNAs in library is {library_pegRNAs}")
    print(f"Number of variants in library is {library_variants}")
    print(f"Number of SNVs in library is {library_SNVs}")
    print(f"Number of MNVs in library is {library_MNVs}")
    print(f"Number of INS in library is {library_INS}")
    print(f"Number of DEL in library is {library_DEL}")

    # Reduce dataframe to variant level
    filter_number = replicate_number -1
    df_variant_score_filtered = df_variant_score.loc[df_variant_score['number_replicates_pegRNA']>filter_number]
    scored_pegRNAs = df_variant_score_filtered["log2_mean_variant_replicates_normalised"].count()
    print(f"Number of scored pegRNAs is {scored_pegRNAs}")
    df_variant_score_filtered = df_variant_score_filtered.drop_duplicates(subset=['var_key']) # important! this dataframe now contains variant level information!

    # Determine: Number of variants scored, number of variants scoring LoF, and LoF threshold (max score of non-LoF variant)
    df_variant_score_filtered = df_variant_score_filtered.pipe(determine_fdr) # Determine fdr
    scored_variants = df_variant_score_filtered["log2_mean_variant_replicates_normalised"].count()
    scored_SNVs = df_variant_score_filtered.loc[df_variant_score_filtered["variant_type"] == "SNP"]["log2_mean_variant_replicates_normalised"].count()
    scored_upstream_SNVs = df_variant_score_filtered.loc[(df_variant_score_filtered["variant_type"] == "SNP") & (df_variant_score_filtered["Consequence"] == "UPSTREAM")]["log2_mean_variant_replicates_normalised"].count()
    LoF_variants = df_variant_score_filtered.loc[df_variant_score_filtered["fdr_0.01"] == True]["log2_mean_variant_replicates_normalised"].count()
    Lof_threshold = df_variant_score_filtered.loc[df_variant_score_filtered["fdr_0.01"] != True]["log2_mean_variant_replicates_normalised"].max()
    print(f"Number of variants scored is {scored_variants}")
    print(f"Number of SNVs scored is {scored_SNVs}")
    print(f"Number of upstream SNVs scored is {scored_upstream_SNVs}")
    print(f"Number of LoF variants is {LoF_variants}")
    print(f"LoF threshold is {Lof_threshold}")

    # Determine Number of all 5'UTR and upstream variants, number of 5'UTR and upstream variants scoring LoF, and ratio of LoF over all upstream variants scored
    df_upstream = df_variant_score_filtered.loc[(df_variant_score_filtered["Consequence"] == "5PRIME_UTR") | (df_variant_score_filtered["Consequence"] == "UPSTREAM") ]
    df_upstream_lof = df_upstream.loc[df_upstream["fdr_0.01"] == True] # important includes 8-bp deletion extending into exon 1, needs to be substracted from ratio
    u_count_total = df_upstream["log2_mean_variant_replicates_normalised"].count() - 1
    u_lof_count = (df_upstream_lof["log2_mean_variant_replicates_normalised"].count()) - 1
    print(f"Number of 5' UTR / upstream variants scored is {u_count_total}")
    print(f"Number of 5' UTR / upstream variants scoring LoF is {u_lof_count}")
    print(f"Ratio of 5' UTR / upstream LoF is {u_lof_count / (u_count_total)}")

    # Determine: CliVar stats, number of variants scored with ClinVar annotations (excluding NaN, and not provided), number of pathogenic, benign and VUS + conflicting scored, number of variants with scores > 2
    df_clinvar = df_variant_score_filtered.dropna(subset = "ClinVar_Significance_Simple")
    df_clinvar = df_clinvar[df_clinvar["ClinVar_Significance_Simple"] != 'not provided']
    clinvar_count = df_clinvar["log2_mean_variant_replicates_normalised"].count()
    path_count_total = df_variant_score_filtered.loc[(df_variant_score_filtered["ClinVar_Significance_Simple"] == "Pathogenic") | (df_variant_score_filtered["ClinVar_Significance_Simple"] == "Likely pathogenic")]["log2_mean_variant_replicates_normalised"].count()
    benign_count_total = df_variant_score_filtered.loc[(df_variant_score_filtered["ClinVar_Significance_Simple"] == "Benign") | (df_variant_score_filtered["ClinVar_Significance_Simple"] == "Likely benign")]["log2_mean_variant_replicates_normalised"].count() 
    VUS_conflicting_count_total = df_variant_score_filtered.loc[(df_variant_score_filtered["ClinVar_Significance_Simple"] == "Uncertain") | (df_variant_score_filtered["ClinVar_Significance_Simple"] == "Conflicting")]["log2_mean_variant_replicates_normalised"].count()
    df_score_over_2 = df_variant_score_filtered.copy().loc[df_variant_score_filtered["log2_mean_variant_replicates_normalised"]>2]
    path_count = df_score_over_2.loc[(df_score_over_2["ClinVar_Significance_Simple"] == "Pathogenic") | (df_score_over_2["ClinVar_Significance_Simple"] == "Likely pathogenic")]["log2_mean_variant_replicates_normalised"].count()
    benign_count = df_score_over_2.loc[(df_score_over_2["ClinVar_Significance_Simple"] == "Benign") | (df_score_over_2["ClinVar_Significance_Simple"] == "Likely benign")]["log2_mean_variant_replicates_normalised"].count()

    print(f"Number of clinvar variants scored is {clinvar_count}")
    print(f"Number of pathogenic variants scored is {path_count_total}")
    print(f"Number of benign variants scored is {benign_count_total}")
    print(f"Number of VUS + conflicting variants scored is {VUS_conflicting_count_total}")
    print(f"Number of pathogenic variants scoring > 2 is {path_count}")
    print(f"Number of benign variants scoring > 2 is {benign_count}")
    print(f"Ratio of pathogenic variants scoring > 2 over all pathogenic varaints is {path_count/path_count_total}")

    # clinvar_stats(df_clinvar)
    # Determine: AUC, confusion matrix, sensitivity and specificity

    # map and dataframe to be used in auc_consequence and cm_consequence calculation
    dict_clinvar_map = {'Pathogenic': 1, 'Likely pathogenic' : 1, 'Benign': 0, 'Likely benign' : 0}
    df_clinvar_auc = df_clinvar.loc[df_clinvar['ClinVar_Significance_Simple'].isin(['Pathogenic','Likely pathogenic', 'Benign', 'Likely benign'])]
    df_clinvar_auc = df_clinvar_auc.assign(expected_lof=df_clinvar_auc['ClinVar_Significance_Simple'].map(dict_clinvar_map))

    # calculate fpr, tpr, tnr, lof threshold
    auc_value_clinvar = roc_auc_score(y_true=df_clinvar_auc['expected_lof'],
                              y_score=df_clinvar_auc["log2_mean_variant_replicates_normalised"])

    df_clinvar_auc['predicted_lof'] = df_clinvar_auc['log2_mean_variant_replicates_normalised'].map(lambda x: 1 if x > Lof_threshold else 0)

    cm = confusion_matrix(y_true=df_clinvar_auc['expected_lof'],
                          y_pred=df_clinvar_auc['predicted_lof'])
    print(cm)

    # calculate sensitivity and specificity (clinsig)
    tn_cl, fp_cl, fn_cl, tp_cl = confusion_matrix(y_true=df_clinvar_auc['expected_lof'],
                          y_pred=df_clinvar_auc['predicted_lof']).ravel()

    tpr_cl = tp_cl/(tp_cl+fn_cl) 
    tnr_cl = tn_cl/(tn_cl+fp_cl)
    
    print(f"ClinVar AUC is {auc_value_clinvar}")
    print(f"sensitivity is {tpr_cl}")
    print(f"specificity is {tnr_cl}")

# ----------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
   main ()