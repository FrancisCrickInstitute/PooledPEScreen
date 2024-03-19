"""
    Project :  Prime editing pilot screens. Analysis script to determine key numbers and statistics for MLH1x10 screen.
    Date : 231107
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

# PATHS

INPUT_DIR_VAR_SCORE = pl.Path(f"/Users/kajbac/code/240116_final_scripts_for_PE_paper/MLH1x10/pegRNA/variant_score/pseudo_pre_freq_{pseudo_pre_freq_filter}/{ST_editing_filter_PEmax}_PEmax_{ST_editing_filter_PEmaxdn}_PEmaxdn")
df_variant_score_path = INPUT_DIR_VAR_SCORE / "data_MLH1x10_variant_score_replicates.csv"

INPUT_DIR_ENDOGENOUS = pl.Path(f"/Users/kajbac/code/240116_final_scripts_for_PE_paper/MLH1x10/endogenous/log2_enriched")
df_endogenous_score_path = INPUT_DIR_ENDOGENOUS / "data_expected_variants_replicates_log2_filtered.csv"

INPUT_DIR_PEGRNA_NO_FILT = pl.Path(f"/Users/kajbac/code/240116_final_scripts_for_PE_paper/MLH1x10/pegRNA/pegRNA_score/pseudo_pre_freq_{pseudo_pre_freq_filter}/0_PEmax_0_PEmaxdn")
df_pegRNA_no_filt_path = INPUT_DIR_PEGRNA_NO_FILT / "data_MLH1x10_pegRNA_score_replicates.csv"


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
    Main function to print key statistics for MLH1 saturation screen.
    """

    df_variant_score = pd.read_csv(
        df_variant_score_path,
        delimiter=",") 
    
    df_endogenous = pd.read_csv(
    df_endogenous_score_path,
    delimiter=",") 

    df_pegRNA_no_filt = pd.read_csv(
        df_pegRNA_no_filt_path,
        delimiter=",") 
    
    # Determine: Number of pegRNAs encoding all variants in library, number of variants (all, SNVs, MNVs), percentage of designed SNVs over all possible SNVs in 200bp region 
    library_pegRNAs = len(df_variant_score["mutant index"]) + 58 # high in negative CTRL
    library_pegRNAs_per_variant = df_variant_score["mutant index"].value_counts()
    library_variants = df_variant_score["mutant index"].nunique()
    library_SNVs = df_variant_score.loc[df_variant_score["Type"] == "SNV"]["mutant index"].nunique()
    library_MNVs = library_variants - library_SNVs
    ratio_designed_vs_possible_SNV = library_SNVs / 600
    print(f"Number of pegRNAs in library is {library_pegRNAs}")
    print(f"Number of pegRNAs per variant in library is {library_pegRNAs_per_variant}")
    print(f"Number of variants in library is {library_variants}")
    print(f"Number of SNVs in library is {library_SNVs}")
    print(f"Number of MNVs in library is {library_MNVs}")
    print(f"Percentage of designed SNVs vs possible SNVs is {ratio_designed_vs_possible_SNV}")

    # Determine: endogenous variant frequencies across conditions
    column_list = ["var_key"]
    condition_list = ["S1", "S1O", "S2", "S2O"]
    for condition in condition_list:
     column_list.append(f"D20_freq_{condition}")
    df_var_freq = df_endogenous[column_list]
    df_var_freq.rename(columns={"D20_freq_S1": "PEmax", "D20_freq_S1O": "PEmax + obn","D20_freq_S2": "PEmax-MLH1dn", "D20_freq_S2O": "PEmax-MLH1dn + obn"},
        inplace=True)
    
    # Mean and median
    PEmax_end_mean = df_var_freq["PEmax"].mean()
    PEmaxO_end_mean = df_var_freq["PEmax + obn"].mean()
    PEmaxdn_end_mean = df_var_freq["PEmax-MLH1dn"].mean()
    PEmaxdnO_end_mean = df_var_freq["PEmax-MLH1dn + obn"].mean()
    PEmax_end_median = df_var_freq["PEmax"].median()
    PEmaxO_end_median = df_var_freq["PEmax + obn"].median()
    PEmaxdn_end_median = df_var_freq["PEmax-MLH1dn"].median()
    PEmaxdnO_end_median = df_var_freq["PEmax-MLH1dn + obn"].median()
    print(f"Mean endogenous frequency in PEmax is {PEmax_end_mean}")
    print(f"Mean endogenous frequency in PEmaxO is {PEmaxO_end_mean}")
    print(f"Mean endogenous frequency in PEmaxdn is {PEmaxdn_end_mean}")
    print(f"Mean endogenous frequency in PEmaxdnO is {PEmaxdnO_end_mean}")
    print(f"Median endogenous frequency in PEmax is {PEmax_end_median}")
    print(f"Median endogenous frequency in PEmaxO is {PEmaxO_end_median}")
    print(f"Median endogenous frequency in PEmaxdn is {PEmaxdn_end_median}")
    print(f"Median endogenous frequency in PEmaxdnO is {PEmaxdnO_end_median}")

    # Determine: ST-editing rates across conditions
    column_list_pegRNA = ["var_key"]
    for condition in D20_list:
     column_list_pegRNA.append(f"{condition}_percentage_editing")
    df_ST_editing = df_pegRNA_no_filt[column_list_pegRNA]
    df_ST_editing.rename(columns={"CK003_percentage_editing": "PEmax", "CK005_percentage_editing": "PEmax + obn","CK004_percentage_editing": "PEmax-MLH1dn", "CK006_percentage_editing": "PEmax-MLH1dn + obn"},
        inplace=True)
    
    # Mean and median
    PEmax_ST_mean = df_ST_editing["PEmax"].mean()
    PEmaxO_ST_mean = df_ST_editing["PEmax + obn"].mean()
    PEmaxdn_ST_mean = df_ST_editing["PEmax-MLH1dn"].mean()
    PEmaxdnO_ST_mean = df_ST_editing["PEmax-MLH1dn + obn"].mean()
    PEmax_ST_median = df_ST_editing["PEmax"].median()
    PEmaxO_ST_median = df_ST_editing["PEmax + obn"].median()
    PEmaxdn_ST_median = df_ST_editing["PEmax-MLH1dn"].median()
    PEmaxdnO_ST_median = df_ST_editing["PEmax-MLH1dn + obn"].median()
    print(f"Mean ST editing in PEmax is {PEmax_ST_mean}")
    print(f"Mean ST editing in PEmaxO is {PEmaxO_ST_mean}")
    print(f"Mean ST editingy in PEmaxdn is {PEmaxdn_ST_mean}")
    print(f"Mean ST editing in PEmaxdnO is {PEmaxdnO_ST_mean}")
    print(f"Median ST editingy in PEmax is {PEmax_ST_median}")
    print(f"Median ST editing in PEmaxO is {PEmaxO_ST_median}")
    print(f"Median ST editing in PEmaxdn is {PEmaxdn_ST_median}")
    print(f"Median ST editing in PEmaxdnO is {PEmaxdnO_ST_median}")

    # Percentage of pegRNAs with ST-editing > 25% in PEmax-MLH1dn-O cells
    ST_editing_all = df_ST_editing["PEmax-MLH1dn + obn"].count()
    ST_editing_over_25 = df_ST_editing.loc[df_ST_editing["PEmax-MLH1dn + obn"] > 25]["PEmax-MLH1dn + obn"].count()
    print(f"pegRNAs with ST-editing rates >25 =  {ST_editing_over_25}")
    print(f"Percentage of pegRNAs with ST-editing > 25% in PEmax-MLH1dn-O cells is {(ST_editing_over_25 / ST_editing_all)*100}")

    # Reduce dataframe to variant level
    filter_number = replicate_number -1
    df_variant_score_filtered = df_variant_score.loc[df_variant_score['number_replicates_pegRNA']>filter_number]
    scored_pegRNAs = df_variant_score_filtered["log2_mean_variant_replicates_normalised"].count()
    print(f"Number of scored pegRNAs is {scored_pegRNAs}")
    df_variant_score_filtered = df_variant_score_filtered.drop_duplicates(subset=['var_key']) # important! this dataframe now contains variant level information!

    # Determine: Number of variants scored and number of variants scoring LoF
    df_variant_score_filtered = df_variant_score_filtered.pipe(determine_fdr) # Determine fdr

    scored_variants = df_variant_score_filtered["log2_mean_variant_replicates_normalised"].count()
    LoF_variants = df_variant_score_filtered.loc[df_variant_score_filtered["fdr_0.01"] == True]["log2_mean_variant_replicates_normalised"].count()
    missense_variants = df_variant_score_filtered.loc[df_variant_score_filtered["Consequence"] == "NON_SYNONYMOUS"]["log2_mean_variant_replicates_normalised"].count()
    print(f"Number of variants scored is {scored_variants}")
    print(f"Number of missense scored is {missense_variants}")
    print(f"Number of LoF variants is {LoF_variants}")

    # Determine: Number of stopGain variants scored
    df_variant_score_filtered_stopgain = df_variant_score_filtered.loc[df_variant_score_filtered["Consequence"] == "STOP_GAINED"]
    df_variant_score_filtered_stopgain_SNV = df_variant_score_filtered_stopgain.loc[df_variant_score_filtered_stopgain["variant_type"] == "SNP"]
    df_variant_score_filtered_stopgain_MNV = df_variant_score_filtered_stopgain.loc[df_variant_score_filtered_stopgain["variant_type"] == "ONP"]
    stopGain_SNVs = df_variant_score_filtered_stopgain_SNV["variant_type"].count() 
    stopGain_MNVs = df_variant_score_filtered_stopgain_MNV["variant_type"].count()
    stopGain_SNV_mean = df_variant_score_filtered_stopgain_SNV["log2_mean_variant_replicates_normalised"].mean()
    stopGain_MNV_mean = df_variant_score_filtered_stopgain_MNV["log2_mean_variant_replicates_normalised"].mean()
    print(f"Number of stopGain SNVs is {stopGain_SNVs}")
    print(f"Number of stopGain MNVs is {stopGain_MNVs}")
    print(f"Mean score of stopGain SNVs is {stopGain_SNV_mean}")
    print(f"Mean score of stopGain MNVs is {stopGain_MNV_mean}")

    # Determine: CliVar stats
    clinvar_count = df_variant_score_filtered.dropna(subset = "ClinVar_Significance_Simple")["log2_mean_variant_replicates_normalised"].count()
    path_count_total = df_variant_score_filtered.loc[(df_variant_score_filtered["ClinVar_Significance_Simple"] == "Pathogenic") | (df_variant_score_filtered["ClinVar_Significance_Simple"] == "Likely pathogenic")]["log2_mean_variant_replicates_normalised"].count()
    benign_count_total = df_variant_score_filtered.loc[(df_variant_score_filtered["ClinVar_Significance_Simple"] == "Benign") | (df_variant_score_filtered["ClinVar_Significance_Simple"] == "Likely benign")]["log2_mean_variant_replicates_normalised"].count() 
    VUS_conflicting_count_total = df_variant_score_filtered.loc[(df_variant_score_filtered["ClinVar_Significance_Simple"] == "Uncertain") | (df_variant_score_filtered["ClinVar_Significance_Simple"] == "Conflicting")]["log2_mean_variant_replicates_normalised"].count()
    df_VUS_conflicting = df_variant_score_filtered.loc[(df_variant_score_filtered["ClinVar_Significance_Simple"] == "Uncertain") | (df_variant_score_filtered["ClinVar_Significance_Simple"] == "Conflicting")]
    VUS_conflicting_LoF_count = df_VUS_conflicting.loc[df_VUS_conflicting["fdr_0.01"] == True]["log2_mean_variant_replicates_normalised"].count()
    df_score_over_2 = df_variant_score_filtered.copy().loc[df_variant_score_filtered["log2_mean_variant_replicates_normalised"]>2]
    path_count = df_score_over_2.loc[(df_score_over_2["ClinVar_Significance_Simple"] == "Pathogenic") | (df_score_over_2["ClinVar_Significance_Simple"] == "Likely pathogenic")]["log2_mean_variant_replicates_normalised"].count()
    benign_count = df_score_over_2.loc[(df_score_over_2["ClinVar_Significance_Simple"] == "Benign") | (df_score_over_2["ClinVar_Significance_Simple"] == "Likely benign")]["log2_mean_variant_replicates_normalised"].count()
    print(f"Number of clinvar variants scored is {clinvar_count}")
    print(f"Number of pathogenic variants scored is {path_count_total}")
    print(f"Number of benign variants scored is {benign_count_total}")
    print(f"Number of VUS + conflicting variants scored is {VUS_conflicting_count_total}")
    print(f"Number of VUS and conflicting variants scoring LoF is {VUS_conflicting_LoF_count}")
    print(f"Percentage of VUS and conflicting variants scoring LoF is {VUS_conflicting_LoF_count/VUS_conflicting_count_total}")
    print(f"Number of pathogenic variants scoring > 2 is {path_count}")
    print(f"Number of benign variants scoring > 2 is {benign_count}")
    print(f"Ratio of pathogenic variants scoring > 2 over all pathogenic varaints is {path_count/path_count_total}")

    # Determine: high-stringency stats
    df_hc = df_variant_score_filtered.loc[(df_variant_score_filtered["hc"] == True)]
    hc_clinvar_count = df_hc.dropna(subset = "ClinVar_Significance_Simple")["log2_mean_variant_replicates_normalised"].count()
    hc_path_count = df_hc.loc[(df_hc["ClinVar_Significance_Simple"] == "Pathogenic") | (df_hc["ClinVar_Significance_Simple"] == "Likely pathogenic")]["log2_mean_variant_replicates_normalised"].count()
    hc_benign_count = df_hc.loc[(df_hc["ClinVar_Significance_Simple"] == "Benign") | (df_hc["ClinVar_Significance_Simple"] == "Likely benign")]["log2_mean_variant_replicates_normalised"].count()
    hc_missense_count = df_hc.loc[df_hc["Consequence"] == "NON_SYNONYMOUS"].count()
    print(f"Number of high-stringency ClinVar variants is {hc_clinvar_count}")
    print(f"Number of high-stringency pathogenic variants is {hc_path_count}")
    print(f"Number of high-stringency benign variants is {hc_benign_count}")
    print(f"Number of high-stringency missense variants is {hc_missense_count}")

    # Determine number of VUS lof in beta sheet
    df_beta_sheet = df_variant_score_filtered.loc[(df_variant_score_filtered["GRCh38startlocation"].isin(range(37017574, 37017600))) & df_variant_score_filtered["fdr_0.01"] == True]
    df_beta_sheet = df_beta_sheet.loc[df_beta_sheet["ClinVar_Significance_Simple"]== "Uncertain"]
    beta_sheet_VUS = df_beta_sheet["ClinVar_Significance_Simple"].count()
    print(f"Number of VUS in betasheet is {beta_sheet_VUS}")

# ----------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
   main ()