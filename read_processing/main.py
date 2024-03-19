"""
[19/03/2024, Michael Herger] main.py 
Description: Main script fastq-to-pegRNA conversion and analysis of pegRNA frequencies and surrogate target editing 
efficiencies
"""

###----- IMPORTS ----------------------------------------------------------------------------------------------------###
import os
import sys
import subprocess
import time
import argparse
import pandas as pd
import math
import numpy as np
from functools import reduce
import matplotlib.pyplot as plt
import seaborn as sns
# import auxiliary functions
import surrogate_target_analysis as pooledpe_sta

###----- PARSING ----------------------------------------------------------------------------------------------------###
parser = argparse.ArgumentParser()
parser.add_argument('-s', '--samples', help='List of sample IDs, separated by space', nargs='*', type=str,
                    required=False)
args = parser.parse_args()

###----- MAIN -------------------------------------------------------------------------------------------------------###
if __name__ == '__main__':
    # ----- Define input parameters -----#
    sBaseDir = '/Users/hergerm/Documents/Data/NGS/230728_NovaSeq/'  # base directory for experiment
    sInputDir = sBaseDir + 'input/'  # input directory (experiment & pegRNA files)
    sPathToSampleDescription = sInputDir + 'sample_description.csv'  # path to sample description file
    sPathToPegRNADataFrame = sInputDir + 'pegRNA_oligos_df_asordered.csv'  # path to csv file of pegRNAs expected in this experiment (allow list)
    sSTFastqDir = sBaseDir + 'fastqs/surrogate_test/'  # directory for processed fastq files (pegRNAs & surrogate targets)
    sSTPreProcessedFastqDir = sSTFastqDir + 'processed/'  
    sPathToCutadaptBash = sSTFastqDir + 'extract_features_from_fastqs.sh'  # path to bash file for pegRNA-ST fastq pre-processing
    sDiMSumOutputDir = sBaseDir + 'input/DiMSum/'  # directory for processed fastq files (endogenous targets)
    sOutputDir = sBaseDir + 'output/' + 'MH_O1/'  # 'MH_%s/' % datetime.date.today().strftime('%y%m%d')  # output directory for analysis

    # ----- Define experiment-specific parameters ----- #
    # specify experiment to analyse: 'MLH1x10', 'MLH1intronic', 'SMARCB1', etc.
    sExperiment = 'SMARCB1'
    #  specify library BC sequences present in pegRNA oligo pool
    dict_LibBCs = {1: 'CTGTATGGCA', 2: 'GTCGGTACTC', 3: 'CAGGAAGTGT',
                   4: 'CGACTGTGAC', 5: 'GCAAGCTAGT', 6: 'GATCTGCACA',
                   7: 'CGCTCACTAG', 8: 'GGTCCATAGT', 9: 'CCTATCACTG'}
    #  read in sample description (IDs should match fastq labels)
    dict_SampleInfo = pd.read_csv(sPathToSampleDescription).set_index('ID')['description'].to_dict()
    #  specify samples to include in surrogate target analysis (can be redefined at each analysis stage)
    ls_STSamples = args.samples
    #  specify pegRNA lookup feature for search against allow list, takes either 'pegRNA-BC' or 'RTT-PBS'
    lookup_feature = 'pegRNA-BC'  # 'RTT-PBS'   # !!! implement across functions !!!
    #  specify samples to include in endogenous target analysis (can be redefined at each analysis stage)
    ls_ETSamples = []

    # ----- Define which stages of pipeline to run ----- #
    #  surrogate target analysis:
    bPreProcessSTFastqs = True  # split original fastqs into separate files for each pegRNA feature
    bGenerateEditSearchTermsAndEditedTargetSeq = True  # updates pegRNA allow list with edit search terms and edited ST
    bCreateVariantCorrespondenceTables = True  # create variant correspondence for stop gains and SNV+
    bProcessSTFastqs = True  # assign pegRNA features from fastq files to each read, filter and count pegRNAs
    bAnalyseReadFiltering = True  # analyse and visualise reads filtered for pegRNA calls
    bIncludeScaffoldAsFeature = False  # treats scaffold as pegRNA feature (e.g. scaffold libraries)
    bAnalyseSTs = True  # calculate editing efficiencies for each pegRNA
    bWriteSTFiles = False  # write ST file for each pegRNA (optional)
    bMergeSampleDFs = True  # merge all sample result dataframes
    bIntegrateExternalData = False  # merge pegRNA data with CADD, ClinVar, and PRIDICT data (optional)

    # ----- Check and set up environment ----- #
    if (sExperiment == 'MLH1x10'):
        sOutputDir = sBaseDir + 'output/MLH1x10/'
        ls_AllowedLibBCs = [8]
        sPathToCADDFile = sInputDir + 'CADD/CADD-GRCh37-v1.6-3-MLH1x10_37058944-37059143_merged.csv'
        sPathToClinVarFile = sInputDir + 'ClinVar/MLH1_clinvar_result.csv'
        if not ls_STSamples:
            ls_STSamples = ['CK001', 'CK002', 'CK003', 'CK004', 'CK005', 'CK006',
                            'CK007', 'CK008', 'CK009', 'CK010', 'CK011', 'CK023']
    elif (sExperiment == 'MLH1intronic'):
        sOutputDir = sBaseDir + 'output/MLH1intronic/'
        ls_AllowedLibBCs = [9]
        sPathToCADDFile = sInputDir + 'CADD/CADD-GRCh37-v1.6-3-MLH1-37034445-37095071_merged.csv'
        sPathToClinVarFile = sInputDir + 'ClinVar/MLH1intronic_clinvar_result_merged.csv'  # merged with output from PEGG!
        if not ls_STSamples:
            ls_STSamples = ['CK012', 'CK013', 'CK014', 'CK015', 'CK016', 'CK017',
                            'CK018', 'CK019', 'CK020', 'CK021', 'CK022', 'CK024']
    elif (sExperiment == 'SMARCB1'):
        sOutputDir = sBaseDir + 'output/SMARCB1/'
        ls_AllowedLibBCs = [1, 2, 4, 5, 6]
        sPathToCADDFile = sInputDir + 'CADD/CADD-GRCh37-v1.6-22-SMARCB1_24175725-24176488_merged.csv'
        sPathToClinVarFile = sInputDir + 'ClinVar/SMARCB1_clinvar_result.csv'
        if not ls_STSamples:
            ls_STSamples = ['MH001', 'MH002', 'MH003', 'MH004', 'MH005', 'MH006', 'MH007', 'MH008', 'MH009', 'MH010',
                            'MH011', 'MH012', 'MH013', 'MH014', 'MH015', 'MH016', 'MH017', 'MH024', 'MH025']
    elif (sExperiment == 'LV465'):  # scaffold comparison
        sOutputDir = sBaseDir + 'output/LV465/'
        ls_AllowedLibBCs = [4]
        if not ls_STSamples:
            ls_STSamples = ['MH022', 'MH023', 'MH028']
        bIncludeScaffoldAsFeature = True
    else:
        sys.exit('Invalid experiment specified: ' + sExperiment)

    if bPreProcessSTFastqs and not os.path.exists(sBaseDir + 'fastqs'):
        sys.exit('Directory of pre-processed fastq files not found at ' + sBaseDir + 'fastqs')
    if not os.path.isfile(sPathToSampleDescription):
        sys.exit('No sample description file found at ' + sPathToSampleDescription)
    if not os.path.exists(sOutputDir):
        os.makedirs(sOutputDir)

    # ----- Run fastq-to-pegRNA conversion ----- #
    if (bGenerateEditSearchTermsAndEditedTargetSeq):  # add search terms and edited target sequence to pegRNA allow List
        print('[pooledpe.py] Generating edited target sequences & search terms')
        df_pegRNAs_AllowList = pooledpe_sta.generate_pegrna_allowlist(sPathToPegRNADataFrame, ls_AllowedLibBCs, dict_LibBCs)
        df_pegRNAs_AllowList.to_csv(sOutputDir + sExperiment + '_pegRNA-allowlist.csv', index=False)

        df_Variant_SearchTerms = pooledpe_sta.generate_edit_search_terms(df_pegRNAs_AllowList)
        df_Variant_SearchTerms.to_csv(sOutputDir + sExperiment + '_variant_search-terms.csv', index=False)
        # Search terms flagged with 'search_terms_pass_qc = False' require manual curation
        # -> make changes in file and save as '..._updated.csv' then proceed with merging

        if (False):  # only set to True once search terms updated
            df_Variant_SearchTerms = pd.read_csv(sOutputDir + sExperiment + '_variant_search-terms_updated.csv')
            df_pegRNAs_AllowList_SearchTerms = pd.merge(df_pegRNAs_AllowList, df_Variant_SearchTerms, on='pegRNA_BC', how='left')
            df_pegRNAs_AllowList_SearchTerms.to_csv(sOutputDir + sExperiment + '_pegRNA-allowlist_search-terms.csv', index=False)
  
    if (bPreProcessSTFastqs):  # pre-process fastqs by splitting into pegRNA elements
        print('[pooledpe.py] Pre-processing fastq files')
        t_start = time.time()
        subprocess.call(sPathToCutadaptBash)
        t_end = time.time()
        print('[pooledpe.py] Pre-processing fastq files (done)')
        print('[pooledpe.py] Completed in ' + '{:.2f}'.format((t_end - t_start) / 60) + ' min')

    if (bProcessSTFastqs):  # assign pegRNA elements to reads and filter out forbidden combinations then determine pegRNA counts
        print('[pooledpe.py] Filter reads with forbidden pegRNA element combinations')
        df_pegRNAs_AllowList = pd.read_csv(sOutputDir + sExperiment + '_pegRNA-allowlist_search-terms.csv')
        ls_FilterLog_Entries = []

        for sample_name in ls_STSamples:
            sSampleOutputDir = sOutputDir + sample_name + '/'
            if not os.path.exists(sOutputDir + sample_name):
                os.makedirs(sOutputDir + sample_name)
            df_ReadFeatures = pooledpe_sta.process_fastq_files(sample_name, sSTPreProcessedFastqDir, bIncludeScaffoldAsFeature)
            df_ReadFeatures.to_csv(sSampleOutputDir + sample_name + '_read-features.csv')

            t_start = time.process_time()
            df_FeatureMatches, df_ReadFeatures_Filtered, df_pegRNA_Counts, dict_FilterLog, df_LibBC_Counts, df_pegRNABC_Counts, df_RTTPBS_Counts = pooledpe_sta.filter_pegrnas(
                sample_name, df_ReadFeatures, df_pegRNAs_AllowList, bIncludeScaffoldAsFeature)
            t_end = time.process_time()
            print('[pooledpe.py] pegRNA filtering for ' + sample_name + ' completed in ' + '{:.2f}'.format(
                (t_end - t_start) / 60) + ' min')

            df_ReadFeatures_Filtered.to_csv(sSampleOutputDir + sample_name + '_read-features_filtered.csv')
            df_pegRNA_Counts.to_csv(sSampleOutputDir + sample_name + '_pegRNA_counts.csv', index=False)
            ls_FilterLog_Entries.append(dict_FilterLog)

        df_FilterLog_Merged = pd.DataFrame(ls_FilterLog_Entries)
        df_FilterLog_Merged.to_csv(sOutputDir + 'pegRNA-filter-log.csv', index=False)

   if (bAnalyseSTs):  # analyse filtered pegRNA-STs
        print('[pooledpe.py] Calculating surrogate target editing rates')
        df_Variant_SearchTerms = pd.read_csv(sOutputDir + sExperiment + '_variant_search-terms_updated.csv', index_col=0)
        df_Variant_SearchTerms = df_Variant_SearchTerms.fillna('X')

        for sample_name in ls_STSamples:
            sSampleOutputDir = sOutputDir + sample_name + '/'
            sSTFilesDir = sSampleOutputDir + sample_name + '_sensors/'
            df_ReadFeatures_Filtered = pd.read_csv(sSampleOutputDir + sample_name + '_read-features_filtered.csv')
            if (bWriteSTFiles):
                if not os.path.exists(sSTFilesDir):
                    os.makedirs(sSTFilesDir)
            t_start = time.process_time()
            df_ST_Editing = pooledpe_sta.calculate_st_editing(sample_name, df_ReadFeatures_Filtered, df_Variant_SearchTerms, bWriteSTFiles, sSTFilesDir)
            t_end = time.process_time()
            print('[pooledpe.py] Calculation of ST editing efficiencies for ' + sample_name + ' completed in ' + '{:.2f}'.format((t_end - t_start) / 60) + ' min')

            # calculate pegRNA (pseudo-)frequencies
            df_pegRNA_Counts = pd.read_csv(sSampleOutputDir + sample_name + '_pegRNA_counts.csv')
            total_counts = df_pegRNA_Counts[sample_name + '_pegRNA_count'].sum()
            df_pegRNA_Counts[sample_name + '_pegRNA_freq'] = df_pegRNA_Counts[sample_name + '_pegRNA_count'].apply(lambda x: x / total_counts)
            df_pegRNA_Counts[sample_name + '_pegRNA_pseudocount'] = df_pegRNA_Counts[sample_name + '_pegRNA_count'] + 1  # add pseudocount of 1
            total_pseudocounts = df_pegRNA_Counts[sample_name + '_pegRNA_pseudocount'].sum()
            df_pegRNA_Counts[sample_name + '_pegRNA_pseudofreq'] = df_pegRNA_Counts[sample_name + '_pegRNA_pseudocount'].apply(lambda x: x / total_pseudocounts)
            
            # merge ST editing efficiencies with pegRNA counts dataframe for final output
            df_ST_Editing = df_ST_Editing.rename(columns={'mut_percentage': sample_name + '_percentage_editing',
                                                          'n_mut': sample_name + '_edited_count'})[[sample_name + '_percentage_editing', sample_name + '_edited_count']]
            total_edited_counts = df_ST_Editing[sample_name + '_edited_count'].sum()
            df_ST_Editing[sample_name + '_edited_freq'] = df_ST_Editing[sample_name + '_edited_count'].apply(lambda x: x / total_edited_counts)
            df_pegRNA_Counts_Efficiencies = df_pegRNA_Counts.merge(df_ST_Editing, left_on='pegRNA_BC', right_index=True, how='left')
            df_pegRNA_Counts_Efficiencies.to_csv(sSampleOutputDir + sample_name + '_pegRNA_counts_efficiencies.csv', index=False)

    if (bAnalyseReadFiltering):
        print('[pooledpe.py] Analysing and visualising fastq to pegRNA processing')
        ls_STSamples = ['CK001', 'CK002', 'CK003','CK004', 'CK005', 'CK006', 'CK007', 'CK008', 'CK009', 'CK010',
                        'CK011', 'CK012', 'CK013', 'CK014', 'CK015', 'CK016', 'CK017', 'CK018', 'CK019', 'CK020',
                        'CK021', 'CK022', 'CK023', 'CK024', 'MH001', 'MH002', 'MH003', 'MH004', 'MH005', 'MH006',
                        'MH007', 'MH008', 'MH009', 'MH010', 'MH011', 'MH012', 'MH013', 'MH014', 'MH015', 'MH016',
                        'MH017', 'MH018', 'MH019', 'MH020', 'MH021', 'MH022', 'MH023', 'MH024', 'MH025', 'MH026',
                        'MH027', 'MH028']
        df_NGSReadProcessing_data = pd.read_csv(sBaseDir + 'output/NGSReadProcessing_mod.csv')
        dict_SampleStatistics_All = {}
        for sample_name in ls_STSamples:
            sSampleOutputDir = sOutputDir + sample_name + '/'
            df_ReadFeatureMatches = pd.read_csv(sSampleOutputDir + sample_name + '_read-feature-matches.csv', index_col=0)
            ls_SumMatches = df_ReadFeatureMatches.sum(axis=0).to_list()

            dict_SampleStatistics = {sample_name: ls_SumMatches[:-1]}
            dict_SampleStatistics_All.update(dict_SampleStatistics)
        df_SampleStatistics = pd.DataFrame.from_dict(dict_SampleStatistics_All, orient='index')
        df_SampleStatistics.columns = [x + '_count' for x in df_ReadFeatureMatches.columns[:-1]]
        df_NGSReadProcessing_data_merged = pd.merge(df_NGSReadProcessing_data, df_SampleStatistics, left_on='ID', right_index=True, how='left')
        df_NGSReadProcessing_data_merged.to_csv(sBaseDir + 'output/NGSReadProcessing_mod.csv', index=False)

        # calculate frequencies
        df_NGSReadProcessing_data = pd.read_csv(sBaseDir + 'output/NGSReadProcessing_merged.csv')
        df_NGSReadProcessing_data['all_features_freq'] = df_NGSReadProcessing_data['all_features_count'].div(df_NGSReadProcessing_data['total_read_count'], axis=0)
        ls_ColNames = ['allowed_pegrna', 'spacer_match', 'sensor_match', 'nonlookup_match']
        ls_ColNames_Count = [x + '_count' for x in ls_ColNames]
        ls_ColNames_Freq = [x + '_freq' for x in ls_ColNames]
        df_NGSReadProcessing_data[ls_ColNames_Freq] = df_NGSReadProcessing_data[ls_ColNames_Count].div(df_NGSReadProcessing_data['all_features_count'], axis=0)
        df_NGSReadProcessing_data.to_csv(sBaseDir + 'output/NGSReadProcessing_merged_freq.csv', index=False)

    if (bMergeSampleDFs):
        print('[pooledpe.py] Merging all sample pegRNA counts and efficiencies into one dataframe')
        df_pegRNA_Data_merged = pd.read_csv(sOutputDir + sExperiment + '_pegRNA-allowlist_search-terms.csv')
        for sample_name in ls_STSamples:
            sSampleOutputDir = sOutputDir + sample_name + '/'
            df_pegRNA_Counts_Efficiencies = pd.read_csv(
                sSampleOutputDir + sample_name + '_pegRNA_counts_efficiencies.csv',
                usecols=['pegRNA_BC', sample_name + '_pegRNA_count',
                         sample_name + '_pegRNA_freq',
                         sample_name + '_pegRNA_pseudocount',
                         sample_name + '_pegRNA_pseudofreq',
                         sample_name + '_edited_count',
                         sample_name + '_edited_freq',
                         sample_name + '_percentage_editing'])
            df_pegRNA_Data_merged = df_pegRNA_Data_merged.merge(df_pegRNA_Counts_Efficiencies, on='pegRNA_BC', how='left')
        df_pegRNA_Data_merged.to_csv(sOutputDir + sExperiment + '_pegRNA_counts_efficiencies_merged.csv', index=False)

    if (bIntegrateExternalData):
        print('[pooledpe.py] Integrating pegRNA data with external databases')
        df_pegRNA_Data_merged = pd.read_csv(sOutputDir + sExperiment + '_pegRNA_counts_efficiencies_merged.csv')
        if os.path.isfile(sPathToCADDFile):  # adding CADD data to pegRNA dataframe
            df_CADD = pd.read_csv(sPathToCADDFile)
            if (sExperiment == 'MLH1x10'):
                df_CADD = df_CADD.loc[df_CADD['GeneName'] == 'MLH1']
                df_pegRNA_Data_merged = pd.merge(df_pegRNA_Data_merged, df_CADD, how='left',
                                                 left_on=['start', 'mut_allele'],
                                                 right_on=['Pos', 'Alt'])
                df_pegRNA_Data_merged.loc[df_pegRNA_Data_merged['variant_type'] == 'ONP', 'Consequence'] = 'STOP_GAINED'
            elif (sExperiment == 'MLH1intronic'):
                df_CADD = df_CADD.loc[df_CADD['GeneName'] == 'MLH1']
                df_pegRNA_Data_merged = pd.merge(df_pegRNA_Data_merged, df_CADD, how='left',
                                                 left_on=['start', 'mut_allele'],
                                                 right_on=['Pos', 'Alt'])
            elif (sExperiment == 'SMARCB1'):
                df_CADD = df_CADD.loc[df_CADD['GeneName'] == 'SMARCB1']
                df_pegRNA_Data_merged = pd.merge(df_pegRNA_Data_merged, df_CADD, how='left',
                                                 left_on=['start', 'mut_allele'],
                                                 right_on=['Pos', 'Alt'])
                df_pegRNA_Data_merged.loc[df_pegRNA_Data_merged['variant_type'] == 'ONP', 'Consequence'] = 'TBD (SNV+)'
                df_pegRNA_Data_merged.loc[df_pegRNA_Data_merged['mut_allele'] == 'TAA', 'Consequence'] = 'STOP_GAINED'
                df_pegRNA_Data_merged.loc[df_pegRNA_Data_merged['mut_allele'] == '-', 'Consequence'] = 'CODON_DELETION'
            df_pegRNA_Data_merged.to_csv(sOutputDir + sExperiment + '_pegRNA_data_merged.csv', index=False)

        if os.path.isfile(sPathToClinVarFile):  # adding ClinVar data to pegRNA dataframe
            df_pegRNA_Data_merged = pd.read_csv(sOutputDir + sExperiment + '_pegRNA_data_merged.csv')
            df_ClinVar = pd.read_csv(sPathToClinVarFile)
            if (sExperiment != 'MLH1intronic'):
                # convert ClinVar variant file to pegRNA dataframe format
                df_ClinVar = df_ClinVar.join(df_ClinVar['Canonical SPDI'].str.split(':', expand=True).iloc[:, [2, 3]])
                df_ClinVar['start'] = df_ClinVar['GRCh37Location'].str.split(' - ', expand=True).iloc[:, 0].astype(int)
                df_ClinVar = df_ClinVar.rename(columns={2: 'ref_allele', 3: 'mut_allele'})
            # extract ClinVar significance and group similar terms
            dict_SimpleClinVarSignificance = {'Benign/Likely benign': 'Likely benign',
                                              'Pathogenic/Likely pathogenic': 'Likely pathogenic',
                                              'Uncertain significance': 'Uncertain',
                                              'Conflicting interpretations of pathogenicity': 'Conflicting'}
            df_ClinVar['ClinVar_Significance'] = df_ClinVar['Clinical significance (Last reviewed)'].str.split('(', expand=True).iloc[:, 0]
            df_ClinVar['ClinVar_Significance_Simple'] = df_ClinVar['ClinVar_Significance'].replace(
                dict_SimpleClinVarSignificance)
            # merge ClinVar data with pegRNA dataframe
            df_pegRNA_Data_merged = pd.merge(df_pegRNA_Data_merged, df_ClinVar, how='left',
                                             on=['start', 'ref_allele', 'mut_allele'])
            df_pegRNA_Data_merged.to_csv(sOutputDir + sExperiment + '_pegRNA_data_merged.csv', index=False)

        if os.path.isfile(sInputDir + 'PRIDICT/pegRNA_oligos_df_asordered_PRIDICTScores.csv'):  # adding PRIDICT scores to pegRNA dataframe
            df_pegRNA_Data_merged = pd.read_csv(sOutputDir + sExperiment + '_pegRNA_data_merged.csv')
            df_PRIDICTScores = pd.read_csv(sInputDir + 'PRIDICT/pegRNA_oligos_df_asordered_PRIDICTScores.csv',
                                           usecols=['pegRNA_BC', 'PRIDICT_editing_Score_deep'])

            df_pegRNA_Data_merged = pd.merge(df_pegRNA_Data_merged, df_PRIDICTScores, on='pegRNA_BC', how='left')
            df_pegRNA_Data_merged.to_csv(sOutputDir + sExperiment + '_pegRNA_data_merged.csv', index=False)

    if (sExperiment == 'LV465'):  # separate processing of scaffold data set
        df_pegRNAs_AllowList = pd.read_csv(sOutputDir + sExperiment + '_pegRNA-allowlist_search-terms.csv')

        # assign scaffold design
        for sample_name in ls_STSamples:
            sSampleOutputDir = sOutputDir + sample_name + '/'
            df_ReadFeatures_Filtered = pd.read_csv(sSampleOutputDir + sample_name + '_read-features_filtered.csv')
            df_ReadFeatures_Filtered_processed, df_Scaffold_Counts = pooledpe_sta.assign_scaffold_design(df_ReadFeatures_Filtered)
            df_ReadFeatures_Filtered_processed.to_csv(sSampleOutputDir + sample_name + '_read-features_filtered_processed.csv', index=False)
            df_pegRNA_Counts = df_ReadFeatures_Filtered_processed.groupby(df_ReadFeatures_Filtered_processed[['pegRNA-BC', 'scaffold_design']].columns.tolist(), as_index=False).size().rename(columns={'pegRNA-BC': 'pegRNA_BC', 'size': 'pegRNA_count'})
            df_pegRNA_Counts = pd.merge(df_pegRNAs_AllowList, df_pegRNA_Counts, on='pegRNA_BC', how='left')
            df_pegRNA_Counts.to_csv(sSampleOutputDir + sample_name + '_pegRNA_counts.csv', index=False)

        # calculate percentage ST editing
        for sample_name in ls_STSamples:
            sSampleOutputDir = sOutputDir + sample_name + '/'
            df_ReadFeatures_Filtered = pd.read_csv(sSampleOutputDir + sample_name + '_read-features_filtered_processed.csv')
            df_Variant_SearchTerms = pd.read_csv(sOutputDir + sExperiment + '_variant_search-terms_updated.csv', index_col=0)
            df_Variant_SearchTerms = df_Variant_SearchTerms.fillna('X')

            df_pegRNA_Data = pd.read_csv(sSampleOutputDir + sample_name + '_pegRNA_counts.csv')
            df_pegRNA_Data['percentage_editing'] = np.nan
            for scaffold in ['original', 'FE', 'a', 'PRIDICT', 't-lock']:
                df_ReadFeatures_Filtered_group = df_ReadFeatures_Filtered.loc[df_ReadFeatures_Filtered['scaffold_design'] == scaffold, :]
                df_ST_Editing = pooledpe_sta.calculate_st_editing(sample_name, df_ReadFeatures_Filtered_group, df_Variant_SearchTerms)
                df_ST_Editing = df_ST_Editing.rename_axis('pegRNA_BC').reset_index()
                df_ST_Editing = df_ST_Editing.rename(columns={'mut_percentage': 'percentage_editing'})
                df_ST_Editing['scaffold_design'] = scaffold
                df_pegRNA_Data = pd.merge(df_pegRNA_Data, df_ST_Editing.loc[:, ['pegRNA_BC', 'scaffold_design', 'percentage_editing']],
                                          on=['pegRNA_BC', 'scaffold_design'],
                                          how='left',
                                          suffixes=('', '_y'))
                df_pegRNA_Data['percentage_editing'] = df_pegRNA_Data['percentage_editing'].fillna(df_pegRNA_Data['percentage_editing_y'])
                df_pegRNA_Data = df_pegRNA_Data.drop(['percentage_editing_y'], axis=1)
            df_pegRNA_Data.to_csv(sSampleOutputDir + sample_name + '_pegRNA_data_mod.csv', index=False)

      
