"""
[]
"""

###----- IMPORTS ----------------------------------------------------------------------------------------------------###
import math
import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import f_oneway
import statsmodels.stats.multicomp as mc

###----- FUNCTIONS --------------------------------------------------------------------------------------------------###
def hamming_distance(s1, s2):
    """ calculates and returns Hamming/Levenshtein/edit distance between two strings provided as input """
    hamming_dist = np.inf
    if len(s1) == len(s2):
        hamming_dist = sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))
    return hamming_dist

def mismatched_positions(mutant_index, start, s1, s2):
    """ determines and returns positions of mismatches between to strings of equal size """
    ls_mismatch_pos = []
    if len(s1) == len(s2):
        for i, nt in enumerate(s1):
            if (nt != s2[i]):
                ls_mismatch_pos.append({'mutant index': mutant_index, 'start': start + i, 'ref_allele': nt, 'mut_allele': s2[i]})
    return ls_mismatch_pos

def annotate_pam_edits_for_snvs(pegrna_row):
    """ xxx """
    variant_type = pegrna_row['variant_type']
    pam_seq = pegrna_row['PAM']
    pam_strand = pegrna_row['PAM strand']
    dist_to_nick = pegrna_row['distance to nick']
    pam_edit = np.nan
    if ((pam_seq in ['AGG', 'CGG', 'GGG', 'TGG']) and (variant_type == 'SNP')):
        if (pam_strand == '-'):  # correct distance to nick if PAM on negative strand (+1)
            dist_to_nick += 1
        pam_edit = True if dist_to_nick in [5, 6] else False
    return pam_edit

def annotate_pam_edits_for_onvs(df_top_pegrna, df_correspondence):
    """ xxx """
    df_pam_edit_snv = df_top_pegrna.loc[df_top_pegrna['variant_type'] == 'SNP', ['mutant index', 'PAM_edit']]
    df_grouped = df_correspondence.groupby('mutant_index_ONV')
    dict_onv_pam_edit_merged = {}
    for name, df_group in df_grouped:
        ls_snv_indexes = list(df_group['mutant_index_SNV'])
        ls_pam_edits_SNVs = list(df_pam_edit_snv.loc[df_pam_edit_snv['mutant index'].isin(ls_snv_indexes), 'PAM_edit'])
        pam_edit_onv = True if ls_pam_edits_SNVs.count(True) > 0 else False
        dict_onv_pam_edit_merged.update({name: pam_edit_onv})
    df_top_pegrna['PAM_edit'] = df_top_pegrna['mutant index'].map(dict_onv_pam_edit_merged).fillna(df_top_pegrna['PAM_edit'])
    return df_top_pegrna

def normalise_scores(df_data, ls_samples, renormalise):
    """ xxx """
    for sample_name in ls_samples:
        # normalise log2FCs to median synonymous (or intronic) scores
        if (renormalise == False):
            col_name = 'log2FC_' + sample_name
        else:
            col_name = 'log2FC_' + sample_name + '_norm'
        median_log2fc_neutral = df_data.loc[df_data['Consequence'] == 'SYNONYMOUS', col_name].median()
        if (np.isnan(median_log2fc_neutral)):
            median_log2fc_neutral = df_data.loc[df_data['Consequence'] == 'INTRONIC', col_name].median()
        df_data['log2FC_' + sample_name + '_norm'] = df_data[col_name] - median_log2fc_neutral
    return df_data

def determine_fdr(df_data, sample_name):
    """ xxx """
    df_data_mod = df_data.dropna(subset=['log2FC_' + sample_name + '_norm'])
    loc, scale = stats.norm.fit(df_data_mod.loc[df_data_mod['Consequence'] == 'SYNONYMOUS', 'log2FC_' + sample_name + '_norm'])
    df_data_mod['p_values'] = stats.norm.cdf(df_data_mod['log2FC_' + sample_name + '_norm'], loc, scale)
    df_data_mod['p_adjust'] = stats.false_discovery_control(df_data_mod['p_values'], method='bh')
    #df_data_mod['fdr_0.01'] = df_data_mod['p_adjust'].map(lambda x: True if x < 0.01 else False)
    #df_data_mod['fdr_0.02'] = df_data_mod['p_adjust'].map(lambda x: True if x < 0.02 else False)
    df_data_mod['fdr_0.05'] = df_data_mod['p_adjust'].map(lambda x: True if x < 0.05 else False)
    return df_data_mod

def mask_st_editing_with_low_pegrna_count(df_pegrna_data, ls_sample_ids, count_threshold=10):
    """ xxx """
    df_pegrna_data_masked = df_pegrna_data
    for sample_id in ls_sample_ids:
        count_colname = sample_id + '_pegRNA_count'
        st_colname = sample_id + '_percentage_editing'
        df_pegrna_data_masked.loc[df_pegrna_data_masked[count_colname] < count_threshold, st_colname] = np.nan
    return df_pegrna_data_masked

def calculate_confidence_interval(row, sample_name, conf_level, lower):
    """ xxx """
    mean = row['log2FC_' + sample_name + '_mean']
    sd = row['log2FC_' + sample_name + '_sd']
    count = row['count']
    ci_adj = np.nan
    if (count >= 1):
        ci_adj = conf_level * sd / math.sqrt(count)
    if lower == True:
        ci = mean - ci_adj
    else:
        ci = mean + ci_adj
    return ci
