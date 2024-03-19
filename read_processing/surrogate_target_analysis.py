"""
[19/03/2024, Michael Herger] surrogate_target_analysis.py 
Description: Functions to facilitate fastq-to-pegRNA conversion and analysis of pegRNA frequencies and surrogate target 
editing efficiencies
"""

###----- IMPORTS ----------------------------------------------------------------------------------------------------###
import os
import subprocess
import gzip
import glob
import regex as re
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq


###----- FUNCTIONS --------------------------------------------------------------------------------------------------###
def generate_pegrna_allowlist(path_to_pegrna_file, list_allowed_lib_bcs, dict_lib_bcs):
    """ function to load reference dataframe of all allowed pegRNAs and transform to dictionaries """
    df_pegrnas_allowlist = pd.read_csv(path_to_pegrna_file)
    df_pegrnas_allowlist = df_pegrnas_allowlist.loc[df_pegrnas_allowlist['Library_ID'].isin(list_allowed_lib_bcs)]
    df_pegrnas_allowlist['Library_BC'] = df_pegrnas_allowlist['Library_ID'].map(dict_lib_bcs)
    return df_pegrnas_allowlist

def add_edit_seq_to_allowlist(path_to_pegrna_file, path_to_search_term_file):
    """ generate and QC search terms for each pegRNA """
    df_pegrnas_allowlist = pd.read_csv(path_to_pegrna_file)
    df_search_terms = pd.read_csv(path_to_search_term_file)
    df_search_terms = df_search_terms.rename(columns={'Unnamed: 0': 'pegRNA_BC'})
    df_pegrnas_allowlist_mod = pd.merge(df_pegrnas_allowlist, df_search_terms, on='pegRNA_BC', how='left')

    ls_edited_target_seqs = []
    ls_qc_pass_search_terms = []
    for indx, row in df_pegrnas_allowlist_mod.iterrows():
        pegrna_bc = row['pegRNA_BC']
        target_seq = row['target sequence']
        ref_search_term = row['ref_search_term']
        edit_search_term = row['mut_search_term'].upper()
        edited_target_seq = target_seq.replace(ref_search_term, edit_search_term)
        ls_edited_target_seqs.append(edited_target_seq)

        # QC of search terms
        qc_pass_search_term = True
        if ((ref_search_term == '') or (edit_search_term == '')):
            qc_pass_search_term = False
            print('No search terms found for pegRNA with BC: ' + pegrna_bc)
        if (edit_search_term in target_seq):
            qc_pass_search_term = False
            print('Edit search term occurs in unedited sensor sequence -> redesign for ' + pegrna_bc)
        if (ref_search_term in edited_target_seq):
            qc_pass_search_term = False
            print('Ref search term occurs in edited sensor sequence -> redesign for ' + pegrna_bc)
        ls_ref_search_term_matches = [m.start() for m in re.finditer(ref_search_term, target_seq)]
        if (len(ls_ref_search_term_matches) != 1):
            qc_pass_search_term = False
            print('Ref search term occurs multiple times in unedited sensor sequence -> redesign for ' + pegrna_bc)
        ls_edit_search_term_matches = [m.start() for m in re.finditer(edit_search_term, edited_target_seq)]
        if (len(ls_edit_search_term_matches) != 1):
            qc_pass_search_term = False
            print('Edit search term occurs multiple times in edited sensor sequence -> redesign for ' + pegrna_bc)
        ls_qc_pass_search_terms.append(qc_pass_search_term)

    df_pegrnas_allowlist_mod['edited target sequence'] = ls_edited_target_seqs
    df_pegrnas_allowlist_mod['search_terms_pass_QC'] = ls_qc_pass_search_terms
    return df_pegrnas_allowlist_mod

def generate_edit_search_terms(df_pegrnas_allowlist, flanking_nts=5, lookup_feature='pegRNA_BC'):
    """ xxx """
    dict_variant_search_terms_merged = {}
    for index, row in df_pegrnas_allowlist.iterrows():
        pegrna_bc = row[lookup_feature]
        variant_type = row['variant_type']
        pam_strand = row['PAM strand']
        pam_seq = row['PAM']
        spacer_seq = row['protospacer']
        dist_nick_to_mut = row['distance to nick']
        ref_seq = row['ref_allele']
        mut_seq = row['mut_allele'].lower()
        target_seq = row['target sequence']
        query_seq = spacer_seq + pam_seq
        if (pam_strand == '+'):
            query_pos = target_seq.find(query_seq)
            #mut_pos = query_pos + len(query_seq) - 6 + dist_nick_to_mut - 1  # 6 nt: 3 nt PAM + 3 nt from PAM to nick
            mut_pos = query_pos + len(query_seq) - len(pam_seq) - 3 + dist_nick_to_mut - 1
        else:
            query_seq = str(Seq(query_seq).reverse_complement()).upper()
            query_pos = target_seq.find(query_seq)
            #mut_pos = query_pos + 6 - dist_nick_to_mut - len(ref_seq)
            mut_pos = query_pos + len(pam_seq) + 3 - dist_nick_to_mut - len(ref_seq)

        if (variant_type == 'DEL'):
            if (mut_seq == '-'):
                mut_seq = ''
        elif (variant_type == 'INS'):
            if (pam_strand == '+'):
                mut_pos = mut_pos - len(mut_seq) + len(ref_seq)
            else:
                mut_pos = mut_pos + len(mut_seq) - len(ref_seq)
        elif (variant_type not in ['SNP', 'ONP', 'DEL', 'INS']):
            raise Exception('Variant type ' + variant_type + ' not supported. Takes values of SNP, ONP, DEL, and INS only.')

        ref_search_term = target_seq[mut_pos - flanking_nts:mut_pos + len(ref_seq) + flanking_nts]
        edit_search_term = target_seq[mut_pos - flanking_nts:mut_pos] + mut_seq + target_seq[mut_pos + len(ref_seq):mut_pos + len(ref_seq) + flanking_nts]
        edited_target_seq = target_seq[0:mut_pos] + mut_seq + target_seq[mut_pos + len(ref_seq):len(target_seq)]

        # QC of search terms
        qc_pass_search_term = search_term_qc(pegrna_bc, target_seq, edited_target_seq, ref_search_term, edit_search_term)
        dict_variant_search_terms = {pegrna_bc: [ref_search_term, edit_search_term, qc_pass_search_term, edited_target_seq]}
        dict_variant_search_terms_merged.update(dict_variant_search_terms)

    df_variant_search_terms_merged = pd.DataFrame.from_dict(dict_variant_search_terms_merged, orient='index',
                                                            columns=['ref_search_term', 'mut_search_term',
                                                                     'search_terms_passed_qc', 'edited target sequence'])
    df_variant_search_terms_merged = df_variant_search_terms_merged.rename_axis('pegRNA_BC').reset_index()
    return df_variant_search_terms_merged

def search_term_qc(pegrna_bc, target_seq, edited_target_seq, ref_search_term, edit_search_term):
    """ performs  quality control for ref and edit search terms """
    edited_target_seq = edited_target_seq.upper()
    edit_search_term = edit_search_term.upper()
    qc_pass_search_term = True
    if ((ref_search_term == '') or (edit_search_term == '')):
        qc_pass_search_term = False
        print('No search terms found for pegRNA with BC: ' + pegrna_bc)
    if (edit_search_term in target_seq):
        qc_pass_search_term = False
        print('Edit search term occurs in unedited sensor sequence -> redesign for ' + pegrna_bc)
    if (ref_search_term in edited_target_seq):
        qc_pass_search_term = False
        print('Ref search term occurs in edited sensor sequence -> redesign for ' + pegrna_bc)
    ls_ref_search_term_matches = [m.start() for m in re.finditer(ref_search_term, target_seq)]
    if (len(ls_ref_search_term_matches) != 1):
        qc_pass_search_term = False
        print('Ref search term occurs multiple times in unedited sensor sequence -> redesign for ' + pegrna_bc)
    ls_edit_search_term_matches = [m.start() for m in re.finditer(edit_search_term, edited_target_seq)]
    if (len(ls_edit_search_term_matches) != 1):
        qc_pass_search_term = False
        print('Edit search term occurs multiple times in edited sensor sequence -> redesign for ' + pegrna_bc)
    return qc_pass_search_term

def process_fastq_files(sample_name, directory_of_preprocessed_fastqs, include_scaffold_as_feature=False):
    """ assigns all retrieved pegRNA elements to each read (assumes pre-processed fastq files to be formatted specifically) """
    dict_read_features = {}
    ls_features = ['protospacer', 'RTT-PBS', 'sensor', 'pegRNA-BC', 'library-BC', 'scaffold']
    if not include_scaffold_as_feature:
        ls_features = ls_features[:-1]
    for feature_name in ls_features:
        print('\tAssigning ' + feature_name + '...')
        path_to_fastq_file = glob.glob(directory_of_preprocessed_fastqs + sample_name + '*' + feature_name + '.fastq.gz')[0]
        with gzip.open(path_to_fastq_file, 'rt') as fastqFile:
            for seq_record in SeqIO.parse(fastqFile, 'fastq'):
                feature_seq = seq_record.seq
                read_identifier = seq_record.id
                if read_identifier in dict_read_features:
                    dict_read_features[read_identifier][feature_name] = feature_seq
                else:
                    dict_read_features[read_identifier] = {feature_name: feature_seq}
    df_read_features = pd.DataFrame.from_dict(dict_read_features, orient='index')
    return df_read_features

def filter_pegrnas(sample_name, df_read_features, df_pegrnas_allowlist, lookup_feature='pegRNA-BC',
                      include_scaffold_as_feature=False, error_tolerant_matching=True):
    """ checks read elements against allowed pegRNAs and discards reads with forbidden element combinations """
    # check if pegRNA with sequenced RTT-PBS/pegRNA-BC matches expected features
    ls_allowed_library_bcs = df_pegrnas_allowlist['Library_BC'].drop_duplicates()
    df_read_features_prefiltered = df_read_features[df_read_features['library-BC'].isin(ls_allowed_library_bcs)]
    ls_features = ['protospacer', 'RTT-PBS', 'sensor', 'pegRNA-BC', 'scaffold']
    if not include_scaffold_as_feature:
        ls_features = ls_features[:-1]
    df_read_features_filtered = df_read_features_prefiltered[ls_features].dropna()
    dict_feature_matches_merged = {}
    ls_indexes_pegrna_match = []
    ls_lookup_sequence_updated = []
    ls_pegrna_counts = [0] * len(df_pegrnas_allowlist)

    for index, row in df_read_features_filtered.iterrows():
        protospacer_read = str(row['protospacer'])
        rttpbs_read = str(row['RTT-PBS'])
        sensor_read = str(row['sensor'])
        pegrnabc_read = str(row['pegRNA-BC'])
        mismatch_threshold = 0.95  # define fraction of nts required to match each allowed feature
        if (lookup_feature == 'RTT-PBS'):
            lookup_sequence = rttpbs_read
            nonlookup_sequence = pegrnabc_read
            lookup_column = df_pegrnas_allowlist['PBS_RTT_5to3']
            nonlookup_column = df_pegrnas_allowlist['pegRNA_BC']
        elif (lookup_feature == 'pegRNA-BC'):
            lookup_sequence = pegrnabc_read
            nonlookup_sequence = rttpbs_read
            lookup_column = df_pegrnas_allowlist['pegRNA_BC']
            nonlookup_column = df_pegrnas_allowlist['PBS_RTT_5to3']
        else:
            raise Exception('Lookup feature "' + lookup_feature + '" not allowed (either "RTT-PBS" or "pegRNA-BC"')
        dict_feature_matches = {'spacer_match': False,
                                'sensor_match': False,
                                'nonlookup_match': False,
                                'pegrna_match': False}
        ls_index = df_pegrnas_allowlist.loc[df_pegrnas_allowlist['pegRNA_BC'] == lookup_sequence].index

        # only proceed if lookup feature is perfectly and uniquely matched to element of allow list
        if (len(ls_index) == 1):
            matching_index = ls_index[0]
            matching_index_loc = df_pegrnas_allowlist.index.get_loc(matching_index)
            lookup_sequence = str(lookup_column.at[matching_index])
            ls_lookup_sequence_updated.append(lookup_sequence)  # update lookup sequence to matching one
            protospacer_wl = str(df_pegrnas_allowlist.loc[matching_index]['protospacer'])
            sensor_wl = str(df_pegrnas_allowlist.loc[matching_index]['target sequence'])
            nonlookup_feature_wl = str(nonlookup_column.at[matching_index])
            dict_feature_matches['spacer_match'] = bool(re.search(r'(%s){s<=%d}'%(protospacer_read, 2), protospacer_wl))
          
            # check 3' of surrogate target for match against unedited and unedited surrogate target
            sensor_match = False
            edited_sensor_wl = str(df_pegrnas_allowlist.loc[matching_index]['edited target sequence'])
            if (bool(re.search(r'(%s){s<=%d}' % (sensor_read[-20:], 2), sensor_wl[-20:]))):
                sensor_match = True
            elif (bool(re.search(r'(%s){s<=%d}' % (sensor_read[-20:], 2), edited_sensor_wl[-20:]))):
                sensor_match = True
            dict_feature_matches['sensor_match'] = sensor_match
            allowed_mismatches = int(np.ceil(len(nonlookup_sequence) * (1 - mismatch_threshold)))
            dict_feature_matches['nonlookup_match'] = bool(re.search(r'(%s){s<=%d}' % (nonlookup_sequence, allowed_mismatches), nonlookup_feature_wl))
            dict_feature_matches['pegrna_match'] = dict_feature_matches['spacer_match'] & \
                                                   dict_feature_matches['sensor_match'] & \
                                                   dict_feature_matches['nonlookup_match']
        else:  # none or multiple pegRNAs with lookup sequence identified
            ls_lookup_sequence_updated.append(np.nan)
        dict_feature_matches_merged.update({index: dict_feature_matches})
        if dict_feature_matches['pegrna_match'] is True:  # read encodes allowed pegRNA
            ls_indexes_pegrna_match.append(index)
            ls_pegrna_counts[matching_index_loc] += 1
    df_read_features_filtered[lookup_feature] = ls_lookup_sequence_updated  # re-assign accepted lookup features
    df_read_features_filtered_final = df_read_features_filtered.loc[ls_indexes_pegrna_match]

    # check distribution of library and pegRNA BCs counts
    library_bc_counts = df_read_features['library-BC'].value_counts(ascending=False)
    pegrna_bc_counts = df_read_features_filtered_final['pegRNA-BC'].value_counts(ascending=False)
    rttpbs_counts = df_read_features_filtered_final['RTT-PBS'].value_counts(ascending=False)

    # determine pegRNA counts
    df_pegrnas_allowlist[sample_name + '_pegRNA_count'] = ls_pegrna_counts

    # return elements passing
    df_feature_matches_merged = pd.DataFrame.from_dict(dict_feature_matches_merged, orient='index')

    return df_feature_matches_merged, df_read_features_filtered_final, df_pegrnas_allowlist, dict_filter_log, library_bc_counts, pegrna_bc_counts, rttpbs_counts

def assign_scaffold_design(df_read_features_filtered):
    """ identify and assign scaffold name to scaffold sequence """
    dict_allowed_scaffolds = {'TAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTG': 'original',
                              'AAGAGCTATGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTG': 'FE',
                              'CAGAGCTATGCTGGAAACAGCATAGCAAGTTGAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTG': 'FEv1',
                              'CAGAGCTATGCTGGAAACAGCATAGCAAGTTGAAATAAGGCTAGTCCGTTATCAACTTGGACTTCGGTCCAAGTG': 'FEv2'}
    df_read_features_filtered['scaffold_design'] = df_read_features_filtered['scaffold'].map(dict_allowed_scaffolds)
    df_read_features_filtered = df_read_features_filtered.dropna()
    scaffold_counts = df_read_features_filtered['scaffold_design'].value_counts(ascending=False)
    return df_read_features_filtered, scaffold_counts

def calculate_st_editing(sample_name, df_read_features_filtered, df_variant_search_terms, write_sensor_files=False, sensor_files_output_dir='', lookup_feature='pegRNA-BC'):
    """ calculate percentage correct ST editing values for each pegRNA """
    df_reads_filtered_grouped = df_read_features_filtered.groupby(lookup_feature)
    dict_st_editing_merged = {}
    for group_name, df_group in df_reads_filtered_grouped:
        file_name = group_name + '_sensor.csv'
        ref_search_term = df_variant_search_terms.loc[group_name]['ref_search_term']
        mut_search_term = df_variant_search_terms.loc[group_name]['mut_search_term'].upper()
        df_group_mod = df_group
        df_group_mod['is_ref'] = df_group['sensor'].apply(lambda x: ref_search_term in x)
        df_group_mod['is_mut'] = df_group['sensor'].apply(lambda x: mut_search_term in x)
        if (write_sensor_files):
            df_group_mod.to_csv(sensor_files_output_dir + file_name, index=False)
        # calculate percentage edited STs
        n_ref_reads = df_group_mod['is_ref'].sum()
        n_mut_reads = df_group_mod['is_mut'].sum()
        n_total_reads = len(df_group_mod)
        n_other_reads = n_total_reads - n_ref_reads - n_mut_reads
        if (n_total_reads > 0):
            mut_percentage = n_mut_reads/n_total_reads*100
        else:
            mut_percentage = np.nan
        dict_st_editing = {group_name: [n_ref_reads, n_mut_reads, n_other_reads, n_total_reads, mut_percentage]}
        dict_st_editing_merged.update(dict_st_editing)
    df_st_editing_merged = pd.DataFrame.from_dict(dict_st_editing_merged, orient='index',
                                                      columns=['n_ref', 'n_mut', 'n_other', 'n_total', 'mut_percentage'])
    return df_st_editing_merged

def write_st_fastq_files(df_pegrna_allowlist, df_read_features, input_fastq_file_path, ls_pegrna_bcs, sensor_files_output_dir=''):
    """ write surrogate target fastq files """
    df_read_features_mod = df_read_features.dropna()
    lines = []
    for pegrna_bc in ls_pegrna_bcs:
        # create file with read identifiers for pegRNA
        sensor_seq = df_pegrna_allowlist.loc[df_pegrna_allowlist['pegRNA_BC'] == pegrna_bc, 'target sequence'].values[0]
        spacer_seq = df_pegrna_allowlist.loc[df_pegrna_allowlist['pegRNA_BC'] == pegrna_bc, 'protospacer'].values[0]
        extension = df_pegrna_allowlist.loc[df_pegrna_allowlist['pegRNA_BC'] == pegrna_bc, 'PBS_RTT_5to3'].values[0]
        ls_read_ids = list(df_read_features_mod.loc[df_read_features_mod['pegRNA-BC'] == pegrna_bc, 'read'])
        filename = sensor_files_output_dir + pegrna_bc + '-' + sensor_seq + '-' + spacer_seq + '.lst'
        with open(filename, 'w') as fp:
            for read_id in ls_read_ids:
                fp.write('%s\n' % read_id)
        # extract STs for pegRNA to fastq file
        output_fastq_path = sensor_files_output_dir + pegrna_bc + '-' + sensor_seq + '-' + spacer_seq + '-' + extension + '.fastq'
        with open(output_fastq_path, 'w') as f:
            subprocess.call(['seqtk', 'subseq', input_fastq_file_path, filename], stdout=f)
        os.remove(filename)
        # create CRISPResso command
        lines.append(generate_crispresso_line(pegrna_bc, sensor_seq, spacer_seq, extension))
    with open(sensor_files_output_dir + 'CRISPResso_STAnalysis.sh', 'w') as f:
        for l in lines:
            f.write(f'{l}\n')
    return 0
