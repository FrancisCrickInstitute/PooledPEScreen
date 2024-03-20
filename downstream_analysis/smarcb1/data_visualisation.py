





def vis_t804n_enrichment_ouabain_selection():
    """ code to reproduce Figure 1D """
    # import data set
    df_edit_data = pd.read_csv('/Users/hergerm/Documents/Data/NGS/230428/ATP1A1/fastqs/grep_results.csv')
    df_edit_data['day'] = df_edit_data['day'].astype(str)
    df_edit_data_filtered = df_edit_data.loc[df_edit_data['is_control'] == False, :]
    df_edit_data_filtered = df_edit_data_filtered.loc[df_edit_data['day'].isin(['4', '11']), :]
    #df_edit_data_filtered = df_edit_data_filtered.loc[df_edit_data['day'].isin(['4', '11', '18', '24']), :]
    df_edit_data_filtered['group'] = df_edit_data_filtered.apply(lambda row: row['day'] + row['ouabain'], axis=1)
    df_edit_data_perfect = df_edit_data_filtered.loc[df_edit_data_filtered['mut_type'] == 'PE', :]
    df_edit_data_partial = df_edit_data_filtered.loc[df_edit_data_filtered['mut_type'] == 'C5T', :]

    fig, ax = plt.subplots(figsize=(1.4, 1.4))
    ax.bar(df_edit_data_perfect['group'], df_edit_data_perfect['percentage'], label='PE (correct)',
           linewidth=0.2, edgecolor='black')
    ax.bar(df_edit_data_partial['group'], df_edit_data_partial['percentage'], bottom=df_edit_data_perfect['percentage'],
           label='PE (silent)',
           linewidth=0.2, edgecolor='black')
    ax.set(xlabel='Day + obn', ylabel='% of Reads')
    ax.spines[['right', 'top']].set_visible(False)
    ax.legend()
    plt.ylim(0, 100)
    plt.tight_layout()
    plt.savefig(sVisOutputDir + 'Figure1X.svg', format='svg')
    plt.close()
    return 0

def vis_temporal_and_coselection_st_editing(df_pegrna_data):
    """ code to reproduce Figure 2F """
    df_pegrna_data_mod = df_pegrna_data
    # average over duplicates
    ls_duplicates = [['MH003', 'MH004'], ['MH005', 'MH006'], ['MH007', 'MH007'], ['MH008', 'MH009'], ['MH010', 'MH011'],
                     ['MH013', 'MH014'], ['MH015', 'MH015']]
    ls_timepoints = ['D4', 'D10', 'D10obn', 'D20', 'D27', 'D34', 'D34obn']
    df_pegrna_data_mod = df_pegrna_data.loc[df_pegrna_data['percentage_editing_control_high'] == False, :]
    ls_timepoint_editing_colnames = [timepoint + '_percentage_editing' for timepoint in ls_timepoints]
    dict_timepoints = {ls_timepoint_editing_colnames[i]: ls_timepoints[i] for i in
                       range(len(ls_timepoint_editing_colnames))}

    for n in range(0, len(ls_duplicates)):
        ls_samples = ls_duplicates[n]
        ls_count_colnames = [sample + '_pegRNA_count' for sample in ls_samples]
        ls_editing_colnames = [sample + '_percentage_editing' for sample in ls_samples]
        # filter out pegRNAs with count < 10
        df_pegrna_data_mod.loc[df_pegrna_data_mod[ls_count_colnames[0]] < 10, ls_editing_colnames[0]] = np.nan
        df_pegrna_data_mod.loc[df_pegrna_data_mod[ls_count_colnames[1]] < 10, ls_editing_colnames[1]] = np.nan
        # average ST editing across duplicates
        df_pegrna_data[ls_timepoint_editing_colnames[n]] = df_pegrna_data_mod[ls_editing_colnames].mean(axis=1, skipna=True)

    # reformat dataframe
    df_pegrna_data = df_pegrna_data.rename(columns=dict_timepoints)
    df_pegrna_data_melt = df_pegrna_data[ls_timepoints].melt(var_name='day', value_name='vals')
    df_pegrna_data_melt['obn'] = df_pegrna_data_melt['day'].apply(lambda x: 'obn' in x)
    df_pegrna_data_melt['day'] = df_pegrna_data_melt['day'].replace({'D10obn': 'D10', 'D34obn': 'D34'})

    # visualise as violin plot with categorical x-axis
    #plt.figure(figsize=(3, 1.4))
    plt.figure(figsize=(2, 2))
    sns_plot = sns.violinplot(data=df_pegrna_data_melt, x='day', y='vals', hue='obn',
                              split=False, inner='quart', fill=True,
                              order=['D4', 'D10', 'D20', 'D27', 'D34'],
                              density_norm='count', width=1, #color=ST_Data_Color,
                              linewidth=0.2, linecolor='k')
    for l in sns_plot.lines:
        l.set_linewidth(0.4)
    for l in sns_plot.lines[1::3]:
        l.set_linestyle('-')
        l.set_linewidth(0.8)
        l.set_solid_capstyle('butt')
    sns_plot.set_ylim(0, 100)
    sns_plot.set(xlabel='Sample', ylabel='% Correct ST Edits')
    sns_plot.spines[['right', 'top']].set_visible(False)
    plt.tight_layout()
    plt.savefig(sVisOutputDir + 'Figure2F.svg', format='svg')
    plt.close()
    return 0


def vis_t804n_pegrnas():
    """ code to reproduce Figure S1 """
    # import data set
    df_edit_data = pd.read_csv('/Users/hergerm/Documents/Code/CRISPResso2/PE Pilot Experiment/CRISPRessoData_processed.csv')
    min_aligned_reads = 0
    df_edit_data_filtered = df_edit_data.loc[(df_edit_data['PE'] == 'PEmax') &
                                             (df_edit_data['Scaffold'] == 'FE') &
                                             (df_edit_data['ngRNA'] == '-') &
                                             (df_edit_data['Amplicon'] == 'x17') &
                                             (df_edit_data['Reads_aligned_all_amplicons'] > min_aligned_reads), :]
    df_edit_data_filtered['PE_percentage'] = df_edit_data_filtered['PE.frequency'] * 100

    plt.figure(figsize=(2, 2))
    sns_plot = sns.barplot(df_edit_data_filtered, x='pegRNA', y='PE_percentage', color='grey',
                           linewidth=0.2, edgecolor='black')
    sns.stripplot(df_edit_data_filtered, x='pegRNA', y='PE_percentage', color='black', #jitter=0.5,
                  size=2.2, linewidth=0.1, edgecolor='black')
    sns_plot.spines[['right', 'top']].set_visible(False)
    plt.tight_layout()
    plt.savefig(sVisOutputDir + 'FigureS1.svg', format='svg')
    plt.close()
    return 0

def vis_read_processing(df_read_processing_data):
    """ code to reproduce Figure S2 """
    ls_samples = ['MH001', 'MH002', 'MH003', 'MH004', 'MH005', 'MH006', 'MH007', 'MH008', 'MH009', 'MH010',
                  'MH011', 'MH012', 'MH013', 'MH014', 'MH015', 'MH016', 'MH017', 'MH024', 'MH025']  # SMARCB1
    ls_samples = ['CK001', 'CK002', 'CK003', 'CK004', 'CK005', 'CK006', 'CK007', 'CK008', 'CK009',
                  'CK010', 'CK011', 'CK023']  # MLH1x10
    ls_samples = ['CK012', 'CK013', 'CK014', 'CK015', 'CK016', 'CK017', 'CK018', 'CK019', 'CK020',
                  'CK021', 'CK022', 'CK024']  # MLH1intronic
    df_read_processing_data = df_read_processing_data.loc[df_read_processing_data['ID'].isin(ls_samples), :]

    # grouped bar plot of matches against allow list (relative to reads with all features)
    ls_col_names = ['spacer_match_freq', 'nonlookup_match_freq', 'sensor_match_freq', 'allowed_pegrna_freq']
    df_read_processing_data_filtered = df_read_processing_data.loc[:, ['ID'] + ls_col_names]
    df_read_processing_data_filtered = pd.melt(df_read_processing_data_filtered, id_vars='ID', value_vars=ls_col_names)
    plt.figure(figsize=(6, 2.2))
    sns_plot = sns.barplot(data=df_read_processing_data_filtered, x='ID', y='value', hue='variable')
    sns_plot.set(xlabel='samples', ylabel='reads passed / reads with all features')
    sns_plot.spines[['right', 'top']].set_visible(False)
    sns.move_legend(sns_plot, 'upper left', bbox_to_anchor=(1, 1))
    plt.ylim(0, 1)
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(sVisOutputDir + 'FigureS2.svg', format='svg')
    plt.close()
    return 0

def smarcb1_ko_indel_frequencies():
    """ code to reproduce Figure S3 """
    # import data set
    df_indel_data = pd.read_csv('/Users/hergerm/Documents/Code/CRISPResso2/221223_SMARCB1 Indel Depletion/IndelFrequencies.csv')
    df_indel_data['Indel.frequency'] = df_indel_data['Indel.frequency'] * 100
    df_indel_data['IF_percentage'] = df_indel_data['Indel.frequency'] * df_indel_data['InFrameFrequency']
    df_indel_data['OoF_percentage'] = df_indel_data['Indel.frequency'] * (1 - df_indel_data['InFrameFrequency'])

    ls_exons = [4, 5, 6]
    for exon in ls_exons:
        df_indel_data_exon = df_indel_data.loc[df_indel_data['Exon'] == exon, :]
        fig, ax = plt.subplots(figsize=(2, 1.8))
        ax.bar(df_indel_data_exon['Day'], df_indel_data_exon['OoF_percentage'], label='Out-Of-Frame Indel Frequency',
               width=2, linewidth=0.2, edgecolor='black')
        ax.bar(df_indel_data_exon['Day'], df_indel_data_exon['IF_percentage'], bottom=df_indel_data_exon['OoF_percentage'],
               label='In-Frame Indel Frequency',
               width=2, linewidth=0.2, edgecolor='black')
        ax.set(xlabel='Day', ylabel='% of sequencing reads')
        ax.spines[['right', 'top']].set_visible(False)
        ax.legend()
        plt.ylim(0, 100)
        plt.tight_layout()
        fig.savefig(sVisOutputDir + 'FigureS3.svg', format='svg')
        plt.close()
    return 0

def vis_pegrna_data_sample_correlation(df_pegrna_data, dict_sample_info, experiment):
    """ code to redproduce Figure S4A and S4B """
    ls_samples = ['MH024', 'MH025', 'MH001', 'MH002', 'MH003', 'MH004', 'MH005', 'MH006', 'MH007', 'MH008', 'MH009', 'MH010', 'MH011',
                      'MH013', 'MH014', 'MH015', 'MH017']

    # pegRNA frequency - calculate Pearson and Spearman correlation
    ls_samples_freq_colnames = [x + '_pegRNA_freq' for x in ls_samples]
    df_freq_correlation_spearman = df_pegrna_data.loc[:, ls_samples_freq_colnames].corr(method='spearman', numeric_only=True)
    df_freq_correlation_pearson = df_pegrna_data.loc[:, ls_samples_freq_colnames].corr(method='pearson', numeric_only=True)
    
    # ST editing - calculate Pearson and Spearman correlation
    df_pegrna_data_filtered = utils.mask_st_editing_with_low_pegrna_count(df_pegrna_data, ls_samples[2:])
    ls_samples_st_colnames = [x + '_percentage_editing' for x in ls_samples]
    df_st_correlation_spearman = df_pegrna_data_filtered.loc[:, ls_samples_st_colnames].corr(method='spearman', numeric_only=True)
    df_st_correlation_pearson = df_pegrna_data_filtered.loc[:, ls_samples_st_colnames].corr(method='pearson', numeric_only=True)

    # visualise as heatmap
    plt.figure(figsize=(3.5, 3.5))
    heatmap_mask = np.triu(np.ones_like(df_freq_correlation_pearson, dtype=bool))
    sns_heatmap = sns.heatmap(df_freq_correlation_pearson, cmap='viridis', square=True,
                              vmin=0, vmax=1, linewidths=1)
                              #annot=True, fmt='.1f',
                              #mask=heatmap_mask)
    plt.title('pegRNA Frequency - Pearson Correlation')  # \n(filtered)
    plt.tight_layout()
    plt.savefig(sVisOutputDir + 'FigureS4A.svg', format='svg')
    plt.close()

    plt.figure(figsize=(3.79, 3.79))
    heatmap_mask = np.triu(np.ones_like(df_st_correlation_pearson, dtype=bool))
    sns_heatmap = sns.heatmap(df_st_correlation_pearson, cmap='viridis', square=True,
                              vmin=0, vmax=1, linewidths=1)
                              #annot=True, fmt='.1f',
                              #mask=heatmap_mask)
    plt.title('% Correct ST Editing - Pearson Correlation')  # \n(filtered)
    plt.tight_layout()
    plt.savefig(sVisOutputDir + 'FigureS4B.svg', format='svg')
    plt.close()
    return 0



