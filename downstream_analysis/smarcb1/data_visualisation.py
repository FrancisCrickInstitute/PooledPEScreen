


###----- IMPORTS ----------------------------------------------------------------------------------------------------###
import pandas as pd
import numpy as np
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

###----- CONSTANTS --------------------------------------------------------------------------------------------------###
sBaseDir = '/Users/hergerm/Documents/Data/NGS/230728_NovaSeq/'
sOutputDir = sBaseDir + 'output/'
sVisOutputDir = sOutputDir + 'figures/'

dict_CADD_ColorMap = {'SYNONYMOUS': '#4e77bb', 'STOP_GAINED': '#e83677', 'CANONICAL_SPLICE': '#f7a83e',
                      'NON_SYNONYMOUS': '#64bb97', 'SPLICE_SITE': '#243672', 'INTRONIC': '#8acdef',
                      '5PRIME_UTR': '#B5A695', '3PRIME_UTR': '#d091bf',
                      'UPSTREAM': '#B5A695', 'DOWNSTREAM': '#B5A695',
                      'TBD (SNV+)': '#a9a9a9', 'STOP_LOST': '#a9a9a9', 'CODON_DELETION': '#a9a9a9'}
dict_ClinVar_ColorMap = {'Conflicting': '#8c6320', 'Uncertain': '#d3a44c', 'Likely pathogenic': '#b97097',
                         'Pathogenic': '#7d2a54', 'Likely benign': '#506af3', 'Benign': '#1432f5'}
ST_Data_Color = '#9ebeff'  # blue
EL_Data_Color = '#9effa9'  # green
BoxPlotProps = {'boxprops': {'facecolor': 'none', 'edgecolor': 'black'},
                'medianprops': {'color': 'black', 'linewidth': 3, 'solid_capstyle': 'butt'},
                'whiskerprops': {'color': 'black'},
                'capprops': {'color': 'black'}}

###----- FUNCTIONS --------------------------------------------------------------------------------------------------###
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
    plt.savefig(sVisOutputDir + 'Figure1D.svg', format='svg')
    plt.close()
    return 0

def vis_scaffold_vs_st_editing(df_pegrna_data, df_pegrna_data_ctrl1, df_pegrna_data_ctrl2):
    """ code to reproduce Figure 1E and 1F """
    ls_scaffolds = ['original', 'FE', 'PRIDICT', 't-lock']
    # identify STs with high background editing
    background_st_editing_threshold = 2
    df_pegrna_data_ctrl1['percentage_editing_control_high'] = df_pegrna_data_ctrl1.filter(items=['percentage_editing']).ge(background_st_editing_threshold).any(axis=1)
    df_pegrna_data_ctrl2['percentage_editing_control_high'] = df_pegrna_data_ctrl2.filter(items=['percentage_editing']).ge(background_st_editing_threshold).any(axis=1)
    df_pegrna_data['percentage_editing_control_high'] = df_pegrna_data_ctrl1['percentage_editing_control_high'] | df_pegrna_data_ctrl2['percentage_editing_control_high']
    df_pegrna_data_filtered = df_pegrna_data.loc[(df_pegrna_data['percentage_editing_control_high'] == False) &
                                                 (df_pegrna_data['pegRNA_count'] >= 10), :]
    df_counts = df_pegrna_data_filtered.groupby('pegRNA_BC').agg({'pegRNA_BC': 'count'})
    ls_pegrna_bc_selection = df_counts.index[df_counts['pegRNA_BC'] == 5].tolist()
    df_pegrna_data_filtered_v2 = df_pegrna_data_filtered.loc[df_pegrna_data_filtered['pegRNA_BC'].isin(ls_pegrna_bc_selection), :]
    df_pegrna_data_filtered_v2.to_csv(sOutputDir + '_SuppTable_ScaffoldComparison', index=False)
    
    # visualise ST editing via boxplot
    plt.figure(figsize=(2, 2.2))
    sns_plot = sns.boxplot(data=df_pegrna_data_filtered_v2, x='scaffold_design', y='percentage_editing',
                           width=0.5, showfliers=False, color='white',  # showfliers=False
                           order=ls_scaffolds, **BoxPlotProps)
    sns_plot.set(title='', xlabel='Scaffold', ylabel='% Correct ST Edits')
    sns_plot.set_xticklabels(sns_plot.get_xticklabels(), rotation=45, horizontalalignment='right')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)
    sns_plot.spines[['right', 'top']].set_visible(False)
    sns_plot.set_ylim(0, 30)
    plt.tight_layout()
    plt.savefig(sVisOutputDir + 'Figure1E.svg', format='svg')
    plt.close()

    # calculate and visualise fold-changes
    ref_scaffold = 'original'
    query_scaffold = 'FE'
    dict_FCs = {}
    for pegrna_bc, df_group in df_pegrna_data_filtered.groupby('pegRNA_BC'):
        fc = np.nan
        l2fc = np.nan
        fe_editing = df_group.loc[df_group['scaffold_design'] == query_scaffold, 'percentage_editing'].values
        og_editing = df_group.loc[df_group['scaffold_design'] == ref_scaffold, 'percentage_editing'].values
        if ((len(fe_editing) == 1) & (len(og_editing) == 1)):
            if ((fe_editing[0] > 0) & (og_editing[0] > 0)):
                fc = fe_editing[0] / og_editing[0]
                l2fc = math.log2(fc)
        dict_FCs.update({pegrna_bc: [fc, l2fc]})
    df_FCs = pd.DataFrame.from_dict(dict_FCs, orient='index')
    df_FCs.columns = ['FC', 'L2FC']
    
    # boxplot of log2FC(FE/original) values
    plt.figure(figsize=(1.0, 1.6))
    sns_plot = sns.boxplot(data=df_FCs, y='L2FC', width=0.5, showfliers=False, color='white', **BoxPlotProps)
    sns.stripplot(data=df_FCs, y='L2FC', alpha=0.2, size=2, linewidth=0.2, edgecolor='black', zorder=0)
    sns_plot.axhline(0, ls='--', c='grey')
    sns_plot.spines[['right', 'top']].set_visible(False)
    sns_plot.set(title='', xlabel='', ylabel='Log2(FE/original)')
    plt.tight_layout()
    plt.savefig(sVisOutputDir + 'Figure1F.svg', format='svg')
    plt.close()
    return 0

def vis_snv_vs_onv_efficiencies(df_pegrna_data, df_correspondence):
    """ code to reproduce Figures 1G and 1H """
    ls_snv_mutant_indexes = list(set(df_correspondence['mutant_index_SNV']))
    ls_onv_mutant_indexes = list(set(df_correspondence['mutant_index_ONV']))
    df_pegrna_data_mod = df_pegrna_data.loc[df_pegrna_data['mutant index'].isin(ls_snv_mutant_indexes + ls_onv_mutant_indexes), :]
    df_pegrna_data_snv = df_pegrna_data.loc[df_pegrna_data['mutant index'].isin(ls_snv_mutant_indexes), :]
    df_pegrna_data_onv = df_pegrna_data.loc[df_pegrna_data['mutant index'].isin(ls_onv_mutant_indexes), :]
    df_top_pegrna_data_mod = df_pegrna_data_mod.sort_values('percentage_editing', ascending=False).drop_duplicates(['mutant index'], keep='first')
    df_top_pegrna_data_snv = df_pegrna_data_snv.sort_values('percentage_editing', ascending=False).drop_duplicates(['mutant index'], keep='first')
    df_top_pegrna_data_onv = df_pegrna_data_onv.sort_values('percentage_editing', ascending=False).drop_duplicates(['mutant index'], keep='first')
    # determine if PAM is edited
    df_top_pegrna_data_mod['PAM_edit'] = df_top_pegrna_data_mod.apply(lambda row: annotate_pam_edits_for_snvs(row), axis=1)
    df_top_pegrna_data_mod = annotate_pam_edits_for_onvs(df_top_pegrna_data_mod, df_correspondence)
    df_pegrna_data_mod['PAM_edit'] = df_pegrna_data_mod.apply(lambda row: annotate_pam_edits_for_snvs(row), axis=1)
    df_pegrna_data_mod = annotate_pam_edits_for_onvs(df_pegrna_data_mod, df_correspondence)
    # add to correspondence table
    ls_onv_pam_edits = []
    ls_snv_pam_edits = []
    for index, row in df_correspondence.iterrows():
        onv_pam_edit = df_top_pegrna_data_mod.loc[df_top_pegrna_data_mod['mutant index'] == row['mutant_index_ONV'], 'PAM_edit']
        snv_pam_edit = df_top_pegrna_data_mod.loc[df_top_pegrna_data_mod['mutant index'] == row['mutant_index_SNV'], 'PAM_edit']
        onv_pam_edit = onv_pam_edit.values[0] if (len(onv_pam_edit) == 1) else np.nan
        snv_pam_edit = snv_pam_edit.values[0] if (len(snv_pam_edit) == 1) else np.nan
        ls_onv_pam_edits.append(onv_pam_edit)
        ls_snv_pam_edits.append(snv_pam_edit)
    df_correspondence['PAM_edit_ONV'] = ls_onv_pam_edits
    df_correspondence['PAM_edit_SNV'] = ls_snv_pam_edits
    df_correspondence['PAM_edit_MNVnotSNV'] = df_correspondence['PAM_edit_ONV'] & ~df_correspondence['PAM_edit_SNV']

    # box plot of variant type and PAM edit vs ST editing
    plt.figure(figsize=(1.5, 1.7))
    sns_plot = sns.boxplot(data=df_top_pegrna_data_mod, x='variant_type', y='percentage_editing',  # df_top_pegrna_data_mod or df_pegrna_data_mod
                           width=0.5, hue='PAM_edit', showfliers=False, color='white', **BoxPlotProps)
    sns.stripplot(data=df_top_pegrna_data_mod, x='variant_type', y='percentage_editing',
                  hue='PAM_edit', jitter=True, dodge=True, zorder=0,
                  alpha=0.2, size=2, color='grey', linewidth=0.2, edgecolor='black')
    sns_plot.set(xlabel='Variant Type', ylabel='% Correct ST Editing (Max)')
    sns_plot.legend([], [], frameon=False)
    sns_plot.spines[['right', 'top']].set_visible(False)
    plt.tight_layout()
    plt.savefig(sVisOutputDir + 'Figure1G.svg', format='svg')
    plt.close()

    # SNV-ONV fold change in ST editing
    ls_onv_snv_ratios = []
    ls_snv_percentage_editing = []
    ls_onv_percentage_editing = []
    for indx, row in df_correspondence.iterrows():
        snv_mut_index = row['mutant_index_SNV']
        onv_mut_index = row['mutant_index_ONV']
        snv_percentage_editing = df_top_pegrna_data_mod.loc[df_top_pegrna_data_mod['mutant index'] == snv_mut_index, 'percentage_editing']
        onv_percentage_editing = df_top_pegrna_data_mod.loc[df_top_pegrna_data_mod['mutant index'] == onv_mut_index, 'percentage_editing']
        onv_snv_ratio = np.nan
        snv_pedited = np.nan
        onv_pedited = np.nan
        if ((len(snv_percentage_editing) == 1) & (len(onv_percentage_editing) == 1)):
            if ((snv_percentage_editing.values[0] > 0) & (onv_percentage_editing.values[0] > 0)):
                onv_snv_ratio = onv_percentage_editing.values[0] / snv_percentage_editing.values[0]
                snv_pedited = snv_percentage_editing.values[0]
                onv_pedited = onv_percentage_editing.values[0]
        ls_onv_snv_ratios.append(onv_snv_ratio)
        ls_snv_percentage_editing.append(snv_pedited)
        ls_onv_percentage_editing.append(onv_pedited)
    df_correspondence['snv_percentage_editing'] = ls_snv_percentage_editing
    df_correspondence['onv_percentage_editing'] = ls_onv_percentage_editing
    df_correspondence['onv_to_snv_ratio'] = ls_onv_snv_ratios
    df_correspondence['log2_onv_to_snv_ratio'] = df_correspondence['onv_to_snv_ratio'].apply(lambda x: np.math.log(x, 2))

    # box plot of log2-FC of ONV/SNV ratio
    plt.figure(figsize=(1.0, 1.6))
    sns_plot = sns.boxplot(data=df_correspondence, y='log2_onv_to_snv_ratio',
                           width=0.5, showfliers=False, color='white', **BoxPlotProps)
    sns.stripplot(data=df_correspondence, y='log2_onv_to_snv_ratio', hue='PAM_edit_MNVnotSNV', zorder=0,
                  alpha=0.5, size=2, linewidth=0.2, edgecolor='black')
    sns_plot.set(xlabel='', ylabel='Log2(ONV/SNV) - % Correct ST Edits')
    sns_plot.axhline(0, ls='--')
    sns_plot.spines[['right', 'top']].set_visible(False)
    plt.tight_layout()
    plt.savefig(sVisOutputDir + 'Figure1H.svg', format='svg')
    plt.close()
    return 0

def vis_pegrna_coverage_per_snv(df_pegrna_data, experiment):
    """ code to reproduce Figure 2C """
    df_pegrna_data_mod = df_pegrna_data.loc[df_pegrna_data['variant_type'] == 'SNP', :]
    # number pegRNAs per SNV
    df_counts = df_pegrna_data_mod.groupby('mutant index').agg({'mutant index': ['first', 'count'], 'start': 'first'})
    df_counts.columns = ['mutant index', 'pegRNA_count', 'start']
    df_averages = df_counts.groupby('start').agg({'start': 'first', 'pegRNA_count': 'mean'})
    df_averages.columns = ['start', 'mean_pegRNA_count']
    # visualise as bar plot
    ls_regions = [[24175725, 24175924], [24176289, 24176488]]  # GRCh37
    fig, (ax1, ax2) = plt.subplots(ncols=2, nrows=1, sharey=True, figsize=(7, 0.3))
    df_plot = df_averages.loc[df_averages['start'].between(ls_regions[0][0], ls_regions[0][1], inclusive='both'), :]
    color_palette = mpl.cm.ScalarMappable(cmap='viridis_r').to_rgba(df_plot['mean_pegRNA_count'])  # viridis_r
    sns_plot = sns.barplot(data=df_plot, x='start', y='mean_pegRNA_count',
                           width=1, linewidth=0,
                           palette=color_palette, ax=ax1)
    sns_plot.set(xlabel=None, ylabel='Mean Number pegRNAs', xticklabels=[])
    sns_plot.tick_params(bottom=False)
    sns_plot.spines[['right', 'top']].set_visible(False)
    df_plot = df_averages.loc[df_averages['start'].between(ls_regions[1][0], ls_regions[1][1], inclusive='both'), :]
    color_palette = mpl.cm.ScalarMappable(cmap='viridis_r').to_rgba(df_plot['mean_pegRNA_count'])
    sns_plot = sns.barplot(data=df_plot, x='start', y='mean_pegRNA_count',
                           width=1, linewidth=0,
                           palette=color_palette, ax=ax2)
    sns_plot.set(xlabel=None, ylabel=None, xticklabels=[])
    sns_plot.tick_params(left=False, bottom=False)
    sns_plot.spines[['left', 'right', 'top']].set_visible(False)
    sns_plot.legend([], [], frameon=False)
    ax1.set_ylim(0, 10)
    ax2.set_ylim(0, 10)
    plt.subplots_adjust(wspace=0.1, hspace=0)
    plt.savefig(sVisOutputDir + 'Figure2C.svg', format='svg')
    plt.close()
    return 0

def vis_positional_st_editing(df_pegrna_data):
    """ code to reproduce Figure 2D and S4C """
    df_pegrna_data_mod = df_pegrna_data.loc[(df_pegrna_data['percentage_editing_control_high'] == False) &
                                            (df_pegrna_data['variant_type'] == 'SNP'), :]
    dict_timepoints = {'D4': ['MH003', 'MH004'],
                       'D10': ['MH005', 'MH006'],
                       'D20': ['MH008', 'MH009'],
                       'D27': ['MH010', 'MH011'],
                       'D34': ['MH013', 'MH014']}
    timepoint = 'D10'
    ls_duplicate = dict_timepoints[timepoint]
    ls_count_colnames = [sample + '_pegRNA_count' for sample in ls_duplicate]
    ls_editing_colnames = [sample + '_percentage_editing' for sample in ls_duplicate]
    ls_region_boundaries = [[24175725, 24175924], [24176289, 24176488]]
    # filter out pegRNAs with count < 10
    df_pegrna_data_mod.loc[(df_pegrna_data_mod[ls_count_colnames[0]] < 10) |
                           (df_pegrna_data_mod[ls_count_colnames[1]] < 10), ls_editing_colnames] = np.nan
    # average ST editing efficiencies across duplicate
    df_pegrna_data_mod[timepoint + '_percentage_editing'] = df_pegrna_data_mod[ls_editing_colnames].mean(axis=1, skipna=True)
    # flag top pegRNA per variant
    df_pegrna_data_mod['top_pegRNA'] = False
    ls_top_pegrna_indexes = df_pegrna_data_mod.groupby('mutant index')[timepoint + '_percentage_editing'].idxmax()
    ls_top_pegrna_indexes = ls_top_pegrna_indexes.dropna()
    df_pegrna_data_mod.loc[ls_top_pegrna_indexes, 'top_pegRNA'] = True
    
    # visualise as scatter plot
    df_data_plot = df_pegrna_data_mod
    #df_data_plot = df_pegrna_data_mod.loc[df_pegrna_data_mod['top_pegRNA'] == True, :]  # -> uncomment for Figure S4C
    fig, (ax1, ax2) = plt.subplots(ncols=2, nrows=1, sharey=True, figsize=(7, 1.5))
    sns_plot = sns.scatterplot(data=df_data_plot, x='start', y=timepoint + '_percentage_editing', hue='mut_allele',
                               ax=ax1, s=5, linewidth=0.2, edgecolor='black', alpha=0.75, palette='Set2')
    sns_plot.set(xlabel='', ylabel='Percentage Correct ST Edits (' + timepoint + ')')
    sns_plot.spines[['right', 'top']].set_visible(False)
    sns_plot.legend([], [], frameon=False)
    sns_plot = sns.scatterplot(data=df_data_plot, x='start', y=timepoint + '_percentage_editing', hue='mut_allele',
                               ax=ax2, s=5, linewidth=0.2, edgecolor='black', alpha=0.75, palette='Set2')
    sns_plot.set(xlabel=None, ylabel=None)
    sns_plot.tick_params(left=False)
    sns_plot.spines[['left', 'right', 'top']].set_visible(False)
    sns_plot.legend([], [], frameon=False)
    #sns_plot.set(xlabel='Position (GRCh37)', ylabel='Percentage Correct ST Edits (' + timepoint + ')')
    #sns_plot.spines[['right', 'top']].set_visible(False)
    ax1.set_xlim(ls_region_boundaries[0][0]-5, ls_region_boundaries[0][1]+5)
    ax2.set_xlim(ls_region_boundaries[1][0]-5, ls_region_boundaries[1][1]+5)
    fig.text(0.5, 0.0, 'Position (GRCh37)', ha='center')
    plt.subplots_adjust(wspace=0.1, hspace=0)
    plt.savefig(sVisOutputDir + 'Figure2D.svg', format='svg')
    plt.close()
    return 0

def vis_features_vs_st_editing(df_pegrna_data):
    """ code to reproduce Figure 2E """
    dict_features = {'percentage_editing': 'percentage_editing', 'PBS GC content': 'PBS GC Content',
                     'PBS length': 'PBS Length', 'RTT length': 'RTT Length',
                     "distance mut to 5' RTT": 'RTT Overhang Length', 'distance to nick': 'Nick to Edit Distance',
                     'CRISPick_score': 'CRISPick Score', 'PRIDICT_editing_Score_deep': 'PRIDICT Score'}
    ls_colnames = ['percentage_editing', 'PBS GC content', 'PBS length', 'RTT length', "distance mut to 5' RTT", 'distance to nick', 'CRISPick_score', 'PRIDICT_editing_Score_deep']
    df_correlation_spearman = df_pegrna_data.loc[:, ls_colnames].corr(method='spearman', numeric_only=True)
    df_correlation_pearson = df_pegrna_data.loc[:, ls_colnames].corr(method='pearson', numeric_only=True)
    df_correlation = df_correlation_spearman
    df_correlation = df_correlation.rename(columns={'percentage_editing': 'correlation'})
    df_correlation = df_correlation.rename_axis('feature').reset_index()
    df_correlation['is_positive'] = df_correlation['correlation'] > 0
    df_correlation = df_correlation.loc[df_correlation['feature'].isin(ls_colnames[1:]), ['feature', 'correlation', 'is_positive']]
    df_correlation.sort_values(by='correlation', axis=0, ascending=False, inplace=True)
    df_correlation.replace({'feature': dict_features}, inplace=True)

    # visualise correlation via barplot
    plt.figure(figsize=(2.3, 1.5))
    sns_plot = sns.barplot(df_correlation, x='correlation', y='feature', hue='is_positive', dodge=False,
                           linewidth=0.2, edgecolor='black')
    sns_plot.spines[['right', 'top']].set_visible(False)
    sns_plot.set(xlabel='Spearman Correlation', ylabel='')
    sns_plot.set_xlim(-0.65, 0.65)
    plt.legend([], [], frameon=False)
    plt.tight_layout()
    plt.savefig(sVisOutputDir + 'Figure2E.svg', format='svg')
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

def vis_el_vs_st_data_correlation(df_pegrna_data, df_el_data):
    """ code to reproduce Figures 2G and S4F """
    ls_variant_types = ['SNP', 'ONP', 'INS', 'DEL'] 
    ls_timepoints = ['D4', 'D10', 'D20', 'D27', 'D34', 'D10obn', 'D34obn', 'D43obn']
    df_pegrna_data = df_pegrna_data.loc[df_pegrna_data['variant_type'].isin(ls_variant_types), :]
    df_el_data = df_el_data.loc[df_el_data['variant_type'].isin(ls_variant_types), :]
    # processing ST data
    df_pegrna_data_mod = df_pegrna_data.loc[df_pegrna_data['percentage_editing_control_high'] == False, :]
    ls_duplicates = [['MH003', 'MH004'], ['MH005', 'MH006'], ['MH008', 'MH009'], ['MH010', 'MH011'], ['MH013', 'MH014'],
                     ['MH007', 'MH007'], ['MH015', 'MH015'], ['MH017', 'MH017'],
                     ['MH012', 'MH012'], ['MH016', 'MH016']]
    ls_timepoint_st_freq_colnames = [timepoint + '_edited_freq' for timepoint in ls_timepoints]
    # convert 'edited_freq' to ratio of edited STs over all STs (rather than only edited STs)
    ls_sample_ids = ['MH003', 'MH004', 'MH005', 'MH006', 'MH007', 'MH008', 'MH009', 'MH010',
                     'MH011', 'MH012', 'MH013', 'MH014', 'MH015', 'MH016', 'MH017']
    for sample_id in ls_sample_ids:
        total_pegrna_counts = df_pegrna_data_mod[sample_id + '_pegRNA_count'].sum()
        df_pegrna_data_mod[sample_id + '_edited_freq'] = df_pegrna_data_mod[sample_id + '_edited_count'].apply(lambda x: x / total_pegrna_counts)
    # average ST editing across duplicates
    for n in range(0, len(ls_duplicates)):
        ls_samples = ls_duplicates[n]
        ls_st_freq_colnames = [sample + '_edited_freq' for sample in ls_samples]
        df_pegrna_data_mod[ls_timepoint_st_freq_colnames[n]] = df_pegrna_data_mod[ls_st_freq_colnames].mean(axis=1, skipna=True)  # True: inclusive
    # collapse to variant level
    df_st_data = df_pegrna_data_mod.groupby(by='mutant index')[ls_timepoint_st_freq_colnames].apply(lambda x: x.sum())
    
    # processing of ET data by averaging over duplicates
    freq_suffix = '_freq_norm'  # '_freq' OR '_freq_norm'
    df_el_data_mod = df_el_data.loc[df_el_data['high_background_frequency'] == False, :]
    ls_duplicates = [['453_D4', '454_D4'], ['453_D10', '454_D10'], ['453_D20', '454_D20'], ['453_D27', '454_D27'], ['453_D34', '454_D34'],
                     ['453_D10obn', '453_D10obn'], ['453_D34obn', '453_D34obn'], ['453_D43obn', '453_D43obn'],
                     ['454_D27tg', '454_D27tg'], ['454_D34tg', '454_D34tg']]
    ls_timepoint_norm_freq_colnames = [timepoint + freq_suffix for timepoint in ls_timepoints]
    for n in range(0, len(ls_duplicates)):
        ls_samples = ls_duplicates[n]
        ls_el_freq_colnames = [sample + freq_suffix for sample in ls_samples]
        df_el_data_mod[ls_timepoint_norm_freq_colnames[n]] = df_el_data_mod[ls_el_freq_colnames].mean(axis=1, skipna=True)  # True: inclusive
    # determine number of variants above background levels
    df_el_data_mod = df_el_data_mod.loc[df_el_data_mod['Consequence'] != 'TBD (SNV+)', :]
    df_el_data_above_bg = (df_el_data_mod == 0).sum(axis=0)

    # correlation analysis
    df_data_merged = pd.merge(df_el_data_mod, df_st_data, left_on='mutant index', right_index=True, how='left')
    ls_all_colnames = ls_timepoint_st_freq_colnames + ls_timepoint_norm_freq_colnames
    # select timepoints to visualise
    ls_timepoints_subset = ['D10', 'D10obn', 'D20', 'D34', 'D34obn']
    ls_selection_colnames = [timepoint + '_edited_freq' for timepoint in ls_timepoints_subset] + \
                            [timepoint + freq_suffix for timepoint in ls_timepoints_subset]
    df_data_merged_subset = df_data_merged.loc[:, ['mutant index'] + ls_selection_colnames]
    df_data_merged_subset.to_csv(sOutputDir + '_SMARCB1_STandET_variant_freq.csv', index=False)
    df_correlation_pearson = df_data_merged.loc[:, ls_selection_colnames].corr(method='pearson')
    df_correlation_spearman = df_data_merged.loc[:, ls_selection_colnames].corr(method='spearman')

    # heatmap of Spearman or Pearson correlation
    df_heatmap_data = df_correlation_spearman
    heatmap_mask = np.triu(np.ones_like(df_heatmap_data, dtype=bool))
    plt.figure(figsize=(2.9, 2.9))
    sns_heatmap = sns.heatmap(df_heatmap_data, cmap='viridis', square=True,
                              vmin=0, vmax=1, linewidths=1,
                              #annot=True, fmt='.1f',
                              mask=heatmap_mask)
    plt.title('Variant Frequencies - Spearman Correlation')  # \n(filtered)
    plt.tight_layout()
    plt.savefig(sVisOutputDir + 'Figure2G.svg', format='svg')
    plt.close()

    # scatter plot
    timepoint = 'D10obn'
    plt.figure(figsize=(2.2, 1.9))
    sns_plot = sns.scatterplot(data=df_data_merged,
                               x=timepoint + freq_suffix,
                               y=timepoint + '_edited_freq',
                               hue='variant_type',
                               s=4, alpha=0.6, linewidth=0.2, edgecolor='black',
                               palette='Set2')
    sns_plot.set(xlabel='EL Variant Frequency (' + timepoint + ')',
                 ylabel='ST Variant Frequency (' + timepoint + ')')
    sns_plot.spines[['right', 'top']].set_visible(False)
    plt.tight_layout()
    plt.savefig(sVisOutputDir + 'FigureS4F.svg', format='svg')
    plt.close()
    return 0

def vis_l2fc_means_by_st_threshold(df_pegrna_data, ls_sample_pairs):
    """ code to reproduce Figure 2H """
    freq_threshold = 0.6e-4
    ls_consequences = ['SYNONYMOUS', 'STOP_GAINED', 'CANONICAL_SPLICE', 'NON_SYNONYMOUS', 'SPLICE_SITE', 'INTRONIC', '3PRIME_UTR']
    df_pegrna_data_mod = df_pegrna_data.loc[df_pegrna_data['percentage_editing_control_high'] == False, :]
    df_pegrna_data_mod = df_pegrna_data_mod.loc[df_pegrna_data['Consequence'].isin(ls_consequences), :]
    # calculate pegRNA log2FCs
    sample_name = ls_sample_pairs[1] + '_' + ls_sample_pairs[0]
    df_pegrna_data_mod['log2FC_' + sample_name] = df_pegrna_data_mod.apply(
        lambda row: math.log(row[ls_sample_pairs[1] + '_pegRNA_pseudofreq'] / row[ls_sample_pairs[0] + '_pegRNA_pseudofreq'], 2), axis=1)
    df_pegrna_data_mod.loc[df_pegrna_data_mod[ls_sample_pairs[0] + '_pegRNA_pseudofreq'] < freq_threshold, 'log2FC_' + sample_name] = np.nan
    # filter pegRNAs at various ST thresholds
    st_step_size = 1.0
    ls_st_thresholds = [*np.arange(0, 100 + st_step_size, st_step_size)]
    ls_col_names = ['Consequence', 'log2FC_' + sample_name + '_mean', 'log2FC_' + sample_name + '_sd', 'count']
    df_pegrna_l2fc_means_merged = pd.DataFrame(columns=ls_col_names)
    df_variant_l2fc_means_merged = pd.DataFrame(columns=ls_col_names)
    for x in range(0, len(ls_st_thresholds)):
        st_threshold = ls_st_thresholds[x]
        df_pegrna_data_filtered = df_pegrna_data_mod.loc[df_pegrna_data_mod[ls_sample_pairs[0] + '_percentage_editing'] >= st_threshold, :]
        # calculate mean pegRNA-log2FCs per consequence
        df_pegrna_l2fc_means = df_pegrna_data_filtered.groupby(by='Consequence').agg({'Consequence': 'first',
                                                                                      'log2FC_' + sample_name: ['mean', 'std', 'count']})
        df_pegrna_l2fc_means.columns = ls_col_names
        df_pegrna_l2fc_means['st_threshold'] = st_threshold
        df_pegrna_l2fc_means_merged = pd.concat([df_pegrna_l2fc_means_merged, df_pegrna_l2fc_means], ignore_index=True)
        # collapse to variant-log2FCs
        df_variant_data_filtered = df_pegrna_data_filtered.groupby(by='mutant index').agg({'mutant index': 'first',
                                                                                           'Consequence': 'first',
                                                                                           'log2FC_' + sample_name: 'mean'})
        df_variant_l2fc_means = df_variant_data_filtered.groupby(by='Consequence').agg({'Consequence': 'first',
                                                                                        'log2FC_' + sample_name: ['mean', 'std', 'count']})
        df_variant_l2fc_means.columns = ls_col_names
        df_variant_l2fc_means['st_threshold'] = st_threshold
        df_variant_l2fc_means_merged = pd.concat([df_variant_l2fc_means_merged, df_variant_l2fc_means], ignore_index=True)
    # background-correction to synonymous distribution (w/o st filter)
    log2fc_syn_bg = df_pegrna_l2fc_means_merged.loc[(df_pegrna_l2fc_means_merged['st_threshold'] == 0) &
                                                    (df_pegrna_l2fc_means_merged['Consequence'] == 'SYNONYMOUS'),
                                                    'log2FC_' + sample_name + '_mean'].values[0]
    df_pegrna_l2fc_means_merged['log2FC_' + sample_name + '_mean'] -= log2fc_syn_bg
    log2fc_syn_bg = df_variant_l2fc_means_merged.loc[(df_variant_l2fc_means_merged['st_threshold'] == 0) &
                                                     (df_variant_l2fc_means_merged['Consequence'] == 'SYNONYMOUS'),
                                                     'log2FC_' + sample_name + '_mean'].values[0]
    df_variant_l2fc_means_merged['log2FC_' + sample_name + '_mean'] -= log2fc_syn_bg
    # calculate confidence interval
    conf_level = 0.95
    df_pegrna_l2fc_means_merged['log2FC_' + sample_name + '_ci_lower'] = df_pegrna_l2fc_means_merged.apply(lambda row: calculate_confidence_interval(row, sample_name, conf_level, lower=True), axis=1)
    df_pegrna_l2fc_means_merged['log2FC_' + sample_name + '_ci_upper'] = df_pegrna_l2fc_means_merged.apply(lambda row: calculate_confidence_interval(row, sample_name, conf_level, lower=False), axis=1)
    df_variant_l2fc_means_merged['log2FC_' + sample_name + '_ci_lower'] = df_variant_l2fc_means_merged.apply(lambda row: calculate_confidence_interval(row, sample_name, conf_level, lower=True), axis=1)
    df_variant_l2fc_means_merged['log2FC_' + sample_name + '_ci_upper'] = df_variant_l2fc_means_merged.apply(lambda row: calculate_confidence_interval(row, sample_name, conf_level, lower=False), axis=1)

    # visualise mean pegRNA scores as line plot
    plt.figure(figsize=(3.2, 2.2))
    sns_plot = sns.lineplot(data=df_pegrna_l2fc_means_merged,
                            x='st_threshold', y='log2FC_' + sample_name + '_mean',
                            hue='Consequence', palette=dict_CADD_ColorMap,
                            linewidth=1)
    for cons in ls_consequences:
        df_conf_int = df_pegrna_l2fc_means_merged.loc[df_pegrna_l2fc_means_merged['Consequence'] == cons, :]
        sns_plot.fill_between(x=df_conf_int['st_threshold'],
                              y1=df_conf_int['log2FC_' + sample_name + '_ci_lower'],
                              y2=df_conf_int['log2FC_' + sample_name + '_ci_upper'],
                              color=dict_CADD_ColorMap[cons], alpha=.15)
    sns_plot.axhline(0, ls='--', color='black', linewidth=1)
    sns_plot.set_xlim(0, 100)
    sns_plot.set(title=sample_name, xlabel='Minimum % Correct ST Edits', ylabel='pegRNA Log2FC')
    sns_plot.set_xticks([*range(0, 101, 10)], labels=map(str, [*range(0, 101, 10)]))
    sns_plot.spines[['top']].set_visible(False)
    ax2 = plt.twinx()
    df_n_pegrnas = df_pegrna_l2fc_means_merged.groupby(['st_threshold']).agg({'count': 'sum'})
    df_n_pegrnas.reset_index(inplace=True)
    sns_plot_2 = sns.lineplot(data=df_n_pegrnas, x='st_threshold', y='count',
                              c='black', linewidth=1, ax=ax2)
    sns_plot_2.spines[['right', 'top']].set_visible(False)
    plt.tight_layout()
    plt.savefig(sVisOutputDir + 'Figure2H.svg', format='svg')
    plt.close()
    return 0

def vis_pegrna_scores_by_consequence(df_data, col_name, plot_title):
    """ code to reproduce Figure 2I """
    ls_consequences = ['SYNONYMOUS', 'STOP_GAINED', 'CANONICAL_SPLICE', 'NON_SYNONYMOUS', 'SPLICE_SITE',
                       'INTRONIC', '3PRIME_UTR']
    df_data = df_data.loc[df_data['Consequence'].isin(ls_consequences), :]

    # variant score vs CADD consequence
    plt.figure(figsize=(4, 3))
    sns_plot = sns.boxplot(data=df_data, x='Consequence', y=col_name,
                           width=0.5, showfliers=False, color='white',
                           order=ls_consequences, **BoxPlotProps)
    sns.stripplot(data=df_data, x='Consequence', y=col_name, hue='Consequence',
                  palette=dict_CADD_ColorMap, order=ls_consequences,
                  alpha=0.5, size=2, linewidth=0.2, edgecolor='black',
                  zorder=0)
    sns_plot.set(title=plot_title, xlabel='Consequence', ylabel='Log2FC')
    sns_plot.spines[['right', 'top']].set_visible(False)
    sns_plot.set_xticklabels(sns_plot.get_xticklabels(), rotation=45, horizontalalignment='right')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)
    plt.tight_layout()
    #plt.show()
    plt.savefig(sVisOutputDir + 'Figure2I.svg', format='svg')
    plt.close()
    return 0

def vis_variant_scores_by_consequence(df_data, col_name, plot_title):
    """ code to reproduce Figure 2J """
    ls_consequences = ['SYNONYMOUS', 'STOP_GAINED', 'CANONICAL_SPLICE', 'NON_SYNONYMOUS', 'SPLICE_SITE',
                       'INTRONIC', '3PRIME_UTR']
    df_data = df_data.loc[df_data['Consequence'].isin(ls_consequences), :]

    # variant score vs CADD consequence
    plt.figure(figsize=(4, 3))
    sns_plot = sns.boxplot(data=df_data, x='Consequence', y=col_name,
                           width=0.5, showfliers=False, color='white',
                           order=ls_consequences, **BoxPlotProps)
    sns.stripplot(data=df_data, x='Consequence', y=col_name, hue='Consequence',
                  palette=dict_CADD_ColorMap, order=ls_consequences,
                  alpha=0.5, size=2, linewidth=0.2, edgecolor='black',
                  zorder=0)
    sns_plot.set(title=plot_title, xlabel='Consequence', ylabel='Log2FC')
    sns_plot.spines[['right', 'top']].set_visible(False)
    sns_plot.set_xticklabels(sns_plot.get_xticklabels(), rotation=45, horizontalalignment='right')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)
    plt.tight_layout()
    plt.savefig(sVisOutputDir + 'Figure2J.svg', format='svg')
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
    df_pegrna_data_filtered = mask_st_editing_with_low_pegrna_count(df_pegrna_data, ls_samples[2:])
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

def vis_pridict_vs_st_correlation(df_pegrna_data):
    """ code to reproduce Figures S4D and S4E """
    # box plot
    plt.figure(figsize=(1.9, 2.1))
    df_pegrna_data['PRIDICT_bin'] = pd.cut(df_pegrna_data['PRIDICT_editing_Score_deep'],
                                           [0, 40, 50, 60, 70, 80, 90, 100],
                                           labels=['0 to 40', '40 to 50', '50 to 60',
                                                   '60 to 70', '70 to 80', '80 to 90',
                                                   '90 to 100'])
    sns_plot = sns.boxplot(data=df_pegrna_data, x='PRIDICT_bin', y='percentage_editing',
                           showfliers=False, color='grey', **BoxPlotProps)
    sns_plot.set(xlabel='PRIDICT score', ylabel='% Correct ST Edits')
    sns_plot.set_xticklabels(sns_plot.get_xticklabels(), rotation=45, horizontalalignment='right')
    sns_plot.spines[['right', 'top']].set_visible(False)
    plt.tight_layout()
    plt.savefig(sVisOutputDir + 'FigureS4D.svg', format='svg')
    plt.close()
    
    # scatter plot
    plt.figure(figsize=(2.2, 1.8))
    df_pegrna_data['PRIDICT_bin_2'] = pd.cut(df_pegrna_data['PRIDICT_editing_Score_deep'],
                                             [0, 70, 100], labels=['0 to 70', '70 to 100'])
    sns_plot = sns.scatterplot(data=df_pegrna_data, x='PRIDICT_editing_Score_deep', y='percentage_editing',
                               hue='PRIDICT_bin_2', s=4, alpha=0.5, linewidth=0.1, edgecolor='black')
    sns_plot.set(xlabel='PRIDICT score', ylabel='% correct ST edits')
    sns_plot.spines[['right', 'top']].set_visible(False)
    sns_plot.set_ylim(0, 100)
    sns_plot.set_xticks([0, 40, 50, 60, 70, 80, 90], ['0', '40', '50', '60', '70', '80', '90'])
    plt.tight_layout()
    plt.savefig(sVisOutputDir + 'FigureS4E.svg', format='svg')
    plt.close()
    return 0

def vis_scores_sample_correlation(df_pegRNA_Data_Log2FCs, df_Variant_Data_Log2FCs, ls_log2fc_colnames):
    """ code to reproduce Figure S5 """
    # pre-filter data on consequences
    ls_consequences = ['SYNONYMOUS', 'STOP_GAINED', 'CANONICAL_SPLICE', 'NON_SYNONYMOUS', 'SPLICE_SITE', 'INTRONIC', '3PRIME_UTR']
    df_pegRNA_Data_Log2FCs = df_pegRNA_Data_Log2FCs.loc[df_pegRNA_Data_Log2FCs['Consequence'].isin(ls_consequences), :]
    df_Variant_Data_Log2FCs = df_Variant_Data_Log2FCs.loc[df_Variant_Data_Log2FCs['Consequence'].isin(ls_consequences), :]

    # replicate correlation plot for pegRNA scores
    plt.figure(figsize=(2, 2))
    sns_plot = sns.pairplot(df_pegRNA_Data_Log2FCs,
                            x_vars=ls_log2fc_colnames,
                            y_vars=ls_log2fc_colnames,
                            height=1.5, aspect=0.62,
                            kind='scatter', diag_kind='hist', plot_kws={'s': 15, 'alpha': 0.75},
                            hue='Consequence', palette=dict_CADD_ColorMap,
                            corner=True)
    sns_plot.set(xlim=(-8, 8), ylim=(-8, 8))
    plt.legend(loc='upper right')
    plt.tight_layout()
    plt.savefig(sVisOutputDir + 'FigureS5A.svg', format='svg')
    plt.close()
    
    # replicate correlation plot for variant scores
    sns_plot = sns.pairplot(df_Variant_Data_Log2FCs,
                            x_vars=ls_log2fc_colnames,
                            y_vars=ls_log2fc_colnames,
                            height=1.5, aspect=0.62,
                            kind='scatter', diag_kind='hist', plot_kws={'s': 15, 'alpha': 0.75},
                            hue='Consequence', palette=dict_CADD_ColorMap,
                            corner=True)
    sns_plot.set(xlim=(-8, 5), ylim=(-8, 5))
    plt.legend(loc='upper right')
    plt.tight_layout()
    plt.savefig(sVisOutputDir + 'FigureS5C.svg', format='svg')
    plt.close()
    return 0

###----- UTILITY FUNCTIONS ------------------------------------------------------------------------------------------###
def annotate_pam_edits_for_snvs(pegrna_row):
    """ function to flag SNVs within PAM """
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
    """ function to flag MNVs with edit within PAM """
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

def calculate_confidence_interval(row, sample_name, conf_level, lower):
    """ function to calculate confindence interval for given input """
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

def mask_st_editing_with_low_pegrna_count(df_pegrna_data, ls_sample_ids, count_threshold=10):
    """ function to mask all ST editing values for pegRNAs with low count """
    df_pegrna_data_masked = df_pegrna_data
    for sample_id in ls_sample_ids:
        count_colname = sample_id + '_pegRNA_count'
        st_colname = sample_id + '_percentage_editing'
        df_pegrna_data_masked.loc[df_pegrna_data_masked[count_colname] < count_threshold, st_colname] = np.nan
    return df_pegrna_data_masked
