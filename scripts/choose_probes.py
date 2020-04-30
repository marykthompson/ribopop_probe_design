'''
choose_probes.py
Choose probes across target sets, as evenly spaced as possible.
Check for heterodimer clashes when adding probes to the set.
'''
import sys
import numpy as np
from Bio import SeqIO
import pandas as pd
import os
import math
import primer3
import matplotlib.pyplot as plt
import logging
from scipy import signal
import probe_site_selection
import probe_helpers

class NotEnoughSpaceException(Exception):
    '''
    Not enough space on target sequence to design requested number of probes.
    Consider designing fewer probes or relaxing your design constraints.
    '''
    pass

def dist_from_neighbors(s):
    '''
    Get simple metric of distance from neighbors.
    A higher number means that neighbors are further away.
    '''
    s2 = s.sort_values(ascending = True)
    idx = s2.index
    a = s2.values
    dist = a[1:] - a[:-1]
    rdist = np.append(dist, dist[-1])
    ldist = np.append(dist[0], dist)
    m = np.mean([rdist, ldist], axis = 0)
    new_s = pd.Series(m, index = idx)
    return new_s

def calc_dimer(df):
    '''
    Calculate the dimer dG between each probe and each other probe. Add the
    min dG to the dataframe, along with the index of the dimer partner.
    '''
    df['index_num'] = df.index
    a = df[['index_num', 'sequence']].to_numpy()
    max_ints = []
    for i in range(0, len(a)):
        l = []
        for j in range(0, len(a)):
            #This includes both homodimers and heterodimers.
            l.append((primer3.calcHeterodimer(a[i][1], a[j][1], mv_conc = 300).dg/1000, a[j][0]))
        maxinteraction = min(l)
        max_ints.append(maxinteraction)
    dimer_dG = pd.DataFrame(max_ints, index = df.index)
    return dimer_dG

def remove_bad_probes(df, dimer_min_dG, all_selected_probes, filter = True):
    '''
    Returns a probeset that doesn't have any pairs with dimer dG less than min_dimer_dG.
    Compare not only to the probes in that target set but also the probes in the other
    target sets that have already been selected.
    If filter = False, simply return the df with the calculated dimer dGs.
    '''
    combo_df = df.append(all_selected_probes)
    combo_df[['dimer_dG','dimer_partner']] = calc_dimer(combo_df)
    #don't modify or drop values from the original df
    df = df.copy()
    df[['dimer_dG', 'dimer_partner']] = combo_df[['dimer_dG','dimer_partner']]
    df.sort_values('dimer_dG', inplace = True)
    dg_col = df.columns.get_loc('dimer_dG')

    if filter == False:
        return df

    while df.iloc[0, df.columns.get_loc('dimer_dG')] < dimer_min_dG:
        #Get a measure of the probes distance from its neighbors
        dist_from_nbrs = dist_from_neighbors(df['start'])
        #check if the top two dG values are the same (means from the same target)
        if df.iloc[0, dg_col] == df.iloc[1, dg_col]:
            indices = df.iloc[[0,1]].index
            #choose between the top 2 rows to find one with the lower distance between neighbors
            lower_dist_idx = dist_from_nbrs.loc[indices].idxmin()
            df.drop(lower_dist_idx, inplace = True)
        #otherwise, drop the probe with the most negative dG
        else:
            df.drop(df.index[0], inplace = True)
        df.drop('dimer_dG', axis = 1, inplace = True)

        #try again with all the problematic probes removed
        combo_df = df.append(all_selected_probes)
        combo_df[['dimer_dG','dimer_partner']] = calc_dimer(combo_df)
        df[['dimer_dG', 'dimer_partner']] = combo_df[['dimer_dG','dimer_partner']]
        df.sort_values('dimer_dG', inplace = True)
    return df

def prune(df, desired_number_probes, target_len, subregions = None):
    '''
    Prune the probes into the desired number of probes per target.
    - find Tm peaks
    - get N evenly spaced probes
    '''

    error_message = '''\
    Not enough space to design requested number of probes.
    Consider designing fewer probes or relaxing your design constraints.'''

    #choose the highest Tm probe at each start site:
    idx = df.groupby(['start'])['Tm'].transform(max) == df['Tm']
    tm_df = df[idx].copy()
    tm_df['unique_id'] = tm_df.index

    if not probe_helpers.range_defined(subregions):
        #make 0-based, half-open to match format of get_subregion_ranges()
        subregions = np.array([[0, target_len - 1]])

    #split the desired number of probes betwen the subregions
    this_probeset_size =  int(math.ceil(desired_number_probes/len(subregions)))
    chosen_probes = []
    sorted_subregions = subregions[np.argsort(subregions[:, 0])]

    #add 1 to the endpts to make half-open
    sorted_subregions[:,1] += 1
    substring = ', '.join(['%s-%s' % (i[0] + 1, i[1]) for i in sorted_subregions])
    logging.info('Choosing probes in subregions %s' % substring)

    for i, subregion in enumerate(subregions):
        #get the mini df that contains the data in the subregion
        sub_df = tm_df[(subregion[0] <= tm_df['target_start']) & (tm_df['target_end'] <= subregion[1])]

        #Add the missing start positons back to but with Tm of 0
        #This way, they will be included in the distance consideration for peak finding
        #but can't be chosen as peaks themselves
        this_distance = 100

        #start earlier because it cannot choose the endpts
        start_range = range(sub_df['start'].min() - 1, sub_df['start'].max()+ 2)
        range_df = pd.DataFrame(start_range, columns = ['start'])
        new_df = pd.merge(range_df[['start']], sub_df[['unique_id', 'Tm', 'start']], 'outer', on = 'start')
        new_df['Tm'].fillna(0, inplace = True)

        #find Tm peaks
        data = new_df['Tm'].values
        maxes, properties = signal.find_peaks(data, distance = this_distance)
        new_df['Tm_peak'] = new_df.index.isin(maxes)
        peak_locs = new_df.loc[new_df['Tm_peak'], 'start'].values
        logging.info('%s Tm peaks found.' % len(peak_locs))
        #get optimal spacing for desired number of probes
        if len(peak_locs) < this_probeset_size:
            logging.info(error_message)
            raise NotEnoughSpaceException(error_message)

        #remove edgecases
        #if only want 1 probe, choose one closest to the middle
        if len(peak_locs) == 1 and this_probeset_size == 1:
            chosen_locs = peak_locs
        elif len(peak_locs) > 1 and this_probeset_size == 1:
            mid_dist = abs((peak_locs - subregion[0])/(subregion[-1] - subregion[0] + 1) - 0.5)
            chosen_locs = [peak_locs[np.argmin(mid_dist)]]
        else:
            chosen_locs = probe_site_selection.choose_combination(peak_locs, this_probeset_size)
        chosen_ids = new_df.loc[new_df['start'].isin(chosen_locs), 'unique_id'].values
        chosen_probes.append(df[df.index.isin(chosen_ids)].copy())

    #combine probes from each subregion into the pruned_df
    pruned_df = pd.concat(chosen_probes)
    return pruned_df

def summarize_results(df, final_df, target_len, outfile):

    '''
    Write the final selected probes, plot the selected regions.
    '''
    fig = plt.figure(figsize = (5, 2.5))
    ax = fig.add_subplot(111)

    grey = '#A5AA99'
    pink = '#CF1C90'

    pre_tm_cols = ['passed_masking', 'passed_sequence', 'passed_structure']
    df['midpt'] = df['target_start'] + (df['length'] - 1)/2
    df.sort_values(by = 'midpt', ascending = True, inplace = True)

    bg = ax.scatter(df['midpt'], df['Tm'], s = 30,
    alpha = 0.3, color = grey, edgecolors = 'none')

    mini_df = df[df.index.isin(final_df.index)].copy()
    selected = ax.scatter(mini_df['midpt'], mini_df['Tm'], s = 30, alpha = 0.3, color = pink, edgecolors = 'none')

    ax.set_xlim(0, target_len)
    ax.set_ylabel('Tm')
    ax.set_xlabel('target position (nt)')

    ax.legend([bg, selected], ['before selection', 'selected'],
           mode = 'expand', fontsize = 8, ncol = 3, bbox_to_anchor=(0., 1.05, 1., .105), loc=3,
           borderaxespad=0., handletextpad=0.1)

    plt.tight_layout()

    plt.savefig(outfile, dpi = 600)

def main(arglist):
    '''
    Pick evenly-space probes corresponding to Tm peaks.
    Remove probes that heterodimerize with other probes.
    '''
    #columns to output after analysis
    col_order = ['sequence', 'target_name',  'target_start',  'target_end', 'length',
    'unique_id', 'Tm', 'GC_content', 'A_content', 'C_content', 'rolling_Tm_quantile_co',
    'hairpin_dG', 'homodimer_dG', 'dimer_dG', 'dimer_partner']

    probe_csvs = snakemake.input['probe_csvs']
    target_fastas = snakemake.input['target_fastas']
    excluded_regions = snakemake.input['excluded_regions']
    logfile = snakemake.params['logfile']
    desired_number_probes = snakemake.params['desired_number_probes']
    target_subregions_consensus = snakemake.params['target_subregions_consensus']
    min_dimer_dG = snakemake.params['min_dimer_dG']
    selected_probes_plots = snakemake.output['plots']
    all_selected_probes_file = snakemake.output['all_selected_probes']

    num_targets = len(target_fastas)
    target_names = [os.path.basename(i).rstrip('.fa') for i in target_fastas]
    target_lens = [len(next(SeqIO.parse(i, 'fasta'))) for i in target_fastas]

    all_selected_probes = pd.DataFrame()

    logging.basicConfig(level=logging.DEBUG, filename = logfile,
    filemode = 'w', format = '%(message)s')

    #stop writing all the font warnings to the log file
    logging.getLogger('matplotlib.font_manager').disabled = True

    for i, target in enumerate(target_names):
        print('target', target)
        logging.info('Target %s: ' % target)

        df = pd.read_csv(probe_csvs[i], index_col = 'unique_id')

        exdf = pd.read_csv(excluded_regions[i])
        #choose the subregion ranges to be used. If there are ranges provided wrt consensus, replace the the calculated ones
        target_subregions = exdf.loc[exdf['region_type'] == 'target_subregions', ['start', 'end']].values
        if not pd.isnull(target_subregions_consensus[i]):
            target_subregions = probe_helpers.get_subregion_ranges(target_subregions_consensus[i])

        logging.info("Starting with %s potential probes." % len(df))

        #Prune the remaining probes to the desired number of probes.
        removed_idx = [1]
        p = 0
        while len(removed_idx) > 0:
            #get evenly spaced probes from the passed ones
            pruned_df = prune(df, desired_number_probes[i], target_lens[i], subregions = target_subregions)
            #drop rows causing the low dimer dG scores
            logging.info('removing probes with dimer clashes')
            screened_df = remove_bad_probes(pruned_df, min_dimer_dG[i], all_selected_probes)
            #get indices of the pruned probes that don't survive heterodimer screening
            removed_idx = pruned_df.index.difference(screened_df.index)
            #drop bad heterodimer probes, will trigger rerun of prune with this probe removed.
            df.drop(removed_idx, inplace = True)
            p += 1

        #Get the properties of the selected probes and write to output file
        final_df = pd.merge(df.reindex(pruned_df.index),
        screened_df[['dimer_dG','dimer_partner']], left_index = True, right_index = True)
        logging.info("%s probes selected." % len(final_df))
        all_selected_probes = all_selected_probes.append(final_df)
        summarize_results(df, final_df, target_lens[i], selected_probes_plots[i])

    #write the combined output file with probes selected for all targets.
    all_selected_probes.sort_values(by = ['target_name', 'target_start'], inplace = True)
    all_selected_probes.reset_index(inplace = True)
    all_selected_probes['probe_num'] = all_selected_probes.index
    all_selected_probes['probe_num'] += 1
    cols = ['probe_num']
    cols.extend(col_order)
    other_cols = df.columns.values
    rule_cols = [i for i in other_cols if i.endswith('rule')]
    cols.extend(rule_cols)
    all_selected_probes[cols].round(2).to_csv(all_selected_probes_file, index = False)

if __name__ == '__main__':
    main(sys.argv[1:])
