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
import scipy.stats

class NotEnoughSpaceException(Exception):

    msg = '''Not enough space on target sequence to design requested number of
          probes. Consider designing fewer probes or relaxing your design
          constraints.'''

    def __init__(self, message=msg):
        self.message = message

    def __str__(self):
        return f'{self.message}'

class AllPotentialProbesRemovedException(Exception):

    msg = '''All potential probes have been removed, for example by setting
          target subregions as too narrow a range. Consider widening the target
          subregion range.'''

    def __init__(self, message=msg):
        self.message = message

    def __str__(self):
        return f'{self.message}'

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
        maxinteraction = sorted(l, key = lambda x: x[0])[0]
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

def thin_bad_peaks(df, bad_peak_positions, window_size = 11):
    '''
    Given the location for the bad peaks, return df with the potential probes
    within given window size removed. This will hopefully allow new peaks without
    dimer clashes to be found.
    window_size: should be an odd number, e.g. if 11, then potential probes from
    5 nt 5' and 3 nt 3' will be removed around the bad peak.
    '''
    halfwin = int(math.floor((window_size - 1)/2))
    starts_to_remove = []
    for i in bad_peak_positions:
        starts_to_remove.extend([*range(i - halfwin, i + halfwin + 1)])

    #set Tm of starts_to_remove to 0 so that they will not be found as peaks
    df.loc[df['start'].isin(starts_to_remove), 'Tm'] = 0
    return df

def find_peaks_and_screen(df, this_distance, min_dimer_dG, all_selected_probes):
    '''
    Find Tm peaks within the set of probes.
    Remove probes that have heterodimer clashes with other probes in the set.
    '''
    data = df['Tm'].values
    maxes, properties = signal.find_peaks(data, distance = this_distance)
    df['Tm_peak'] = df.index.isin(maxes)
    #screen the peaks to remove ones that cause dimers with the others:
    screened_df = remove_bad_probes(df[df['Tm_peak']], min_dimer_dG, all_selected_probes)
    return df, screened_df

def peaks_testval(peak_locs, range, nprobes, half_probelen):
    '''
    This is a measure of peak distribution evenness.
    One metric is the number of nonzero bins in a range.
    Another metric, currently in use, is the largest gap between peak starts.
    '''
    mid_locs = peak_locs + half_probelen

    #This metric returns the number of nonzero bins in range
    #nonzero = np.count_nonzero(np.histogram(mid_locs, range = range, bins = nprobes)[0])
    #return nonzero

    #special cases for 1 and 2 probes desired:
    if nprobes == 1:
        #we want to minize the distance of a probe to the midpt
        mid = math.ceil(range[1] - range[0]/2)
        dist = [abs(i - mid) for i in mid_locs]
        return min(dist)

    elif nprobes == 2:
        #we want to minimize the distance to the ends of the target
        end_dist_sum = sum([mid_locs[0] - range[0], range[1] - mid_locs[-1]])
        return end_dist_sum

    else:
        #This metric returns the largest gap between probe mid pts
        dist = [t - s for s, t in zip(mid_locs, mid_locs[1:])]
        end_dist = [mid_locs[0] - range[0], range[1] - mid_locs[-1]]
        dist.extend(end_dist)
        return max(dist)

def prune(df, desired_number_probes, target_len, min_dimer_dG, all_selected_probes, subregions = None):
    '''
    Prune the probes into the desired number of probes per target.
    - find Tm peaks
    - get N evenly spaced probes
    '''

    #choose the highest Tm probe at each start site:
    idx = df.groupby(['start'])['Tm'].transform(max) == df['Tm']
    tm_df = df[idx].copy()
    tm_df['unique_id'] = tm_df.index

    if not probe_helpers.range_defined(subregions):
        #make 0-based, closed to match format of get_subregion_ranges()
        subregions = np.array([[0, target_len - 1]])

    #split the desired number of probes between the subregions
    this_probeset_size =  int(math.ceil(desired_number_probes/len(subregions)))
    chosen_probes = []
    sorted_subregions = subregions[np.argsort(subregions[:, 0])]
    #add 1 to the endpts to make half-open
    sorted_subregions[:,1] += 1
    substring = ', '.join(['%s-%s' % (i[0] + 1, i[1]) for i in sorted_subregions])
    logging.info('Choosing probes in subregions %s' % substring)

    for i, subregion in enumerate(subregions):
        #get the mini df that contains the data in the subregion
        sub_df = tm_df[(subregion[0] <= tm_df['start']) & (tm_df['end'] <= subregion[1] + 1)]
        if sub_df.empty:
            error = AllPotentialProbesRemovedException()
            logging.info(error.message)
            raise error

        #Add the missing start positons back to but with Tm of 0
        #This way, they will be included in the distance consideration for peak finding
        #but can't be chosen as peaks themselves
        #add a point before and after endpts so endpts can be chosen as peaks
        start_range = range(sub_df['start'].min() - 1, sub_df['start'].max()+ 2)
        if this_probeset_size < 3:
            this_distance = 20
        else:
            this_distance = math.ceil((start_range[-1] - start_range[0] - 2)/(this_probeset_size*3))
            if this_distance < 20:
                this_distance = 20
        range_df = pd.DataFrame(start_range, columns = ['start'])
        new_df = pd.merge(range_df[['start']], sub_df[['unique_id', 'Tm', 'start', 'length', 'sequence']], 'outer', on = 'start')
        new_df['Tm'].fillna(0, inplace = True)
        half_probelen = int(round(new_df['length'].mean()/2, 0))

        #Find peaks and screen for dimers. Drop peaks which cause dimers and repeat.
        #Choose the set of peaks with best testval (i.e. min distance between neighboring peaks).
        new_df, screened_df = find_peaks_and_screen(new_df, this_distance, min_dimer_dG, all_selected_probes)

        #1) Get original and screened peaks
        prescreened = new_df[new_df['Tm_peak']]['start'].values
        screened = screened_df['start'].values
        screened = np.sort(screened)
        screened_attempts = [screened]
        p = 1
        num_tests = 5
        while p < num_tests:
            #find peaks that were screened out in the dimer removal step:
            bad_peaks = set(prescreened).difference(set(screened))
            #remove probes from around the bad peaks and find new peaks:
            new_df = thin_bad_peaks(new_df, bad_peaks)
            new_df, screened_df = find_peaks_and_screen(new_df, this_distance, min_dimer_dG, all_selected_probes)
            prescreened = new_df[new_df['Tm_peak']]['start'].values
            screened = screened_df['start'].values
            screened = np.sort(screened)
            screened_attempts.append(screened)
            p += 1

        allvals = [peaks_testval(i, (subregion[0], subregion[1]), this_probeset_size, half_probelen) for i in screened_attempts]
        #get index max if the value is good, like nonzero bins
        #best_attempt = allvals.index(max(allvals))
        #get index min if the value is bad, like for largest distance between points
        best_attempt = allvals.index(min(allvals))
        peak_locs = screened_attempts[best_attempt]
        logging.info('%s Tm peaks found.' % len(peak_locs))
        #get optimal spacing for desired number of probes
        if len(peak_locs) < this_probeset_size:
            error = NotEnoughSpaceException()
            logging.info(error.message)
            raise error

        #remove edgecases
        #if only want 1 probe, choose one closest to the middle
        if len(peak_locs) == 1 and this_probeset_size == 1:
            chosen_locs = peak_locs
        elif this_probeset_size == 1:
            mid_dist = abs((peak_locs + half_probelen - subregion[0])/(subregion[-1] - subregion[0] + 1) - 0.5)
            chosen_locs = [peak_locs[np.argmin(mid_dist)]]
        elif this_probeset_size == 2:
            #for two probes, choose ones closest to the ends
            chosen_locs = [peak_locs[0], peak_locs[-1]]
        else:
            chosen_locs = probe_site_selection.choose_combination(peak_locs, this_probeset_size)
        chosen_ids = new_df.loc[new_df['start'].isin(chosen_locs), 'unique_id'].values
        chosen_df = df[df.index.isin(chosen_ids)].copy()
        #After the final choice, append to chosen probes and all_selected_probes
        chosen_probes.append(chosen_df)
        all_selected_probes = all_selected_probes.append(chosen_df)

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

    target_names = [os.path.basename(i).split('.fa')[0] for i in target_fastas]
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

        #Get evenly spaced probes from the passed ones, also screening for heterodimers
        pruned_df = prune(df, desired_number_probes[i], target_lens[i], min_dimer_dG[i], all_selected_probes, subregions = target_subregions)
        logging.info("%s probes selected." % len(pruned_df))
        all_selected_probes = all_selected_probes.append(pruned_df)
        summarize_results(df, pruned_df, target_lens[i], selected_probes_plots[i])

    #write the combined output file with probes selected for all targets.
    #add the dimer dG and dimer partner for all probes:
    all_selected_probes[['dimer_dG', 'dimer_partner']] = calc_dimer(all_selected_probes)
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
