'''
190430 MKT
Design probes to bind a set of targets.
'''
import sys
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import os
import math
import primer3
import subprocess
import argparse
import matplotlib.pyplot as plt
from Bio.SeqUtils import MeltingTemp as mt
import itertools
from collections import defaultdict
from collections.abc import Iterable
import logging
from scipy import signal
import probe_site_selection

class NotEnoughSpaceException(Exception):
    '''
    Not enough space on target sequence to design requested number of probes.
    Consider designing fewer probes or relaxing your design constraints.
    '''
    pass

def range_defined(ranges):
    '''
    Check whether ranges is defined. Let's us differentiate between nan/None and
    a potential list of ranges.
    '''
    if isinstance(ranges, Iterable):
        return(not pd.isnull(ranges).all())
    else:
        return(not pd.isnull(ranges))

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
    df[['dimer_dG', 'dimer_partner']] = combo_df[['dimer_dG','dimer_partner']]
    df.sort_values('dimer_dG', inplace = True)

    if filter == False:
        return df
    while df.iloc[0, df.columns.get_loc('dimer_dG')] < dimer_min_dG:
        #choose between the top 2 rows to find one with the lower Tm
        lower_tm_idx = df.iloc[[0, 1], df.columns.get_loc('Tm')].idxmin()
        df.drop(lower_tm_idx, inplace = True)
        df.drop('dimer_dG', axis = 1, inplace = True)

        #try again with all the problematic probes removed
        combo_df = df.append(all_selected_probes)
        combo_df[['dimer_dG','dimer_partner']] = calc_dimer(combo_df)
        df[['dimer_dG', 'dimer_partner']] = combo_df[['dimer_dG','dimer_partner']]
        df.sort_values('dimer_dG', inplace = True)
    return df

def Tm_from_position(seq, length, salt_conc):
    '''
    Calculate Tm of the probe of indicated length starting at each possible position.
    Use Biopython's MeltingTemp functions:
    Use probe conc (dnac1) = 250 nM (typical concentration for FISH experiments).
    Use template conc (dnac1) = 0 to assume that probe is in excess.
    This procedure goes from left to right along the target (RNA) sequence.
    Reverse complement to get the probe (DNA) sequence.
    '''
    probe_list = []
    target_len = len(seq)
    for i in range(0, target_len - length + 1):
        target_seq = seq[i: i + length]
        probe_seq = str(target_seq.reverse_complement())
        if 'N' in probe_seq:
            #Add nan for probes containing Ns
            probe_list.append([i, i + length, length, np.nan, probe_seq])
        else:
            #using the R_DNA NN table requires the RNA sequence
            Tm =  mt.Tm_NN(target_seq.transcribe(), nn_table = mt.R_DNA_NN1, Na = salt_conc, saltcorr = 4, dnac1 = 250, dnac2 = 0)
            probe_list.append([i, i + length, length, Tm, probe_seq])

    probe_df = pd.DataFrame(probe_list, columns = ['start', 'end', 'length', 'Tm', 'sequence'])
    #adjust to 1-based, inclusive indices for output
    probe_df['target_start'] = probe_df['start'] + 1
    probe_df['target_end'] = probe_df['end']
    return probe_df

def choose_nonoverlapping(probe_df, nt_spacing = 2):
    '''
    Choose largest set of non-overlapping probes with the greedy interval scheduling algorithm.
    Add nt_spacing to the finish time to add spacing between probes
    If nt_spacing is a list, then it will use the spacing specified at each index
    If nt_spacing is an integer, it will make a list with this integer the length of probe_df
    Initialize with the first acceptable probe, not subject to the gap.
    If probe_df is empty return an error
    '''
    error_message = 'Not enough space to choose probes. Consider relaxing your design constraints.'

    if len(probe_df) < 1:
        raise NotEnoughSpaceException(error_message)

    if type(nt_spacing) == int:
        nt_spacing = [nt_spacing]*len(probe_df)
        #nt spacing is irrelevant to the first chosen probe
        nt_spacing[0] = 0

    #sort probes by finish time
    sorted_df = probe_df.sort_values('end')
    sorted_finish_times = sorted_df['end'].values
    start_times = sorted_df['start'].values

    chosen_indices = []
    #set prev_finish_time to - nt spacing will allow selection of a probe starting at 0
    prev_finish_time = - 1
    for i in range(0, len(sorted_finish_times)):
        #if the start time of the next element is later than previous finish time, add
        if start_times[i] > prev_finish_time + nt_spacing[i]:
            chosen_indices.append(i)
            prev_finish_time = sorted_finish_times[i]

    chosen_probes_df = sorted_df.iloc[chosen_indices]
    return chosen_probes_df

def trim_to_size(probe_df, desired_number_probes, gap_size_1, gap_size_2, overshot_probe_num):
    '''
    We know that neither gap_size_1 nor gap_size_2 produced the exact right number of probes.
    We know that gap size 2 produced too many (gap size too small) and that gap size 1 produced too few (too big).
    Therefore, run interval overlap with a mix of gap sizes until desired probe size is reached.
    overshot_probe_num tells us how many probes are in the too large set, to initialize the gap size vector
    '''
    #arrange indices to go approximately by quarters through the target space
    num_groups =  int(math.ceil(desired_number_probes/4))
    ll = []
    for i in range(0, num_groups):
        ll.append([j for j in range(i, overshot_probe_num, num_groups)])

    indices = [item for sublist in ll for item in sublist]
    #nt_spacing is a list that gives spacing for each probe
    attempt_num = 0
    nt_spacing = [gap_size_2]*overshot_probe_num
    while attempt_num < len(nt_spacing):
        nt_spacing[indices[attempt_num]] = gap_size_1
        suggested_probe_df = choose_nonoverlapping(probe_df, nt_spacing = nt_spacing)
        if len(suggested_probe_df) == desired_number_probes:
            return suggested_probe_df
        attempt_num += 1

    #if you get all the way back to the gap_size_1 and still not the right number,
    #return the probes as is.
    return suggested_probe_df

def melting_temp(seq, Na_conc = 300):
    seq = Seq(seq)
    targeted_RNA = seq.reverse_complement()
    Tm =  mt.Tm_NN(targeted_RNA, nn_table = mt.R_DNA_NN1, Na = Na_conc, saltcorr = 4, dnac1 = 250, dnac2 = 0)
    return Tm

def get_subregion_ranges(ranges_string):
    '''
    Convert the provided ranges string into an array, 0-based, closed interval.
    The reason for making it closed is to match the alignment-derived intervals.
    e.g. "1-2000,2300-3400" -> [[0, 1999], [2299, 3399]]
    '''

    subregion_ranges = []
    subregion_list = [i for i in ranges_string.replace(', ',',').split(',')]
    for i in subregion_list:
        start, end = [int(j.strip()) for j in i.split('-')]
        #convert to 0-based, half-open
        start -= 1
        end -= 1
        subregion_ranges.append([start, end])
    subregion_ranges = np.array(subregion_ranges)
    return subregion_ranges

class ProbeSet(object):
    '''
    Set of probes covering one target.
    '''
    def __init__(self, Na_conc, outdir, target_name = ''):
        '''
        Initiate probe set with Na concentration.
        '''
        self.Na_conc = Na_conc
        self.outdir = outdir
        self.target_name = target_name

    def scan_sequence(self, infasta, min_probe_len, max_probe_len):
        '''
        Initiate probe set using min and max lengths of probes.
        '''
        self.target_seq = SeqIO.read(infasta, 'fasta').upper()
        self.target_len = len(self.target_seq)
        self.min_probe_len = min_probe_len
        self.max_probe_len = max_probe_len

        lens_to_test = range(min_probe_len, max_probe_len + 1)
        all_probes = []
        for probe_len in lens_to_test:
            these_probes = Tm_from_position(self.target_seq.seq, probe_len, self.Na_conc)
            all_probes.append(these_probes)
        self.probe_df = pd.concat(all_probes)
        #do not change the index in any of the subsequent steps -- used as probe identifier
        self.probe_df.reset_index(drop = True, inplace = True)
        self.probe_df.index.name = 'unique_id'
        self.probe_df['target_name'] = self.target_name

    def load_premade_probes(self, oligo_file):
        '''
        Load user-input probes in csv format.
        Useful for calculating physical properties of probes designed elsewhere.
        Infile must be a csv file with a column labeled 'sequence'.
        '''
        oligo_df = pd.read_csv(oligo_file)
        oligo_df['Tm'] = oligo_df['sequence'].apply(melting_temp, Na_conc = self.Na_conc)
        oligo_df['length'] = oligo_df['sequence'].apply(lambda x: len(x))
        self.passed_df = oligo_df

    def sequence_composition_filter(self, rules, filter = True):
        '''
        Remove sequences with homopolymeric repeats and other issues described in the rules.
        If filter == True, then it will replace the probe_df with filtered probes only.
        Set to False if you want to score the probes but not filter.
        '''
        self.passed_df['first_half'] = self.passed_df['length'].apply(lambda x: math.floor(x/2))
        self.passed_df['GC_content'] = self.passed_df['sequence'].apply(lambda x: (x.count('G') + x.count('C'))/len(x))
        self.passed_df['A_content'] = self.passed_df['sequence'].apply(lambda x: x.count('A')/len(x))
        self.passed_df['C_content'] = self.passed_df['sequence'].apply(lambda x: x.count('C')/len(x))
        #Set to True if the probe has this characteristic
        self.passed_df['GC_content_rule'] = self.passed_df.apply(lambda x: 0.4 > x['GC_content'] or x['GC_content'] > 0.6, axis = 1)
        self.passed_df['A_composition_rule'] = self.passed_df.apply(lambda x: x['A_content'] > 0.28, axis = 1)
        self.passed_df['C_composition_rule'] = self.passed_df.apply(lambda x: 0.22 > x['C_content'] or x['C_content'] > 0.28, axis = 1)
        self.passed_df['4xA_stack_rule'] = self.passed_df['sequence'].apply(lambda x: 'AAAA' in x)
        self.passed_df['4xC_stack_rule'] = self.passed_df.apply(lambda x: 'CCCC' in x['sequence'][0:x['first_half']], axis = 1)
        self.passed_df['earlyCs_rule'] = self.passed_df.apply(lambda x: np.any([(x['sequence'][i:i+6].count('C') >= 4) for i in range(0, x['first_half'] - 5)]), axis = 1)
        self.passed_df['any5_rule'] = self.passed_df['sequence'].apply(lambda x: np.any([N*5 in x for N in ['A','T','C','G']]))
        self.passed_df['passed_sequence'] = self.passed_df.apply(lambda row: not row[rules].any(), axis = 1)
        if filter == True:
            self.passed_df = self.passed_df[self.passed_df['passed_sequence']].copy()

    def structure_filter(self, hairpin_min, filter = True):
        '''
        Use primer3 to calculate energy of hairpin structure.
        https://libnano.github.io/primer3-py/quickstart.html#thermodynamic-analysis
        '''
        self.passed_df['hairpin_dG'] = self.passed_df['sequence'].apply(lambda x: primer3.calcHairpin(x, mv_conc = self.Na_conc).dg/1000)
        self.passed_df['passed_structure'] = self.passed_df['hairpin_dG'] >= hairpin_min

        if filter == True:
            self.passed_df = self.passed_df[self.passed_df['passed_structure']].copy()

    def masking_filter(self, masking_file):
        '''
        Take in file with masked starts by length and remove probes at those positions.
        '''
        mask_df = pd.read_csv(masking_file)
        mask_df.rename(columns = {i: int(i.split('_')[-1]) for i in mask_df.columns}, inplace = True)
        mask_series = mask_df.to_dict(orient = 'series')
        #remove the nans and convert indices to integers
        mask_dict = {k: mask_series[k].dropna().astype(int).values for k in mask_series}
        #Check if potential probe start is masked:
        self.passed_df['passed_masking'] = self.passed_df.apply(lambda x: (x['start'] not in mask_dict[x['length']] if x['length'] in mask_dict else True), axis = 1)
        self.passed_df = self.passed_df[self.passed_df['passed_masking']].copy()

    def excluded_nt_filter(self, excluded = None, excluded_consensus = None):
        '''
        Remove probes in user-provided excluded regions, for example regions of high
        tertiary structure. Can be from excluded (calculated in pipeline relative to
        one transcript) or in excluded_conserved (relative to the pre-calculated
        consensus sequence)

        '''
        ll = []
        for ranges in [excluded, excluded_consensus]:
            if range_defined(ranges):
                #range is already 0-based, but closed interval.
                region_dict = defaultdict(set)
                for i in ranges:
                    logging.info('Excluded nts %s-%s in consensus sequence.' % (i[0]+1, i[1]+1))
                    for j in range(self.min_probe_len, self.max_probe_len + 1):
                        #+1 because we want reads that include that region
                        #mask any of the start positions that would overlap the excluded region
                        masked_start = i[0] - j + 1
                        if masked_start < 0:
                            masked_start = 0
                        #endpt is excluded. +1 to cover it.
                        masked_end = i[1] + 1
                        new_range = set(range(masked_start, masked_end))
                        region_dict[j].update(new_range)

            else:
                region_dict = {}
            ll.append(region_dict)

        #combine the excluded dicts
        bad_start_dict = {k: set() for k in set().union(*ll)}
        for l in bad_start_dict:
            bad_start_dict[l] = set().union(*[i[l] for i in ll if l in i])

        #Check if potential probe start is masked:
        ##self.passed_df['passed_excluded'] = self.passed_df.apply(lambda x: (x['start']) not in bad_start_dict[x['length']], axis = 1)
        self.passed_df['passed_excluded'] = self.passed_df.apply(lambda x: (x['start'] not in bad_start_dict[x['length']] if x['length'] in bad_start_dict else True), axis = 1)
        self.passed_df = self.passed_df[self.passed_df['passed_excluded']].copy()

    def quantile_filter(self, window_size, Tm_quantile):
        '''
        Calculate a rolling Tm quantile and determine if the probe passes the
        quantile threshold. The quantile is calculated using the physical window size,
        e.g. 200 nt. However, only the remaining probes passing filters will have
        their Tm used for calculating the quantiles.
        '''
        #Calculate rolling Tm quantile and find high Tm probes.
        quantile_df = self.probe_df.sort_values('start', ascending = True).copy()
        failed_idx = quantile_df.index.difference(self.passed_df.index)
        #set the failed indices Tm to np.nan so won't be used for the rolling quantile
        quantile_df.loc[failed_idx, 'Tm'] = np.nan
        #window size is adjusted here to be scaled to the number of probes per position
        window_size = window_size * (self.max_probe_len - self.min_probe_len + 1)
        #Calculate the Tm of the provided quantile at each position
        #Only the probes that passed the filters are used for the Tm quantile calculations
        quantile_df['rolling_quantile_co'] = quantile_df['Tm'].rolling(
        window_size, center = True, min_periods = 1).quantile(Tm_quantile)

        quantile_df['passed_Tm_quantile'] = quantile_df['Tm'] > quantile_df['rolling_quantile_co']
        self.passed_df['passed_Tm_quantile'] = quantile_df['passed_Tm_quantile']
        self.passed_df = self.passed_df[self.passed_df['passed_Tm_quantile']].copy()

    def Tm_window_filter(self, min_tm, max_tm):
        '''
        Remove probes with Tm that fall outside the user-specified Tm window.
        '''
        sorted_probes_df = self.passed_df.sort_values('Tm')
        passed_indices = np.searchsorted(sorted_probes_df['Tm'].values, [min_tm, max_tm])
        passed_probes_df = sorted_probes_df.iloc[passed_indices[0]: passed_indices[1]].copy()
        passed_probes_df['passed_Tm_window'] = True
        self.passed_df = passed_probes_df.copy()

    def prune(self, desired_number_probes, quantile_filter = True, subregions = None):
        '''
        Prune the probes into the desired number of probes per target.
        - find Tm peaks
        - get N evenly spaced probes
        - if quantile_filter = True, require the peaks to also pass the Tm quantile threshold
        '''

        error_message = '''\
        Not enough space to design requested number of probes.
        Consider designing fewer probes or relaxing your design constraints.'''

        #choose the highest Tm probe at each start site:
        idx = self.passed_df.groupby(['start'])['Tm'].transform(max) == self.passed_df['Tm']
        tm_df = self.passed_df[idx].copy()
        #tm_df['passed_quantile'] = tm_df.index.isin(self.hitm_df[self.hitm_df['high_Tm']].index)
        tm_df['unique_id'] = tm_df.index

        if not range_defined(subregions):
            subregions = np.array([[1, self.target_len - 1]])

        #split the desired number of probes betwen the subregions
        this_probeset_size =  int(math.ceil(desired_number_probes/len(subregions)))
        chosen_probes = []
        sorted_subregions = subregions[np.argsort(subregions[:, 0])]

        #add 1 to the endpts to make half-open
        sorted_subregions[:,1] += 1
        substring = ', '.join(['%s-%s' % (i[0], i[1]) for i in sorted_subregions])
        logging.info('Choosing probes in subregions %s' % substring)

        for i, subregion in enumerate(subregions):

            #get the mini df that contains the data in the subregion
            sub_df = tm_df[(subregion[0] <= tm_df['target_start']) & (tm_df['target_end'] <= subregion[1])]

            #Add the missing start positons back to but with Tm of 0
            #This way, they will be included in the distance consideration for peak finding
            #but can't be chosen as peaks themselves
            this_distance = 50

            #start earlier because it cannot choose the endpts
            start_range = range(sub_df['start'].min() - 1, sub_df['start'].max()+ 2)
            range_df = pd.DataFrame(start_range, columns = ['start'])
            new_df = pd.merge(range_df[['start']], sub_df[['unique_id', 'Tm', 'start']], 'outer', on = 'start')

            ##new_df = pd.merge(range_df[['start']], sub_df[['unique_id', 'Tm', 'start', 'passed_quantile']], 'outer', on = 'start')
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
            #edge case where the number of peaks is two and the number of probes is 1, choose one closer to the middle
            if len(peak_locs) == 1 and this_probeset_size == 1:
                chosen_locs = peak_locs
            elif len(peak_locs) == 2 and this_probeset_size == 1:
                mid_dist = abs(peak_locs/(subregion[-1] - subregion[0] + 1) - 0.5)
                chosen_locs = [peak_locs[np.argmin(mid_dist)]]
            else:
                chosen_locs = probe_site_selection.choose_combination(peak_locs, this_probeset_size)
            chosen_ids = new_df.loc[new_df['start'].isin(chosen_locs), 'unique_id'].astype('int')
            sub_df = self.passed_df[self.passed_df.index.isin(chosen_ids)].copy()

            chosen_probes.append(sub_df)

        #combine probes from each subregion into the pruned_df
        self.pruned_df = pd.concat(chosen_probes)

    def plot_selected_probes(self):
        '''
        Plot the Tm by position of the target and the locations of the selected
        probes. First plot the Tms from the predfs -- including all Tms for probes without Ns.
        Next, plot the Tms from the probes passing all the filters, including all
        Tms for the probes that passed the masking, sequence, and structure filters.
        Finally, plot the positions of the selected probes.
        '''
        fig = plt.figure(figsize = (5, 2.5))
        ax = fig.add_subplot(111)

        grey = '#A5AA99'
        purple = '#4b4b8f'
        pink = '#CF1C90'

        pre_tm_cols = ['passed_masking', 'passed_sequence', 'passed_structure']

        self.pre_tm_filter_df['midpt'] = self.pre_tm_filter_df['target_start'] + (self.pre_tm_filter_df['length'] - 1)/2
        self.pre_tm_filter_df.sort_values(by = 'midpt', ascending = True, inplace = True)

        bg = ax.scatter(self.pre_tm_filter_df['midpt'], self.pre_tm_filter_df['Tm'], s = 30,
        alpha = 0.3, color = purple, edgecolors = 'none')

        mini_df = self.pre_tm_filter_df[self.pre_tm_filter_df.index.isin(self.final_df.index)].copy()
        selected = ax.scatter(mini_df['midpt'], mini_df['Tm'], s = 30, alpha = 0.3, color = pink, edgecolors = 'none')

        #minidf = self.probe_df.loc[self.probe_df.index.isin(self.hitm_df.index)].copy()

        #minidf = self.probe_df.loc[self.probe_df.index.isin(self.hitm_df.index)].copy()
        #pre = ax.scatter(minidf['midpt'], minidf['Tm'], s = 30, alpha = 0.3,
        #color = pink, edgecolors = 'none')

        for p in self.final_df.itertuples():
            probe = ax.axvspan(p.target_start - 1, p.target_end, alpha=0.5, ymin = 0.0,
            ymax = 1, color = grey, linewidth = 0)

        ax.set_xlim(0, self.target_len)
        #ax.set_xlim(self.probe_df['target_start'].min(), self.probe_df['target_end'].max())

        ax.set_ylabel('Tm')
        ax.set_xlabel('target position (nt)')

        ax.legend([bg, selected], ['before selection', 'selected'],
               mode = 'expand', fontsize = 8, ncol = 3, bbox_to_anchor=(0., 1.05, 1., .105), loc=3,
               borderaxespad=0., handletextpad=0.1)

        plt.tight_layout()

        plt.savefig('%s.png' % os.path.join(self.outdir, 'selected_probes_%s' % self.target_name), dpi = 600)

    def summarize_results(self, col_order):
        '''
        Write the final selected probes, plot the selected regions.
        '''
        self.plot_selected_probes()
        self.final_df.reset_index()[col_order].to_csv(os.path.join(self.outdir, 'selected_probes_%s.csv' % self.target_name), index = False)

def main(arglist):
    '''
    Run the probe design pipeline or calculate the properties of the provided oligos.
    '''
    #The sequence composition rules
    default_rules = ['GC_content_rule', '4xA_stack_rule', '4xC_stack_rule', 'any5_rule']
    possible_rules = ['GC_content_rule', 'A_composition_rule', 'C_composition_rule', '4xA_stack_rule', '4xC_stack_rule', 'earlyCs_rule', 'any5_rule']

    parser = argparse.ArgumentParser()
    #This set of args is generally provided by snakemake
    parser.add_argument('-target_names', nargs = '+', help = 'names of the targets, otherwise will infer from the names of the fasta files')
    parser.add_argument('-target_fastas', nargs = '+', help = 'fasta file of target sequence with any excluded nts (i.e. non-conserved in set of transcript targets) marked by Ns')
    parser.add_argument('-min_probe_length', nargs = '+', type = int, default = [26])
    parser.add_argument('-max_probe_length', nargs = '+', type = int, default = [35])
    parser.add_argument('-min_Tm', nargs = '+', type = int, default = [0])
    parser.add_argument('-max_Tm', nargs = '+', type = int, default = [10000])
    parser.add_argument('-outdir', help = 'name of output directory')
    parser.add_argument('-desired_number_probes', nargs = '+', type = int, default = [10], help = 'input an integer number and the probe set will be pruned to this number of probes.' )
    parser.add_argument('-Tm_quantile', nargs = '+', type = float, default = [0.9], help = 'Tm of probe must be above this quantile to be selected')
    parser.add_argument('-Tm_window_size', nargs = '+', type = int, default = [200], help ='# nt to include in calculating the Tm quantile')
    parser.add_argument('-masked_nts', nargs = '+', default = [None], help = 'csv file containing masked starting nts by length')
    parser.add_argument('-excluded_regions_consensus', nargs = '+', default = [None], help = 'Regions to avoid placing probes in. Specify as a list of 1-based closed interval regions.')
    parser.add_argument('-excluded_regions', nargs = '+', default = [None], help = 'Regions to avoid placing probes in, calculated from one transcript in each target. Not provided directly.')
    parser.add_argument('-target_subregions_consensus', nargs = '+', default = [None], help = 'Regions to split the number of probes bewteen. Specify as a list of 1-based closed interval regions.')
    parser.add_argument('-min_hairpin_dG', nargs = '+', default = [-3], help = 'deltaG values for probes must be > than this value')
    parser.add_argument('-min_dimer_dG', nargs = '+', default = [-10], help = 'deltaG values for probes must be > than this value')

    #These args are generally taken from the defaults but they could be overriden by snakemake
    #These could be in a list/series (specified at the per target level)
    parser.add_argument('-sequence_filter_rules', nargs = '+', default = [default_rules], help = 'remove probes not passing these, choose from: %s' % ', '.join(possible_rules))
    #These will not be in a list
    parser.add_argument('-Na_conc', default = 300, help = 'Na+ concentration of hybridization in mM')
    parser.add_argument('-premade_probes', help = 'csv file containing sequences of probes already designed elsewhere')
    parser.add_argument('--analyze_oligos', action = 'store_true', help = 'input a csv file containing sequences to calculate Tm and dG values, dont run the rest')
    parser.add_argument('--design_probes', action = 'store_true', help = 'run probe design')
    parser.add_argument('--quick_test', action = 'store_true', help = 'run with only a subset of probes to check output')

    args, unknown = parser.parse_known_args()
    argdict = vars(args)

    #columns to output after analysis
    col_order = ['sequence', 'target_name',  'target_start',  'target_end', 'length',\
    'unique_id', 'Tm', 'GC_content', 'A_content', 'C_content', 'dimer_dG',\
    'dimer_partner', 'GC_content_rule', 'A_composition_rule', 'C_composition_rule',\
    '4xA_stack_rule', '4xC_stack_rule', 'earlyCs_rule', 'any5_rule', 'passed_masking',\
    'passed_sequence', 'passed_structure', 'passed_Tm_window', 'passed_Tm_quantile']

    #A hack to get the argparse defaults overriden by snakemake args, if provided
    if 'snakemake' in globals():
        #If the script is called by snakemake, reassign args from input, output, and params, where specified
        toset = defaultdict(list)
        for arg in argdict:
            if hasattr(snakemake.input, arg):
                toset['input'].append(arg)
            if hasattr(snakemake.params, arg):
                toset['params'].append(arg)
            if hasattr(snakemake.output, arg):
                toset['output'].append(arg)

        for block in toset:
            for arg in toset[block]:
                setattr(args, arg, getattr(getattr(snakemake, block), arg))

    #Duplicate any default arguments to match the number of targets.
    target_param_l = ['min_probe_length', 'max_probe_length', 'min_Tm', 'max_Tm', 'desired_number_probes',
    'Tm_quantile', 'Tm_window_size', 'masked_nts', 'excluded_regions', 'excluded_regions_consensus',
    'target_subregions_consensus', 'sequence_filter_rules',
    'min_hairpin_dG', 'min_dimer_dG']

    if args.analyze_oligos:
        probeset = ProbeSet(args.Na_conc, args.outdir, args.target_names)
        probeset.load_premade_probes(oligo_df)
        probeset.sequence_composition_filter(args.sequence_filter_rules, filter = False)
        probeset.structure_filter(args.min_hairpin_dG, filter = False)
        probeset.passed_df = remove_bad_probes(probeset.passed_df, args.min_dimer_dG, filter = False)
        probeset.passed_df.sort_values('target_start', inplace = True)
        probeset.passed_df.to_csv(os.path.join(args.outdir, 'oligo_properties.csv'), index = False)

    elif args.design_probes:
        num_targets = len(args.target_fastas)

        for arg in target_param_l:
            if len(argdict[arg]) < num_targets:
                setattr(args, arg, argdict[arg]*num_targets)

        os.makedirs(args.outdir, exist_ok = True)

        if args.target_names is None:
            args.target_names = [os.path.basename(i).split('.')[0] for i in args.target_fasta]

        all_selected_probes = pd.DataFrame()
        starting_index = 0

        logging.basicConfig(level=logging.DEBUG, filename = os.path.join(args.outdir,
        'logfile.txt'), filemode = 'w', format = '%(message)s')
        #stop writing all the font warnings to the log file
        logging.getLogger('matplotlib.font_manager').disabled = True
        for i, target in enumerate(args.target_names):
            assert args.min_Tm[i] < args.max_Tm[i]
            assert args.min_probe_length[i] <= args.max_probe_length[i]
            target_outdir = os.path.join(args.outdir, target)
            os.makedirs(target_outdir, exist_ok = True)

            logging.info('Target %s: ' % target)

            exdf = pd.read_csv(args.excluded_regions[i])
            #choose the subregion ranges to be used. If there are ranges provided wrt consensus, replace the the calculated ones
            target_subregions = exdf.loc[exdf['region_type'] == 'target_subregions', ['start', 'end']].values
            if not pd.isnull(args.target_subregions_consensus[i]):
                target_subregions = get_subregion_ranges(args.target_subregions_consensus[i])
            if not pd.isnull(args.excluded_regions_consensus[i]):
                args.excluded_regions_consensus[i] = get_subregion_ranges(args.excluded_regions_consensus[i])


            probeset = ProbeSet(args.Na_conc, target_outdir, target)
            probeset.scan_sequence(args.target_fastas[i], args.min_probe_length[i], args.max_probe_length[i])
            probeset.probe_df.index += starting_index
            probeset.passed_df = probeset.probe_df.copy()

            #update the starting index
            starting_index += len(probeset.probe_df)

            if args.quick_test:
                #Just for debugging, remove all but a subset of probes
                l = np.arange(0, len(probeset.probe_df), 10)
                probeset.probe_df = probeset.probe_df.iloc[l]
                probeset.passed_df = probeset.probe_df.copy()

            logging.info("Starting with %s potential probes." % len(probeset.passed_df))
            print('filtering out masked nts')
            if args.masked_nts:
                probeset.masking_filter(args.masked_nts[i])
                logging.info("%s potential probes remaning after masked nt filter." % len(probeset.passed_df))

            print('filtering out excluded regions')
            ex_calc = exdf.loc[exdf['region_type'] == 'excluded_regions', ['start', 'end']].values
            probeset.excluded_nt_filter(excluded = ex_calc, excluded_consensus = args.excluded_regions_consensus[i])
            logging.info("%s potential probes remaning after masked nt filter." % len(probeset.passed_df))

            print('filtering by sequence')
            probeset.sequence_composition_filter(args.sequence_filter_rules[i])
            logging.info("%s potential probes remaning after sequence composition filter." % len(probeset.passed_df))

            print('filtering by structure')
            probeset.structure_filter(args.min_hairpin_dG[i])
            logging.info("%s potential probes remaning after structure filter." % len(probeset.passed_df))
            probeset.pre_tm_filter_df = probeset.passed_df.copy()
            probeset.passed_df.to_csv(os.path.join(target_outdir, 'before_tm_filter_%s.csv' % target))

            print('getting probes in selected Tm quantile')
            probeset.quantile_filter(args.Tm_window_size[i], args.Tm_quantile[i])
            logging.info("%s potential probes passed quantile threshold." % len(probeset.passed_df))

            print('getting probes in selected Tm range')
            probeset.Tm_window_filter(args.min_Tm[i], args.max_Tm[i])
            logging.info("%s potential probes in selected Tm range." % len(probeset.passed_df))

            probeset.passed_df.to_csv(os.path.join(target_outdir, 'before_pruning_%s.csv' % target))

            #Prune the remaining probes to the desired number of probes.
            removed_idx = [1]
            p = 0
            while len(removed_idx) > 0:
                #get evenly spaced probes from the passed ones
                probeset.prune(args.desired_number_probes[i], subregions = target_subregions)
                #drop rows causing the low dimer dG scores
                print('removing probes with dimer clashes')
                screened_df = remove_bad_probes(probeset.pruned_df, args.min_dimer_dG[i], all_selected_probes)
                #get indices of the pruned probes that don't survive heterodimer screening
                removed_idx = probeset.pruned_df.index.difference(screened_df.index)
                #drop bad heterodimer probes, will trigger rerun of prune with this probe removed.
                probeset.passed_df.drop(removed_idx, inplace = True)
                p += 1

            #Get the properties of the selected probes and write to output file
            final_probe_df = pd.merge(probeset.passed_df.reindex(probeset.pruned_df.index),
            screened_df[['dimer_dG','dimer_partner']], left_index = True, right_index = True)
            #probeset.final_df = final_probe_df.reset_index(drop = True)[col_order]
            probeset.final_df = final_probe_df
            logging.info("%s probes selected." % len(probeset.final_df))
            all_selected_probes = all_selected_probes.append(probeset.final_df)
            probeset.summarize_results(col_order)

        #write the combined output file with probes selected for all targets.
        all_selected_probes.sort_values(by = ['target_name', 'target_start'], inplace = True)
        all_selected_probes.reset_index(inplace = True)
        all_selected_probes['probe_num'] = all_selected_probes.index
        all_selected_probes['probe_num'] += 1
        cols = ['probe_num']
        cols.extend(col_order)
        all_selected_probes[cols].to_csv(os.path.join(args.outdir, 'all_selected_probes.csv'), index = False)
    return args
if __name__ == '__main__':
    main(sys.argv[1:])
