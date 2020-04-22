'''
screen_kmers.py
Generate kmers of the target sequence and screen out ones that don't pass
the sequence composition rules, that overlap Ns, and that don't fall within
specified Tm range.
'''
from Bio import SeqIO
import sys
import os
import pandas as pd
import numpy as np
from Bio.SeqUtils import MeltingTemp as mt
import math
import primer3
import probe_helpers
import argparse
from collections import defaultdict
import logging
os.environ['NUMEXPR_NUM_THREADS'] = '8'

def Tm_window_filter(df, min_tm, max_tm):
    '''Select probes with Tms between min and max'''

    sorted_df = df.sort_values('Tm')
    passed_indices = np.searchsorted(sorted_df['Tm'].values, [min_tm, max_tm])
    passed_df = sorted_df.iloc[passed_indices[0]: passed_indices[1]].copy()
    return passed_df

def sequence_composition_filter(df, min_gc, max_gc, rules, filter = True):
    '''
    Remove sequences with homopolymeric repeats and other issues described in the rules.
    If filter == True, then it will replace the probe_df with filtered probes only.
    Set to False if you want to score the probes but not filter.
    '''
    df['first_half'] = df['length'].apply(lambda x: math.floor(x/2))
    df['GC_content'] = df['sequence'].apply(lambda x: (x.count('G') + x.count('C'))/len(x))
    df['A_content'] = df['sequence'].apply(lambda x: x.count('A')/len(x))
    df['C_content'] = df['sequence'].apply(lambda x: x.count('C')/len(x))
    #Set to True if the probe has this characteristic
    df['GC_content_rule'] = df.apply(lambda x: min_gc > x['GC_content'] or x['GC_content'] > max_gc, axis = 1)
    df['A_composition_rule'] = df.apply(lambda x: x['A_content'] > 0.28, axis = 1)
    df['C_composition_rule'] = df.apply(lambda x: 0.22 > x['C_content'] or x['C_content'] > 0.28, axis = 1)
    df['4xA_stack_rule'] = df['sequence'].apply(lambda x: 'AAAA' in x)
    df['4xC_stack_rule'] = df.apply(lambda x: 'CCCC' in x['sequence'][0:x['first_half']], axis = 1)
    df['earlyCs_rule'] = df.apply(lambda x: np.any([(x['sequence'][i:i+6].count('C') >= 4) for i in range(0, x['first_half'] - 5)]), axis = 1)
    df['any5_rule'] = df['sequence'].apply(lambda x: np.any([N*5 in x for N in ['A','T','C','G']]))
    df['passed_sequence'] = df.apply(lambda row: not row[rules].any(), axis = 1)
    if filter == True:
        df = df[df['passed_sequence']].copy()
    return df

def scan_sequence(target_name, seq, min_len, max_len, Na_conc):
    '''
    Return probes from sequence for each length.
    '''
    lens_to_test = range(min_len, max_len + 1)
    all_probes = []
    for probe_len in lens_to_test:
        these_probes = Tm_from_position(seq, probe_len, Na_conc)
        all_probes.append(these_probes)

    df = pd.concat(all_probes)
    df['target_name'] = target_name
    #do not change the index in any of the subsequent steps -- used as probe identifier
    df.reset_index(drop = True, inplace = True)
    df['unique_id'] = df.index.map(lambda x: '%s_%s' % (target_name, x))
    df.set_index('unique_id', inplace = True)
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
            probe_list.append([i, i + length, length, np.nan, str(target_seq), probe_seq])
        else:
            #using the R_DNA NN table requires the RNA sequence
            Tm =  mt.Tm_NN(target_seq.transcribe(), nn_table = mt.R_DNA_NN1, Na = salt_conc, saltcorr = 4, dnac1 = 250, dnac2 = 0)
            probe_list.append([i, i + length, length, Tm, str(target_seq), probe_seq])

    probe_df = pd.DataFrame(probe_list, columns = ['start', 'end', 'length', 'Tm', 'target_sequence', 'sequence'])
    #adjust to 1-based, inclusive indices for output
    probe_df['target_start'] = probe_df['start'] + 1
    probe_df['target_end'] = probe_df['end']
    return probe_df

def quantile_filter(df, original_df, window_size, Tm_quantile, min_probe_len, max_probe_len):
    '''
    Calculate a rolling Tm quantile and determine if the probe passes the
    quantile threshold. Set the Tms of the failed probes to NaN so that
    the quantile will be calculated relative to original position in the sequence
    but the failed probes will not be included in the calculation.
    '''
    quantile_df = original_df.sort_values('start', ascending = True)
    failed_idx = quantile_df.index.difference(df.index)
    quantile_df.loc[failed_idx, 'Tm'] = np.nan
    #window size is adjusted here to be scaled to the number of probes per position
    window_size = window_size * (max_probe_len - min_probe_len + 1)
    #Calculate the Tm of the provided quantile at each position
    quantile_df['rolling_Tm_quantile_co'] = quantile_df['Tm'].rolling(
    window_size, center = True, min_periods = 1).quantile(Tm_quantile)
    quantile_df['passed_Tm_quantile'] = quantile_df['Tm'] > quantile_df['rolling_Tm_quantile_co']
    df['rolling_Tm_quantile_co'] = quantile_df['rolling_Tm_quantile_co']
    df['passed_Tm_quantile'] = quantile_df['passed_Tm_quantile']
    return df[df['passed_Tm_quantile']]

def structure_filter(df, hairpin_min, dimer_min, Na_conc, filter = True):
    '''
    Use primer3 to calculate energy of hairpin structure.
    https://libnano.github.io/primer3-py/quickstart.html#thermodynamic-analysis
    '''
    df['hairpin_dG'] = df['sequence'].apply(lambda x: primer3.calcHairpin(x, mv_conc = Na_conc).dg/1000)
    df['homodimer_dG'] = df['sequence'].apply(lambda x: primer3.calcHomodimer(x, mv_conc = Na_conc).dg/1000)

    df['passed_structure'] = (df['hairpin_dG'] >= hairpin_min) & (df['homodimer_dG'] >= dimer_min)

    if filter == True:
        df = df[df['passed_structure']].copy()
    return df

def excluded_nt_filter(df, min_probe_len, max_probe_len, excluded = None, excluded_consensus = None):
    '''
    Remove probes in user-provided excluded regions, for example regions of high
    tertiary structure. Can be from excluded (calculated in pipeline relative to
    one transcript) or in excluded_conserved (relative to the pre-calculated
    consensus sequence)
    '''
    ll = []
    for ranges in [excluded, excluded_consensus]:
        if probe_helpers.range_defined(ranges):
            #range is already 0-based, but closed interval.
            region_dict = defaultdict(set)
            for i in ranges:
                logging.info('Excluded nts %s-%s in consensus sequence.' % (i[0]+1, i[1]+1))

                for j in range(min_probe_len, max_probe_len + 1):
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
    df['passed_excluded'] = df.apply(lambda x: (x['start'] not in bad_start_dict[x['length']] if x['length'] in bad_start_dict else True), axis = 1)
    return df[df['passed_excluded']]

def main(arglist):
    #The sequence composition rules
    default_rules = ['GC_content_rule', '4xA_stack_rule', '4xC_stack_rule', 'any5_rule']
    possible_rules = ['GC_content_rule', 'A_composition_rule', 'C_composition_rule', '4xA_stack_rule', '4xC_stack_rule', 'earlyCs_rule', 'any5_rule']

    parser = argparse.ArgumentParser()
    #This set of args is generally provided by snakemake
    parser.add_argument('-target_fasta', help = 'fasta file of target sequence with any excluded nts (i.e. non-conserved in set of transcript targets) marked by Ns')
    parser.add_argument('-min_probe_len', type = int, default = 26)
    parser.add_argument('-max_probe_len', type = int, default = 35)
    parser.add_argument('-min_Tm', type = int, default = 0)
    parser.add_argument('-max_Tm', type = int, default = 10000)
    parser.add_argument('-min_gc', type = float, default = 0.4)
    parser.add_argument('-max_gc', type = float, default = 0.6)
    parser.add_argument('-outdir', help = 'name of output directory')
    parser.add_argument('-Tm_quantile', type = float, default = 0.9, help = 'Tm of probe must be above this quantile to be selected')
    parser.add_argument('-Tm_window_size', type = int, default = 200, help ='# nt to include in calculating the Tm quantile')
    parser.add_argument('-excluded_regions_consensus', help = 'Regions to avoid placing probes in. Specify as a list of 1-based closed interval regions.')
    parser.add_argument('-excluded_regions', help = 'Regions to avoid placing probes in, calculated from one transcript in each target. Not provided directly.')
    parser.add_argument('-min_hairpin_dG', default = -3, help = 'deltaG values for probes must be > than this value')
    parser.add_argument('-min_dimer_dG', default = -10, help = 'deltaG values for probes must be > than this value')
    parser.add_argument('-probe_csv', default = 'potential_probes.csv')
    parser.add_argument('-probe_fa', default = 'potential_probes.fa')
    parser.add_argument('-logfile', default = 'logfile.txt')
    #These args are generally taken from the defaults but they could be overriden by snakemake
    #These could be in a list/series (specified at the per target level)
    parser.add_argument('-sequence_filter_rules', nargs = '+', default = default_rules, help = 'remove probes not passing these, choose from: %s' % ', '.join(possible_rules))
    #These will not be in a list
    parser.add_argument('-Na_conc', default = 300, help = 'Na+ concentration of hybridization in mM')
    parser.add_argument('-premade_probes', help = 'csv file containing sequences of probes already designed elsewhere')
    parser.add_argument('--analyze_oligos', action = 'store_true', help = 'input a csv file containing sequences to calculate Tm and dG values, dont run the rest')
    parser.add_argument('--design_probes', action = 'store_true', help = 'run probe design')
    parser.add_argument('--quick_test', action = 'store_true', help = 'run with only a subset of probes to check output')

    args, unknown = parser.parse_known_args()
    if 'snakemake' in globals():
        args = probe_helpers.set_snake_args(args, snakemake)

    target_name = os.path.basename(args.target_fasta).rstrip('.fa')

    logging.basicConfig(level=logging.DEBUG, filename = args.logfile, filemode = 'w', format = '%(message)s')
    logging.info('Target %s: ' % target_name)

    exdf = pd.read_csv(args.excluded_regions)

    if not pd.isnull(args.excluded_regions_consensus):
        args.excluded_regions_consensus = probe_helpers.get_subregion_ranges(args.excluded_regions_consensus)

    target_seq = next(SeqIO.parse(args.target_fasta, 'fasta')).seq.upper()

    kmer_df = scan_sequence(target_name, target_seq, args.min_probe_len, args.max_probe_len, args.Na_conc)

    #keep kmer_df as the original one and df as the one to be filtered.
    df = kmer_df.copy()
    df.dropna(subset = ['Tm'], inplace = True)

    logging.info('Starting with %s potential probes.' % len(df))

    print('removing probes in excluded regions')
    ex_calc = exdf.loc[exdf['region_type'] == 'excluded_regions', ['start', 'end']].values
    df = excluded_nt_filter(df, args.min_probe_len, args.max_probe_len, excluded = ex_calc, excluded_consensus = args.excluded_regions_consensus)
    logging.info('%s potential probes remaning after excluding excluded regions.' % len(df))

    print('filtering by structure filter')
    df = structure_filter(df, args.min_hairpin_dG, args.min_dimer_dG, args.Na_conc, filter = True)
    logging.info('%s potential probes remaning after structure filter.' % len(df))

    print('filtering by sequence composition')
    df = sequence_composition_filter(df, args.min_gc, args.max_gc, args.sequence_filter_rules, filter = True)
    logging.info('%s potential probes remaning after sequence composition filter.' % len(df))

    #only perform quantile filtering on the remaining probes
    print('filtering by Tm quantile')
    df = quantile_filter(df, kmer_df, args.Tm_window_size, args.Tm_quantile,
    args.min_probe_len, args.max_probe_len)
    logging.info('%s potential probes remaning after Tm quantile filter.' % len(df))

    #further filter out ones not in specified Tm range.
    print('filtering by Tm window')
    df = Tm_window_filter(df, args.min_Tm, args.max_Tm)
    logging.info('%s potential probes remaning after Tm window filter.' % len(df))

    #write a csv file for later to take the probe properties from.
    #write a kmer fasta file to use for the blast check.
    df.to_csv(args.probe_csv)
    with open(args.probe_fa, 'w') as f:
        for i in df.itertuples():
            f.write('>%s\n%s\n' % (i.Index, i.target_sequence))

if __name__ == '__main__':
    main(sys.argv[1:])
