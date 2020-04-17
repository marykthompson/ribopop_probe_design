'''
Given a multi fasta file, use mafft to create the
multiple sequence alignment for the given organism.
Then parse the alignments to create a consensus sequence.
'''

import pandas as pd
import os
import alignment_funcs

def adjust_indices(id, regions, alignment, name, region_type):
    '''
    Convert regions [[start1, end1], [start2, end2], [], ...]
    to a dataframe and annotate with the sequence name and the region type
    '''
    adj_regions = []
    for i in regions:
        new_start = alignment_funcs.column_from_residue_number(alignment, id, i[0])
        new_end = alignment_funcs.column_from_residue_number(alignment, id, i[1])
        adj_regions.append([new_start, new_end])

    df = pd.DataFrame(adj_regions, columns = ['start', 'end'])
    df['ID'] = name
    df['region_type'] = region_type
    return df

def main(arglist):
    fasta = snakemake.input['fasta']
    outfile = snakemake.output['outfasta']
    excluded1 = snakemake.output['excluded1']
    target_df = snakemake.params['target_df']
    name = snakemake.params['name']

    alignment = alignment_funcs.write_alignment(fasta, name, outfile)
    adj_dfs = []
    for rec in alignment:
        excluded_region_string = target_df.loc[target_df['ID'] == rec.id, 'excluded_regions'].values[0]
        if not pd.isnull(excluded_region_string):
            excluded_regions = alignment_funcs.get_subregion_ranges(excluded_region_string)
            ex_df = adjust_indices(rec.id, excluded_regions, alignment, name, 'excluded_regions')
            adj_dfs.append(ex_df)

        target_region_string = target_df.loc[target_df['ID'] == rec.id, 'target_subregions'].values[0]
        if not pd.isnull(target_region_string):
            target_regions = alignment_funcs.get_subregion_ranges(target_region_string)
            t_df = adjust_indices(rec.id, target_regions, alignment, name, 'target_subregions')
            adj_dfs.append(t_df)

    if adj_dfs == []:
        adj_df = pd.DataFrame(columns = ['ID', 'start', 'end', 'region_type'])
    else:
        adj_df = pd.concat(adj_dfs)

    adj_df.to_csv(excluded1, index = False)

if __name__ == '__main__':
    main(sys.argv[1:])
