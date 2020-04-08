'''
Given a multi fasta file, use mafft to create the
multiple sequence alignment for the given organism.
Then parse the alignments to create a consensus sequence.
'''

import pandas as pd
import os
import alignment_funcs

def main(arglist):
    fasta = snakemake.input['fasta']
    outfile = snakemake.output['outfasta']
    excluded1 = snakemake.output['excluded1']
    target_df = snakemake.params['target_df']
    name = snakemake.params['name']

    alignment = alignment_funcs.write_alignment(fasta, name, outfile)
    adj_regions = []
    for rec in alignment:
        excluded_region_string = target_df.loc[target_df['ID'] == rec.id, 'excluded_regions'].values[0]
        if not pd.isnull(excluded_region_string):
            excluded_regions = alignment_funcs.get_subregion_ranges(excluded_region_string)
            for i in excluded_regions:
                new_start = alignment_funcs.column_from_residue_number(alignment, rec.id, i[0])
                new_end = alignment_funcs.column_from_residue_number(alignment, rec.id, i[1])
                adj_regions.append([new_start, new_end])

    new_df = pd.DataFrame(adj_regions, columns = ['start', 'end'])
    new_df['ID'] = name
    new_df.to_csv(excluded1, index = False)

if __name__ == '__main__':
    main(sys.argv[1:])
