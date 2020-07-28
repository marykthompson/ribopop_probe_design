'''
Given multiple fasta files (corresponding to different organisms),
use mafft to create the multiple sequence alignment for the given target.
Then parse the alignments to create a consensus sequence.
'''

import pandas as pd
import os
import alignment_funcs
from Bio import SeqIO

def convert_indices(x, alignment =  None, col = None):
    '''
    Call column_from_residue_number to add the new index to the df
    '''
    new_index = alignment_funcs.column_from_residue_number(alignment, x['ID'], x[col])
    return new_index

def main(arglist):
    fastas = snakemake.input['fastas']
    outfile = snakemake.output['outfasta']
    excluded2 = snakemake.output['excluded2']
    excluded1_files = snakemake.input['excluded_regions_files']
    name = snakemake.params['name']

    #combine fastas to single file
    temp_fasta = 'temp_multi_%s.fa' % name
    record_list = []
    with open(temp_fasta, "w") as g:
        for i in fastas:
            records = SeqIO.parse(i, "fasta")
            for j in records:
                record_list.append(j)
        SeqIO.write(record_list, temp_fasta, "fasta")

    alignment = alignment_funcs.write_alignment(temp_fasta, name, outfile)
    os.remove(temp_fasta)

    ex_df = pd.concat([pd.read_csv(i) for i in excluded1_files])
    if not ex_df.empty:
        ex_df['new_start'] = ex_df.apply(convert_indices, alignment = alignment, col = 'start', axis = 1)
        ex_df['new_end'] = ex_df.apply(convert_indices, alignment = alignment, col = 'end', axis = 1)
        ex_df.drop(['start', 'end'], axis = 1, inplace = True)
        ex_df['ID'] = name
        ex_df.rename(columns = {'new_start':'start', 'new_end':'end'}, inplace = True)
    ex_df.to_csv(excluded2, index = False)

if __name__ == '__main__':
    main(sys.argv[1:])
