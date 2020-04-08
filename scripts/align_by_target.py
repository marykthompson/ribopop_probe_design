'''
Given multiple fasta files (corresponding to different organisms),
use mafft to create the multiple sequence alignment for the given target.
Then parse the alignments to create a consensus sequence.
'''

import pandas as pd
import os
import alignment_funcs
from Bio import SeqIO

def main(arglist):
    fastas = snakemake.input['fastas']
    outfile = snakemake.output['outfasta']
    excluded2 = snakemake.output['excluded2']
    excluded1_files = snakemake.input['excluded_regions_files']
    name = snakemake.params['name']

    #combine fastas to single file
    temp_fasta = 'temp_multi.fa'
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
    adj_regions = []
    for rec in alignment:
        if rec.id in ex_df['ID'].values:
            excluded_regions = ex_df.loc[ex_df['ID'] == rec.id, ['start', 'end']].values
            for i in excluded_regions:
                new_start = alignment_funcs.column_from_residue_number(alignment, rec.id, i[0])
                new_end = alignment_funcs.column_from_residue_number(alignment, rec.id, i[1])
                adj_regions.append([new_start, new_end])

    new_df = pd.DataFrame(adj_regions, columns = ['start', 'end'])
    new_df['ID'] = name
    new_df.to_csv(excluded2, index = False)

if __name__ == '__main__':
    main(sys.argv[1:])
