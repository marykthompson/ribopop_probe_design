#190526 MKT
#Extract the target sequences from the ncrna file or the provided file from the targets.csv file
from Bio import SeqIO
import pandas as pd
import os

def collect_seqs_by_id(infasta, ids):
    '''
    Go through the provided fasta and return records matching the provided ids.
    '''
    records = []
    with open(infasta, 'r') as f:
        for record in SeqIO.parse(f, 'fasta'):
            this_id = record.id.replace('|', ' ').split(' ')[0]
            if this_id in ids:
                records.append(record)
    return records

def main(arglist):
    infasta = snakemake.input['fasta_file']
    target_df = snakemake.params['sub_target_df']
    outfile = snakemake.output[0]

    #Extract targets from the given file or the ncrna.fa file
    #If there is no file given, assume we want the ncrna.fa file
    target_df['file'] = target_df['file'].fillna(infasta)

    #split these into set of ids and fasta
    records = []
    for i in target_df['file'].unique():
        ids = target_df.loc[target_df['file'] == i, 'ID'].values
        records.append(collect_seqs_by_id(i, ids))

    with open(outfile, 'w') as out:
        for record in records:
            SeqIO.write(record, out, 'fasta')

if __name__ == '__main__':
    main(sys.argv[1:])
