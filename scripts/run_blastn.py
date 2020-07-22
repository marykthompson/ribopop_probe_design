#190526 MKT
#this script is going to be run from blastn env
#will build blast index if required with --build_index option
#otherwise will blast the sequences and make the accepted homology bed file
import argparse
import subprocess
import os
import sys
import pandas as pd
from collections import defaultdict
from pathlib import Path
import probe_helpers

def build_blast_index(fasta_file, title):
    print('building index')
    subprocess.run(['makeblastdb', '-in', fasta_file, '-parse_seqids', '-title', title, '-dbtype', 'nucl'], shell = False)
    #create file to show snakemake that the index is completed
    touch_file = os.path.join(os.getcwd(), '%s_completed.txt' % title)
    Path(touch_file).touch()

def blast_kmers(subject_fasta, query_fasta, outfile, min_bitscore = 30, evalue = 50):
    '''
    Blast the kmers to genome or cDNA collection with blastn.
    '''
    cmd = ['blastn', '-task', 'blastn-short', '-dust', 'no', '-soft_masking',
    'false', '-db', subject_fasta, '-query', query_fasta, '-outfmt', '10', '-evalue', str(evalue), '-out', outfile]
    subprocess.run(cmd, shell = False)
    df = pd.read_csv(outfile, names = ['qseqid', 'sseqid', 'pident', 'length',
    'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])
    aln_df = df[df['bitscore'] >= min_bitscore].copy()
    #write alignments passing score cutoff
    aln_df.to_csv(outfile, index = False)

def blast_txts(subject_fasta, query_fasta, outfile):
    cmd = ['blastn', '-db', subject_fasta, '-query', query_fasta, '-outfmt', '10', '-out', outfile]
    subprocess.run(cmd, shell = False)
    #add header
    df = pd.read_csv(outfile, names = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])
    df.to_csv(outfile, index = False)

def main(arglist):
    parser = argparse.ArgumentParser()
    parser.add_argument('-fasta', help = 'name of fasta file for building index')
    parser.add_argument('-outname', help = 'name of output database')
    parser.add_argument('--build_index', action = 'store_true', help = 'set if you want to build the blast index')
    parser.add_argument('--blast_kmers', action = 'store_true', help = 'run short blast on potential probes')
    parser.add_argument('--blast_txts', action = 'store_true', help = 'run blast longer sequences')
    parser.add_argument('-subject_fasta', help = 'name of fasta file to be searched')
    parser.add_argument('-query_fasta', help = 'name of query fasta file')
    parser.add_argument('-min_bitscore', type = float, default = 30, help = 'only write alignments with >= than this score')
    parser.add_argument('-evalue', type = float, default = 50, help = 'use this as the E-value cutoff for the blast search.')
    parser.add_argument('-outfile', help ='name for csv output file')
    args = parser.parse_args()

    #A hack to get the argparse defaults overriden by snakemake args, if provided
    #snakemake is a global variable
    if 'snakemake' in globals():
        args = probe_helpers.set_snake_args(args, snakemake)

    #convert in case snakemake is passing a named list from pandas dataframe
    args.subject_fasta = str(args.subject_fasta)
    if args.build_index == True:
        build_blast_index(args.fasta, args.outname)
    elif args.blast_txts:
        blast_txts(args.subject_fasta, args.query_fasta, args.outfile)
    elif args.blast_kmers:
        blast_kmers(args.subject_fasta, args.query_fasta, args.outfile, min_bitscore = args.min_bitscore, evalue = args.evalue)

if __name__ == '__main__':
    main(sys.argv[1:])

#the subject needs to be the name of the fasta, the query is the .fa file
