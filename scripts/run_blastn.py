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

def build_blast_index(fasta_file, title):
    print('building index')
    subprocess.run(['makeblastdb', '-in', fasta_file, '-parse_seqids', '-title', title, '-dbtype', 'nucl'], shell = False)
    #create file to show snakemake that the index is completed
    touch_file = os.path.join(os.getcwd(), '%s_completed.txt' % title)
    Path(touch_file).touch()

def run_blast(subject_fasta, query_fasta, outname):
    temp_csv = '{outname}.csv'.format(outname = outname)
    cmd = ' '.join(['blastn', '-db', subject_fasta, '-query', query_fasta, '-outfmt', '10', '-out', temp_csv])
    subprocess.check_call(cmd, shell = True)
    df = pd.read_csv(temp_csv, names = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])
    df.to_csv(outname)
    os.remove(temp_csv)

def main(arglist):
    parser = argparse.ArgumentParser()
    parser.add_argument('-fasta', help = 'name of fasta file for building index')
    parser.add_argument('--build_index', action = 'store_true', help = 'set if you want to build the blast index')
    parser.add_argument('-subject_fasta', help = 'name of fasta file to be searched')
    parser.add_argument('-query_fasta', help = 'name of query fasta file')
    parser.add_argument('-outname', help ='prefix for output file')
    args = parser.parse_args()

    #A hack to get the argparse defaults overriden by snakemake args, if provided
    #snakemake is a global variable
    if 'snakemake' in globals():
        #If the script is called by snakemake, reassign args from input, output, and params, where specified
        toset = defaultdict(list)
        argdict = vars(args)
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

    if args.build_index == True:
        build_blast_index(args.fasta, args.outname)
    else:
        run_blast(args.subject_fasta, args.query_fasta, args.outname)

if __name__ == '__main__':
    main(sys.argv[1:])

#the subject needs to be the name of the fasta, the query is the .fa file
