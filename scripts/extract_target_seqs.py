#190526 MKT
#Extract the target sequences from the ncrna file
from Bio import SeqIO
import pandas as pd
import os
#import snakemake
#snake_dir = snakemake.workflow.basedir
#print('snakedir', snake_dir)

def target_fasta_provided():
    '''Figure out if there is a target fasta provided already'''
    target_df = pd.read_csv(os.path.join(snakemake.config['parameter_dir'], snakemake.config['targets']))
    file = target_df.loc[((target_df['target'] == snakemake.wildcards.target) & (target_df['organism'] == snakemake.wildcards.org)), 'file']
    if file.isnull().all():
        print('no files given')
        return None
    else:
        target_fasta = file.values[0]
        return target_fasta

def collect_seqs_by_id(infasta, ids, outfile):
    '''
    Filter fasta, writing only ones matching the ids to the outfile
    '''
    with open(infasta, "rU") as f:
        with open(outfile, "w") as out:
            for record in SeqIO.parse(f, "fasta"):
                this_id = record.id.split('|')[0]
                if this_id in ids:
                    SeqIO.write(record, out, "fasta")
                else:
                    continue

def main(arglist):
    #print(dir(snakemake))
    infasta = snakemake.input['fasta_file']
    ids = snakemake.params['ids']
    target_file = snakemake.params['target_file']
    snakedir = snakemake.params['snakedir']
    outfile = snakemake.output[0]

    print('snakedir', snakedir)
    print('cwd', os.getcwd())
    print('outfile', outfile)

    target_fasta = target_fasta_provided()
    if target_fasta:
        os.rename(os.path.join(os.getcwd(),target_fasta), os.path.join(os.getcwd(),outfile))
    else:
        collect_seqs_by_id(infasta, ids, outfile)

if __name__ == '__main__':
    main(sys.argv[1:])
'''
FileNotFoundError: [Errno 2] No such file or directory: '/Users/maryk.thompson/Desktop/Davislab/OneDriveNexus365/2.25_rRNA_sub_2b/scripts/ribopop_pipeline/testyeast3/myfastas/Scer-18S.fa' -> '/Users/maryk.thompson/Desktop/Davislab/OneDriveNexus365/2.25_rRNA_sub_2b/scripts/ribopop_pipeline/target_sequences/Scer-18S.fa'
'''
