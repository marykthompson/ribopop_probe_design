#Check a probeset against the genome and transcriptome
#Useful for example to see if genes that are changing between subtracted
#and unsubtracted have more homology to the probes than other genes.

import pandas as pd
import os
snakedir = workflow.basedir
configfile: 'config_bifid_homology.yml'

outdir = config['outname']

probe_df = pd.read_csv(config['probe_file']).set_index('name', drop = False)

rule all:
    input:
        genome_alnmt_file = os.path.join(outdir, 'genome_offtarget.csv'),
        txt_alnmt_file = os.path.join(outdir, 'txt_offtarget.csv')

rule blast_probes:
    output:
        genome_alnmt_file = os.path.join(outdir, 'genome_offtarget.csv'),
        txt_alnmt_file = os.path.join(outdir, 'txt_offtarget.csv')
    input:
        probe_file = config['probe_file']
    params:
        genome_fasta = config['genome_fasta'],
        txt_fasta = config['txt_fasta'],
        genome_homology = config['genome_homology_file'],
        txt_homology = config['txt_homology_file'],
        gtf = config['gtf_file'],
        db_file = config['db_file'],
        evalue = config['evalue'],
        outdir = outdir
    conda:
        'envs/probe_design.yaml'
    script:
        'scripts/blast_probes.py'
