#Snakefile to design oligos for Ribopop
import pandas as pd
import os
#import pdb; pdb.set_trace()

snakedir = workflow.basedir
configfile: 'config.yml'

ann_df = pd.read_csv(os.path.join(config['parameter_dir'], config['seqs_and_annotations'])).set_index('organism', drop = False)
target_df = pd.read_csv(os.path.join(config['parameter_dir'], config['targets']))
targets = target_df['target'].unique()
orgs = target_df['organism'].unique()
targets.sort()
orgs.sort()
param_df = pd.read_csv(os.path.join(config['parameter_dir'], config['params'])).set_index('target', drop = False)

def get_targets(wildcards):
    '''
    Get the target IDs from the target.csv file.
    '''
    return target_df.loc[((target_df['target'] == wildcards.target) & (target_df['organism'] == wildcards.org)), 'ID'].values

def get_masked_files(wildcards):
    '''
    Return all the _masked.csv files for a target for the choose_probes rule.
    Adding the function instead of expand() makes it not break on cases where
    some targets are not provided for all organisms.
    If offtarget_screening is False, return [] and it will not be performed.
    '''
    if config['offtarget_screening']:
        masked_nts = expand('offtarget_filtering/masked/{target}_masked.csv', target = targets)
        return masked_nts
    else:
        return []

rule all:
    input:
        #"target_sequences/aln_by_org/Scer/18S.fa",
        ##expand('target_sequences/consensus/{target}.fa', target = targets)
        #"target_sequences/aln_by_org/Spom/18S.fa",
        expand('probe_design/{target}/selected_probes_{target}.csv', target = targets)

rule download_ncrna:
    params:
        outdir = 'offtarget_filtering/{org}/',
        to_download = {'ncrna':'temp_ncrna.fa'}
    output:
        ncrna_file = 'offtarget_filtering/{org}/temp_ncrna.fa'
    script:
        "scripts/download_seqs.py"

rule extract_targets:
    input:
        fasta_file = 'offtarget_filtering/{org}/temp_ncrna.fa'
    params:
        target_file = os.path.join(config['parameter_dir'], config['targets']),
        ids = get_targets,
        snakedir = snakedir
    conda:
        'envs/probe_design.yaml'
    output:
        "target_sequences/original/{org}/{target}.fa"
    script:
        'scripts/extract_target_seqs.py'

#first convert the excluded regions from specific target to the organism consensus
rule make_alignment_by_organism:
    input:
        fasta = 'target_sequences/original/{org}/{target}.fa'
    params:
        target_df = target_df,
        name = '{org}-{target}'
    output:
        outfasta = 'target_sequences/aln_by_org/{org}/{target}.fa',
        excluded1 = 'target_sequences/{org}-{target}.excluded1'
    conda:
        'envs/probe_design.yaml'
    script:
        'scripts/align_by_organism.py'

rule make_alignment_by_family:
    input:
        fastas = expand("target_sequences/aln_by_org/{org}/{{target}}.fa", org = orgs),
        excluded_regions_files = expand('target_sequences/{org}-{{target}}.excluded1', org = orgs)
    params:
        name = "{target}"
    output:
        outfasta = "target_sequences/consensus/{target}.fa",
        excluded2 = 'target_sequences/{target}.excluded2'
    conda:
        'envs/probe_design.yaml'
    script:
        'scripts/align_by_target.py'

rule download_seqs:
    output:
        ncrna_file = 'offtarget_filtering/{org}/ncrna.fa',
        cdna_file = 'offtarget_filtering/{org}/cdna.fa',
        gtf_file = 'offtarget_filtering/{org}/ann.gtf',
        genome_file = 'offtarget_filtering/{org}/genome.fa'
    params:
        outdir = 'offtarget_filtering/{org}/',
        to_download = {'genome':'genome.fa', 'cdna': 'cdna.fa', 'ncrna':'ncrna.fa', 'gtf': 'ann.gtf'}
    script:
        'scripts/download_seqs.py'

rule build_index_genome:
    input:
        fasta = 'offtarget_filtering/{org}/genome.fa'
    params:
        build_index = True,
        outname = 'genometest'
    output:
        an_index_file = 'offtarget_filtering/{org}/genome.fa.nhr'
    conda:
        'envs/probe_design.yaml'
    script:
        'scripts/run_blastn.py'

rule concatenate_transcripts:
    input:
        cdna_fasta = 'offtarget_filtering/{org}/cdna.fa',
        ncrna_fasta = 'offtarget_filtering/{org}/ncrna.fa'
    output:
        txt_fasta = 'offtarget_filtering/{org}/txts.fa'
    shell:
        'cat {input.cdna_fasta} {input.ncrna_fasta} > {output.txt_fasta}'

rule build_index_transcripts:
    input:
        fasta = 'offtarget_filtering/{org}/txts.fa'
    params:
        build_index = True,
        outname = 'txttest'
    output:
        an_index_file = 'offtarget_filtering/{org}/txts.fa.nhr'
    conda:
        'envs/probe_design.yaml'
    script:
        'scripts/run_blastn.py'

rule make_rRNA_homology_genome:
    input:
        an_index_file = 'offtarget_filtering/{org}/genome.fa.nhr',
        subject_fasta = 'offtarget_filtering/{org}/genome.fa',
        query_fasta = 'target_sequences/consensus/{target}.fa'
    output:
        outname = 'offtarget_filtering/{org}/{target}-homology_genome_blast.csv'
    conda:
        'envs/probe_design.yaml'
    script:
        'scripts/run_blastn.py'

rule make_rRNA_homology_txts:
    input:
        an_index_file = 'offtarget_filtering/{org}/txts.fa.nhr',
        subject_fasta = 'offtarget_filtering/{org}/txts.fa',
        query_fasta = 'target_sequences/consensus/{target}.fa'
    output:
        outname = 'offtarget_filtering/{org}/{target}-homology_txts_blast.csv'
    conda:
        'envs/probe_design.yaml'
    script:
        'scripts/run_blastn.py'
#an_index_file is a way to link the rules. Blast indexing does not rename the file
#you still just give it the fasta file and it looks for the related files in the same directory
rule mask_nts:
    input:
        some_index_files = expand(['offtarget_filtering/{org}/genome.fa.nhr', 'offtarget_filtering/{org}/txts.fa.nhr'], org = orgs),
        target_fasta = 'target_sequences/consensus/{target}.fa',
        genome_fasta = expand('offtarget_filtering/{org}/genome.fa', org = orgs),
        txt_fasta = expand('offtarget_filtering/{org}/txts.fa', org = orgs),
        genome_homology_file = expand('offtarget_filtering/{org}/{{target}}-homology_genome_blast.csv', org = orgs),
        txt_homology_file = expand('offtarget_filtering/{org}/{{target}}-homology_txts_blast.csv', org = orgs),
        gtf = expand('offtarget_filtering/{org}/ann.gtf', org = orgs)
    params:
        orgs = orgs,
        min_probe_length = lambda wildcards: param_df.loc[wildcards.target, 'min_probe_length'],
        max_probe_length = lambda wildcards: param_df.loc[wildcards.target, 'max_probe_length'],
        min_bitscore = config['min_bitscore'],
        outdir = 'offtarget_filtering/masked/{target}'
    output:
        outfile = "offtarget_filtering/masked/{target}_masked.csv" #doesn't work only uses one target
    conda:
        'envs/probe_design.yaml'
    script:
        'scripts/mask_gene.py'

rule choose_probes:
    input:
        target_fastas = expand('target_sequences/consensus/{target}.fa', target = targets),
        masked_nts = get_masked_files,
        excluded_regions = expand('target_sequences/{target}.excluded2', target = targets)
    params:
        target_names = targets,
        min_probe_length = param_df.loc[targets, 'min_probe_length'],
        max_probe_length = param_df.loc[targets, 'max_probe_length'],
        min_Tm = param_df.loc[targets, 'min_Tm'],
        max_Tm = param_df.loc[targets, 'max_Tm'],
        Tm_quantile = param_df.loc[targets, 'Tm_quantile'],
        desired_number_probes = param_df.loc[targets, 'max_number_probes'],
        Tm_window_size = param_df.loc[targets, 'Tm_window_size'],
        min_hairpin_dG = param_df.loc[targets, 'min_hairpin_dG'],
        min_dimer_dG = param_df.loc[targets, 'min_dimer_dG'],
        target_subregions_consensus = param_df.loc[targets, 'target_subregions_consensus'],
        excluded_regions_consensus = param_df.loc[targets, 'excluded_regions_consensus'],
        #excluded_regions = get_excluded_regions,
        masked_nts = lambda wildcards, input: input.masked_nts if input.masked_nts != [] else None,
        quick_test = True,
        design_probes = True,
        outdir = 'probe_design/'
    output:
        expand('probe_design/{target}/selected_probes_{target}.csv', target = targets)
    conda:
        'envs/probe_design.yaml'
    script:
        'scripts/choose_fish_probes.py'
