#use_prebuilt_indices.smk
#run this version if you want to design another set of probes using a blast database that is already built.
import pandas as pd
import os

snakedir = workflow.basedir
thisdir = os.getcwd()

configfile: 'config.yml'
ann_df = pd.read_csv(os.path.join(config['parameter_dir'], config['indices'])).set_index(['organism', 'target'], drop = False)
target_df = pd.read_csv(os.path.join(config['parameter_dir'], config['targets']))
param_df = pd.read_csv(os.path.join(config['parameter_dir'], config['params'])).set_index('target', drop = False)

targets = target_df['target'].unique()
orgs = target_df['organism'].unique()
index_orgs = ann_df['organism'].unique()

targets.sort()
orgs.sort()
index_orgs.sort()

#softlink the target_homology files if they are specified in the index file so that they won't
#be recreated for other organisms if they are not being included in the target consensus
for i in ann_df.index:
    homol_file = ann_df.loc[i, 'target_homology']
    target = ann_df.loc[i, 'target']
    org = ann_df.loc[i, 'organism']
    if not pd.isnull(homol_file):
        new_path = os.path.join(thisdir, f'offtarget_filtering/{org}/{target}-homology_target_blast.csv')
        os.makedirs(os.path.dirname(new_path), exist_ok = True)
        if os.path.exists(new_path):
            os.remove(new_path)
        os.symlink(homol_file, new_path)

def get_targets(wildcards):
    '''
    Get the target IDs from the target.csv file.
    '''
    return target_df.loc[((target_df['target'] == wildcards.target) & (target_df['organism'] == wildcards.org)), 'ID'].values

rule all:
    input:
        all_selected_probes = 'probe_design/all_selected_probes.csv',
        plots = expand('probe_design/{target}/selected_probes_{target}.png', target = targets)

#get target sequences out of the provided file in targets.csv
#note that for this version, you need to provide the filename in the target sequences file.
rule extract_targets:
    params:
        sub_target_df = lambda wildcards: target_df[(target_df['target'] == wildcards.target) & (target_df['organism'] == wildcards.org)]
    output:
        'target_sequences/original/{org}/{target}.fa'
    script:
        'scripts/extract_target_seqs.py'

#align targets by organism, then target. Convert exlcuded region coordinates to alignment coordinates.
rule make_alignment_by_organism:
    input:
        fasta = 'target_sequences/original/{org}/{target}.fa'
    params:
        target_df = target_df,
        name = '{org}-{target}'
    output:
        outfasta = 'target_sequences/aln_by_org/{org}/{target}.fa',
        excluded1 = temp('target_sequences/{org}-{target}.excluded1')
    script:
        'scripts/align_by_organism.py'

rule make_alignment_by_target:
    input:
        fastas = expand('target_sequences/aln_by_org/{org}/{{target}}.fa', org = orgs),
        excluded_regions_files = expand('target_sequences/{org}-{{target}}.excluded1', org = orgs)
    params:
        name = '{target}'
    output:
        outfasta = 'target_sequences/consensus/{target}.fa',
        excluded2 = 'target_sequences/{target}.excluded.csv'
    script:
        'scripts/align_by_target.py'

#mask simple repeats so they won't end up in the probes
rule tantan_mask:
    input:
        'target_sequences/consensus/{target}.fa'
    output:
        'target_sequences/masked/{target}.fa'
    shell:
        'tantan -x N {input} > {output}'

#find regions in transcript db that are homologous to the target --i.e. rRNAs
rule make_rRNA_homology_target:
    input:
        subject_fasta = lambda wildcards: ann_df.loc[[wildcards.target, wildcards.index_org], 'blast_db'],
        query_fasta = 'target_sequences/aln_by_org/{index_org}/{target}.fa'
    params:
        blast_txts = True
    output:
        outfile = 'offtarget_filtering/{index_org}/{target}-homology_target_blast.csv'
    script:
        'scripts/run_blastn.py'

#find alignments of potential probes with transcripts
rule make_rRNA_homology_kmers:
    input:
        subject_fasta = lambda wildcards: ann_df.loc[[wildcards.target, wildcards.index_org], 'blast_db'],
        query_fasta = 'probe_design/{target}/potential_probes.fa'
    params:
        blast_kmers = True,
        min_bitscore = config['min_bitscore'],
        evalue = config['evalue']
    output:
        outfile = 'offtarget_filtering/{index_org}/{target}-homology_kmers_blast.csv'
    script:
        'scripts/run_blastn.py'

#remove potential probes that overlap with non-rRNA transcripts
rule remove_homologous_kmers:
    input:
        probe_csv = 'probe_design/{target}/potential_probes.csv',
        kmer_homology_files = expand('offtarget_filtering/{index_org}/{{target}}-homology_kmers_blast.csv', index_org = index_orgs),
        target_homology_files = expand('offtarget_filtering/{index_org}/{{target}}-homology_target_blast.csv', index_org = index_orgs)
    output:
        filtered_probe_csv = 'probe_design/{target}/potential_probes_filt.csv'
    script:
        'scripts/filter_probes.py'

#make potential probes and screen by provided constraints
rule screen_kmers:
    input:
        target_fasta = 'target_sequences/masked/{target}.fa',
        excluded_regions = 'target_sequences/{target}.excluded.csv'
    params:
        min_probe_len = lambda wildcards: param_df.loc[wildcards.target, 'min_probe_length'],
        max_probe_len = lambda wildcards: param_df.loc[wildcards.target, 'max_probe_length'],
        min_Tm = lambda wildcards: param_df.loc[wildcards.target, 'min_Tm'],
        max_Tm = lambda wildcards: param_df.loc[wildcards.target, 'max_Tm'],
        min_gc = lambda wildcards: param_df.loc[wildcards.target, 'min_GC'],
        max_gc = lambda wildcards: param_df.loc[wildcards.target, 'max_GC'],
        Tm_quantile = lambda wildcards: param_df.loc[wildcards.target, 'Tm_quantile'],
        Tm_window_size = lambda wildcards: param_df.loc[wildcards.target, 'Tm_window_size'],
        min_hairpin_dG = lambda wildcards: param_df.loc[wildcards.target, 'min_hairpin_dG'],
        min_dimer_dG = lambda wildcards: param_df.loc[wildcards.target, 'min_dimer_dG'],
        excluded_regions_consensus = lambda wildcards: param_df.loc[wildcards.target, 'excluded_regions_consensus'],
        logfile = 'probe_design/{target}/log.txt'
    output:
        probe_csv = 'probe_design/{target}/potential_probes.csv',
        probe_fa = 'probe_design/{target}/potential_probes.fa'
    script:
        'scripts/screen_kmers.py'

#choose N number of probes from potential probes, as evenly spaced as possible.
rule choose_probes:
    input:
        probe_csvs = expand('probe_design/{target}/potential_probes_filt.csv', target = targets),
        target_fastas = expand('target_sequences/consensus/{target}.fa', target = targets),
        excluded_regions = expand('target_sequences/{target}.excluded.csv', target = targets)
    params:
        logfile = 'probe_design/log.txt',
        desired_number_probes = param_df.loc[targets, 'number_probes'],
        target_subregions_consensus = param_df.loc[targets, 'target_subregions_consensus'],
        min_dimer_dG = param_df.loc[targets, 'min_dimer_dG']
    output:
        plots = expand('probe_design/{target}/selected_probes_{target}.png', target = targets),
        all_selected_probes = 'probe_design/all_selected_probes.csv'
    script:
        'scripts/choose_probes.py'
