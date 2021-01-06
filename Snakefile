#Snakefile to design oligos for Ribopop
import pandas as pd
import os

snakedir = workflow.basedir
configfile: 'config.yml'

ann_df = pd.read_csv(os.path.join(config['parameter_dir'], config['seqs_and_annotations'])).set_index('organism', drop = False)
target_df = pd.read_csv(os.path.join(config['parameter_dir'], config['targets']))
targets = target_df['target'].unique()
orgs = target_df['organism'].unique()
param_df = pd.read_csv(os.path.join(config['parameter_dir'], config['params'])).set_index('target', drop = False)

rule all:
    input:
        all_selected_probes = 'probe_design/all_selected_probes.csv',
        plots = expand('probe_design/{target}/selected_probes_{target}.png', target = targets)

#get target sequences out of the ncrna file or other provided file
rule extract_targets:
    input:
        fasta_file = 'offtarget_filtering/{org}/ncrna.fa'
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

#download seqs for the offtarget screening
rule download_seqs:
    output:
        ncrna_file = temp('offtarget_filtering/{org}/ncrna.fa'),
        cdna_file = temp('offtarget_filtering/{org}/cdna.fa'),
        gtf_file = temp('offtarget_filtering/{org}/ann.gtf'),
        genome_file = temp('offtarget_filtering/{org}/genome.fa')
    params:
        outdir = 'offtarget_filtering/{org}/',
        to_download = {'genome':'genome.fa', 'cdna': 'cdna.fa', 'ncrna':'ncrna.fa', 'gtf': 'ann.gtf'}
    script:
        'scripts/download_seqs.py'

#extract intronic sequences from the genome to include in the blast database
rule extract_intronic_seqs:
    input:
        genome_fasta = 'offtarget_filtering/{org}/genome.fa',
        gtf_file = 'offtarget_filtering/{org}/ann.gtf'
    params:
        flanking_nt = 40
    output:
        outfasta = temp('offtarget_filtering/{org}/introns.fa')
    script:
        'scripts/extract_intronic_seqs.py'

#make a combined file of mRNAs, ncRNAs, and intronic sequences
rule concatenate_transcripts:
    input:
        cdna_fasta = 'offtarget_filtering/{org}/cdna.fa',
        ncrna_fasta = 'offtarget_filtering/{org}/ncrna.fa',
        intron_fasta = 'offtarget_filtering/{org}/introns.fa'
    output:
        txt_fasta = 'offtarget_filtering/{org}/txts.fa'
    shell:
        'cat {input.cdna_fasta} {input.ncrna_fasta} {input.intron_fasta} > {output.txt_fasta}'

#build transcript blast index
rule build_index_transcripts:
    input:
        fasta = 'offtarget_filtering/{org}/txts.fa'
    params:
        build_index = True,
        outname = 'offtarget_filtering/{org}/txts'
    output:
        an_index_file = 'offtarget_filtering/{org}/txts_completed.txt'
    script:
        'scripts/run_blastn.py'

#find regions in transcript db that are homologous to the target --i.e. rRNAs
rule make_rRNA_homology_target:
    input:
        an_index_file = 'offtarget_filtering/{org}/txts_completed.txt',
        subject_fasta = 'offtarget_filtering/{org}/txts.fa',
        query_fasta = 'target_sequences/aln_by_org/{org}/{target}.fa'
    params:
        blast_txts = True
    output:
        outfile = 'offtarget_filtering/{org}/{target}-homology_target_blast.csv'
    script:
        'scripts/run_blastn.py'

#find alignments of potential probes with transcripts
rule make_rRNA_homology_kmers:
    input:
        an_index_file = 'offtarget_filtering/{org}/txts_completed.txt',
        subject_fasta = 'offtarget_filtering/{org}/txts.fa',
        query_fasta = 'probe_design/{target}/potential_probes.fa'
    params:
        blast_kmers = True,
        min_bitscore = config['min_bitscore'],
        evalue = config['evalue']
    output:
        outfile = 'offtarget_filtering/{org}/{target}-homology_kmers_blast.csv'
    script:
        'scripts/run_blastn.py'

#remove potential probes that overlap with non-rRNA transcripts
rule remove_homologous_kmers:
    input:
        probe_csv = 'probe_design/{target}/potential_probes.csv',
        kmer_homology_files = expand('offtarget_filtering/{org}/{{target}}-homology_kmers_blast.csv', org = orgs),
        target_homology_files = expand('offtarget_filtering/{org}/{{target}}-homology_target_blast.csv', org = orgs)
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
