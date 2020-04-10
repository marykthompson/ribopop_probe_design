'''
mask_gene.py
Purpose:
To determine which sites in the target sequence should be avoided due to
specificity issues in the genome.
- Read in target sequence and break up into kmers
- Align kmers against the genome and transcriptome
- Screen aligned kmers against previously detected regions of rRNA homology
and discard those alignments (keeping those kmers as probe candidates).
- Screen genome-aligned kmers against gene regions and discard alignments
falling outside annotated gene regions
- Write a file containining the start position of the remaining alignments
(Probes overlapping these sites will be discarded in the probe design step).
'''

import os
import subprocess
import argparse
import sys
from Bio import SeqIO
import HTSeq
from collections import defaultdict
import pandas as pd

def create_genomic_interval(start, end, chrom):
    '''
    Given 1-based blast start, end and chromosome, return an HTSeq genomic interval.
    '''
    #set start < end (blast returns start > end if maps to - strand)
    if start > end:
        start2 = end
        end2 = start
        strand = '-'
    else:
        start2 = start
        end2 = end
        strand = '+'

    start2 = start2 - 1
    iv = HTSeq.GenomicInterval(chrom, start2, end2, strand)
    return iv

def create_alnmt(start, end, chrom, name):
    '''
    Given 1-based blast start, end and chromosome, return an HTSeq Alignment object
    '''
    iv = create_genomic_interval(start, end, chrom)
    #note that name is string of 0-based start of alignment
    r = HTSeq.SequenceWithQualities(b"", str(name), b"")
    alnmt = HTSeq.Alignment(r, iv)
    return alnmt

def generate_kmers(fastafile, klength, outdir, outname):
    '''
    Generate kmers from the target sequence.
    Note that there will sometimes be Ns in the kmers from the consensus sequence.
    These will not be output because they will be screened out anyway in a downstream
    probe design step.
    '''
    print('using klength %s' % klength)
    target = str(next(SeqIO.parse(fastafile, 'fasta')).seq)
    kmers = []
    for i in range(0, len(target) - klength + 1):
        kmer = target[i: i + klength]
        if 'N' in kmer:
            continue
        else:
            kmers.append((i, kmer))

    #write kmer fasta
    kmer_file = os.path.join(outdir, '%s_kmers.fa' % outname)
    with open(kmer_file, 'w') as g:
        for pos, kmer in kmers:
            g.write('>%s\n%s\n' % (pos, kmer))
    return kmer_file

def blast_kmers(kmer_file, fasta, outdir, outname, min_bitscore = 30, evalue = 10):
    '''
    Blast the kmers to genome or cDNA collection with blastn.
    '''
    outfile = '{outprefix}.csv'.format(outprefix = os.path.join(outdir, outname))
    cmd = ' '.join(['blastn', '-task', 'blastn-short', '-dust', 'no', '-soft_masking',
    'false', '-db', fasta, '-query', kmer_file, '-outfmt', '10', '-evalue', str(evalue), '-out', outfile])
    subprocess.check_call(cmd, shell = True)
    df = pd.read_csv(outfile, names = ['qseqid', 'sseqid', 'pident', 'length',
    'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])
    aln_df = df[df['bitscore'] >= min_bitscore].copy()

    #write alignments passing score cutoff
    aln_df.to_csv(outfile)
    return aln_df

def overlapper(blast_df, regions, discard_minus_strand = False, mode = 'keep'):
    '''
    Given alignments and genomic array, return aligments that do or do not overlap the array
    mode = keep -> return alignments that overlap the array (use for genes)
    mode = discard -> discard alignments that overlap the array (use for rRNA)
    The set of alignments will be used to create the masked target sequence.
    '''
    indices = []
    for row in blast_df.itertuples():
        if (discard_minus_strand == True) and (row.send < row.sstart):
            continue

        alnmt = create_alnmt(row.sstart, row.send, row.sseqid, row.qseqid)

        iset = None
        for iv2, step_set in regions[alnmt.iv].steps():
            if iset is None:
                iset = step_set.copy()
            else:
                iset.intersection_update(step_set)

        #if aligned to at least 1 feature, keep alignment
        blast_df.at[row.Index, 'gene'] = ','.join(iset)

        if len(iset) > 0:
            if mode == 'keep':
                indices.append(row.Index)

        #otherwise, iset = 0, doesn't align to rRNA, so keep that position masked if mode == discard
        else:
            if mode == 'discard':
                indices.append(row.Index)

    filtered_df = blast_df.reindex(indices)
    return filtered_df

def ol_alnmts(blast_df, homol_csv = None, gtf = None, discard_minus_strand = False):
    '''
    This function has two main uses. Either:
    1) discard genomic alignments if they overlap with regions homologous to rRNA
    and then screen these for alignment to genes in the gtf. If they align, mask
    this start site.
    2) discard transcriptomic alignments if they overlap with rRNA transcript
    homology, otherwise mask these start sites.
    The discard_minus_strand option is for alignment to the transcripts.
    The order of the alignment calling is important in this case. When you pass it
    the rRNA homology csv, it only returns kmers that do not align to those regions.
    At the next step it screens these passed alignments against the gtf to determine
    if they fall in gene regions.
    '''
    #remove alignments that overlap rRNA homology regions
    if homol_csv:
        homol_df = pd.read_csv(homol_csv)
        regions = [create_genomic_interval(x, y, z) for x, y, z in zip(homol_df['sstart'], homol_df['send'], homol_df['sseqid'])]
        rRNA_genes = HTSeq.GenomicArrayOfSets("auto", stranded = True)
        for r in regions:
            rRNA_genes[r] += 'rRNA'

        blast_df = overlapper(blast_df, rRNA_genes, discard_minus_strand = discard_minus_strand, mode = 'discard')

    #now check which of the remaining alignments fall within gene regions.
    if gtf:
        gtf_file = HTSeq.GFF_Reader(gtf, end_included = True)
        genes = HTSeq.GenomicArrayOfSets( "auto", stranded = True)

        for feature in gtf_file:
            if feature.type == "gene":
                genes[feature.iv] += feature.name

        blast_df = overlapper(blast_df, genes, discard_minus_strand = discard_minus_strand, mode = 'keep')

    masked_indices = set(blast_df['qseqid'].tolist())
    return masked_indices, blast_df

def main(arglist):
    parser = argparse.ArgumentParser()
    parser.add_argument('-target_fasta', help = 'fasta file containing target sequence(s)')
    parser.add_argument('-min_probe_length', type = int, help = 'min probe length to use for mask generation')
    parser.add_argument('-max_probe_length', type = int, help = 'max probe length to use for mask generation')
    parser.add_argument('-orgs', nargs = '+', help = 'list of organism names in same oreder as genome fastas')
    parser.add_argument('-outdir')
    parser.add_argument('-outfile', help = 'the name of the masked_nts file, i.e. 25S_masked.csv')
    parser.add_argument('-genome_fasta', nargs = '+', help = 'fasta files for the genome')
    parser.add_argument('-txt_fasta', nargs = '+', help = 'fasta files for the cDNA')
    parser.add_argument('-genome_homology_file', nargs = '+', help = 'csv files(s) containing accepted homology regions from blast')
    parser.add_argument('-txt_homology_file', nargs = '+', help = 'csv files(s) containing accepted homology regions from blast')
    parser.add_argument('-gtf', nargs = '+', help = 'gtf files, expecting a feature "gene" to be present, based on Ensembl format')
    parser.add_argument('-min_bitscore', type = float, help = 'bitscore of blast hit must be at least this high to be screened out of potential probes')
    args = parser.parse_args(arglist)

    #A hack to get the argparse defaults overriden by snakemake args, if provided
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

    os.makedirs(args.outdir, exist_ok = True)
    kmer_outdir = os.path.join(args.outdir, 'kmers')
    os.makedirs(kmer_outdir, exist_ok = True)
    kmer_file_dict = {}

    for probe_len in range(args.min_probe_length, args.max_probe_length + 1):
        kmer_file_dict[probe_len] = generate_kmers(args.target_fasta, probe_len, kmer_outdir, 'probelen_%s' % probe_len)

    #store a list of dictionaries with {length:{masked indices}}
    masked_start_l = []
    for i in range(0, len(args.genome_fasta)):
        aln_check_dir = os.path.join(args.outdir, 'aln_files', args.orgs[i])
        os.makedirs(aln_check_dir, exist_ok = True)
        for probe_len in range(args.min_probe_length, args.max_probe_length + 1):
            #note that I'm requiring it to be on the same strand to be filtered out by alignment
            blast_df_genome = blast_kmers(kmer_file_dict[probe_len], args.genome_fasta[i], aln_check_dir, 'probelen_%snt_genome' % probe_len, min_bitscore = args.min_bitscore)
            masked_indices_genome, filtered_df_genome = ol_alnmts(blast_df_genome, homol_csv = args.genome_homology_file[i], gtf = args.gtf[i])
            masked_start_l.append({probe_len: masked_indices_genome})

            #add masked starts based on kmer alignments to transcriptome
            blast_df_txt = blast_kmers(kmer_file_dict[probe_len], args.txt_fasta[i], aln_check_dir, 'probelen_%snt_cdna' % probe_len, min_bitscore = args.min_bitscore)
            masked_indices_txt, filtered_df_txt = ol_alnmts(blast_df_txt, homol_csv = args.txt_homology_file[i], discard_minus_strand = True)
            masked_start_l.append({probe_len: masked_indices_txt})

            #print('not in cdna', masked_indices_genome.difference(masked_indices_txt))
            #print('not in genome', masked_indices_txt.difference(masked_indices_genome))

    #get the union of all masked starts from all probe lengths and all species
    all_masked_start_dict = {k: set() for k in set().union(*masked_start_l)}
    for l in all_masked_start_dict:
        all_masked_start_dict[l] = set().union(*[i[l] for i in masked_start_l if l in i])

    masked_df = pd.DataFrame.from_dict(all_masked_start_dict, orient = 'index')
    transposed_df = masked_df.transpose()
    cols = sorted(transposed_df.columns)
    new_cols = {i: 'len_%s' % i for i in cols}
    transposed_df = transposed_df[cols]
    transposed_df.rename(columns = new_cols, inplace = True)
    transposed_df.to_csv(args.outfile, index = False)

if __name__ == "__main__":
    main(sys.argv[1:])
