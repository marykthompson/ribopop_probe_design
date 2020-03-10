#Blast provided probes against genome and transcriptome.
import pandas as pd
from mask_gene import blast_kmers
from mask_gene import ol_alnmts
from Bio.Seq import Seq
import gffutils
import os

probe_df = pd.read_csv(snakemake.input['probe_file']).set_index('name', drop = False)
probe_df['target_sequence'] = probe_df['sequence'].apply(lambda x: str(Seq(x).reverse_complement()))

genome_fasta = snakemake.params['genome_fasta']
txt_fasta = snakemake.params['txt_fasta']
genome_homology = snakemake.params['genome_homology']
txt_homology = snakemake.params['txt_homology']
gtf = snakemake.params['gtf']
outdir = snakemake.params['outdir']
genome_alnmt_file = snakemake.output['genome_alnmt_file']
txt_alnmt_file = snakemake.output['txt_alnmt_file']
db_file = snakemake.params['db_file']
evalue = str(snakemake.params['evalue'])

def write_target_file(outname, aln_df, db, convert_to_gene = False):
    if convert_to_gene == True:
        aln_df['target'] = aln_df['sseqid'].apply(lambda x: next(db.parents(x, featuretype = 'gene')).id)
    aln_df.to_csv(outname, index = False)

with open('temp.fa', 'w') as g:
    for row in probe_df.itertuples():
        g.write('>%s\n%s\n' % (row.name, row.target_sequence))

#evalue 500 seems to do the trick of getting all the 12mers aligned
blast_df_genome = blast_kmers('temp.fa', genome_fasta, outdir, 'genome_homology', evalue = '500', min_bitscore = 0)
blast_df_txt = blast_kmers('temp.fa', txt_fasta, outdir, 'txt_homology', evalue = '500', min_bitscore = 0)

masked_indices_genome, alnmts_genome = ol_alnmts(blast_df_genome, homol_csv = genome_homology, gtf = gtf)
masked_indices_txt, alnmts_txt = ol_alnmts(blast_df_txt, homol_csv = txt_homology, discard_minus_strand = True)

db = gffutils.FeatureDB(db_file)

write_target_file(genome_alnmt_file, alnmts_genome, db)
write_target_file(txt_alnmt_file, alnmts_txt, db, convert_to_gene = True)

os.remove('temp.fa')
os.remove(os.path.join(outdir, 'genome_homology.csv'))
os.remove(os.path.join(outdir, 'txt_homology.csv'))
