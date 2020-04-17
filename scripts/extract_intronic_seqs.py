'''
extract_intronic_seqs.py

Use the gtf file and the genomic fasta to generate a fasta file with the
intronic sequences + flanking nt. Because the goal is to keep the file size
down, only include unique intronic sequences (i.e. group by gene and don't
write a different intron for constituitive introns of different transcript variants).
'''

import sys
import gffutils
import numpy as np
from pyfaidx import Fasta

def main(arglist):
    genome_fasta = snakemake.input['genome_fasta']
    gtf_file = snakemake.input['gtf_file']
    db_file = '%s_db' % gtf_file
    outfasta = snakemake.output['outfasta']
    flanking_nt = int(snakemake.params['flanking_nt'])

    genome = Fasta(genome_fasta)

    #create_unique necessary for yeast gtf file. Some genes and transcripts have same ID
    db = gffutils.create_db(gtf_file, db_file, force = True,
    merge_strategy = 'create_unique', disable_infer_genes = True,
    disable_infer_transcripts = True)
    introns = list(db.create_introns())
    db.update(introns, merge_strategy = 'create_unique',
    disable_infer_genes = True, disable_infer_transcripts = True, make_backup = False)

    all_genes = db.features_of_type('gene')
    with open(outfasta, 'w') as g:
        for gene in all_genes:
            strand = db[gene].strand
            chrom = db[gene].chrom
            gene_id = db[gene].id
            introns = [i for i in db.children(gene, featuretype = 'intron')]
            a = np.array([[i.start - flanking_nt - 1, i.end + flanking_nt] for i in introns])
            if a.size == 0:
                continue
            #get unique intronic intervals to write to file
            unq = np.unique(a, axis=0)
            j = 1
            for intron in unq:
                name = '%s_I%s' % (gene_id, j)
                if strand == '+':
                    seq = genome[chrom][intron[0]: intron[1]].seq
                else:
                    seq = genome[chrom][intron[0]: intron[1]].complement.seq[::-1]
                g.write('>%s\n%s\n' % (name, seq))
                j += 1

if __name__ == '__main__':
    main(sys.argv[1:])
