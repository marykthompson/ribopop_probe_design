'''
Given a single fasta file or list of fasta files, use mafft to create the
multiple sequence alignment. Then parse the alignments to create a consensus sequence.
'''
from Bio import SeqIO
import subprocess
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import os

def write_alignment(fasta, outname):
    '''
    Make an alignment and consensus sequence for seqs in the provided fasta
    If you provide the outdir as the family outdir, then this will be distinguished from
    the species specific consensus by being one directory up.
    '''

    name = snakemake.params['name']

    #Check if only 1 sequence in the fasta file. If so, then just write this sequence
    with open(fasta, 'r') as f:
        num_seqs = len([1 for line in f if line.startswith(">")])

    if num_seqs == 1:
        record = SeqIO.parse(fasta, 'fasta')
        SeqIO.write(record, outname, 'fasta')

    else:
        subprocess.run('mafft --auto --clustalout %s > %s.clustal' % (fasta, outname), shell = True)
        align = AlignIO.read('%s.clustal' % outname, 'clustal')
        info = AlignInfo.SummaryInfo(align)
        #Get the sequence with alignment gaps or mismatches turned into Ns, but with the ends filled in with any sequence that covers it
        # ungapped cons will not take gaps into acount, we can get the start and end sequence from here after applying the gapped_cons
        ungapped_cons = SeqRecord(info.dumb_consensus(threshold = 1, ambiguous = 'N'), id = '%s_consensus' % name, description = "")

        #'require multiple' means it will place an N if only one sequence unaligned
        # this is the desired behavior for internal gaps but not for the end gaps -- because some of the sequences are incomplete
        gapped_cons = SeqRecord(info.gap_consensus(threshold = 1, ambiguous = 'N', require_multiple = 1), id = '%s_consensus' % name, description = "")

        gapped_len = len(gapped_cons.seq)
        #get the indices where the gapped alignment should be used
        left_strip = gapped_cons.seq.lstrip('N')
        left_start = gapped_len - len(left_strip)

        right_strip = gapped_cons.seq.rstrip('N')
        right_start = len(right_strip)

        left_seq = ungapped_cons.seq[0:left_start]
        right_seq = ungapped_cons.seq[right_start:]
        middle_seq = gapped_cons.seq[left_start: right_start]

        whole_seq = left_seq + middle_seq + right_seq
        record = SeqRecord(whole_seq, id = '%s_cons' % name, description = '')
        SeqIO.write(record, outname, 'fasta')

def main(arglist):
    #I think this is a more robust way to check if we are dealing with a single file or list because snakemake converts the type of things
    outfile = snakemake.output['outfasta']
    fasta_list = snakemake.input['fastas']
    #if it's a list of files, the first item would always have more than 1 character
    if len(fasta_list) > 1:
    #check if there is one fasta file or >1 fasta file and combine to make a temporary file
        print('more than one fasta!')
        temp_fasta = 'temp_multi.fa'
        record_list = []
        with open(temp_fasta, "w") as g:
            for i in fasta_list:
                records = SeqIO.parse(i, "fasta")
                for j in records:
                    record_list.append(j)

            SeqIO.write(record_list, temp_fasta, "fasta")
        this_fasta = temp_fasta
    else:
        this_fasta = fasta_list[0]

    write_alignment(this_fasta, outfile)
    if os.path.exists(temp_fasta):
        os.remove(temp_fasta)

if __name__ == '__main__':
    main(sys.argv[1:])
