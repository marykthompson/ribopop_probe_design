from Bio import SeqIO
import subprocess
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def get_subregion_ranges(ranges_string):
    '''
    Convert the provided ranges string into a list of integer lists,
    Subtract 0 from both to make python based, but closed.
    If you get the index of end + 1, this might be on the other side of a gap.
    e.g. "1-2000,2300-3400" -> [(0, 1999), (2299, 3399)]
    '''
    subregion_ranges = []
    #subregion_list = [i for i in ranges_string.split(',')]
    subregion_list = [i for i in ranges_string.replace(', ',',').split(',')]
    for i in subregion_list:
        start, end = [int(j.strip()) for j in i.split('-')]
        #convert to 0-based, closed
        start -= 1
        end -= 1
        subregion_ranges.append([start, end])
    return subregion_ranges

def column_from_residue_number(aln, id, res_no):
    '''
    Get the position in the alignment corresponding to a given position in a
    specific sequence.
    #https://www.biostars.org/p/93805/
    '''
    rec = next((r for r in aln if r.id == id), None)
    j = 0
    for i, res in enumerate(rec.seq):
        if res!='-':
            if j==res_no:
                return i
            j+=1

def write_alignment(fasta, name, outname):
    '''
    Make an alignment and consensus sequence for seqs in the provided fasta
    If you provide the outdir as the family outdir, then this will be distinguished from
    the species specific consensus by being one directory up.
    '''
    #Check if only 1 sequence in the fasta file. If so, then just write this sequence
    with open(fasta, 'r') as f:
        num_seqs = len([1 for line in f if line.startswith(">")])

    if num_seqs == 1:
        record = SeqIO.parse(fasta, 'fasta')
        SeqIO.write(record, outname, 'fasta')
    else:
        subprocess.run('mafft --auto --clustalout %s > %s.clustal' % (fasta, outname), shell = True)
        alignment = AlignIO.read('%s.clustal' % outname, 'clustal')
        info = AlignInfo.SummaryInfo(alignment)
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
        record = SeqRecord(whole_seq, id = name, description = '')
        SeqIO.write(record, outname, 'fasta')
        return alignment

def convert_indices(x, alignment =  None, col = None):
    '''
    Call column_from_residue_number to add the new index to the df
    '''
    new_index = column_from_residue_number(alignment, x['ID'], x[col])
    return new_index
