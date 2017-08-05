"""
This modules contains functions used to process paired-end cells.

"""

import sys
import os
import subprocess
import datetime
from Bio import SeqIO
import pysam
import operator
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

from utils import get_c_info


def analyze_chain(fasta_dict, vdj_dict, output, bam, unmapped, id_name_dict, bases, chain, strand, low_q, top, by_exp,
                  read_overlap):
    junction_segs = make_junction_file(bam, chain, output, bases, vdj_dict, fasta_dict, id_name_dict, top, by_exp,
                                       read_overlap)
    un_dict = write_reads_file(bam, unmapped, junction_segs, output, vdj_dict, chain, strand, low_q)
    return un_dict


def make_junction_file(bam, chain, output, bases, vdj_dict, fasta_dict, id_name_dict, top, by_exp, read_overlap):
    mapped_file = pysam.AlignmentFile(bam, "rb")
    if chain == 'A':
        vdj_chain_dict = vdj_dict['Alpha']
        out_name = output + '.alpha.junctions.txt'
    elif chain == 'B':
        vdj_chain_dict = vdj_dict['Beta']
        out_name = output + '.beta.junctions.txt'
    else:
        sys.stderr.write(
            str(datetime.datetime.now()) + ' Error! chain parameter for function analyzeChain can only be A or B\n')
        sys.stderr.flush()
    j_segs = vdj_chain_dict['J']
    v_segs = vdj_chain_dict['V']
    (c_seq, c_id, c_name) = get_c_info(vdj_chain_dict['C'][0], id_name_dict, fasta_dict)
    vj_segs = []
    for x in j_segs:
        vj_segs.append(x)
    for y in v_segs:
        vj_segs.append(y)
    vj_reads = dict()
    (vj_reads, vj_counts) = load_reads_to_dict(vj_segs, mapped_file, vj_reads, read_overlap)
    junction_segs = write_junctions(vj_reads, out_name, bases, fasta_dict, id_name_dict, c_seq, top, vj_counts, by_exp)
    if len(junction_segs) == 0:
        sys.stdout.write(
            str(datetime.datetime.now()) + ' Did not find any V-J reads, searching for V-C and J-C reads:\n')
        sys.stdout.flush()
        c_reads = dict()
        (c_reads, c_counts_dict) = load_reads_to_dict(vdj_chain_dict['C'], mapped_file, c_reads, read_overlap)
        junction_segs = write_junctions_with_c(vj_reads, out_name, bases, fasta_dict, id_name_dict, c_reads)
    mapped_file.close()
    return junction_segs


def write_junctions(vj_reads, out_name, bases, fasta_dict, id_name_dict, c_seq, top, vj_counts_dict, by_exp):
    out = open(out_name, 'w')
    f_arr = []
    pair_count_dict = dict()
    for seg in vj_reads:
        if id_name_dict[seg].find('J') != -1:
            if len(vj_reads[seg]) > 0:
                for s_seg in vj_reads:
                    if ((id_name_dict[s_seg].find('V') != -1) & (len(vj_reads[s_seg]) > 0)):
                        if (len([val for val in vj_reads[seg]['first'] if val in vj_reads[s_seg]['second']]) > 0) | \
                                (len([val for val in vj_reads[seg]['second'] if val in vj_reads[s_seg]['first']]) > 0):
                            if seg not in f_arr:
                                f_arr.append(seg)
                            if s_seg not in f_arr:
                                f_arr.append(s_seg)
                            v_seq = fasta_dict[s_seg]
                            j_seq = fasta_dict[seg]
                            len_seg = min(len(v_seq), len(j_seq))
                            if bases != -10:
                                if len_seg < bases:
                                    if bases > len(v_seq):
                                        sys.stdout.write(str(
                                            datetime.datetime.now()) + ' Bases parameter is bigger than the length of '
                                                                       'the V segment, taking the length' \
                                                                       'of the V/J segment instead, which is: ' + str(
                                            len_seg) + '\n')
                                        sys.stdout.flush()
                                    else:
                                        sys.stdout.write(str(
                                            datetime.datetime.now()) + ' Bases parameter is bigger than the length of '
                                                                       'the J segment, appending the C segment to the '
                                                                       'J segment\n')
                                        sys.stdout.flush()
                                        j_seq = j_seq + c_seq
                                        len_seg = bases
                                else:
                                    len_seg = bases
                            j_trim = j_seq[:len_seg]
                            v_trim = v_seq[-1 * len_seg:]
                            junc = v_trim + j_trim
                            record_name = s_seg + '.' + seg + '(' + id_name_dict[s_seg] + '-' + id_name_dict[seg] + ')'
                            record = SeqRecord(Seq(junc, IUPAC.ambiguous_dna), id=record_name, description='')
                            cur_cont = vj_counts_dict[seg] + vj_counts_dict[s_seg]
                            pair_count_dict[record] = cur_cont
    sorted_pairs = sorted(pair_count_dict.items(), key=operator.itemgetter(1), reverse=True)
    if ((top == -1) | (top > len(sorted_pairs))):
        for rec, count in sorted_pairs:
            SeqIO.write(rec, out, 'fasta')
    else:
        if not by_exp:
            for i in range(0, top):
                SeqIO.write(sorted_pairs[i][0], out, 'fasta')
        else:
            wrote = 1
            SeqIO.write(sorted_pairs[0][0], out, 'fasta')
            cur_count = sorted_pairs[0][1]
            wrote_second = False
            for i in range(1, len(sorted_pairs)):
                if sorted_pairs[i][1] == cur_count:
                    if not wrote_second:
                        wrote_second = True
                        SeqIO.write(sorted_pairs[i][0], out, 'fasta')
                        wrote += 1
                else:
                    cur_count = sorted_pairs[i][1]
                    wrote_second = False
                    SeqIO.write(sorted_pairs[i][0], out, 'fasta')
                    wrote += 1
                if wrote == top:
                    break

    out.close()
    return f_arr


# Load all the reads from a list of segments into a dictionary.
# INPUT: segsDict: A dict, where the key is a segment, and the value is an array of bed entries of this segment
#        mappedFile: Bam file of the mapped reaeds
#       readDict: A dictionary. The keys are the segment name, the values is a dictionary 'first':[] and 'second':[]
#                 where 'first' array holds the query name of R1's that overlap this segment, and 'second' holds the
#                 query name of R2's that overlap the segment.
def load_reads_to_dict(segs_dict, mapped_file, read_dict, read_overlap):
    count_dict = dict()
    for seg in segs_dict:
        l_arr = seg.strip('\n').split('\t')
        seg_name = l_arr[3]
        read_dict[seg_name] = {'first': [], 'second': []}
        l_arr = seg.strip('\n').split('\t')
        chr = l_arr[0]
        start = int(l_arr[1])
        end = int(l_arr[2])
        reads_iter = mapped_file.fetch(chr, start - 1, end + 1)
        read_counter = 0
        for read in reads_iter:
            overlap = read.get_overlap(start - 1, end + 1)
            if (end - start) < read_overlap:
                read_overlap = end - start - 15
            if overlap >= read_overlap:
                curr_name = read.query_name
                if read.is_read1:
                    if curr_name not in read_dict[seg_name]['first']:
                        read_dict[seg_name]['first'].append(curr_name)
                        read_counter += 1
                elif read.is_read2:
                    if curr_name not in read_dict[seg_name]['second']:
                        read_dict[seg_name]['second'].append(curr_name)
                        read_counter += 1
        count_dict[seg_name] = read_counter
    return (read_dict, count_dict)


# Similar to "writeJunctionsWithC, only that instead of looking for V-J paired-reads, it looks for
# V-C and J-C paired-reads
# INPUT:
#       vjReads - reads dict of the V and J segments created by loadReadsToDict
#       outName - output name for junction file
#       bases - number of bases to take from V and J for the junction
#       fastaDict
#       idNameDict
#       cReads - reads dict of the C segments created by loadReadsToDict
# OUTPUT:
#        fArr - the V and J segments for which we found a junction for
def write_junctions_with_c(vj_reads, out_name, bases, fasta_dict, id_name_dict, c_reads):
    out = open(out_name, 'w')
    f_arr = []
    v_arr = []
    j_arr = []
    for seg in vj_reads:
        if len(vj_reads[seg]) > 0:
            for c_seg in c_reads:
                if (len(c_reads[c_seg]) > 0):
                    if (len([val for val in vj_reads[seg]['first'] if val in c_reads[c_seg]['second']]) > 0) | \
                            (len([val for val in vj_reads[seg]['second'] if val in c_reads[c_seg]['first']]) > 0):
                        if id_name_dict[seg].find('J') != -1:
                            if seg not in j_arr:
                                j_arr.append(seg)
                        elif id_name_dict[seg].find('V') != -1:
                            if seg not in v_arr:
                                v_arr.append(seg)
                        f_arr.append(seg)
    for v_seg in v_arr:
        for j_seg in j_arr:
            v_seq_fa = fasta_dict[v_seg]
            j_seq_fa = fasta_dict[j_seg]
            len_seg = min(len(v_seq_fa), len(j_seq_fa))
            if bases != -10:
                if len_seg < bases:
                    sys.stdout.write(str(
                        datetime.datetime.now()) + ' Bases parameter is bigger than the length of the V or J segment, '
                                                   'taking the length' \
                                                   'of the V/J segment instead, which is: ' + str(len_seg) + '\n')
                    sys.stdout.flush()
                else:
                    len_seg = bases
            j_trim = j_seq_fa[:len_seg]
            v_trim = v_seq_fa[-1 * len_seg:]
            junc = v_trim + j_trim
            record_name = v_seg + '.' + j_seg + '(' + id_name_dict[v_seg] + '-' + id_name_dict[j_seg] + ')'
            record = SeqRecord(Seq(junc, IUPAC.ambiguous_dna), id=record_name, description='')
            SeqIO.write(record, out, 'fasta')
    out.close()
    return f_arr


def write_reads_file(bam, unmapped, junction_segs, output, vdj_dict, chain, strand, low_q):
    if chain == 'A':
        vdj_chain_dict = vdj_dict['Alpha']
        out_reads = output + '.alpha.mapped.and.unmapped.fa'
        paired_reads_1 = output + '.alpha.R1.fa'
        paired_reads_2 = output + '.alpha.R2.fa'
    elif chain == 'B':
        vdj_chain_dict = vdj_dict['Beta']
        out_reads = output + '.beta.mapped.and.unmapped.fa'
        paired_reads_1 = output + '.beta.R1.fa'
        paired_reads_2 = output + '.beta.R2.fa'
    else:
        sys.stderr.write(
            str(datetime.datetime.now()) + ' Error! chain parameter for function analyzeChain can only be A or B\n')
        sys.stderr.flush()
    out = open(out_reads, 'w')
    const_dict = vdj_chain_dict['C']
    # This dict for unmapped reads has reads that should be rev.comp in the revcomp arr, otherwise in id.
    # For every read, the value is a tuple - the first value is first/second, to make sure there are no errors.
    # The second value is id/revcomp, to see what sequence should be written.
    # Note: This classification is about what should happen to the unmapped reads, not how their paired maaped
    # reads were read.
    unmapped_dict = dict()
    seq_dict = dict()
    aligned_dict = dict()
    mapped_pairs_dict = dict()
    low_q_dict = dict()
    for seg in const_dict:
        (unmapped_dict, aligned_dict, seq_dict, mapped_pairs_dict, low_q_dict) = add_reads_to_dict(unmapped_dict, seg,
                                                                                                   bam, out,
                                                                                                   False, aligned_dict,
                                                                                                   seq_dict,
                                                                                                   strand, 'C',
                                                                                                   mapped_pairs_dict,
                                                                                                   low_q_dict)
    v_segs = vdj_chain_dict['V']
    for v_seg in v_segs:
        v_seg_name = v_seg.strip('\n').split('\t')[3]
        if v_seg_name in junction_segs:
            (unmapped_dict, aligned_dict, seq_dict, mapped_pairs_dict, low_q_dict) = add_reads_to_dict(unmapped_dict,
                                                                                                       v_seg, bam,
                                                                                                       out, True,
                                                                                                       aligned_dict,
                                                                                                       seq_dict, strand,
                                                                                                       'V',
                                                                                                       mapped_pairs_dict,
                                                                                                       low_q_dict)
    jSegs = vdj_chain_dict['J']
    for jSeg in jSegs:
        jSegName = jSeg.strip('\n').split('\t')[3]
        if jSegName in junction_segs:
            (unmapped_dict, aligned_dict, seq_dict, mapped_pairs_dict, low_q_dict) = add_reads_to_dict(unmapped_dict,
                                                                                                       jSeg, bam,
                                                                                                       out, True,
                                                                                                       aligned_dict,
                                                                                                       seq_dict, strand,
                                                                                                       'J',
                                                                                                       mapped_pairs_dict,
                                                                                                       low_q_dict)
    un_dict = dict()
    (seq_dict, un_dict) = write_unmapped_reads(unmapped_dict, out, unmapped, seq_dict, un_dict, aligned_dict,
                                               low_q_dict, low_q)
    seq_dict = add_mapped_pairs_to_seq_dict(seq_dict, bam, out, low_q, aligned_dict)
    write_seq_dict(seq_dict, paired_reads_1, paired_reads_2)
    out.close()
    return un_dict


# Aligned dict - all the reads (with _1/_2) that were already written to the mapped.unmapped.fa file
def add_reads_to_dict(unmapped_dict, seg_bed, bam, out, mapped_read, aligned_dict, seq_dict, strand, seg_type,
                      mapped_pairs_dict,
                      low_q_dict):
    bed_arr = seg_bed.strip('\n').split('\t')
    chr = bed_arr[0]
    start = int(bed_arr[1])
    end = int(bed_arr[2])
    mapped_file = pysam.AlignmentFile(bam, "rb")
    reads_iter = mapped_file.fetch(chr, start - 1, end + 1)
    for read in reads_iter:
        curr_name = read.query_name
        if curr_name not in seq_dict:
            seq_dict[curr_name] = ['0', '1']
        if read.is_read1:
            pair_read = 'second'
            read_name = curr_name + '\\1'
            pair_name = curr_name + '\\2'
            seq_pos = 0
        elif read.is_read2:
            pair_read = 'first'
            read_name = curr_name + '\\2'
            pair_name = curr_name + '\\1'
            seq_pos = 1
        else:
            sys.stderr.write(str(datetime.datetime.now()) + ' Error! Read is not read1 and not read2\n')
            sys.stderr.flush()
        curr_seq = Seq(read.query_sequence, IUPAC.ambiguous_dna)
        if read.is_reverse:
            pair_or = 'id'
            read_strand = 'minus'
            curr_seq = curr_seq.reverse_complement()
        else:
            read_strand = 'plus'
            pair_or = 'rev'
        seq_dict[curr_name][seq_pos] = curr_seq
        if read.mate_is_unmapped:
            take_pair = to_take_pair(seg_type, strand, read_strand)
            if take_pair == False:
                low_q_dict[pair_name] = '1'
            # takePair = True
            # if takePair:
            if curr_name in unmapped_dict:
                if unmapped_dict[curr_name] != (pair_read, pair_or):
                    sys.stderr.write(str(
                        datetime.datetime.now()) + ' Error! Read %s has more than one unmppaed mate with differnet '
                                                   'strand/mate\n' % curr_name)
                    sys.stderr.flush()
            unmapped_dict[curr_name] = (pair_read, pair_or)
        if mapped_read == True:
            if read_name in aligned_dict:
                if aligned_dict[read_name] != read.query_sequence:
                    sys.stderr.write(str(
                        datetime.datetime.now()) + ' Error! Read %s has two instances but different seuqences\n' %
                                     read.query_name)
                    sys.stderr.flush()
            else:
                aligned_dict[read_name] = read.query_sequence
                record = SeqRecord(Seq(read.query_sequence, IUPAC.ambiguous_dna), id=read_name, description='')
                SeqIO.write(record, out, 'fasta')
    mapped_file.close()
    # print "Removed counter: " + str(cc)
    return (unmapped_dict, aligned_dict, seq_dict, mapped_pairs_dict, low_q_dict)


def write_seq_dict(seq_dict, r1, r2):
    r1f = open(r1, 'w')
    r2f = open(r2, 'w')
    for seq in seq_dict:
        if ((seq_dict[seq][0] != '0') & (seq_dict[seq][1] != '1')):
            seq1 = seq
            seq2 = seq
            rec1 = SeqRecord(seq_dict[seq][0], id=seq1, description='')
            rec2 = SeqRecord(seq_dict[seq][1], id=seq2, description='')
            SeqIO.write(rec1, r1f, 'fasta')
            SeqIO.write(rec2, r2f, 'fasta')
        else:
            sys.stderr.write(str(datetime.datetime.now()) + ' The read %s has only one mate found, ignoring it\n' % seq)
            sys.stderr.flush()
    r1f.close()
    r2f.close()


### For minus, add an unmapped pair if the current mate is: 1. V and Plus 2. J/C and minus
### For plus, exactly the opposite

def to_take_pair(seg_type, strand, read_strand):
    if strand == 'none':
        return True
    if ((read_strand != 'minus') & (read_strand != 'plus')):
        sys.stderr.write(str(datetime.datetime.now()) + ' Error! Read strand should be plus or minus only\n')
        sys.stderr.flush()
        return True
    if ((seg_type == 'C') | (seg_type == 'J')):
        if strand == 'minus':
            if read_strand == 'minus':
                return True
            else:
                return False
        else:
            if read_strand == 'minus':
                return False
            else:
                return True
    else:
        if strand == 'minus':
            if read_strand == 'minus':
                return False
            else:
                return True
        else:
            if read_strand == 'minus':
                return True
            else:
                return False
    return True


def write_unmapped_reads(unmapped_dict, out, unmapped, seq_dict, un_dict, aligned_dict, low_q_dict, low_q):
    f = pysam.AlignmentFile(unmapped, "rb")
    reads_iter = f.fetch(until_eof=True)
    for read in reads_iter:
        name = read.query_name
        if name in unmapped_dict:
            c_name = name
            un_dict_name = name
            (strand, ori) = unmapped_dict[name]
            if (((strand == 'first') & (read.is_read2)) | ((strand == 'second') & (read.is_read1))):
                sys.stderr.write(str(
                    datetime.datetime.now()) + ' Error! unmapped read is inconsistent regarding first/second read\n')
                sys.stderr.flush()
            else:
                if strand == 'first':
                    name += '\\1'
                    un_dict_name += '_1'
                else:
                    name += '\\2'
                    un_dict_name += '_2'
                if un_dict_name in un_dict:
                    sys.stderr.write(str(
                        datetime.datetime.now()) + ' Error! unmapped read %s appear twice in unmapped bam file\n' %
                                     c_name)
                    sys.stderr.flush()
                un_dict[un_dict_name] = '1'
                q_seq = Seq(read.query_sequence, IUPAC.ambiguous_dna)
                if ori == 'rev':
                    q_seq = q_seq.reverse_complement()
                if name in aligned_dict:
                    if aligned_dict[name] != str(q_seq):
                        sys.stderr.write(str(
                            datetime.datetime.now()) + ' Error! unmapped read %s appear twice in alignedDict with '
                                                       'differnet seqs\n' % name)
                        sys.stderr.flush()
                else:
                    if ((name not in low_q_dict) | ((name in low_q_dict) & (not low_q))):
                        aligned_dict[name] = str(q_seq)
                        record = SeqRecord(q_seq, id=name, description='')
                        SeqIO.write(record, out, 'fasta')
                if ori == 'rev':
                    q_seq = q_seq.reverse_complement()
            if c_name not in seq_dict:
                sys.stderr.write(str(
                    datetime.datetime.now()) + ' Error! unmapped read is in unmappedDict but not in seqDict %s\n' %
                                 read.query_name)
                sys.stderr.flush()
            else:
                if strand == 'first':
                    seq_dict[c_name][0] = q_seq
                else:
                    seq_dict[c_name][1] = q_seq
    f.close()
    return (seq_dict, un_dict)


def add_mapped_pairs_to_seq_dict(seq_dict, bam, out, low_q, aligned_dict):
    first_dict = dict()
    second_dict = dict()
    for name in seq_dict:
        if seq_dict[name][0] == '0':
            first_dict[name] = '1'
        if seq_dict[name][1] == '1':
            second_dict[name] = '1'
        if ((seq_dict[name][1] == '1') & (seq_dict[name][0] == '0')):
            sys.stderr.write(str(datetime.datetime.now()) + ' Error! empty record insdie seqDict\n')
            sys.stderr.flush()
    f = pysam.AlignmentFile(bam, "rb")
    reads_iter = f.fetch()
    for read in reads_iter:
        name = read.query_name
        pos = -1
        if ((read.query_name in first_dict) & (read.is_read1)):
            pos = 0
            check = '0'
        elif ((read.query_name in second_dict) & (read.is_read2)):
            pos = 1
            check = '1'
        if pos != -1:
            qSeq = Seq(read.query_sequence, IUPAC.ambiguous_dna)
            if read.is_read1:
                curr_name = name + '\\1'
            else:
                curr_name = name + '\\2'
            if low_q:
                if read.mate_is_reverse:
                    if read.is_reverse:
                        r_seq = qSeq.reverse_complement()
                    else:
                        r_seq = qSeq
                else:
                    if not read.is_reverse:
                        r_seq = qSeq.reverse_complement()
                    else:
                        r_seq = qSeq
                if curr_name not in aligned_dict:
                    aligned_dict[curr_name] = str(r_seq)
                    record = SeqRecord(r_seq, id=curr_name, description='')
                    SeqIO.write(record, out, 'fasta')
                else:
                    if str(aligned_dict[curr_name]) != str(r_seq):
                        sys.stderr.write(str(
                            datetime.datetime.now()) + ' Error! read %s has two different sequences in alignedDict\n'
                                         % curr_name)
                        sys.stderr.flush()
            if read.is_reverse:
                qSeq = qSeq.reverse_complement()
            if seq_dict[name][pos] != check:
                if str(qSeq) != str(seq_dict[name][pos]):
                    sys.stderr.write(str(
                        datetime.datetime.now()) + ' Error! read %s has two different mapped sequences not in the V/J '
                                                   'region\n' % name)
                    sys.stderr.flush()
            seq_dict[name][pos] = qSeq
    f.close()
    return seq_dict


def run_rsem(out_dir, rsem, bowtie2, full_tcr_file_alpha, full_tcr_file_beta, output, samtools):
    if samtools != '':
        if samtools[-1] != '/':
            rsem += '/'
    rsem_ind_dir = out_dir + 'rsem_ind'
    if os.path.exists(rsem_ind_dir) == False:
        os.makedirs(rsem_ind_dir)
    if rsem != '':
        if rsem[-1] != '/':
            rsem += '/'
    if bowtie2 != '':
        if bowtie2[-1] != '/':
            bowtie2 += '/'
    if os.path.exists(full_tcr_file_alpha):
        if bowtie2 != '':
            subprocess.call([rsem + 'rsem-prepare-reference', '--bowtie2', '--bowtie2-path', bowtie2,
                             '-q', full_tcr_file_alpha, rsem_ind_dir + '/VDJ.alpha.seq'])
            subprocess.call([rsem + 'rsem-calculate-expression', '--no-qualities', '-q',
                             '--bowtie2', '--bowtie2-path', bowtie2, '--bowtie2-mismatch-rate', '0.0', '--paired-end',
                             output + '.alpha.R1.fa',
                             output + '.alpha.R2.fa', rsem_ind_dir + '/VDJ.alpha.seq', output + '.alpha.rsem.out'])
        else:
            subprocess.call([rsem + 'rsem-prepare-reference', '--bowtie2',
                             '-q', full_tcr_file_alpha, rsem_ind_dir + '/VDJ.alpha.seq'])
            subprocess.call([rsem + 'rsem-calculate-expression', '--no-qualities', '-q',
                             '--bowtie2', '--bowtie2-mismatch-rate', '0.0', '--paired-end', output + '.alpha.R1.fa',
                             output + '.alpha.R2.fa', rsem_ind_dir + '/VDJ.alpha.seq', output + '.alpha.rsem.out'])
        unsorted_bam = output + '.alpha.rsem.out.transcript.bam'
        if not os.path.exists(unsorted_bam):
            print "RSEM did not produce any transcript alignment files for alpha chain, please check the -rsem " \
                  "parameter"
        else:
            sorted_bam = output + '.alpha.rsem.out.transcript.sorted.bam'
            if not os.path.exists(sorted_bam):
                subprocess.call([samtools + 'samtools', 'sort', '-o', sorted_bam, unsorted_bam])
                subprocess.call([samtools + 'samtools', 'index', sorted_bam])

    else:
        sys.stdout.write(
            str(datetime.datetime.now()) + " Did not reconstruct any alpha chains, not running RSEM on alpha\n")
        sys.stdout.flush()
    if os.path.exists(full_tcr_file_beta):
        if bowtie2 != '':
            subprocess.call([rsem + 'rsem-prepare-reference', '--bowtie2', '--bowtie2-path', bowtie2,
                             '-q', full_tcr_file_beta, rsem_ind_dir + '/VDJ.beta.seq'])
            subprocess.call([rsem + 'rsem-calculate-expression', '--no-qualities', '-q', '--bowtie2', '--bowtie2-path',
                             bowtie2, '--bowtie2-mismatch-rate', '0.0', '--paired-end', output + '.beta.R1.fa',
                             output + '.beta.R2.fa',
                             rsem_ind_dir + '/VDJ.beta.seq', output + '.beta.rsem.out'])
        else:
            subprocess.call([rsem + 'rsem-prepare-reference', '--bowtie2',
                             '-q', full_tcr_file_beta, rsem_ind_dir + '/VDJ.beta.seq'])
            subprocess.call([rsem + 'rsem-calculate-expression', '--no-qualities', '-q', '--bowtie2',
                             '--bowtie2-mismatch-rate', '0.0', '--paired-end', output + '.beta.R1.fa',
                             output + '.beta.R2.fa',
                             rsem_ind_dir + '/VDJ.beta.seq', output + '.beta.rsem.out'])
        unsorted_bam = output + '.beta.rsem.out.transcript.bam'
        if not os.path.exists(unsorted_bam):
            print "RSEM did not produce any transcript alignment files for beta chain, please check the -rsem parameter"
        else:
            sorted_bam = output + '.beta.rsem.out.transcript.sorted.bam'
            if not os.path.exists(sorted_bam):
                subprocess.call([samtools + 'samtools', 'sort', '-o', sorted_bam, unsorted_bam])
                subprocess.call([samtools + 'samtools', 'index', sorted_bam])
    else:
        sys.stdout.write(
            str(datetime.datetime.now()) + " Did not reconstruct any beta chains, not running RSEM on beta\n")
        sys.stdout.flush()
