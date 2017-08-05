"""
This modules contains functions used to process single-end cells.

"""

import sys
import os
import subprocess
import datetime
from Bio import SeqIO
import pysam
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC


def analyze_chain_single_end(fastq, trimmomatic, trans_ind, bowtie2, id_name_dict, output, fasta_dict, bases):
    mapped_reads_dict_alpha = dict()
    mapped_reads_dict_beta = dict()
    len_arr = [999, 50, 25]
    for currLen in len_arr:
        if currLen == 999:
            steps = ['none']
        else:
            steps = ['left', 'right']
        for side in steps:
            if side == 'left':
                crop = 'CROP:' + str(currLen)
            elif side == 'right':
                crop = 'HEADCROP:' + str(currLen)
            else:
                crop = ''
                # TODO: make sure we delete those files
            trim_fq = fastq + '.' + str(currLen) + '.' + str(side) + '.trimmed.fq'
            # TODO: use bowtie trimmer instead
            # TODO: make sure about minus strand alignment
            if crop == '':
                subprocess.call(
                    ['java', '-jar', trimmomatic, 'SE', '-phred33', fastq, trim_fq, 'LEADING:15', 'TRAILING:15',
                     'MINLEN:20'])
            else:
                subprocess.call(
                    ['java', '-jar', trimmomatic, 'SE', '-phred33', fastq, trim_fq, 'LEADING:15', 'TRAILING:15', crop,
                     'MINLEN:20'])
            sam_f = trim_fq + '.sam'
            if bowtie2 != '':
                if bowtie2.endswith('/'):
                    bowtie_call = bowtie2 + 'bowtie2'
                else:
                    bowtie_call = bowtie2 + '/bowtie2'
            else:
                bowtie_call = 'bowtie2'
            subprocess.call([bowtie_call, '-q --phred33  --score-min L,0,0', '-x', trans_ind, '-U', trim_fq, '-S', sam_f])
            if os.path.isfile(sam_f):
                mapped_reads_dict_alpha = find_reads_and_segments(sam_f, mapped_reads_dict_alpha, id_name_dict, 'A')
                mapped_reads_dict_beta = find_reads_and_segments(sam_f, mapped_reads_dict_beta, id_name_dict, 'B')
    alpha_out = output + '.alpha.junctions.txt'
    alpha_out_reads = output + '.alpha.mapped.and.unmapped.fa'
    beta_out_reads = output + '.beta.mapped.and.unmapped.fa'
    beta_out = output + '.beta.junctions.txt'
    write_junction_file_se(mapped_reads_dict_alpha, id_name_dict, alpha_out, fasta_dict, bases, 'alpha')
    write_junction_file_se(mapped_reads_dict_beta, id_name_dict, beta_out, fasta_dict, bases, 'beta')
    write_reads_file_se(mapped_reads_dict_alpha, alpha_out_reads, fastq)
    write_reads_file_se(mapped_reads_dict_beta, beta_out_reads, fastq)
    sys.exit(1)


def find_reads_and_segments(sam_f, mapped_reads_dict, id_name_dict, chain):
    sam_file = pysam.AlignmentFile(sam_f, 'r')
    reads_iter = sam_file.fetch(until_eof=True)
    for read in reads_iter:
        if read.is_unmapped == False:
            seg = sam_file.getrname(read.reference_id)
            if seg in id_name_dict:
                if id_name_dict[seg].find(chain) != -1:
                    read_name = read.query_name
                    if read_name not in mapped_reads_dict:
                        mapped_reads_dict[read_name] = []
                    if seg not in mapped_reads_dict[read_name]:
                        mapped_reads_dict[read_name].append(seg)
    sam_file.close()
    return mapped_reads_dict


def write_reads_file_se(mapped_reads_dict, out_reads, fastq):
    if fastq.endswith('.gz'):
        subprocess.call(['gunzip', fastq])
        new_fq = fastq.replace('.gz', '')
    else:
        new_fq = fastq
    out = open(out_reads, 'w')
    fqF = open(new_fq, 'rU')
    for record in SeqIO.parse(fqF, 'fastq'):
        if record.id in mapped_reads_dict:
            new_rec = SeqRecord(record.seq, id=record.id, description='')
            SeqIO.write(new_rec, out, 'fasta')
    out.close()
    fqF.close()
    if fastq.endswith('.gz'):
        subprocess.call(['gzip', new_fq])


def write_junction_file_se(mapped_reads_dict, id_name_dict, output, fasta_dict, bases, chain):
    out = open(output, 'w')
    v_segs = []
    j_segs = []
    c_segs = []
    for read in mapped_reads_dict:
        for seg in mapped_reads_dict[read]:
            if id_name_dict[seg].find('V') != -1:
                if seg not in v_segs:
                    v_segs.append(seg)
            elif id_name_dict[seg].find('J') != -1:
                if seg not in j_segs:
                    j_segs.append(seg)
            elif id_name_dict[seg].find('C') != -1:
                if seg not in c_segs:
                    c_segs.append(seg)
            else:
                print "Error! not V/J/C in fasta dict"
    if len(v_segs) == 0:
        print "Did not find any V segments for " + chain + " chain"
    else:
        if len(c_segs) == 0:
            print "Did not find any C segments for " + chain + " chain"
            c_segs = ['NA']
        if len(j_segs) == 0:
            print "Did not find any J segments for " + chain + " chain"
            j_segs = ['NA']
        for v_seg in v_segs:
            for j_seg in j_segs:
                for c_seg in c_segs:
                    add_segment_to_junction_file_se(v_seg, j_seg, c_seg, out, fasta_dict, bases, id_name_dict)
    out.close()


def add_segment_to_junction_file_se(v_seg, j_seg, c_seg, out, fasta_dict, bases, id_name_dict):
    v_seq = fasta_dict[v_seg]
    if j_seg != 'NA':
        j_name = id_name_dict[j_seg]
        j_seq = fasta_dict[j_seg]
    else:
        j_seq = ''
        j_name = 'NoJ'
    if c_seg != 'NA':
        c_name = id_name_dict[c_seg]
        c_seq = fasta_dict[c_seg]
    else:
        c_name = 'NoC'
        c_seq = ''
    jc_seq = j_seq + c_seq
    len_seg = min(len(v_seq), len(jc_seq))
    if bases != -10:
        if len_seg < bases:
            sys.stdout.write(str(
                datetime.datetime.now()) + ' Bases parameter is bigger than the length of the V or J segment, '
                                           'taking the length'
                                           'of the V/J segment instead, which is: ' + str(len_seg) + '\n')
            sys.stdout.flush()
        else:
            len_seg = bases
    j_trim = jc_seq[:len_seg]
    v_trim = v_seq[-1 * len_seg:]
    junc = v_trim + j_trim
    record_name = v_seg + '.' + j_seg + '.' + c_seg + '(' + id_name_dict[v_seg] + '-' + j_name + '-' + c_name + ')'
    record = SeqRecord(Seq(junc, IUPAC.ambiguous_dna), id=record_name, description='')
    SeqIO.write(record, out, 'fasta')
