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


def analyze_chain_single_end(fastaDict, vdjDict, output, bam, unmapped, idNameDict, bases, strand,
                          lowQ, bowtie2, refInd, trim):
    if bowtie2 != '':
        if bowtie2.endswith('/'):
            bowtieCall = bowtie2 + 'bowtie2'
        else:
            bowtieCall = bowtie2 + '/bowtie2'
    else:
        bowtieCall = 'bowtie2'
     
    fastq1 = output + ".unmapped.fq"
    fastq2 = output + ".sorted.fq"
    
    if not os.path.isfile(fastq1):
        fastq1_file = open(fastq1, "w")
        subprocess.call(["samtools", "bam2fq", unmapped], stdout=fastq1_file)
        fastq1_file.close()

    if not os.path.isfile(fastq2):
        fastq2_file = open(fastq2, "w")
        subprocess.call(["samtools", "bam2fq", bam], stdout=fastq2_file)
        fastq2_file.close()

    mappedReadsDictAlpha = dict()
    mappedReadsDictBeta = dict()
    alphaOut = output + '.alpha.junctions.txt'
    alphaOutReads = output + '.alpha.mapped.and.unmapped.fa'
    betaOutReads = output + '.beta.mapped.and.unmapped.fa'
    betaOut = output + '.beta.junctions.txt'

    temp1 = output + '.temp1.sam'
    temp2 = output + '.temp2.sam'
    temp3 = output + '.temp3.sam'
    bam_trimmed = output + '.trimmed.bam'
    bam_filtered = output + '.filtered.bam'
    bam_sorted = output + '.sorted.bam'

    subprocess.call([bowtieCall ,'-q --phred33', '-x', refInd, '-U', fastq2, '-S', temp1, "--trim3", str(trim)])
    subprocess.call([bowtieCall ,'-q --phred33', '-x', refInd, '-U', fastq2, '-S', temp2, "--trim5", str(trim)])
    subprocess.call([bowtieCall ,'-q --phred33', '-x', refInd, '-U', fastq2, '-S', temp3, "--trim5", str(trim//2), "--trim3", str(trim//2)])
    subprocess.call(["samtools", "merge", bam_trimmed, temp1, temp2, temp3])

    bam_filtered_file = open(bam_filtered, "w")
    subprocess.call(["samtools", "view", "-b", "-F", "4", bam_trimmed], stdout=bam_filtered_file)
    bam_filtered_file.close()

    bam_sorted_file = open(bam_sorted, "w")
    subprocess.call(["samtools", "sort", bam_filtered], stdout=bam_sorted_file)
    bam_sorted_file.close()

    subprocess.call(["samtools", "index", bam_sorted])

    mappedReadsDictAlpha = findReadsAndSegments(bam_sorted, mappedReadsDictAlpha, idNameDict, 'A')
    mappedReadsDictBeta = findReadsAndSegments(bam_sorted, mappedReadsDictBeta, idNameDict, 'B')
    vSegsA, jSegsA, cSegsA = write_junction_file_se(mappedReadsDictAlpha, idNameDict, alphaOut, fastaDict, bases, 'alpha')
    vSegsB, jSegsB, cSegsB = write_junction_file_se(mappedReadsDictBeta, idNameDict, betaOut, fastaDict, bases, 'beta')

    allSegsA = vSegsA + jSegsA + cSegsA
    allSegsB = vSegsB + jSegsB + cSegsB

    mappedReadsDictAlpha = {}
    mappedReadsDictBeta = {}
    mappedReadsDictAlpha = findReadsAndSegments2(bam_sorted, mappedReadsDictAlpha, idNameDict, 
        'A', allSegsA)
    mappedReadsDictBeta = findReadsAndSegments2(bam_sorted, mappedReadsDictBeta, idNameDict, 
        'B', allSegsB)

    remove_file(temp1, temp2, temp3, bam_trimmed, bam_filtered, bam_sorted)
     
    subprocess.call([bowtieCall ,'-q --phred33  --score-min L,0,0', '-x', refInd, '-U', fastq1, 
        '-S', temp1, "--trim3", str(trim)])
    subprocess.call([bowtieCall ,'-q --phred33  --score-min L,0,0', '-x', refInd, '-U', fastq1, 
        '-S', temp2, "--trim5", str(trim)])
    subprocess.call([bowtieCall ,'-q --phred33  --score-min L,0,0', '-x', refInd, '-U', fastq1, 
        '-S', temp3, "--trim5", str(trim//2), "--trim3", str(trim//2)])    
    subprocess.call(["samtools", "merge", bam_trimmed, temp1, temp2, temp3])
    remove_file(temp1, temp2, temp3)

    bam_filtered_file = open(bam_filtered, "w")
    subprocess.call(["samtools", "view", "-b", "-F", "4", bam_trimmed], stdout=bam_filtered_file)
    bam_filtered_file.close()
 
    bam_sorted_file = open(bam_sorted, "w")
    subprocess.call(["samtools", "sort", bam_filtered], stdout=bam_sorted_file)
    bam_sorted_file.close()

    subprocess.call(["samtools", "index", bam_sorted])
    mappedReadsDictAlpha = findReadsAndSegments2(bam_sorted, mappedReadsDictAlpha, idNameDict, 
        'A', allSegsA)
    mappedReadsDictBeta = findReadsAndSegments2(bam_sorted, mappedReadsDictBeta, idNameDict, 
        'B', allSegsB)

    writeReadsFileSE(mappedReadsDictAlpha, alphaOutReads, fastq1, fastq2)
    writeReadsFileSE(mappedReadsDictBeta, betaOutReads, fastq1, fastq2)


def find_reads_and_segments(bam_f, mapped_reads_dict, id_name_dict, chain):
    sam_file = pysam.AlignmentFile(bam_f, 'rb')
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


def write_reads_file_se(mapped_reads_dict, out_reads, fastq, fastq_2):
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
    fqF.close()
    fqF2 = open(fastq_2, 'rU')
    for record in SeqIO.parse(fqF2, 'fastq'):
        if record.id in mapped_reads_dict:
            new_rec = SeqRecord(record.seq, id = record.id, description = '')
            SeqIO.write(new_rec, out,'fasta')
    fqF2.close()
    out.close()
    if fastq.endswith('.gz'):
        subprocess.call(['gzip', new_fq])


def write_junction_file_se(mapped_reads_dict, id_name_dict, output, fasta_dict, bases, chain):
    out = open(output, 'w')
    v_segs = {}
    j_segs = {}
    c_segs = {}
    for read in mapped_reads_dict:
        for seg in mapped_reads_dict[read]:
            if id_name_dict[seg].find('V') != -1:
                if seg not in v_segs:
                    v_segs[seg] = 1
                else:
                    v_segs[seg] = v_segs[seg] + 1
            elif id_name_dict[seg].find('J') != -1:
                if seg not in j_segs:
                    j_segs[seg] = 1
                else:
                    j_segs[seg] = j_segs[seg] + 1
            elif id_name_dict[seg].find('C') != -1:
                if seg not in c_segs:
                    c_segs[seg] = 1
                else:
                    c_segs[seg] = c_segs[seg] + 1                 
            else:
                print "Error! not V/J/C in fasta dict"

    v_segs = sorted(v_segs.items(), key=lambda x: -x[1])[:2]
    v_segs = [key for key,val in v_segs]
    j_segs = sorted(j_segs.items(), key=lambda x: -x[1])[:5]
    j_segs = [key for key,val in j_segs]
    c_segs = sorted(c_segs.items(), key=lambda x: -x[1])[:1]
    c_segs = [key for key,val in c_segs]

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


def write_unmapped_reads_to_dict_se(unmapped):
    un_dict = {}
    f = pysam.AlignmentFile(unmapped, "rb")
    readsIter = f.fetch(until_eof = True)
    for read in readsIter:
        name = read.query_name
        un_dict[name] = '1'
    return un_dict


def runRsemSE(outDir, rsem, bowtie2, fullTcrFileAlpha, fullTcrFileBeta, output, samtools):
    if samtools != '':
        if samtools[-1] != '/':
            rsem += '/'
    rsemIndDir = outDir + 'rsem_ind'
    if os.path.exists(rsemIndDir) == False:
        os.makedirs(rsemIndDir)
    if rsem != '':
        if rsem[-1] != '/':
            rsem += '/'
    if bowtie2 != '':
        if bowtie2[-1] != '/':
            bowtie2 += '/'
    if os.path.exists(fullTcrFileAlpha):
        if bowtie2 != '':
            subprocess.call([rsem + 'rsem-prepare-reference' , '--bowtie2', '--bowtie2-path', bowtie2 ,
                             '-q', fullTcrFileAlpha, rsemIndDir + '/VDJ.alpha.seq'])
            subprocess.call([rsem + 'rsem-calculate-expression', '--no-qualities', '-q',
                             '--bowtie2', '--bowtie2-path',bowtie2, '--bowtie2-mismatch-rate', '0.0' , output + '.alpha.mapped.and.unmapped.fa',
                             rsemIndDir + '/VDJ.alpha.seq', output + '.alpha.rsem.out'])
        else:
            subprocess.call([rsem + 'rsem-prepare-reference' , '--bowtie2',
                             '-q', fullTcrFileAlpha, rsemIndDir + '/VDJ.alpha.seq'])
            subprocess.call([rsem + 'rsem-calculate-expression', '--no-qualities', '-q',
                             '--bowtie2', '--bowtie2-mismatch-rate', '0.0', output + '.alpha.mapped.and.unmapped.fa',
                             rsemIndDir + '/VDJ.alpha.seq', output + '.alpha.rsem.out'])
        unsortedBam = output + '.alpha.rsem.out.transcript.bam'
        if not os.path.exists(unsortedBam):
            print "RSEM did not produce any transcript alignment files for alpha chain, please check the -rsem parameter"
        else:
            sortedBam = output + '.alpha.rsem.out.transcript.sorted.bam'
            if not os.path.exists(sortedBam):
                subprocess.call([samtools + 'samtools', 'sort','-o',sortedBam, unsortedBam])
                subprocess.call([samtools + 'samtools', 'index', sortedBam])

    else:
        sys.stdout.write(str(datetime.datetime.now()) + " Did not reconstruct any alpha chains, not running RSEM on alpha\n")
        sys.stdout.flush()
    if os.path.exists(fullTcrFileBeta):
        if bowtie2 != '':
            subprocess.call([rsem + 'rsem-prepare-reference' , '--bowtie2', '--bowtie2-path', bowtie2 ,
                             '-q', fullTcrFileBeta, rsemIndDir + '/VDJ.beta.seq'])
            subprocess.call([rsem + 'rsem-calculate-expression', '--no-qualities', '-q', '--bowtie2', '--bowtie2-path',
                             bowtie2, '--bowtie2-mismatch-rate', '0.0', output + '.beta.mapped.and.unmapped.fa',
                             rsemIndDir + '/VDJ.beta.seq', output + '.beta.rsem.out'])
        else:
            subprocess.call([rsem + 'rsem-prepare-reference' , '--bowtie2',
                             '-q', fullTcrFileBeta, rsemIndDir + '/VDJ.beta.seq'])
            subprocess.call([rsem + 'rsem-calculate-expression', '--no-qualities', '-q', '--bowtie2',
                              '--bowtie2-mismatch-rate', '0.0', output + '.beta.mapped.and.unmapped.fa',
                             rsemIndDir + '/VDJ.beta.seq', output + '.beta.rsem.out'])
        unsortedBam = output + '.beta.rsem.out.transcript.bam'
        if not os.path.exists(unsortedBam):
            print "RSEM did not produce any transcript alignment files for beta chain, please check the -rsem parameter"
        else:
            sortedBam = output + '.beta.rsem.out.transcript.sorted.bam'
            if not os.path.exists(sortedBam):
                subprocess.call([samtools + 'samtools', 'sort','-o',sortedBam, unsortedBam])
                subprocess.call([samtools + 'samtools', 'index', sortedBam])
    else:
        sys.stdout.write(str(datetime.datetime.now()) + " Did not reconstruct any beta chains, not running RSEM on beta\n")
        sys.stdout.flush()

