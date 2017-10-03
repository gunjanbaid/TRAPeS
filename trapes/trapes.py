#!/usr/bin/env python

"""
This module contains functions to perform argument parsing, argument checking, and processing of each cell in the
input data.

"""

import sys
import os
import argparse
import datetime
import utils
import write_output_files
import process_single_cell


def check_parameters(genome, strand, single_cell, path, sum_f):
    if ((genome != 'hg38') & (genome != 'mm10') & (genome != 'hg19') & (genome != 'mm10_ncbi')):
        sys.exit("-genome only accept one of the following: mm10, mm10_ncbi, hg38, hg19")
    if strand.lower() not in ['none', 'minus', 'plus']:
        sys.exit("-strand should be one of: none, minus, plus")
    if not single_cell:
        if path == '':
            sys.exit("when running on multiple cells you must include the -path parameter")
        if sum_f == '':
            sys.exit("when running on multiple cells you must include the -sumF parameter")
        if not os.path.isdir(path):
            sys.exit("%s path does not exists. Please check your -path parameter and run again" % path)


# def runTCRpipe(fasta, bed, output, bam, unmapped, mapping, bases, strand, reconstruction, aaF , numIterations,
# thresholdScore, minOverlap,
# rsem, bowtie2, singleCell, path, subpath, sumF, lowQ, singleEnd, fastq, trimmomatic, transInd):
def run_tcr_pipe(genome, output, bam, unmapped, bases, strand, num_iterations, threshold_score, min_overlap, rsem,
                 bowtie2, single_cell, path, sum_f, low_q, samtools, top, by_exp, read_overlap, one_side):
    check_parameters(genome, strand, single_cell, path, sum_f)
    if single_cell == True:
        # TODO: Fix this, won't work for SE
        # runSingleCell(fasta, bed, output, bam, unmapped, mapping, bases, strand, reconstruction, aaF ,
        # numIterations, thresholdScore, minOverlap,
        #          rsem, bowtie2, lowQ, singleEnd, fastq, trimmomatic, transInd)
        sys.exit(0)
    if path == './':
        path = os.getcwd()
    if not path.endswith('/'):
        path = path + '/'
    final_stat_dict = dict()
    tcr_fout = open(sum_f + '.TCRs.txt', 'w')
    opened = False
    for cell_folder in os.listdir(path):
        full_path = path + cell_folder + '/'
        if ((os.path.exists(full_path)) & (os.path.isdir(full_path))):
            sys.stdout.write(str(datetime.datetime.now()) + " Working on: " + cell_folder + '\n')
            sys.stdout.flush()
            (found, nbam, nunmapped, noutput) = utils.format_files(full_path, bam, unmapped, output)
            if not found:
                sys.stderr.write(str(datetime.datetime.now()) + " There is not a bam or unmapped file in "
                                                                "this folder, moving to the next folder\n")
                sys.stderr.flush()
            else:
                DATA_PATH = pkg_resources.resource_filename('trapes', 'data/')
                reconstruction = pkg_resources.resource_filename('trapes', 'vdj.alignment')
                if genome == 'hg38':
                    fasta = DATA_PATH + 'hg38/hg38.TCR.fa'
                    bed = DATA_PATH + 'hg38/hg38.TCR.bed'
                    mapping = DATA_PATH + 'hg38/hg38.id.name.mapping.TCR.txt'
                    aa_f = DATA_PATH + 'hg38/hg38.TCR.conserved.AA.txt'
                    ref_ind = DATA_PATH + 'hg38/index/hg38'
                if genome == 'mm10':
                    fasta = DATA_PATH + 'mm10/mm10.TCR.fa'
                    bed = DATA_PATH + 'mm10/mm10.TCR.bed'
                    mapping = DATA_PATH + 'mm10/mm10.gene.id.mapping.TCR.txt'
                    aa_f = DATA_PATH + 'mm10/mm10.conserved.AA.txt'
                    ref_ind = DATA_PATH + 'mm10/index/mm10'
                if genome == 'mm10_ncbi':
                    fasta = DATA_PATH + 'mm10_ncbi/mm10.TCR.fa'
                    bed = DATA_PATH + 'mm10_ncbi/mm10.TCR.bed'
                    mapping = DATA_PATH + 'mm10_ncbi/mm10.gene.id.mapping.TCR.txt'
                    aa_f = DATA_PATH + 'mm10_ncbi/mm10.conserved.AA.txt'
                    ref_ind = DATA_PATH + 'mm10_ncbi/index/mm10'
                if genome == 'hg19':
                    fasta = DATA_PATH + 'hg19/hg19.TCR.fa'
                    bed = DATA_PATH + 'hg19/hg19.TCR.bed'
                    mapping = DATA_PATH + 'hg19/hg19.gene.id.mapping.TCR.txt'
                    aa_f = DATA_PATH + 'hg19/hg19.conserved.AA.txt'
                    ref_ind = DATA_PATH + 'hg19/index/hg19'

                process_single_cell.run_single_cell(fasta, bed, noutput, nbam, nunmapped, mapping, bases, strand,
                                                    reconstruction, aa_f,
                                                    num_iterations, threshold_score,
                                                    min_overlap, rsem, bowtie2, low_q, samtools, top, by_exp,
                                                    read_overlap, one_side)
                opened = write_output_files.add_cell_to_tcr_sum(cell_folder, noutput, opened, tcr_fout)
                final_stat_dict = write_output_files.add_to_stat_dict(noutput, cell_folder, final_stat_dict)
    sum_f_out = open(sum_f + '.summary.txt', 'w')
    sum_f_out.write(']ample\talpha\tbeta\n')
    for cell in sorted(final_stat_dict):
        f_out = cell + '\t' + final_stat_dict[cell]['alpha'] + '\t' + final_stat_dict[cell]['beta'] + '\n'
        sum_f_out.write(f_out)
    sum_f_out.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-genome', '-g', '-G', help='Alignment genome. Currently supported: mm10, mm10_ncbi and hg38',
                        required=True)
    parser.add_argument('-single_cell', help='add if you are only running on a single cell. If so,'
                                             'it will ignore -path and -subpath arguments', action='store_true')
    parser.add_argument('-low_q', help='add if you want to add \"low quality\" reads as input to the reconstruction '
                                       'algorithm', action='store_true')
    parser.add_argument('-one_side', help='add if you want to observe reconstrctuion only from the V side',
                        action='store_true')
    parser.add_argument('-path', '-p', '-P', help='The path for the data directory. Assumes that every subdirectory'
                                                  'is a single cell', default='')
    parser.add_argument('-sum_f', help='prefix for summary outputs', default='')
    parser.add_argument('-bowtie2', '-bw', '-BW', help='Path to bowtie2. If not used assumes that bowtie2 is in the'
                                                       'default path', default='')
    parser.add_argument('-rsem', '-RSEM', help='Path to rsem. If not used assumes that rsem is in the'
                                               'default path', default='')
    parser.add_argument('-strand',
                        help='Strand of the right most read in genomic coordinates. Options are: [minus, plus, '
                             'none]. Defualt is minus', default='minus')
    parser.add_argument('-output', '-out', '-o', '-O', help='output prefix, relative to /path/singleCellFolder',
                        required=True)
    parser.add_argument('-bam',
                        help='Input bam alignment file, relative to /path/singleCellFolder/ if working on multiple '
                             'files',
                        default='./tophat_output/picard_output/sorted.bam')
    parser.add_argument('-unmapped', '-u', '-U',
                        help='bam file of the unmapped reads, relative to /path/singleCellFolder/',
                        default='./tophat_output/unmapped.bam')
    parser.add_argument('-bases', '-b', '-B',
                        help='Number of bases to take from each V and J segments, default is min(len(V), len(J) ',
                        type=int, default=-10)
    parser.add_argument('-iterations', '-iter', '-i', '-I', help='Number of iterations for the reconstruction'
                                                                 'algorithm, default is 20', type=int, default=20)
    parser.add_argument('-samtools', help='Path to samtools. If not used assumes that samtools is in the default path',
                        default='')
    parser.add_argument('-score', '-sc', '-SC', help='Alignment score threshold. Default is 15', type=int, default=15)
    parser.add_argument('-top', '-t', '-T', help='Take only the top x combination of V and J, based on the sum '
                                                 'number of reads that map to both. Default is to take all', type=int,
                        default=-1)
    parser.add_argument('-read_overlap', '-ro', '-readoverlap',
                        help='Add a read to list of mapped reads only if it maps at least X bases'
                             'to the V/J/C segment. Default is 1', type=int, default=1)
    parser.add_argument('-by_exp',
                        help='if using the Top option, add this tag if you want to take only two chains from each' \
                             'read count, until top is reached', action='store_true')

    parser.add_argument('-overlap', '-ol', '-OL', help='Number of minimum bases that overlaps V and J ends,'
                                                       'default is 10', type=int, default=10)
    args = parser.parse_args()
    run_tcr_pipe(args.genome, args.output, args.bam, args.unmapped, args.bases, args.strand,
                 args.iterations, args.score, args.overlap, args.rsem, args.bowtie2,
                 args.single_cell, args.path, args.sum_f, args.low_q, args.samtools, args.top, args.by_exp,
                 args.read_overlap,
                 args.one_side)


if __name__ == '__main__':
    main()
