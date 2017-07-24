#!/usr/bin/env python

"""
This module contains functions to perform argument parsing, argument checking, and processing of each cell in the
input data.

"""

import sys
import os
import argparse
import datetime

from utils import formatFiles
from write_output_files import addCellToTCRsum, addToStatDict
from process_single_cell import runSingleCell


def checkParameters(genome, strand, singleCell, path, sumF):
    if ((genome != 'hg38') & (genome != 'mm10') & (genome != 'hg19') & (genome != 'mm10_ncbi')):
        sys.exit("-genome only accept one of the following: mm10, mm10_ncbi, hg38, hg19")
    if strand.lower() not in ['none', 'minus', 'plus']:
        sys.exit("-strand should be one of: none, minus, plus")
    if not singleCell:
        if path == '':
            sys.exit("when running on multiple cells you must include the -path parameter")
        if sumF == '':
            sys.exit("when running on multiple cells you must include the -sumF parameter")
        if not os.path.isdir(path):
            sys.exit("%s path does not exists. Please check your -path parameter and run again" % path)


# def runTCRpipe(fasta, bed, output, bam, unmapped, mapping, bases, strand, reconstruction, aaF , numIterations,
# thresholdScore, minOverlap,
# rsem, bowtie2, singleCell, path, subpath, sumF, lowQ, singleEnd, fastq, trimmomatic, transInd):
def runTCRpipe(genome, output, bam, unmapped, bases, strand, numIterations, thresholdScore, minOverlap, rsem, bowtie2,
               singleCell, path, sumF, lowQ, samtools, top, byExp, readOverlap, oneSide):
    checkParameters(genome, strand, singleCell, path, sumF)
    if singleCell == True:
        # TODO: Fix this, won't work for SE
        # runSingleCell(fasta, bed, output, bam, unmapped, mapping, bases, strand, reconstruction, aaF ,
        # numIterations, thresholdScore, minOverlap,
        #          rsem, bowtie2, lowQ, singleEnd, fastq, trimmomatic, transInd)
        sys.exit(0)
    if path == './':
        path = os.getcwd()
    if not path.endswith('/'):
        path = path + '/'
    finalStatDict = dict()
    tcrFout = open(sumF + '.TCRs.txt', 'w')
    opened = False
    for cellFolder in os.listdir(path):
        fullPath = path + cellFolder + '/'
        if ((os.path.exists(fullPath)) & (os.path.isdir(fullPath))):
            sys.stdout.write(str(datetime.datetime.now()) + " Working on: " + cellFolder + '\n')
            sys.stdout.flush()
            (found, nbam, nunmapped, noutput) = formatFiles(fullPath, bam, unmapped, output)
            if not found:
                sys.stderr.write(str(datetime.datetime.now()) + " There is not a bam or unmapped file in "
                                                                "this folder, moving to the next folder\n")
                sys.stderr.flush()
            else:
                currFolder = os.path.abspath(os.path.dirname(sys.argv[0])) + '/'
                reconstruction = currFolder + '/vdj.alignment'
                if genome == 'hg38':
                    fasta = currFolder + 'Data/hg38/hg38.TCR.fa'
                    bed = currFolder + 'Data/hg38/hg38.TCR.bed'
                    mapping = currFolder + 'Data/hg38/hg38.id.name.mapping.TCR.txt'
                    aaF = currFolder + 'Data/hg38/hg38.TCR.conserved.AA.txt'
                if genome == 'mm10':
                    fasta = currFolder + 'Data/mm10/mm10.TCR.fa'
                    bed = currFolder + 'Data/mm10/mm10.TCR.bed'
                    mapping = currFolder + 'Data/mm10/mm10.gene.id.mapping.TCR.txt'
                    aaF = currFolder + 'Data/mm10/mm10.conserved.AA.txt'
                if genome == 'mm10_ncbi':
                    fasta = currFolder + 'Data/mm10_ncbi/mm10.TCR.fa'
                    bed = currFolder + 'Data/mm10_ncbi/mm10.TCR.bed'
                    mapping = currFolder + 'Data/mm10_ncbi/mm10.gene.id.mapping.TCR.txt'
                    aaF = currFolder + 'Data/mm10_ncbi/mm10.conserved.AA.txt'
                if genome == 'hg19':
                    fasta = currFolder + 'Data/hg19/hg19.TCR.fa'
                    bed = currFolder + 'Data/hg19/hg19.TCR.bed'
                    mapping = currFolder + 'Data/hg19/hg19.gene.id.mapping.TCR.txt'
                    aaF = currFolder + 'Data/hg19/hg19.conserved.AA.txt'

                runSingleCell(fasta, bed, noutput, nbam, nunmapped, mapping, bases, strand, reconstruction, aaF,
                              numIterations, thresholdScore,
                              minOverlap, rsem, bowtie2, lowQ, samtools, top, byExp, readOverlap, oneSide)
                opened = addCellToTCRsum(cellFolder, noutput, opened, tcrFout)
                finalStatDict = addToStatDict(noutput, cellFolder, finalStatDict)
    sumFout = open(sumF + '.summary.txt', 'w')
    sumFout.write('sample\talpha\tbeta\n')
    for cell in sorted(finalStatDict):
        fout = cell + '\t' + finalStatDict[cell]['alpha'] + '\t' + finalStatDict[cell]['beta'] + '\n'
        sumFout.write(fout)
    sumFout.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-genome', '-g', '-G', help='Alignment genome. Currently supported: mm10, mm10_ncbi and hg38',
                        required=True)
    parser.add_argument('-singleCell', help='add if you are only running on a single cell. If so,'
                                            'it will ignore -path and -subpath arguments', action='store_true')
    parser.add_argument('-lowQ', help='add if you want to add \"low quality\" reads as input to the reconstruction '
                                      'algorithm', action='store_true')
    parser.add_argument('-oneSide', help='add if you want to observe reconstrctuion only from the V side',
                        action='store_true')
    parser.add_argument('-path', '-p', '-P', help='The path for the data directory. Assumes that every subdirectory'
                                                  'is a single cell', default='')
    parser.add_argument('-sumF', help='prefix for summary outputs', default='')
    parser.add_argument('-bowtie2', '-bw', '-BW', help='Path to bowtie2. If not used assumes that bowtie2 is in the'
                                                       'default path', default='')
    parser.add_argument('-rsem', '-RSEM', help='Path to rsem. If not used assumes that rsem is in the'
                                               'default path', default='/data/yosef/users/safik/bin/rsem-1.2.21/')
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
    parser.add_argument('-readOverlap', '-ro', '-readoverlap',
                        help='Add a read to list of mapped reads only if it maps at least X bases'
                             'to the V/J/C segment. Default is 1', type=int, default=1)
    parser.add_argument('-byExp',
                        help='if using the Top option, add this tag if you want to take only two chains from each' \
                             'read count, until top is reached', action='store_true')

    parser.add_argument('-overlap', '-ol', '-OL', help='Number of minimum bases that overlaps V and J ends,'
                                                       'default is 10', type=int, default=10)
    args = parser.parse_args()
    runTCRpipe(args.genome, args.output, args.bam, args.unmapped, args.bases, args.strand,
               args.iterations, args.score, args.overlap, args.rsem, args.bowtie2,
               args.singleCell, args.path, args.sumF, args.lowQ, args.samtools, args.top, args.byExp, args.readOverlap,
               args.oneSide)
