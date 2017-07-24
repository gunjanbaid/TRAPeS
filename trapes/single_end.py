#! /usr/bin/env python

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


def analyzeChainSingleEnd(fastq, trimmomatic, transInd, bowtie2, idNameDict, output, fastaDict, bases):
    mappedReadsDictAlpha = dict()
    mappedReadsDictBeta = dict()
    lenArr = [999, 50, 25]
    for currLen in lenArr:
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
            trimFq = fastq + '.' + str(currLen) + '.' + str(side) + '.trimmed.fq'
            # TODO: use bowtie trimmer instead
            # TODO: make sure about minus strand alignment
            if crop == '':
                subprocess.call(
                    ['java', '-jar', trimmomatic, 'SE', '-phred33', fastq, trimFq, 'LEADING:15', 'TRAILING:15',
                     'MINLEN:20'])
            else:
                subprocess.call(
                    ['java', '-jar', trimmomatic, 'SE', '-phred33', fastq, trimFq, 'LEADING:15', 'TRAILING:15', crop,
                     'MINLEN:20'])
            samF = trimFq + '.sam'
            if bowtie2 != '':
                if bowtie2.endswith('/'):
                    bowtieCall = bowtie2 + 'bowtie2'
                else:
                    bowtieCall = bowtie2 + '/bowtie2'
            else:
                bowtieCall = 'bowtie2'
            subprocess.call([bowtieCall, '-q --phred33  --score-min L,0,0', '-x', transInd, '-U', trimFq, '-S', samF])
            if os.path.isfile(samF):
                mappedReadsDictAlpha = findReadsAndSegments(samF, mappedReadsDictAlpha, idNameDict, 'A')
                mappedReadsDictBeta = findReadsAndSegments(samF, mappedReadsDictBeta, idNameDict, 'B')
    alphaOut = output + '.alpha.junctions.txt'
    alphaOutReads = output + '.alpha.mapped.and.unmapped.fa'
    betaOutReads = output + '.beta.mapped.and.unmapped.fa'
    betaOut = output + '.beta.junctions.txt'
    writeJunctionFileSE(mappedReadsDictAlpha, idNameDict, alphaOut, fastaDict, bases, 'alpha')
    writeJunctionFileSE(mappedReadsDictBeta, idNameDict, betaOut, fastaDict, bases, 'beta')
    writeReadsFileSE(mappedReadsDictAlpha, alphaOutReads, fastq)
    writeReadsFileSE(mappedReadsDictBeta, betaOutReads, fastq)
    sys.exit(1)


def findReadsAndSegments(samF, mappedReadsDict, idNameDict, chain):
    samFile = pysam.AlignmentFile(samF, 'r')
    readsIter = samFile.fetch(until_eof=True)
    for read in readsIter:
        if read.is_unmapped == False:
            seg = samFile.getrname(read.reference_id)
            if seg in idNameDict:
                if idNameDict[seg].find(chain) != -1:
                    readName = read.query_name
                    if readName not in mappedReadsDict:
                        mappedReadsDict[readName] = []
                    if seg not in mappedReadsDict[readName]:
                        mappedReadsDict[readName].append(seg)
    samFile.close()
    return mappedReadsDict


def writeReadsFileSE(mappedReadsDict, outReads, fastq):
    if fastq.endswith('.gz'):
        subprocess.call(['gunzip', fastq])
        newFq = fastq.replace('.gz', '')
    else:
        newFq = fastq
    out = open(outReads, 'w')
    fqF = open(newFq, 'rU')
    for record in SeqIO.parse(fqF, 'fastq'):
        if record.id in mappedReadsDict:
            newRec = SeqRecord(record.seq, id=record.id, description='')
            SeqIO.write(newRec, out, 'fasta')
    out.close()
    fqF.close()
    if fastq.endswith('.gz'):
        subprocess.call(['gzip', newFq])


def writeJunctionFileSE(mappedReadsDict, idNameDict, output, fastaDict, bases, chain):
    out = open(output, 'w')
    vSegs = []
    jSegs = []
    cSegs = []
    for read in mappedReadsDict:
        for seg in mappedReadsDict[read]:
            if idNameDict[seg].find('V') != -1:
                if seg not in vSegs:
                    vSegs.append(seg)
            elif idNameDict[seg].find('J') != -1:
                if seg not in jSegs:
                    jSegs.append(seg)
            elif idNameDict[seg].find('C') != -1:
                if seg not in cSegs:
                    cSegs.append(seg)
            else:
                print "Error! not V/J/C in fasta dict"
    if len(vSegs) == 0:
        print "Did not find any V segments for " + chain + " chain"
    else:
        if len(cSegs) == 0:
            print "Did not find any C segments for " + chain + " chain"
            cSegs = ['NA']
        if len(jSegs) == 0:
            print "Did not find any J segments for " + chain + " chain"
            jSegs = ['NA']
        for vSeg in vSegs:
            for jSeg in jSegs:
                for cSeg in cSegs:
                    addSegmentToJunctionFileSE(vSeg, jSeg, cSeg, out, fastaDict, bases, idNameDict)
    out.close()


def addSegmentToJunctionFileSE(vSeg, jSeg, cSeg, out, fastaDict, bases, idNameDict):
    vSeq = fastaDict[vSeg]
    if jSeg != 'NA':
        jName = idNameDict[jSeg]
        jSeq = fastaDict[jSeg]
    else:
        jSeq = ''
        jName = 'NoJ'
    if cSeg != 'NA':
        cName = idNameDict[cSeg]
        cSeq = fastaDict[cSeg]
    else:
        cName = 'NoC'
        cSeq = ''
    jcSeq = jSeq + cSeq
    lenSeg = min(len(vSeq), len(jcSeq))
    if bases != -10:
        if lenSeg < bases:
            sys.stdout.write(str(
                datetime.datetime.now()) + ' Bases parameter is bigger than the length of the V or J segment, '
                                           'taking the length'
                                           'of the V/J segment instead, which is: ' + str(lenSeg) + '\n')
            sys.stdout.flush()
        else:
            lenSeg = bases
    jTrim = jcSeq[:lenSeg]
    vTrim = vSeq[-1 * lenSeg:]
    junc = vTrim + jTrim
    recordName = vSeg + '.' + jSeg + '.' + cSeg + '(' + idNameDict[vSeg] + '-' + jName + '-' + cName + ')'
    record = SeqRecord(Seq(junc, IUPAC.ambiguous_dna), id=recordName, description='')
    SeqIO.write(record, out, 'fasta')
