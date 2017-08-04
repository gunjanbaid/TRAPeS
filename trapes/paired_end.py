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


def analyzeChain(fastaDict, vdjDict, output, bam, unmapped, idNameDict, bases, chain, strand, lowQ, top, byExp,
                 readOverlap):
    junctionSegs = makeJunctionFile(bam, chain, output, bases, vdjDict, fastaDict, idNameDict, top, byExp, readOverlap)
    unDict = writeReadsFile(bam, unmapped, junctionSegs, output, vdjDict, chain, strand, lowQ)
    return unDict


def makeJunctionFile(bam, chain, output, bases, vdjDict, fastaDict, idNameDict, top, byExp, readOverlap):
    mappedFile = pysam.AlignmentFile(bam, "rb")
    if chain == 'A':
        vdjChainDict = vdjDict['Alpha']
        outName = output + '.alpha.junctions.txt'
    elif chain == 'B':
        vdjChainDict = vdjDict['Beta']
        outName = output + '.beta.junctions.txt'
    else:
        sys.stderr.write(
            str(datetime.datetime.now()) + ' Error! chain parameter for function analyzeChain can only be A or B\n')
        sys.stderr.flush()
    jSegs = vdjChainDict['J']
    vSegs = vdjChainDict['V']
    (cSeq, cId, cName) = getCInfo(vdjChainDict['C'][0], idNameDict, fastaDict)
    vjSegs = []
    for x in jSegs:
        vjSegs.append(x)
    for y in vSegs:
        vjSegs.append(y)
    vjReads = dict()
    (vjReads, vjCounts) = loadReadsToDict(vjSegs, mappedFile, vjReads, readOverlap)
    junctionSegs = writeJunctions(vjReads, outName, bases, fastaDict, idNameDict, cSeq, top, vjCounts, byExp)
    if len(junctionSegs) == 0:
        sys.stdout.write(
            str(datetime.datetime.now()) + ' Did not find any V-J reads, searching for V-C and J-C reads:\n')
        sys.stdout.flush()
        cReads = dict()
        (cReads, cCountsDict) = loadReadsToDict(vdjChainDict['C'], mappedFile, cReads, readOverlap)
        junctionSegs = writeJunctionsWithC(vjReads, outName, bases, fastaDict, idNameDict, cReads)
    mappedFile.close()
    return junctionSegs


def writeJunctions(vjReads, outName, bases, fastaDict, idNameDict, cSeq, top, vjCountsDict, byExp):
    out = open(outName, 'w')
    fArr = []
    pairCountDict = dict()
    for seg in vjReads:
        if idNameDict[seg].find('J') != -1:
            if len(vjReads[seg]) > 0:
                for sSeg in vjReads:
                    if ((idNameDict[sSeg].find('V') != -1) & (len(vjReads[sSeg]) > 0)):
                        if (len([val for val in vjReads[seg]['first'] if val in vjReads[sSeg]['second']]) > 0) | \
                                (len([val for val in vjReads[seg]['second'] if val in vjReads[sSeg]['first']]) > 0):
                            if seg not in fArr:
                                fArr.append(seg)
                            if sSeg not in fArr:
                                fArr.append(sSeg)
                            vSeq = fastaDict[sSeg]
                            jSeq = fastaDict[seg]
                            lenSeg = min(len(vSeq), len(jSeq))
                            if bases != -10:
                                if lenSeg < bases:
                                    if bases > len(vSeq):
                                        sys.stdout.write(str(
                                            datetime.datetime.now()) + ' Bases parameter is bigger than the length of '
                                                                       'the V segment, taking the length' \
                                                                       'of the V/J segment instead, which is: ' + str(
                                            lenSeg) + '\n')
                                        sys.stdout.flush()
                                    else:
                                        sys.stdout.write(str(
                                            datetime.datetime.now()) + ' Bases parameter is bigger than the length of '
                                                                       'the J segment, appending the C segment to the '
                                                                       'J segment\n')
                                        sys.stdout.flush()
                                        jSeq = jSeq + cSeq
                                        lenSeg = bases
                                else:
                                    lenSeg = bases
                            jTrim = jSeq[:lenSeg]
                            vTrim = vSeq[-1 * lenSeg:]
                            junc = vTrim + jTrim
                            recordName = sSeg + '.' + seg + '(' + idNameDict[sSeg] + '-' + idNameDict[seg] + ')'
                            record = SeqRecord(Seq(junc, IUPAC.ambiguous_dna), id=recordName, description='')
                            curCont = vjCountsDict[seg] + vjCountsDict[sSeg]
                            pairCountDict[record] = curCont
    sorted_pairs = sorted(pairCountDict.items(), key=operator.itemgetter(1), reverse=True)
    if ((top == -1) | (top > len(sorted_pairs))):
        for rec in pairCountDict:
            SeqIO.write(rec, out, 'fasta')
    else:
        if not byExp:
            for i in range(0, top):
                SeqIO.write(sorted_pairs[i][0], out, 'fasta')
        else:
            wrote = 1
            SeqIO.write(sorted_pairs[0][0], out, 'fasta')
            curCount = sorted_pairs[0][1]
            wroteSecond = False
            for i in range(1, len(sorted_pairs)):
                if sorted_pairs[i][1] == curCount:
                    if not wroteSecond:
                        wroteSecond = True
                        SeqIO.write(sorted_pairs[i][0], out, 'fasta')
                        wrote += 1
                else:
                    curCount = sorted_pairs[i][1]
                    wroteSecond = False
                    SeqIO.write(sorted_pairs[i][0], out, 'fasta')
                    wrote += 1
                if wrote == top:
                    break

    out.close()
    return fArr


# Load all the reads from a list of segments into a dictionary.
# INPUT: segsDict: A dict, where the key is a segment, and the value is an array of bed entries of this segment
#        mappedFile: Bam file of the mapped reaeds
#       readDict: A dictionary. The keys are the segment name, the values is a dictionary 'first':[] and 'second':[]
#                 where 'first' array holds the query name of R1's that overlap this segment, and 'second' holds the
#                 query name of R2's that overlap the segment.
def loadReadsToDict(segsDict, mappedFile, readDict, readOverlap):
    countDict = dict()
    for seg in segsDict:
        lArr = seg.strip('\n').split('\t')
        segName = lArr[3]
        readDict[segName] = {'first': [], 'second': []}
        lArr = seg.strip('\n').split('\t')
        chr = lArr[0]
        start = int(lArr[1])
        end = int(lArr[2])
        readsIter = mappedFile.fetch(chr, start - 1, end + 1)
        readCounter = 0
        for read in readsIter:
            overlap = read.get_overlap(start - 1, end + 1)
            if (end - start) < readOverlap:
                readOverlap = end - start - 15
            if overlap >= readOverlap:
                currName = read.query_name
                if read.is_read1:
                    if currName not in readDict[segName]['first']:
                        readDict[segName]['first'].append(currName)
                        readCounter += 1
                elif read.is_read2:
                    if currName not in readDict[segName]['second']:
                        readDict[segName]['second'].append(currName)
                        readCounter += 1
        countDict[segName] = readCounter
    return (readDict, countDict)


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
def writeJunctionsWithC(vjReads, outName, bases, fastaDict, idNameDict, cReads):
    out = open(outName, 'w')
    fArr = []
    vArr = []
    jArr = []
    for seg in vjReads:
        if len(vjReads[seg]) > 0:
            for cSeg in cReads:
                if (len(cReads[cSeg]) > 0):
                    if (len([val for val in vjReads[seg]['first'] if val in cReads[cSeg]['second']]) > 0) | \
                            (len([val for val in vjReads[seg]['second'] if val in cReads[cSeg]['first']]) > 0):
                        if idNameDict[seg].find('J') != -1:
                            if seg not in jArr:
                                jArr.append(seg)
                        elif idNameDict[seg].find('V') != -1:
                            if seg not in vArr:
                                vArr.append(seg)
                        fArr.append(seg)
    for vSeg in vArr:
        for jSeg in jArr:
            vSeqFa = fastaDict[vSeg]
            jSeqFa = fastaDict[jSeg]
            lenSeg = min(len(vSeqFa), len(jSeqFa))
            if bases != -10:
                if lenSeg < bases:
                    sys.stdout.write(str(
                        datetime.datetime.now()) + ' Bases parameter is bigger than the length of the V or J segment, '
                                                   'taking the length' \
                                                   'of the V/J segment instead, which is: ' + str(lenSeg) + '\n')
                    sys.stdout.flush()
                else:
                    lenSeg = bases
            jTrim = jSeqFa[:lenSeg]
            vTrim = vSeqFa[-1 * lenSeg:]
            junc = vTrim + jTrim
            recordName = vSeg + '.' + jSeg + '(' + idNameDict[vSeg] + '-' + idNameDict[jSeg] + ')'
            record = SeqRecord(Seq(junc, IUPAC.ambiguous_dna), id=recordName, description='')
            SeqIO.write(record, out, 'fasta')
    out.close()
    return fArr


def writeReadsFile(bam, unmapped, junctionSegs, output, vdjDict, chain, strand, lowQ):
    if chain == 'A':
        vdjChainDict = vdjDict['Alpha']
        outReads = output + '.alpha.mapped.and.unmapped.fa'
        pairedReads1 = output + '.alpha.R1.fa'
        pairedReads2 = output + '.alpha.R2.fa'
    elif chain == 'B':
        vdjChainDict = vdjDict['Beta']
        outReads = output + '.beta.mapped.and.unmapped.fa'
        pairedReads1 = output + '.beta.R1.fa'
        pairedReads2 = output + '.beta.R2.fa'
    else:
        sys.stderr.write(
            str(datetime.datetime.now()) + ' Error! chain parameter for function analyzeChain can only be A or B\n')
        sys.stderr.flush()
    out = open(outReads, 'w')
    constDict = vdjChainDict['C']
    # This dict for unmapped reads has reads that should be rev.comp in the revcomp arr, otherwise in id.
    # For every read, the value is a tuple - the first value is first/second, to make sure there are no errors.
    # The second value is id/revcomp, to see what sequence should be written.
    # Note: This classification is about what should happen to the unmapped reads, not how their paired maaped
    # reads were read.
    unmappedDict = dict()
    seqDict = dict()
    alignedDict = dict()
    mappedPairsDict = dict()
    lowQDict = dict()
    for seg in constDict:
        (unmappedDict, alignedDict, seqDict, mappedPairsDict, lowQDict) = addReadsToDict(unmappedDict, seg, bam, out,
                                                                                         False, alignedDict, seqDict,
                                                                                         strand, 'C', mappedPairsDict,
                                                                                         lowQDict)
    vSegs = vdjChainDict['V']
    for vSeg in vSegs:
        vSegName = vSeg.strip('\n').split('\t')[3]
        if vSegName in junctionSegs:
            (unmappedDict, alignedDict, seqDict, mappedPairsDict, lowQDict) = addReadsToDict(unmappedDict, vSeg, bam,
                                                                                             out, True, alignedDict,
                                                                                             seqDict, strand, 'V',
                                                                                             mappedPairsDict, lowQDict)
    jSegs = vdjChainDict['J']
    for jSeg in jSegs:
        jSegName = jSeg.strip('\n').split('\t')[3]
        if jSegName in junctionSegs:
            (unmappedDict, alignedDict, seqDict, mappedPairsDict, lowQDict) = addReadsToDict(unmappedDict, jSeg, bam,
                                                                                             out, True, alignedDict,
                                                                                             seqDict, strand, 'J',
                                                                                             mappedPairsDict, lowQDict)
    unDict = dict()
    (seqDict, unDict) = writeUnmappedReads(unmappedDict, out, unmapped, seqDict, unDict, alignedDict, lowQDict, lowQ)
    seqDict = addMappedPairsToSeqDict(seqDict, bam, out, lowQ, alignedDict)
    writeSeqDict(seqDict, pairedReads1, pairedReads2)
    out.close()
    return unDict


# Aligned dict - all the reads (with _1/_2) that were already written to the mapped.unmapped.fa file
def addReadsToDict(unmappedDict, segBed, bam, out, mappedRead, alignedDict, seqDict, strand, segType, mappedPairsDict,
                   lowQDict):
    bedArr = segBed.strip('\n').split('\t')
    chr = bedArr[0]
    start = int(bedArr[1])
    end = int(bedArr[2])
    mappedFile = pysam.AlignmentFile(bam, "rb")
    readsIter = mappedFile.fetch(chr, start - 1, end + 1)
    for read in readsIter:
        currName = read.query_name
        if currName not in seqDict:
            seqDict[currName] = ['0', '1']
        if read.is_read1:
            pairRead = 'second'
            readName = currName + '\\1'
            pairName = currName + '\\2'
            seqPos = 0
        elif read.is_read2:
            pairRead = 'first'
            readName = currName + '\\2'
            pairName = currName + '\\1'
            seqPos = 1
        else:
            sys.stderr.write(str(datetime.datetime.now()) + ' Error! Read is not read1 and not read2\n')
            sys.stderr.flush()
        currSeq = Seq(read.query_sequence, IUPAC.ambiguous_dna)
        if read.is_reverse:
            pairOr = 'id'
            readStrand = 'minus'
            currSeq = currSeq.reverse_complement()
        else:
            readStrand = 'plus'
            pairOr = 'rev'
        seqDict[currName][seqPos] = currSeq
        if read.mate_is_unmapped:
            takePair = toTakePair(segType, strand, readStrand)
            if takePair == False:
                lowQDict[pairName] = '1'
            # takePair = True
            # if takePair:
            if currName in unmappedDict:
                if unmappedDict[currName] != (pairRead, pairOr):
                    sys.stderr.write(str(
                        datetime.datetime.now()) + ' Error! Read %s has more than one unmppaed mate with differnet '
                                                   'strand/mate\n' % currName)
                    sys.stderr.flush()
            unmappedDict[currName] = (pairRead, pairOr)
        if mappedRead == True:
            if readName in alignedDict:
                if alignedDict[readName] != read.query_sequence:
                    sys.stderr.write(str(
                        datetime.datetime.now()) + ' Error! Read %s has two instances but different seuqences\n' %
                                     read.query_name)
                    sys.stderr.flush()
            else:
                alignedDict[readName] = read.query_sequence
                record = SeqRecord(Seq(read.query_sequence, IUPAC.ambiguous_dna), id=readName, description='')
                SeqIO.write(record, out, 'fasta')
    mappedFile.close()
    # print "Removed counter: " + str(cc)
    return (unmappedDict, alignedDict, seqDict, mappedPairsDict, lowQDict)


def writeSeqDict(seqDict, r1, r2):
    r1f = open(r1, 'w')
    r2f = open(r2, 'w')
    for seq in seqDict:
        if ((seqDict[seq][0] != '0') & (seqDict[seq][1] != '1')):
            seq1 = seq
            seq2 = seq
            rec1 = SeqRecord(seqDict[seq][0], id=seq1, description='')
            rec2 = SeqRecord(seqDict[seq][1], id=seq2, description='')
            SeqIO.write(rec1, r1f, 'fasta')
            SeqIO.write(rec2, r2f, 'fasta')
        else:
            sys.stderr.write(str(datetime.datetime.now()) + ' The read %s has only one mate found, ignoring it\n' % seq)
            sys.stderr.flush()
    r1f.close()
    r2f.close()


### For minus, add an unmapped pair if the current mate is: 1. V and Plus 2. J/C and minus
### For plus, exactly the opposite

def toTakePair(segType, strand, readStrand):
    if strand == 'none':
        return True
    if ((readStrand != 'minus') & (readStrand != 'plus')):
        sys.stderr.write(str(datetime.datetime.now()) + ' Error! Read strand should be plus or minus only\n')
        sys.stderr.flush()
        return True
    if ((segType == 'C') | (segType == 'J')):
        if strand == 'minus':
            if readStrand == 'minus':
                return True
            else:
                return False
        else:
            if readStrand == 'minus':
                return False
            else:
                return True
    else:
        if strand == 'minus':
            if readStrand == 'minus':
                return False
            else:
                return True
        else:
            if readStrand == 'minus':
                return True
            else:
                return False
    return True


def writeUnmappedReads(unmappedDict, out, unmapped, seqDict, unDict, alignedDict, lowQDict, lowQ):
    f = pysam.AlignmentFile(unmapped, "rb")
    readsIter = f.fetch(until_eof=True)
    for read in readsIter:
        name = read.query_name
        if name in unmappedDict:
            cName = name
            unDictName = name
            (strand, ori) = unmappedDict[name]
            if (((strand == 'first') & (read.is_read2)) | ((strand == 'second') & (read.is_read1))):
                sys.stderr.write(str(
                    datetime.datetime.now()) + ' Error! unmapped read is inconsistent regarding first/second read\n')
                sys.stderr.flush()
            else:
                if strand == 'first':
                    name += '\\1'
                    unDictName += '_1'
                else:
                    name += '\\2'
                    unDictName += '_2'
                if unDictName in unDict:
                    sys.stderr.write(str(
                        datetime.datetime.now()) + ' Error! unmapped read %s appear twice in unmapped bam file\n' %
                                     cName)
                    sys.stderr.flush()
                unDict[unDictName] = '1'
                qSeq = Seq(read.query_sequence, IUPAC.ambiguous_dna)
                if ori == 'rev':
                    qSeq = qSeq.reverse_complement()
                if name in alignedDict:
                    if alignedDict[name] != str(qSeq):
                        sys.stderr.write(str(
                            datetime.datetime.now()) + ' Error! unmapped read %s appear twice in alignedDict with '
                                                       'differnet seqs\n' % name)
                        sys.stderr.flush()
                else:
                    if ((name not in lowQDict) | ((name in lowQDict) & (not lowQ))):
                        alignedDict[name] = str(qSeq)
                        record = SeqRecord(qSeq, id=name, description='')
                        SeqIO.write(record, out, 'fasta')
                if ori == 'rev':
                    qSeq = qSeq.reverse_complement()
            if cName not in seqDict:
                sys.stderr.write(str(
                    datetime.datetime.now()) + ' Error! unmapped read is in unmappedDict but not in seqDict %s\n' %
                                 read.query_name)
                sys.stderr.flush()
            else:
                if strand == 'first':
                    seqDict[cName][0] = qSeq
                else:
                    seqDict[cName][1] = qSeq
    f.close()
    return (seqDict, unDict)


def addMappedPairsToSeqDict(seqDict, bam, out, lowQ, alignedDict):
    firstDict = dict()
    secondDict = dict()
    for name in seqDict:
        if seqDict[name][0] == '0':
            firstDict[name] = '1'
        if seqDict[name][1] == '1':
            secondDict[name] = '1'
        if ((seqDict[name][1] == '1') & (seqDict[name][0] == '0')):
            sys.stderr.write(str(datetime.datetime.now()) + ' Error! empty record insdie seqDict\n')
            sys.stderr.flush()
    f = pysam.AlignmentFile(bam, "rb")
    readsIter = f.fetch()
    for read in readsIter:
        name = read.query_name
        pos = -1
        if ((read.query_name in firstDict) & (read.is_read1)):
            pos = 0
            check = '0'
        elif ((read.query_name in secondDict) & (read.is_read2)):
            pos = 1
            check = '1'
        if pos != -1:
            qSeq = Seq(read.query_sequence, IUPAC.ambiguous_dna)
            if read.is_read1:
                currName = name + '\\1'
            else:
                currName = name + '\\2'
            if lowQ:
                if read.mate_is_reverse:
                    if read.is_reverse:
                        rSeq = qSeq.reverse_complement()
                    else:
                        rSeq = qSeq
                else:
                    if not read.is_reverse:
                        rSeq = qSeq.reverse_complement()
                    else:
                        rSeq = qSeq
                if currName not in alignedDict:
                    alignedDict[currName] = str(rSeq)
                    record = SeqRecord(rSeq, id=currName, description='')
                    SeqIO.write(record, out, 'fasta')
                else:
                    if str(alignedDict[currName]) != str(rSeq):
                        sys.stderr.write(str(
                            datetime.datetime.now()) + ' Error! read %s has two different sequences in alignedDict\n'
                                         % currName)
                        sys.stderr.flush()
            if read.is_reverse:
                qSeq = qSeq.reverse_complement()
            if seqDict[name][pos] != check:
                if str(qSeq) != str(seqDict[name][pos]):
                    sys.stderr.write(str(
                        datetime.datetime.now()) + ' Error! read %s has two different mapped sequences not in the V/J '
                                                   'region\n' % name)
                    sys.stderr.flush()
            seqDict[name][pos] = qSeq
    f.close()
    return seqDict


def runRsem(outDir, rsem, bowtie2, fullTcrFileAlpha, fullTcrFileBeta, output, samtools):
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
            subprocess.call([rsem + 'rsem-prepare-reference', '--bowtie2', '--bowtie2-path', bowtie2,
                             '-q', fullTcrFileAlpha, rsemIndDir + '/VDJ.alpha.seq'])
            subprocess.call([rsem + 'rsem-calculate-expression', '--no-qualities', '-q',
                             '--bowtie2', '--bowtie2-path', bowtie2, '--bowtie2-mismatch-rate', '0.0', '--paired-end',
                             output + '.alpha.R1.fa',
                             output + '.alpha.R2.fa', rsemIndDir + '/VDJ.alpha.seq', output + '.alpha.rsem.out'])
        else:
            subprocess.call([rsem + 'rsem-prepare-reference', '--bowtie2',
                             '-q', fullTcrFileAlpha, rsemIndDir + '/VDJ.alpha.seq'])
            subprocess.call([rsem + 'rsem-calculate-expression', '--no-qualities', '-q',
                             '--bowtie2', '--bowtie2-mismatch-rate', '0.0', '--paired-end', output + '.alpha.R1.fa',
                             output + '.alpha.R2.fa', rsemIndDir + '/VDJ.alpha.seq', output + '.alpha.rsem.out'])
        unsortedBam = output + '.alpha.rsem.out.transcript.bam'
        if not os.path.exists(unsortedBam):
            print "RSEM did not produce any transcript alignment files for alpha chain, please check the -rsem " \
                  "parameter"
        else:
            sortedBam = output + '.alpha.rsem.out.transcript.sorted.bam'
            if not os.path.exists(sortedBam):
                subprocess.call([samtools + 'samtools', 'sort', '-o', sortedBam, unsortedBam])
                subprocess.call([samtools + 'samtools', 'index', sortedBam])

    else:
        sys.stdout.write(
            str(datetime.datetime.now()) + " Did not reconstruct any alpha chains, not running RSEM on alpha\n")
        sys.stdout.flush()
    if os.path.exists(fullTcrFileBeta):
        if bowtie2 != '':
            subprocess.call([rsem + 'rsem-prepare-reference', '--bowtie2', '--bowtie2-path', bowtie2,
                             '-q', fullTcrFileBeta, rsemIndDir + '/VDJ.beta.seq'])
            subprocess.call([rsem + 'rsem-calculate-expression', '--no-qualities', '-q', '--bowtie2', '--bowtie2-path',
                             bowtie2, '--bowtie2-mismatch-rate', '0.0', '--paired-end', output + '.beta.R1.fa',
                             output + '.beta.R2.fa',
                             rsemIndDir + '/VDJ.beta.seq', output + '.beta.rsem.out'])
        else:
            subprocess.call([rsem + 'rsem-prepare-reference', '--bowtie2',
                             '-q', fullTcrFileBeta, rsemIndDir + '/VDJ.beta.seq'])
            subprocess.call([rsem + 'rsem-calculate-expression', '--no-qualities', '-q', '--bowtie2',
                             '--bowtie2-mismatch-rate', '0.0', '--paired-end', output + '.beta.R1.fa',
                             output + '.beta.R2.fa',
                             rsemIndDir + '/VDJ.beta.seq', output + '.beta.rsem.out'])
        unsortedBam = output + '.beta.rsem.out.transcript.bam'
        if not os.path.exists(unsortedBam):
            print "RSEM did not produce any transcript alignment files for beta chain, please check the -rsem parameter"
        else:
            sortedBam = output + '.beta.rsem.out.transcript.sorted.bam'
            if not os.path.exists(sortedBam):
                subprocess.call([samtools + 'samtools', 'sort', '-o', sortedBam, unsortedBam])
                subprocess.call([samtools + 'samtools', 'index', sortedBam])
    else:
        sys.stdout.write(
            str(datetime.datetime.now()) + " Did not reconstruct any beta chains, not running RSEM on beta\n")
        sys.stdout.flush()
