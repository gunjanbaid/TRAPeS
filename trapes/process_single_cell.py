"""
This module contains functions used to process each individual cell.

"""

import sys
import os
import subprocess
import datetime
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

from paired_end import analyze_chain, run_rsem
from write_output_files import makeSingleCellOutputFile
from utils import get_c_info

def runSingleCell(fasta, bed, output, bam, unmapped, mapping, bases, strand, reconstruction, aaF, numIterations,
                  thresholdScore, minOverlap, rsem, bowtie2, lowQ, samtools, top, byExp, readOverlap, oneSide):
    idNameDict = makeIdNameDict(mapping)
    fastaDict = makeFastaDict(fasta)
    vdjDict = makeVDJBedDict(bed, idNameDict)
    sys.stdout.write(str(datetime.datetime.now()) + " Pre-processing alpha chain\n")
    sys.stdout.flush()
    unDictAlpha = analyze_chain(fastaDict, vdjDict, output, bam, unmapped, idNameDict, bases, 'A', strand, lowQ, top,
                                byExp, readOverlap)
    sys.stdout.write(str(datetime.datetime.now()) + " Pre-processing beta chain\n")
    sys.stdout.flush()
    unDictBeta = analyze_chain(fastaDict, vdjDict, output, bam, unmapped, idNameDict, bases, 'B', strand, lowQ, top,
                               byExp, readOverlap)
    sys.stdout.write(str(datetime.datetime.now()) + " Reconstructing beta chains\n")
    sys.stdout.flush()
    subprocess.call([reconstruction, output + '.beta.mapped.and.unmapped.fa', output + '.beta.junctions.txt',
                     output + '.reconstructed.junctions.beta.fa', str(numIterations), str(thresholdScore),
                     str(minOverlap)])
    sys.stdout.write(str(datetime.datetime.now()) + " Reconstructing alpha chains\n")
    sys.stdout.flush()
    subprocess.call([reconstruction, output + '.alpha.mapped.and.unmapped.fa', output + '.alpha.junctions.txt',
                     output + '.reconstructed.junctions.alpha.fa', str(numIterations), str(thresholdScore),
                     str(minOverlap)])
    sys.stdout.write(str(datetime.datetime.now()) + " Creating full TCR sequencing\n")
    sys.stdout.flush()
    fullTcrFileAlpha = output + '.alpha.full.TCRs.fa'
    tcrF = output + '.reconstructed.junctions.alpha.fa'
    (cSeq, cName, cId) = get_c_info(vdjDict['Alpha']['C'][0], idNameDict, fastaDict)
    createTCRFullOutput(fastaDict, tcrF, fullTcrFileAlpha, bases, idNameDict, cSeq, cName, cId, oneSide)
    fullTcrFileBeta = output + '.beta.full.TCRs.fa'
    tcrF = output + '.reconstructed.junctions.beta.fa'
    (cSeq, cName, cId) = get_c_info(vdjDict['Beta']['C'][0], idNameDict, fastaDict)
    createTCRFullOutput(fastaDict, tcrF, fullTcrFileBeta, bases, idNameDict, cSeq, cName, cId, oneSide)
    sys.stdout.write(str(datetime.datetime.now()) + " Running RSEM to quantify expression of all possible isoforms\n")
    sys.stdout.flush()
    outDirInd = output.rfind('/')
    if outDirInd != -1:
        outDir = output[:outDirInd + 1]
    else:
        outDir = os.getcwd()
    run_rsem(outDir, rsem, bowtie2, fullTcrFileAlpha, fullTcrFileBeta, output, samtools)
    pickFinalIsoforms(fullTcrFileAlpha, fullTcrFileBeta, output)
    bestAlpha = output + '.alpha.full.TCRs.bestIso.fa'
    bestBeta = output + '.beta.full.TCRs.bestIso.fa'
    sys.stdout.write(str(datetime.datetime.now()) + " Finding productive CDR3\n")
    sys.stdout.flush()
    aaDict = makeAADict(aaF)
    if os.path.isfile(bestAlpha):
        fDictAlpha = findCDR3(bestAlpha, aaDict, fastaDict)
    else:
        fDictAlpha = dict()
    if os.path.isfile(bestBeta):
        fDictBeta = findCDR3(bestBeta, aaDict, fastaDict)
    else:
        fDictBeta = dict()
    betaRsemOut = output + '.beta.rsem.out.genes.results'
    alphaRsemOut = output + '.alpha.rsem.out.genes.results'
    alphaBam = output + '.alpha.rsem.out.transcript.sorted.bam'
    betaBam = output + '.beta.rsem.out.transcript.sorted.bam'
    sys.stdout.write(str(datetime.datetime.now()) + " Writing results to summary file\n")
    sys.stdout.flush()
    makeSingleCellOutputFile(fDictAlpha, fDictBeta, output, betaRsemOut, alphaRsemOut, alphaBam, betaBam, fastaDict,
                             unDictAlpha, unDictBeta, idNameDict)


# Creates a dictionary of ENSEMBL ID -> Gene name
def makeIdNameDict(mapping):
    f = open(mapping, 'r')
    fDict = dict()
    linesArr = f.read().split('\n')
    f.close()
    for line in linesArr:
        lineArr = line.split('\t')
        id = lineArr[0]
        name = lineArr[1]
        ind = name.find('Gene:')
        if ind != -1:
            name = name[ind + len('Gene:'):]
        if id in fDict:
            sys.stderr.write(str(datetime.datetime.now()) + ' Error! %s appear twice in mapping file\n' % id)
            sys.stderr.flush()
        fDict[id] = name
    return fDict


# Creates a dictionary of ENSEMBL ID -> fasta sequence
def makeFastaDict(fasta):
    inF = open(fasta, 'rU')
    fastaDict = dict()
    for record in SeqIO.parse(inF, 'fasta'):
        fastaDict[record.id] = str(record.seq)
    inF.close()
    return fastaDict


# Create a dict {'Alpha':{'C':[bed],'V':[bed],'J':[bed]}, 'Beta':{'C':[],'V':[],'J':[]}}
def makeVDJBedDict(bed, idNameDict):
    fDict = {'Alpha': {'C': [], 'V': [], 'J': []}, 'Beta': {'C': [], 'V': [], 'J': []}}
    f = open(bed, 'r')
    l = f.readline()
    while l != '':
        lArr = l.strip('\n').split('\t')
        gID = lArr[3]
        gName = idNameDict[gID]
        chain = ''
        if (gName.startswith('TRA')):
            chain = 'Alpha'
        elif (gName.startswith('TRB')):
            chain = 'Beta'
        else:
            sys.stderr.write(
                str(datetime.datetime.now()) + ' Error! %s name is not alpha or beta chain, ignoring it\n' % gName)
            sys.stderr.flush()
        if gName.find('C') != -1:
            fDict[chain]['C'].append(l)
        elif gName.find('V') != -1:
            fDict[chain]['V'].append(l)
        elif gName.find('J') != -1:
            fDict[chain]['J'].append(l)
        l = f.readline()
    f.close()
    return fDict


def makeAADict(aaF):
    fDict = dict()
    f = open(aaF, 'r')
    l = f.readline()
    while l != '':
        lArr = l.strip('\n').split('\t')
        if lArr[0] in fDict:
            sys.stderr.write(str(datetime.datetime.now()) + ' Warning! %s appear twice in AA file\n' % lArr[0])
            sys.stderr.flush()
        fDict[lArr[0]] = lArr[1]
        l = f.readline()
    f.close()
    return fDict

def createTCRFullOutput(fastaDict, tcr, outName, bases, mapDict, cSeq, cName, cId, oneSide):
    tcrF = open(tcr, 'rU')
    found = False
    ffound = False
    recNameArr = []
    for tcrRecord in SeqIO.parse(tcrF, 'fasta'):
        addedC = False
        tcrSeq = str(tcrRecord.seq)
        if tcrSeq.find('NNNNN') == -1:
            if ffound == False:
                ffound = True
                outF = open(outName, 'w')
            idArr = tcrRecord.id.split('.')
            vEns = idArr[0]
            jEns = idArr[1].split('(')[0]
            vSeq = fastaDict[vEns]
            jSeq = fastaDict[jEns]
            recNameArr = writeRecord(tcrRecord, tcrSeq, addedC, vEns, jEns, vSeq, jSeq, mapDict, bases, cSeq, cId,
                                     cName, outF, fastaDict, recNameArr)
        elif oneSide:
            curSeq = tcrSeq.split('NNNN')[0]
            jSeg = findBestJforSeq(curSeq, fastaDict, mapDict)
            if jSeg != 'NA':
                if ffound == False:
                    ffound = True
                    outF = open(outName, 'w')
                idArr = tcrRecord.id.split('.')
                vEns = idArr[0]
                vSeq = fastaDict[vEns]
                for jEns in jSeg:
                    jSeq = fastaDict[jEns]
                    recNameArr = writeRecord(tcrRecord, curSeq, addedC, vEns, jEns, vSeq, jSeq, mapDict, bases, cSeq,
                                             cId, cName, outF, fastaDict, recNameArr)
    tcrF.close()
    if found == True:
        outF.close()


def findBestJforSeq(curSeq, fastaDict, idNameDict):
    jArrOld = findJsPerLen(curSeq, fastaDict, idNameDict, 20)
    if len(jArrOld) == 0:
        return 'NA'
    for x in range(21, len(curSeq)):
        newArr = findJsPerLen(curSeq, fastaDict, idNameDict, x)
        if len(newArr) == 0:
            return jArrOld
        else:
            jArrOld = newArr
    print 'Found a full J segment as the V/J junction, ignoring this reconstruction'
    return 'NA'


def findJsPerLen(curSeq, fastaDict, idNameDict, trim):
    fArr = []
    for seq in fastaDict:
        if idNameDict[seq].find('J') != -1:
            jSeq = fastaDict[seq]
            lenJ = len(jSeq)
            for i in range(0, lenJ):
                if ((i + trim) <= lenJ):
                    trimJ = jSeq[i:i + trim]
                    if curSeq.find(trimJ) != -1:
                        if seq not in fArr:
                            fArr.append(seq)
                            break
    return fArr


def writeRecord(tcrRecord, tcrSeq, addedC, vEns, jEns, vSeq, jSeq, mapDict, bases, cSeq, cId, cName, outF, fastaDict,
                recNameArr):
    vSeqTrim = ''
    jSeqTrim = ''
    if bases == -10:
        bases = min(len(vSeq), len(jSeq))
    elif bases > len(jSeq):
        jSeq = jSeq + cSeq
        addedC = True
    found = False
    for i in reversed(range(20, bases)):
        juncStart = tcrSeq[:i]
        vInd = vSeq.find(juncStart)
        if (vInd != -1):
            found = True
            vSeqTrim = vSeq[:vInd]
            break
    if found == False:
        vSeqTrim = vSeq[:-bases]
    found = False
    for j in reversed(range(20, bases)):
        juncEnd = tcrSeq[-j:]
        jInd = jSeq.find(juncEnd)
        if (jInd != -1):
            found = True
            jSeqTrim = jSeq[jInd + j:]
            break
    if found == False:
        jSeqTrim = jSeq[bases:]
    # Add TRBC or TRAC
    cArr = []
    if (str(tcrRecord.id).find('TRB') != -1):
        for ens in mapDict:
            if mapDict[ens].find('TRBC') != -1:
                cArr.append(ens)
    elif (str(tcrRecord.id).find('TRA') != -1):
        for ens in mapDict:
            if mapDict[ens].find('TRAC') != -1:
                cArr.append(ens)
    else:
        sys.stderr.write(str(datetime.datetime.now()) + " Error! no TRBC or TRAC\n")
        sys.stderr.flush()
    if not addedC:
        for ens in cArr:
            cSeq = fastaDict[ens]
            newSeq = vSeqTrim + tcrSeq + jSeqTrim + cSeq
            newId = mapDict[vEns] + '.' + mapDict[jEns] + '.' + mapDict[ens] + '.' + vEns + '.' + jEns + '.' + ens
            while newId in recNameArr:
                newId += '_2'
            recNameArr.append(newId)
            record = SeqRecord(Seq(newSeq, IUPAC.ambiguous_dna), id=newId, description='')
            SeqIO.write(record, outF, 'fasta')
    else:
        newSeq = vSeqTrim + tcrSeq + jSeqTrim
        newId = mapDict[vEns] + '.' + mapDict[jEns] + '.' + cName + '.' + vEns + '.' + jEns + '.' + cId
        while newId in recNameArr:
            newId += '_2'
        recNameArr.append(newId)
        record = SeqRecord(Seq(newSeq, IUPAC.ambiguous_dna), id=newId, description='')
        SeqIO.write(record, outF, 'fasta')
    return recNameArr


def pickFinalIsoforms(fullTcrFileAlpha, fullTcrFileBeta, output):
    pickFinalIsoformChain(fullTcrFileAlpha, output + '.alpha.full.TCRs.bestIso.fa',
                          output + '.alpha.rsem.out.genes.results')
    pickFinalIsoformChain(fullTcrFileBeta, output + '.beta.full.TCRs.bestIso.fa',
                          output + '.beta.rsem.out.genes.results')


def pickFinalIsoformChain(fullTCRfa, newFasta, rsemF):
    if os.path.isfile(fullTCRfa):
        f = open(fullTCRfa, 'rU')
        outFa = open(newFasta, 'w')
        fastaDict = dict()
        byVJDict = dict()
        for record in SeqIO.parse(f, 'fasta'):
            if record.id in fastaDict:
                # record.id = record.id + '_2'
                sys.stderr.write(
                    str(datetime.datetime.now()) + 'error! same name for two fasta entries %s\n' % record.id)
                sys.stderr.flush()
            fastaDict[record.id] = record.seq
            onlyVJrec = str(record.id)
            idArr = onlyVJrec.strip('\n').split('.')
            vjStr = idArr[0] + '.' + idArr[1]
            if vjStr not in byVJDict:
                byVJDict[vjStr] = []
            byVJDict[vjStr].append(record.id)
        for vjStr in byVJDict:

            if len(byVJDict[vjStr]) == 1:
                cId = byVJDict[vjStr][0]
                cSeq = fastaDict[cId]
                newRec = SeqRecord(cSeq, id=cId, description='')
                SeqIO.write(newRec, outFa, 'fasta')
            else:
                # print vjStr
                # print byVJDict[vjStr]
                bestId = findBestC(byVJDict[vjStr], rsemF)
                # print "best: " + bestId
                bSeq = fastaDict[bestId]
                newRec = SeqRecord(bSeq, id=bestId, description='')
                SeqIO.write(newRec, outFa, 'fasta')
        outFa.close()
        f.close()


def findBestC(vjArr, rsemF):
    if (os.path.exists(rsemF)):
        f = open(rsemF, 'r')
        f.readline()
        l = f.readline()
        bestSeq = 'name'
        maxCount = 0.0
        while l != '':
            lArr = l.strip('\n').split('\t')
            if lArr[0] in vjArr:
                currCount = float(lArr[4])
                if currCount > maxCount:
                    bestSeq = lArr[0]
                    maxCount = currCount
            l = f.readline()
        f.close()
        if bestSeq == 'name':
            return vjArr[0]
        return bestSeq
    else:
        return vjArr[0]


def findCDR3(fasta, aaDict, vdjFaDict):
    f = open(fasta, 'rU')
    fDict = dict()
    for record in SeqIO.parse(f, 'fasta'):
        if record.id in fDict:
            sys.stderr.write(str(datetime.datetime.now()) + ' Error! same name for two fasta entries %s\n' % record.id)
            sys.stderr.flush()
        else:
            idArr = record.id.split('.')
            vSeg = idArr[0]
            jSeg = idArr[1]
            if ((vSeg in aaDict) & (jSeg in aaDict)):
                currDict = findVandJaaMap(aaDict[vSeg], aaDict[jSeg], record.seq)
            else:
                if vSeg in aaDict:
                    newVseg = aaDict[vSeg]
                else:
                    vId = idArr[3]
                    currSeq = vdjFaDict[vId]
                    newVseg = getBestVaa(Seq(currSeq))
                if jSeg in aaDict:
                    newJseg = aaDict[jSeg]
                else:
                    jId = idArr[4]
                    currSeq = vdjFaDict[jId]
                    newJseg = getBestJaa(Seq(currSeq))
                currDict = findVandJaaMap(newVseg, newJseg, record.seq)
            fDict[record.id] = currDict
    f.close()
    return fDict


def getNTseq(fullSeq):
    mod = len(fullSeq) % 3
    if mod != 0:
        fSeq = fullSeq[:-mod].translate()
    else:
        fSeq = fullSeq.translate()
    return fSeq


def findVandJaaMap(vSeg, jSeg, fullSeq):
    fDict = dict()
    firstSeq = getNTseq(fullSeq)
    secondSeq = getNTseq(fullSeq[1:])
    thirdSeq = getNTseq(fullSeq[2:])
    ntArr = [fullSeq, fullSeq[1:], fullSeq[2:]]
    aaSeqsArr = [firstSeq, secondSeq, thirdSeq]
    cdrArr = []
    posArr = []
    fPharr = []
    for aaSeq in aaSeqsArr:
        (cdr, pos, curPh) = getCDR3(aaSeq, vSeg, jSeg)
        cdrArr.append(cdr)
        posArr.append(pos)
        fPharr.append(curPh)
    maxLen = 0
    bestCDR = ''
    bestSeq = ''
    hasStop = False
    bestPos = -1
    bestCDRnt = ''
    foundGood = False
    vPos = -1
    jPos = -1
    fPh = False
    for i in range(0, 3):
        if posArr[i] != -1:
            if ((cdrArr[i] != 'Only J') & (cdrArr[i] != 'Only V')):
                if len(cdrArr[i]) > maxLen:
                    if cdrArr[i].find('*') == -1:
                        foundGood = True
                        bestCDR = cdrArr[i]
                        bestPos = posArr[i]
                        maxLen = len(cdrArr[i])
                        bestSeq = ntArr[i]
                        fPh = fPharr[i]
                    else:
                        if maxLen == 0:
                            foundGood = True
                            bestPos = posArr[i]
                            bestCDR = cdrArr[i]
                            maxLen = len(cdrArr[i])
                            bestSeq = ntArr[i]
                            hasStop = True
                            fPh = fPharr[i]
                else:
                    if hasStop == True:
                        if cdrArr[i].find('*') == -1:
                            foundGood = True
                            bestPos = posArr[i]
                            hasStop = False
                            bestCDR = cdrArr[i]
                            maxLen = len(cdrArr[i])
                            bestSeq = ntArr[i]
                            fPh = fPharr[i]
            else:
                if not foundGood:
                    fPh = fPharr[i]
                    if (cdrArr[i] == 'Only J'):
                        jPos = posArr[i] - i
                    elif (cdrArr[i] == 'Only V'):
                        vPos = posArr[i] - i
    if ((vPos != -1) & (jPos != -1) & (not foundGood)):
        bestCDRnt = fullSeq[3 * vPos:3 * jPos]
        bestCDR = 'NA'
    elif bestPos != -1:
        bestCDRnt = bestSeq[3 * bestPos: 3 * bestPos + 3 * len(bestCDR)]
    if bestCDR.find('*') != -1:
        stat = 'Unproductive - stop codon'
    elif fPh:
        stat = 'Productive (no 118 PHE found)'
    else:
        stat = 'Productive'
    if maxLen == 0:
        if (('Only J' in cdrArr) & ('Only V' in cdrArr)):
            stat = 'Unproductive - Frame shift'
        else:
            if (('Only J' not in cdrArr) & ('Only V' in cdrArr)):
                stat = 'Unproductive - found only V segment'
            elif (('Only J' in cdrArr) & ('Only V' not in cdrArr)):
                stat = 'Unproductive - found only J segment'
            elif (('Only J' not in cdrArr) & ('Only V' not in cdrArr)):
                stat = 'Unproductive - didn\'t find V and J segment'
            else:
                stat = 'Unproductive'
            bestCDR = 'NA'
            bestCDRnt = 'NA'
    fDict['stat'] = stat
    fDict['CDR3 AA'] = bestCDR
    fDict['CDR3 NT'] = bestCDRnt
    fDict['Full Seq'] = fullSeq
    return fDict


def getCDR3(aaSeq, vSeq, jSeq):
    minDist = 14
    pos = -1
    for i in range(0, len(aaSeq) - len(vSeq) + 1):
        subAA = aaSeq[i:i + len(vSeq)]
        if len(subAA) != len(vSeq):
            sys.stderr.write(str(datetime.datetime.now()) + ' Error! Wrong sub length\n')
            sys.stderr.flush()
        dist = 0
        for k in range(0, len(vSeq)):
            if vSeq[k] != subAA[k]:
                dist += 1
        if ((dist < minDist) & (subAA.endswith('C'))):
            minDist = dist
            pos = i + len(vSeq)
    jPos = -1
    minDistJ = 4
    curPh = False
    for j in range(pos + 1, len(aaSeq) - len(jSeq) + 1):
        subAA = aaSeq[j: j + len(jSeq)]
        if len(subAA) != len(jSeq):
            sys.stderr.write(str(datetime.datetime.now()) + ' Error! Wrong subj length\n')
            sys.stderr.flush()
        dist = 0
        for m in range(0, len(jSeq)):
            if jSeq[m] != subAA[m]:
                dist += 1
        if (dist <= minDistJ):
            if isLegal(subAA):
                jPos = j
                minDistJ = dist
                curPh = False
            else:
                if dist < minDistJ:
                    curPh = True
                    jPos = j
                    minDistJ = dist
    if pos == -1:
        if jPos != -1:
            return ('Only J', jPos, curPh)
        else:
            return ('No V/J found', -1, curPh)
    else:
        if jPos == -1:
            return ('Only V', pos, curPh)
    return (aaSeq[pos:jPos], pos, curPh)


# Checks that the conserved amino acids remain
def isLegal(subAA):
    if (len(subAA) < 4):
        return False
    if (subAA[0] == 'F'):
        if ((subAA[1] == 'G') | (subAA[3] == 'G')):
            return True
    if ((subAA[1] == 'G') & (subAA[3] == 'G')):
        return True
    if ((subAA[0] == 'G') & (subAA[2] == 'G')):
        return True
    if subAA.startswith('GR'):
        return True
    if subAA.startswith('SR'):
        return True
    return False


def getBestJaa(currSeq):
    firstSeq = getNTseq(currSeq)
    secondSeq = getNTseq(currSeq[1:])
    thirdSeq = getNTseq(currSeq[2:])
    pos = 10
    seq = ''
    found = False
    for s in [firstSeq, secondSeq, thirdSeq]:
        tempSeq = s[:8]
        indF = tempSeq.find('F')
        if indF != -1:
            if ((indF + 3) <= len(s)):
                if s[indF + 3] == 'G':
                    found = True
                    if indF < pos:
                        pos = indF
                        seq = s[indF:]
        indG = tempSeq.find('G')
        if (indG != -1):
            if found == False:
                if ((indG + 2) <= len(s)):
                    if s[indG + 2] == 'G':
                        if indG < pos:
                            found = True
                            seq = s[indG:]
                            pos = indG
    if ((found == False) & (indF != -1)):
        seq = s[indF:]
    if seq != '':
        return seq
    else:
        return firstSeq


def getBestVaa(currSeq):
    firstSeq = getNTseq(currSeq)
    secondSeq = getNTseq(currSeq[1:])
    thirdSeq = getNTseq(currSeq[2:])
    pos = 10
    seq = ''
    for s in [firstSeq, secondSeq, thirdSeq]:
        # print "S: " + s
        tempSeq = s[-8:]
        # print "tempSeq: " + tempSeq
        ind = tempSeq.find('C')
        stopInd = tempSeq.find('*')
        # print "Ind: " + str(ind)
        if ((ind != -1) & (stopInd == -1)):
            # print "inside the ind"
            if ind < pos:
                goodRF = isGoodRF(s)
                if goodRF:
                    pos = ind
                    seq = s[:-8 + ind + 1]
    if seq != '':
        return seq
    else:
        return firstSeq


def isGoodRF(s):
    mInd = s.find('M')
    if mInd == -1:
        return False
    stopInd = s.find('*')
    if stopInd == -1:
        return True
    stopIndNext = s[stopInd + 1:].find('*')
    while stopIndNext != -1:
        stopInd = stopIndNext + stopInd + 1
        stopIndNext = s[stopInd + 1:].find('*')
        mInd = s[stopInd + 1:].find('M')
        mInd = mInd + stopInd + 1
    if mInd != -1:
        return True
    else:
        return False
