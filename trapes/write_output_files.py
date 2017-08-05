"""
This module contains functions used to write output files containing summary statistics for individual cells and
summary statistics for all cells.

"""

import os


def makeSingleCellOutputFile(alphaDict, betaDict, output, betaRsem, alphaRsem, alphaBam, betaBam, fastaDict,
                             unDictAlpha, unDictBeta, idNameDict):
    outF = open(output + '.summary.txt', 'w')
    outF.write(
        'Chain\tStatus\tRank of TCR\tV\tJ\tC\tCDR3 NT\tCDR3 AA\t#reads in TCR\t#reads in CDR3\t#reads in V\t#reads in '
        'J\t#reads in C\t%unmapped reads used in the reconstruction\t# unmapped reads used in the '
        'reconstruction\t%unmapped reads in CDR3\t#unmapped reads in CDR3\tV ID\tJ ID\tC ID\n')
    if (len(alphaDict) > 0):
        writeChain(outF, 'alpha', alphaDict, alphaRsem, alphaBam, fastaDict, unDictAlpha, output, idNameDict)
    if (len(betaDict) > 0):
        writeChain(outF, 'beta', betaDict, betaRsem, betaBam, fastaDict, unDictBeta, output, idNameDict)
    outF.close()


def writeChain(outF, chain, cdrDict, rsemF, bamF, fastaDict, unDict, output, idNameDict):
    writtenArr = []
    if os.path.exists(rsemF):
        noRsem = False
        (rsemDict, unRsemDict) = makeRsemDict(rsemF, cdrDict)
    else:
        noRsem = True
    for tcr in cdrDict:
        jStart = -1
        cdrInd = -1
        cInd = -1
        if cdrDict[tcr]['stat'] == 'Productive':
            isProd = True
        else:
            isProd = False
        fLine = chain + '\t' + cdrDict[tcr]['stat'] + '\t'
        if noRsem:
            rank = 'NA'
        else:
            rank = getRank(tcr, rsemDict, unRsemDict, isProd, noRsem)
        fLine += str(rank) + '\t'
        nameArr = tcr.split('.')
        fLine += nameArr[0] + '\t' + nameArr[1] + '\t' + nameArr[2] + '\t'
        fLine += cdrDict[tcr]['CDR3 NT'] + '\t' + cdrDict[tcr]['CDR3 AA'] + '\t'
        fullSeq = cdrDict[tcr]['Full Seq'].upper()
        if not noRsem:
            totalCount = findCountsInRegion(bamF, 0, len(fullSeq), tcr)
            fLine += str(totalCount) + '\t'
            cName = nameArr[5]
            while cName.endswith('_2'):
                cName = cName[:-2]
            cSeq = fastaDict[cName].upper()
            cInd = fullSeq.find(cSeq)
            if cInd == -1:
                sys.stderr.write(
                    str(datetime.datetime.now()) + 'Error! could not find C segment sequence in the full sequence\n')
                sys.stderr.flush()
                cCounts = 'NA'
            else:
                cCounts = findCountsInRegion(bamF, cInd, len(fullSeq), tcr)
            if cdrDict[tcr]['CDR3 NT'] != 'NA':
                cdrInd = fullSeq.find(cdrDict[tcr]['CDR3 NT'].upper())
            else:
                cdrInd = -1
            if ((cdrInd == -1) & (cdrDict[tcr]['CDR3 NT'] != 'NA')):
                sys.stderr.write(
                    str(datetime.datetime.now()) + ' Error! Cound not find CDR3 NT sequence in the full sequence\n')
                sys.stderr.flush()
            if cdrInd != -1:
                cdrCounts = findCountsInRegion(bamF, cdrInd, cdrInd + len(cdrDict[tcr]['CDR3 NT']), tcr)
                jStart = cdrInd + len(cdrDict[tcr]['CDR3 NT'])
                if cInd != -1:
                    jCounts = findCountsInRegion(bamF, jStart, cInd, tcr)
                else:
                    jCounts = findCountsInRegion(bamF, jStart, jStart + 50, tcr)
                vCounts = findCountsInRegion(bamF, 0, cdrInd, tcr)
                fLine += str(cdrCounts) + '\t' + str(vCounts) + '\t' + str(jCounts) + '\t' + str(cCounts) + '\t'
            else:
                fLine += 'NA\tNA\tNA\t' + str(cCounts) + '\t'
            vId = nameArr[3]
            jId = nameArr[4]
            cId = nameArr[5]
            if cdrDict[tcr]['CDR3 NT'] != 'NA':
                (unDictRatioCDR, unCDRcount) = getUnDictRatio(bamF, cdrInd, cdrInd + len(cdrDict[tcr]['CDR3 NT']), tcr,
                                                              unDict)
                (unDictRatioALL, unAllcount) = getUnDictRatio(bamF, 0, len(fullSeq), tcr, unDict)
                fLine += str(unDictRatioALL) + '\t' + str(unAllcount) + '\t' + str(unDictRatioCDR) + '\t' + str(
                    unCDRcount) + '\t'
            else:
                fLine += 'NA\tNA\tNA\tNA\t'
            writtenArr.append(vId)
            writtenArr.append(jId)
            fLine += vId + '\t' + jId + '\t' + cId + '\n'
        else:
            fLine += 'NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t' + nameArr[3] + '\t' + nameArr[4] + '\t' + nameArr[
                5] + '\n'

            # print fLine

        outF.write(str(fLine))
    writeFailedReconstructions(outF, chain, writtenArr, output, idNameDict, fastaDict)


def makeRsemDict(rsemF, cdrDict):
    fDict = dict()
    unDict = dict()
    f = open(rsemF, 'r')
    f.readline()
    l = f.readline()
    while l != '':
        lArr = l.strip('\n').split('\t')
        name = lArr[1]
        if name in cdrDict:
            if cdrDict[name]['stat'] == 'Productive':
                fDict[name] = float(lArr[4])
            unDict[name] = float(lArr[4])
        l = f.readline()
    f.close()
    return (fDict, unDict)


def getRank(tcr, rsemDict, unRsemDict, isProd, noRsem):
    if isProd:
        currDict = rsemDict
    else:
        currDict = unRsemDict
    if not noRsem:
        currCount = currDict[tcr]
        rank = 1
        for rec in currDict:
            if rec != tcr:
                if unRsemDict[rec] > currCount:
                    rank += 1
        return rank
    else:
        return 'NA'


def findCountsInRegion(bamF, start, end, tcr):
    readsArr = []
    mappedFile = pysam.AlignmentFile(bamF, "rb")
    readsIter = mappedFile.fetch(tcr, start, end)
    for read in readsIter:
        if read.is_read1:
            newName = read.query_name + '_1'
        else:
            newName = read.query_name + '_2'
        if newName not in readsArr:
            readsArr.append(newName)
    mappedFile.close()
    counts = len(readsArr)
    return counts


def getUnDictRatio(bamF, start, end, tcr, unDict):
    unMappedCount = 0
    usedArr = []
    mappedFile = pysam.AlignmentFile(bamF, "rb")
    readsIter = mappedFile.fetch(tcr, start, end)
    for read in readsIter:
        if read.is_read1:
            newName = read.query_name + '_1'
        else:
            newName = read.query_name + '_2'
        if newName not in usedArr:
            usedArr.append(newName)
            if newName in unDict:
                unMappedCount += 1
    mappedFile.close()
    return (float(float(unMappedCount) / len(unDict)), unMappedCount)


def writeFailedReconstructions(outF, chain, writtenArr, output, idNameDict, fastaDict):
    recF = output + '.reconstructed.junctions.' + chain + '.fa'
    if os.path.isfile(recF):
        f = open(recF, 'rU')
        segDict = dict()
        for tcrRecord in SeqIO.parse(f, 'fasta'):
            tcrSeq = str(tcrRecord.seq)
            if tcrSeq.find('NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN') != -1:
                status = 'Failed reconstruction - reached maximum number of iterations'
                segDict = addSegmentsToDict(segDict, status, writtenArr, tcrRecord, idNameDict, fastaDict)
            elif tcrSeq.find('NNNN') != -1:
                status = 'Failed reconstruction - V and J segment do not overlap'
                segDict = addSegmentsToDict(segDict, status, writtenArr, tcrRecord, idNameDict, fastaDict)
        f.close()
        if len(segDict) > 0:
            writeSegDict(segDict, outF, chain)


def addSegmentsToDict(segDict, status, writtenArr, tcrRecord, idNameDict, fastaDict):
    head = tcrRecord.id
    headArr = head.split('.')
    vId = headArr[0]
    jId = headArr[1].split('(')[0]
    currSeqArr = tcrRecord.seq.split('N')
    vSeq = currSeqArr[0]
    jSeq = currSeqArr[-1]
    minLen = min(len(fastaDict[vId]), len(fastaDict[jId]))
    tupArr = [(vId, vSeq), (jId, jSeq)]
    for i in range(0, len(tupArr)):
        (id, seq) = tupArr[i]
        if id not in writtenArr:
            if id in segDict:
                if i == 0:
                    if idNameDict[jId] not in segDict[id]['pairs']:
                        segDict[id]['pairs'].append(idNameDict[jId])
                    if str(segDict[id]['seq'][-20:]) != str(seq[-20:]):
                        sys.stderr.write(str(
                            datetime.datetime.now()) + ' Error! reconstructed two different sequences from the same '
                                                       'V-segment %s\n' % id)
                        sys.stderr.flush()
                else:
                    if idNameDict[vId] not in segDict[id]['pairs']:
                        segDict[id]['pairs'].append(idNameDict[vId])
                    if str(segDict[id]['seq'][:20]) != str(seq[:20]):
                        sys.stderr.write(str(
                            datetime.datetime.now()) + ' Error! reconstructed two different sequences from the same '
                                                       'J-segment %s\n' % id)
                        sys.stderr.flush()
            else:
                segDict[id] = dict()
                segDict[id]['status'] = status
                segDict[id]['seq'] = seq
                segDict[id]['len'] = len(seq) - minLen
                segDict[id]['pairs'] = []

                if i == 0:
                    segDict[id]['type'] = 'V'
                    segDict[id]['pairs'].append(idNameDict[jId])
                else:
                    segDict[id]['type'] = 'J'
                    segDict[id]['pairs'].append(idNameDict[vId])
                segDict[id]['name'] = idNameDict[id]

    return segDict


def writeSegDict(segDict, outF, chain):
    for seg in segDict:
        currDict = segDict[seg]
        pairs = ''
        for pair in currDict['pairs']:
            pairs += pair + '.'
        pairs = pairs[:-1]
        if currDict['len'] > 0:
            fLine = chain + '\t' + currDict['status'] + '\t'
            rank = findCurrRank(segDict, seg, currDict['len'])
            fLine += str(rank) + '\t'
            if currDict['type'] == 'V':
                fLine += currDict['name'] + '\t' + 'paired with: ' + pairs + '\t'
            else:
                fLine += 'paired with: ' + pairs + '\t' + currDict['name'] + '\t'
            fLine += 'NA\t' + currDict['seq'] + '\tNA\tNA\tNA\t'
            fLine += 'NA\tNA\tNA\tNA\tNA\tNA\tNA\t'
            if currDict['type'] == 'V':
                fLine += seg + '\tNA\tNA\n'
            else:
                fLine += 'NA\t' + seg + '\tNA\n'
                # print fLine
            outF.write(str(fLine))


def findCurrRank(segDict, seg, currLen):
    rank = 1
    for s in segDict:
        if s != seg:
            if segDict[s]['len'] > currLen:
                rank += 1
    return rank


def addCellToTCRsum(cellFolder, noutput, opened, tcrFout):
    if os.path.isfile(noutput + '.summary.txt'):
        currOut = open(noutput + '.summary.txt', 'r')
        if not opened:
            opened = True
            head = currOut.readline()
            head = 'cell\t' + head
            tcrFout.write(head)
        else:
            currOut.readline()
        l = currOut.readline()
        while l != '':
            newL = cellFolder + '\t' + l
            tcrFout.write(newL)
            l = currOut.readline()
        currOut.close()
    return opened


def addToStatDict(noutput, cellFolder, finalStatDict):
    if cellFolder in finalStatDict:
        print "Error! %s appear more than once in final stat dictionary" % cellFolder
    finalStatDict[cellFolder] = {'alpha': 'Failed - found V and J segments but wasn\'t able to extend them',
                                 'beta': 'Failed - found V and J segments but wasn\'t able to extend them'}
    if os.path.isfile(noutput + '.summary.txt'):
        currOut = open(noutput + '.summary.txt', 'r')
        msgA = 'None'
        msgB = 'None'
        currOut.readline()
        l = currOut.readline()
        while l != '':
            lArr = l.strip('\n').split('\t')
            chain = lArr[0]
            stat = lArr[1]
            if stat == 'Productive':
                if chain == 'alpha':
                    msgA = 'Productive'
                else:
                    msgB = 'Productive'
            elif stat == 'Productive (no 118 PHE found)':
                if chain == 'alpha':
                    msgA = 'Productive (no 118 PHE found)'
                else:
                    msgB = 'Productive (no 118 PHE found)'
            elif stat.startswith('Unproductive'):
                if chain == 'alpha':
                    if msgA != 'Productive':
                        msgA = 'Unproductive'
                else:
                    if msgB != 'Productive':
                        msgB = 'Unproductive'
            elif stat.startswith('Failed reconstruction'):
                if stat == 'Failed reconstruction - reached maximum number of iterations':
                    if chain == 'alpha':
                        if msgA == 'None':
                            msgA = 'Failed - reconstruction didn\'t converge'
                    else:
                        if msgB == 'None':
                            msgB = 'Failed - reconstruction didn\'t converge'
                elif stat == 'Failed reconstruction - V and J segment do not overlap':
                    if chain == 'alpha':
                        if msgA == 'None':
                            msgA = 'Failed - V and J reconstruction don\'t overlap'
                    else:
                        if msgB == 'None':
                            msgB = 'Failed - V and J reconstruction don\'t overlap'
            l = currOut.readline()
        currOut.close()
        if msgA == 'None':
            alphaJunc = noutput + '.alpha.junctions.txt'
            if (os.path.isfile(alphaJunc) == True):
                if os.stat(alphaJunc).st_size == 0:
                    msgA = 'Failed - didn\'t find any V and J segments in original mapping'
                else:
                    msgA = 'Failed - found V and J segments but wasn\'t able to extend them'
            else:
                msgA = 'Failed - didn\'t find any V and J segments in original mapping'
        if msgB == 'None':
            betaJunc = noutput + '.beta.junctions.txt'
            if (os.path.isfile(betaJunc) == True):
                if os.stat(betaJunc).st_size == 0:
                    msgB = 'Failed - didn\'t find any V and J segments in original mapping'
                else:
                    msgB = 'Failed - found V and J segments but wasn\'t able to extend them'
            else:
                msgB = 'Failed - didn\'t find any V and J segments in original mapping'

    else:
        betaJunc = noutput + '.beta.junctions.txt'
        alphaJunc = noutput + '.alpha.junctions.txt'
        if os.path.isfile(betaJunc) == True:
            if os.stat(betaJunc).st_size == 0:
                msgB = 'Failed - didn\'t find any V and J segments in original mapping'
        else:
            msgB = 'Failed - didn\'t find any V and J segments in original mapping'
        if (os.path.isfile(alphaJunc) == True):
            if os.stat(alphaJunc).st_size == 0:
                msgA = 'Failed - didn\'t find any V and J segments in original mapping'
        else:
            msgA = 'Failed - didn\'t find any V and J segments in original mapping'
    finalStatDict[cellFolder]['alpha'] = msgA
    finalStatDict[cellFolder]['beta'] = msgB
    return finalStatDict
