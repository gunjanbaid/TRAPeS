#! /usr/bin/env python

"""
This module contains utility functions for various tasks.

"""


def formatFiles(fullPath, bam, unmapped, output):
    found = True
    nbam = fullPath + bam
    if bam.startswith('/'):
        nbam = fullPath + bam[1:]
    elif bam.startswith('./'):
        nbam = fullPath + bam[2:]
    nunmapped = fullPath + unmapped
    if unmapped.startswith('/'):
        nunmapped = fullPath + unmapped[1:]
    if unmapped.startswith('./'):
        nunmapped = fullPath + unmapped[2:]
    if ((os.path.isfile(nunmapped)) & (os.path.isfile(nbam))):
        noutput = makeOutputDir(output, fullPath)
    else:
        noutput = output
        found = False
    return (found, nbam, nunmapped, noutput)


def makeOutputDir(output, fullPath):
    noutput = output
    if output.startswith('/'):
        noutput = output[1:]
    if output.startswith('./'):
        noutput = output[2:]
    if output.endswith('/'):
        noutput = output[:-1]
    if output.find('/') != -1:
        outArr = noutput.split('/')
        currPath = fullPath
        for i in range(0, len(outArr) - 1):
            currPath = currPath + outArr[i] + '/'
            if not os.path.exists(currPath):
                os.makedirs(currPath)
    noutput = fullPath + noutput
    return noutput


def findSeqAndLengthOfAA(aaSeq):
    fLen = 0
    fSeq = ''
    startArr = []
    stopArr = []
    startM = aaSeq.find('M')
    while startM != -1:
        startArr.append(startM)
        startM = aaSeq.find('M', startM + 1)
    stopPos = aaSeq.find('*')
    while stopPos != -1:
        stopArr.append(stopPos)
        stopPos = aaSeq.find('*', stopPos + 1)

    if ((len(startArr) == 0) | (len(stopArr) == 0)):
        return (fSeq, fLen)
    for stP in startArr:
        currStop = findStop(stP, stopArr)
        if currStop == -1:
            return (fSeq, fLen)
        else:
            currLen = currStop - stP
            if currLen >= fLen:
                fLen = currLen
                fSeq = aaSeq[stP:currStop]
    return fSeq


def findStop(stP, stopArr):
    for x in stopArr:
        if x > stP:
            return x
    return -1
