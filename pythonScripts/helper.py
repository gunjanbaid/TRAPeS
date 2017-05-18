#! /bin/env python
import sys
import os
import argparse
import subprocess
import datetime
from Bio import SeqIO
import pysam
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

def check_mapped(file):
	total = 0
	mate_unmapped = 0
	mate_mapped = 0
	mappedFile = pysam.AlignmentFile(file, "rb")
	readsIter = mappedFile.fetch(until_eof = True)
	for read in readsIter:
		if not read.is_unmapped:
			print("ahhhh")
		total += 1
		if read.mate_is_unmapped:
			mate_unmapped += 1
		else:
			mate_mapped += 1
	print("mate unmapped for {} reads".format(mate_unmapped))
	print("mate mapped for {} reads".format(mate_mapped))
	print("total is " + str(total))


def check_read_both(file1, file2):
	file_1_reads = []
	in_both = 0

	mappedFile = pysam.AlignmentFile(file1, "rb")
	readsIter = mappedFile.fetch(until_eof = True)	
	for read in readsIter:
		if read.mate_is_unmapped:
			file_1_reads.append(read)

	mappedFile2 = pysam.AlignmentFile(file2, "rb")
	readsIter2 = mappedFile.fetch(until_eof = True)	
	for read in readsIter2:
		if read in file_1_reads:
			in_both += 1
	print("{} reads with an unmapped mate in both files".format(in_both))


def get_pairs(file1, file2):
	qnames = {}
	qname_count = 0
	count = 0
	mappedFile = pysam.AlignmentFile(file1, "rb")
	readsIter = mappedFile.fetch(until_eof = True)	
	for read in readsIter:
		if read.query_name not in qnames:
			qnames[read.query_name] = 1
		else:
			qnames[read.query_name] = qnames[read.query_name] + 1

	mappedFile = pysam.AlignmentFile(file2, "rb")
	readsIter = mappedFile.fetch(until_eof = True)	
	for read in readsIter:
		if read.query_name in qnames:
			qnames[read.query_name] = qnames[read.query_name] + 1

	print(qnames.values())

# read will have \1 or \2, so tophat will get rid of tophat
# so these should be indentical 

def remove_poly_a(filename, output):
    f = pysam.AlignmentFile(filename, "rb")
    readsIter = f.fetch(until_eof=True)
    cleanedReads = pysam.AlignmentFile(output, "wb", template=f)
    tt = "T"*30
    aa = "A"*30
    for read in readsIter:
    	if aa in read.query_sequence or tt in read.query_sequence:
    		print(read.query_sequence)
        else:
            cleanedReads.write(read)
    f.close()
    cleanedReads.close()

if __name__ == '__main__':
        filename = sys.argv[1]
	output = sys.argv[2]
	remove_poly_a(filename, output)

