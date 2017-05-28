#!/bin/bash

source /home/eecs/gunjan_baid/.conda/envs/my_root/bin/activate my_root

data_loc=/data/yosef2/Published_Data/TraCeR/proc_data_100bp/day0
DIRS=`ls -l $data_loc | egrep '^d' | awk '{print $9}'`
OUT_DIR=$1

rm -r -f /home/eecs/gunjan_baid/trapes/${OUT_DIR}/
mkdir /home/eecs/gunjan_baid/trapes/${OUT_DIR}/
mkdir /home/eecs/gunjan_baid/trapes/${OUT_DIR}/cells
cp -r cleanProcessedData/* /home/eecs/gunjan_baid/trapes/${OUT_DIR}/cells
cp /home/eecs/gunjan_baid/trapes/bashScripts/run-clean-data.bash /home/eecs/gunjan_baid/trapes/${OUT_DIR}/
cp /home/eecs/gunjan_baid/trapes/trapes.py /home/eecs/gunjan_baid/trapes/${OUT_DIR}/
cp /home/eecs/gunjan_baid/trapes/vdj.alignment /home/eecs/gunjan_baid/trapes/${OUT_DIR}/
cp -r /home/eecs/gunjan_baid/trapes/Data /home/eecs/gunjan_baid/trapes/${OUT_DIR}

python /home/eecs/gunjan_baid/trapes/${OUT_DIR}/trapes.py -path /home/eecs/gunjan_baid/trapes/${OUT_DIR}/cells/ -bam sorted.bam -unmapped unmapped.bam -output output -sumF /home/eecs/gunjan_baid/trapes/${OUT_DIR}/summary -genome mm10_ncbi -trim 30 -score 35 -lowQ
