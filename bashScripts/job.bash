#!/bin/bash
#PBS -l walltime=48:00:00

source /home/eecs/gunjan_baid/.conda/envs/my_root/bin/activate my_root

nohup python /home/eecs/gunjan_baid/trapes/trapes.py -path /home/eecs/gunjan_baid/trapes/cells/ -bam sorted.bam -unmapped unmapped.bam -output output -sumF /home/eecs/gunjan_baid/trapes/cells/summary -genome mm10_ncbi -score 45 -trim 40 -bases 200

