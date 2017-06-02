#!/bin/bash

source /home/eecs/gunjan_baid/.conda/envs/my_root/bin/activate my_root

data_loc=/data/yosef2/Published_Data/TraCeR/proc_data_100bp/day0
#data_loc=/data/yosef2/Published_Data/TraCeR/proc_data_100bp/day49
DIRS=`ls -l $data_loc | egrep '^d' | awk '{print $9}'`
OUT_DIR=$1

rm -r -f /home/eecs/gunjan_baid/trapes/${OUT_DIR}/
mkdir /home/eecs/gunjan_baid/trapes/${OUT_DIR}/
mkdir /home/eecs/gunjan_baid/trapes/${OUT_DIR}/cells
cp /home/eecs/gunjan_baid/trapes/bashScripts/run-all-simple.bash /home/eecs/gunjan_baid/trapes/${OUT_DIR}/
cp /home/eecs/gunjan_baid/trapes/trapes.py /home/eecs/gunjan_baid/trapes/${OUT_DIR}/
cp /home/eecs/gunjan_baid/trapes/vdj.alignment /home/eecs/gunjan_baid/trapes/${OUT_DIR}/
cp -r /home/eecs/gunjan_baid/trapes/Data /home/eecs/gunjan_baid/trapes/${OUT_DIR}

for CELL_DIR in $DIRS
do
        echo ${CELL_DIR}
        mkdir /home/eecs/gunjan_baid/trapes/${OUT_DIR}/cells/${CELL_DIR}
        old_path="/data/yosef2/Published_Data/TraCeR/proc_data_100bp/day0/${CELL_DIR}"
        new_path="/home/eecs/gunjan_baid/trapes/${OUT_DIR}/cells/${CELL_DIR}"
        # add regions for the chromosomes
        samtools view -h ${old_path}/tophat_output/picard_output/sorted.bam NC_000080.6 NC_000072.6 > ${new_path}/sorted.bam
        samtools view -h ${old_path}/tophat_output/unmapped.bam > ${new_path}/unmapped.bam
        samtools bam2fq ${new_path}/sorted.bam > ${new_path}/output.sorted.fq
        samtools bam2fq ${new_path}/unmapped.bam > ${new_path}/output.unmapped.fq        
        samtools sort ${new_path}/sorted.bam -o ${new_path}/sorted.bam
        samtools index ${new_path}/sorted.bam
        if [ ${CELL_DIR} == "ERR1146641" ]; then
               break   
        fi
done

python /home/eecs/gunjan_baid/trapes/${OUT_DIR}/trapes.py -path /home/eecs/gunjan_baid/trapes/${OUT_DIR}/cells/ -bam sorted.bam -unmapped unmapped.bam -output output -sumF /home/eecs/gunjan_baid/trapes/${OUT_DIR}/summary -genome mm10_ncbi -trim 70 -score 75 -lowQ -bases 200
