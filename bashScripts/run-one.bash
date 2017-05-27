#!/bin/bash

source /home/eecs/gunjan_baid/.conda/envs/my_root/bin/activate my_root

data_loc=/data/yosef2/Published_Data/TraCeR/proc_data_100bp/day0
DIRS=`ls -l $data_loc | egrep '^d' | awk '{print $9}'`
CELL_DIR=$1
OUT_DIR=$2

rm -r -f /home/eecs/gunjan_baid/trapes/${OUT_DIR}/
mkdir /home/eecs/gunjan_baid/trapes/${OUT_DIR}/

echo ${CELL_DIR}
mkdir ${OUT_DIR}/${CELL_DIR}
old_path="/data/yosef2/Published_Data/TraCeR/proc_data_100bp/day0/${CELL_DIR}"
new_path="/home/eecs/gunjan_baid/trapes/${OUT_DIR}/${CELL_DIR}"
# add regions for the chromosomes
samtools view -h -s 0.5 ${old_path}/tophat_output/picard_output/sorted.bam NC_000080.6 NC_000072.6 > ${new_path}/sorted.bam
samtools view -h -s 0.5 ${old_path}/tophat_output/unmapped.bam > ${new_path}/unmapped.bam
samtools bam2fq ${new_path}/sorted.bam > ${new_path}/output.sorted.fq
python pythonScripts/helper.py ${new_path}/sorted.bam ${new_path}/sortedTEMP.bam
python pythonScripts/helper.py ${new_path}/unmapped.bam ${new_path}/unmappedTEMP.bam
rm ${new_path}/sorted.bam ${new_path}/unmapped.bam
mv ${new_path}/sortedTEMP.bam ${new_path}/sorted.bam
mv ${new_path}/unmappedTEMP.bam ${new_path}/unmapped.bam 
bowtie2 -q --phred33 -x Data/mm10_ncbi/index/mm10 -U ${new_path}/output.sorted.fq -S ${new_path}/temp1.sam
# trim sorted reads
bowtie2 -q --phred33 -x Data/mm10_ncbi/index/mm10 -U ${new_path}/output.sorted.fq -S ${new_path}/temp2.sam --trim3 65
bowtie2 -q --phred33 -x Data/mm10_ncbi/index/mm10 -U ${new_path}/output.sorted.fq -S ${new_path}/temp3.sam --trim5 65
bowtie2 -q --phred33 -x Data/mm10_ncbi/index/mm10 -U ${new_path}/output.sorted.fq -S ${new_path}/temp4.sam --trim3 32 --trim5 32
samtools merge ${new_path}/all.bam ${new_path}/temp1.sam ${new_path}/temp2.sam ${new_path}/temp3.sam ${new_path}/temp4.sam
rm ${new_path}/sorted.bam ${new_path}/temp1.sam ${new_path}/temp2.sam ${new_path}/temp3.sam ${new_path}/temp4.sam
samtools view -h -b -F 4 ${new_path}/all.bam | samtools sort > ${new_path}/sorted.bam
# samtools view -h -b -f 4 ${new_path}/all.bam | samtools sort > ${new_path}/unmapped.bam
samtools index ${new_path}/sorted.bam
rm ${new_path}/output.sorted.fq
samtools bam2fq ${new_path}/sorted.bam > ${new_path}/output.sorted.fq
samtools bam2fq ${new_path}/unmapped.bam > ${new_path}/output.unmapped.fq
#if [ ${CELL_DIR} == "ERR1146638" ]; then
#	break	
#fi

python /home/eecs/gunjan_baid/trapes/trapes.py -path /home/eecs/gunjan_baid/trapes/${OUT_DIR}/ -bam sorted.bam -unmapped unmapped.bam -output output -sumF /home/eecs/gunjan_baid/trapes/${OUT_DIR}/summary -genome mm10_ncbi -trim 30 -score 28 -lowQ
