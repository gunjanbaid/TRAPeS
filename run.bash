#!/usr/bin/env bash

rm -r -f cells/
mkdir cells/

data_loc=/data/yosef2/Published_Data/TraCeR/proc_data_100bp/day0 
DIRS=`ls -l $data_loc | egrep '^d' | awk '{print $9}'`

for DIR in $DIRS
do
	echo ${DIR}
	mkdir cells/${DIR}
	old_path="/data/yosef2/Published_Data/TraCeR/proc_data_100bp/day0/${DIR}"
	new_path="cells/${DIR}"
	samtools view -h -b -s 0.5 ${old_path}/tophat_output/picard_output/sorted.bam  > ${new_path}/sorted.bam
	samtools view -h -b -s 0.5 ${old_path}/tophat_output/unmapped.bam  > ${new_path}/unmapped.bam
        # cp ${old_path}/tophat_output/picard_output/sorted.bam ${new_path}/sorted.bam
        # cp ${old_path}/tophat_output/unmapped.bam ${new_path}/unmapped.bam
	samtools bam2fq ${new_path}/sorted.bam > ${new_path}/sorted.fq
        samtools bam2fq ${new_path}/unmapped.bam > ${new_path}/unmapped.fq
	bowtie2 -q --phred33 -x Data/mm10_ncbi/index/mm10 -U ${new_path}/sorted.fq -S ${new_path}/sorted.sam
        bowtie2 -q --phred33 -x Data/mm10_ncbi/index/mm10 -U ${new_path}/unmapped.fq -S ${new_path}/unmapped.sam
        rm ${new_path}/sorted.bam ${new_path}/unmapped.bam        
        samtools view -h -b ${new_path}/sorted.sam | samtools sort > ${new_path}/sorted.bam
        samtools view -h -b ${new_path}/unmapped.sam | samtools sort > ${new_path}/unmapped.bam
        samtools index ${new_path}/sorted.bam
        if [ "${DIR}" == "ERR1146638" ]; then
        	break
	fi
done

python trapes.py -path cells/ -bam sorted.bam -unmapped unmapped.bam -output output -sumF cells/ -genome mm10_ncbi -score 33 -trim 40

