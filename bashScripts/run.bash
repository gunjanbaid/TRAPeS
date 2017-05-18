#!/bin/bash
#PBS -l walltime=48:00:00

source /home/eecs/gunjan_baid/.conda/envs/my_root/bin/activate my_root

rm -r -f /home/eecs/gunjan_baid/trapes/cells/
mkdir /home/eecs/gunjan_baid/trapes/cells/

data_loc=/data/yosef2/Published_Data/TraCeR/proc_data_100bp/day0 
DIRS=`ls -l $data_loc | egrep '^d' | awk '{print $9}'`

for DIR in $DIRS
do
	echo ${DIR}
	mkdir cells/${DIR}
	old_path="/data/yosef2/Published_Data/TraCeR/proc_data_100bp/day0/${DIR}"
	new_path="/home/eecs/gunjan_baid/trapes/cells/${DIR}"
        # add regions for the chromosomes
        samtools view -h ${old_path}/tophat_output/picard_output/sorted.bam NC_000080.6 NC_000072.6 > ${new_path}/sorted.bam
        cp ${old_path}/tophat_output/unmapped.bam ${new_path}/unmapped.bam
        samtools bam2fq ${new_path}/sorted.bam > ${new_path}/output.sorted.fq
        samtools bam2fq ${new_path}/unmapped.bam > ${new_path}/output.unmapped.fq
	bowtie2 -q --phred33 -x Data/mm10_ncbi/index/mm10 -U ${new_path}/output.sorted.fq -S ${new_path}/temp1.sam
        # trim sorted reads
        bowtie2 -q --phred33 -x Data/mm10_ncbi/index/mm10 -U ${new_path}/output.sorted.fq -S ${new_path}/temp2.sam --trim3 80
        bowtie2 -q --phred33 -x Data/mm10_ncbi/index/mm10 -U ${new_path}/output.sorted.fq -S ${new_path}/temp3.sam --trim5 80
        bowtie2 -q --phred33 -x Data/mm10_ncbi/index/mm10 -U ${new_path}/output.sorted.fq -S ${new_path}/temp4.sam --trim3 40 --trim5 40
        #bowtie2 -q --phred33 -x Data/mm10_ncbi/index/mm10 -U ${new_path}/output.unmapped.fq -S ${new_path}/temp2.sam --trim3 40
        #bowtie2 -q --phred33 -x Data/mm10_ncbi/index/mm10 -U ${new_path}/output.unmapped.fq -S ${new_path}/temp3.sam --trim5 40
        samtools merge ${new_path}/all.bam ${new_path}/temp1.sam ${new_path}/temp2.sam ${new_path}/temp3.sam ${new_path}/temp4.sam
        rm ${new_path}/sorted.bam ${new_path}/temp1.sam ${new_path}/temp2.sam ${new_path}/temp3.sam ${new_path}/temp4.sam
	samtools view -h -b -F 4 ${new_path}/all.bam | samtools sort > ${new_path}/sorted.bam
	# samtools view -h -b -f 4 ${new_path}/all.bam | samtools sort > ${new_path}/unmapped.bam
        samtools index ${new_path}/sorted.bam
        #if [ ${DIR} == "ERR1146638" ]; then
	#	break	
	#fi
done

python /home/eecs/gunjan_baid/trapes/trapes.py -path /home/eecs/gunjan_baid/trapes/cells/ -bam sorted.bam -unmapped unmapped.bam -output output -sumF /home/eecs/gunjan_baid/trapes/cells/summary -genome mm10_ncbi -trim 50 


