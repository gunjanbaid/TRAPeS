#!/usr/bin/env bash

for d in paired_100_cells/* ; do
    # echo "$d"
    # mkdir $d/temp
    # mv $d/* $d/temp
    # ls $d/temp
    rm $d/*
    rm -r $d/rsem*
    samtools view -h -b -s 0.1 $d/temp/sorted.bam > $d/sorted.bam
    samtools view -h -b -s 0.1 $d/temp/unmapped.bam > $d/unmapped.bam
    samtools index $d/sorted.bam
done

python trapes.py -path paired_100_cells/ -bam sorted.bam -unmapped unmapped.bam -output output -sumF paired_100_cells/output -genome mm10_ncbi -score 25 -lowQ


for d in paired_100_cells/* ; do
    # echo "$d"
    # mkdir $d/temp
    # mv $d/* $d/temp
    # ls $d/temp
    rm $d/*
    rm -r $d/rsem*
    samtools view -h -b -s 0.1 $d/temp/sorted.bam > $d/sorted.bam
    samtools view -h -b -s 0.1 $d/temp/unmapped.bam > $d/unmapped.bam
    samtools index $d/sorted.bam
done

python trapes.py -path paired_100_cells/ -bam sorted.bam -unmapped unmapped.bam -output output -sumF paired_100_cells/output -genome mm10_ncbi -score 35 -lowQ


for d in paired_100_cells/* ; do
    # echo "$d"
    # mkdir $d/temp
    # mv $d/* $d/temp
    # ls $d/temp
    rm $d/*
    rm -r $d/rsem*
    samtools view -h -b -s 0.1 $d/temp/sorted.bam > $d/sorted.bam
    samtools view -h -b -s 0.1 $d/temp/unmapped.bam > $d/unmapped.bam
    samtools index $d/sorted.bam
done

python trapes.py -path paired_100_cells/ -bam sorted.bam -unmapped unmapped.bam -output output -sumF paired_100_cells/output -genome mm10_ncbi -score 45 -lowQ


for d in paired_100_cells/* ; do
    # echo "$d"
    # mkdir $d/temp
    # mv $d/* $d/temp
    # ls $d/temp
    rm $d/*
    rm -r $d/rsem*
    samtools view -h -b -s 0.1 $d/temp/sorted.bam > $d/sorted.bam
    samtools view -h -b -s 0.1 $d/temp/unmapped.bam > $d/unmapped.bam
    samtools index $d/sorted.bam
done

python trapes.py -path paired_100_cells/ -bam sorted.bam -unmapped unmapped.bam -output output -sumF paired_100_cells/output -genome mm10_ncbi -score 55 -lowQ

