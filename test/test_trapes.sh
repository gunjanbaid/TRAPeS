#!/usr/bin/env sh

# set up test files, create two copies of data, one for testing and one for
rm -r -f test_output
mkdir test_output
cp -r test_data test_output/test_data_old
cp -r test_data test_output/test_data_new

# test paired-end data
python trapes.py -path test_output/test_data_old/paired_end -bam sorted.bam -unmapped unmapped.bam -output \
output -sumF test_output/summary -genome hg38 > test_output/output_old.out

python ../trapes/trapes.py -path test_output/test_data_new/paired_end -bam sorted.bam -unmapped unmapped.bam -output \
output -sum_f test_output/summary -genome hg38 > test_output/output_new.out

FILES="output.alpha.junctions.txt
output.alpha.mapped.and.unmapped.fa
output.beta.junctions.txt
output.beta.mapped.and.unmapped.fa
output.reconstructed.junctions.alpha.fa
output.reconstructed.junctions.beta.fa"

echo "Comparing old and new single end output"

for f in $FILES
do
	diff test_output/test_data_old/paired_end/cell1/$f test_output/test_data_new/paired_end/cell1/$f
	diff test_output/test_data_old/paired_end/cell2/$f test_output/test_data_new/paired_end/cell2/$f
done
