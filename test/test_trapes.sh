#!/usr/bin/env sh

# set up test files, create two copies of data, one for testing and one for
rm -r -f test_output
mkdir test_output
cp -r test_data test_output/test_data_old
cp -r test_data test_output/test_data_new

# test single-end data

# run last working version of trapes
./trapes -path test_output/test_data_old/single_end -bam sorted.bam -unmapped unmapped.bam -output output -sumF \
summary

# run updated trapes
../trapes/trapes.py -path test_output/test_data_new/single_end -bam sorted.bam -unmapped unmapped.bam -output output \
 -sumF summary

# diff single-end outputs
echo "Single-end Diffs"
echo "Differences between old and new summary.TCR.txt"
diff test_output/test_data_old/single_end/summary.TCRs.txt test_output/test_data_new/single_end/summary.TCRs.txt
echo "Differences between old and new summary.summary.txt"
diff test_output/test_data_old/single_end/summary.summary.txt test_output/test_data_new/single_end/summary.summary.txt


# test paired-end data
./trapes -path test_output/test_data_old/single_end -bam sorted.bam -unmapped unmapped.bam -output output -sumF summary
../trapes/trapes.py -path test_output/test_data_new/single_end -bam sorted.bam -unmapped unmapped.bam -output output -sumF summary

# diff paired-end outputs
echo "Paired-end Diffs"
echo "Differences between old and new summary.TCR.txt"
diff test_output/test_data_old/paired_end/summary.TCRs.txt test_output/test_data_new/single_end/summary.TCRs.txt
echo "Differences between old and new summary.summary.txt"
diff test_output/test_data_old/paired_end/summary.summary.txt test_output/test_data_new/single_end/summary.summary.txt
