#!/bin/sh

cd $1
rm *TCRs.txt *summary.txt
cd cell1 && rm -r -f rsem* *output* *fq *sam  *temp* && cd ..
cd cell2 && rm -r -f rsem* *output* *fq *sam  *temp*