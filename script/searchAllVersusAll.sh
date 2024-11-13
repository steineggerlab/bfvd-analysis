#!/bin/bash

#SBATCH -J BFVD-allVall
#SBATCH -p compute
#SBATCH -t 0
#SBATCH -o log/sbatch_bfvd_allVall.out
#SBATCH -w super002
#SBATCH -c 64

BFVD=/fast/databases/foldseek/bfvd_logan/bfvd
OUT=distribute/all_v_all

foldseek search $BFVD $BFVD $OUT/bfvd_all_versus_all $SCRATCH -s 9.5 -e 0.01
foldseek convertalis $BFVD $BFVD $OUT/bfvd_all_versus_all $OUT/bfvd_all_versus_all.tsv
awk -F"\t" 'BEGIN {OFS="\t"} {print $1,$2,$11}' $OUT/bfvd_all_versus_all.tsv > $OUT/bfvd_all_vs_all-queryId_targetId_evalues.tsv

mv $SCRATCH $OUT/tmp
