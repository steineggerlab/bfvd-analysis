#!/bin/bash


#BFVD=<FOLDSEEK_DB>/bfvd_logan/bfvd
BFVD=<FOLDSEEK_DB>/bfvd/bfvd
OUT=distribute/all_v_all

foldseek search $BFVD $BFVD $OUT/bfvd_all_versus_all $SCRATCH -s 9.5 -e 0.01
foldseek convertalis $BFVD $BFVD $OUT/bfvd_all_versus_all $OUT/bfvd_all_versus_all.tsv
awk -F"\t" 'BEGIN {OFS="\t"} {print $1,$2,$11}' $OUT/bfvd_all_versus_all.tsv > $OUT/bfvd_all_vs_all-queryId_targetId_evalues.tsv

mv $SCRATCH $OUT/tmp
