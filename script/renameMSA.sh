#!/bin/bash

SEQFILE=/home/seamustard52/bfvd-analysis/proteome_cover/concat_processed.fasta
MSADIR=/home/seamustard52/bfvd-analysis/proteome_cover/msa/473655

awk -F">" -v dir=$MSADIR 'BEGINFILE {idx=0} $0~/^>/ {cmd="mv "dir"/"idx".a3m "dir"/"$2".a3m";system(cmd) ;idx++}' $SEQFILE
