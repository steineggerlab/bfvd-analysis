#!/bin/bash

#SBATCH -J proteome
#SBATCH -p compute
#SBATCH -t 0
#SBATCH -o log/sbatch_proteomesearch_logan.out
#SBATCH -w hulk
#SBATCH -c 64

PROTEOME=/home/seamustard52/bfvd-analysis/proteome_cover/prediction
#BFVD=/fast/databases/foldseek/bfvd/bfvd
BFVD=/fast/databases/foldseek/bfvd_logan/bfvd
OUT=/home/seamustard52/bfvd-analysis/proteome_cover

foldseek easy-search $PROTEOME/*unrelaxed_rank_001*.pdb $BFVD $OUT/proteome_bfvd_logan.tsv $SCRATCH --format-mode 0 --format-output "query,target,evalue,pident,fident,bits,qcov,tcov,qtmscore,ttmscore" --remove-tmp-files false

mv $SCRATCH $OUT/tmp_logan

