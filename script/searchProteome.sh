#!/bin/bash


PROTEOME=<BFVD_ANALYSIS>/proteome_cover/prediction
#BFVD=<FOLDSEEK_DB>/bfvd/bfvd
BFVD=<FOLDSEEK_DB>/bfvd_logan/bfvd
OUT=<BFVD_ANALYSIS>/proteome_cover

foldseek easy-search $PROTEOME/*unrelaxed_rank_001*.pdb $BFVD $OUT/proteome_bfvd_logan.tsv $SCRATCH --format-mode 0 --format-output "query,target,evalue,pident,fident,bits,qcov,tcov,qtmscore,ttmscore" --remove-tmp-files false

mv $SCRATCH $OUT/tmp_logan

