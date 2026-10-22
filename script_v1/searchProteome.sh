#!/bin/bash


PROTEOME=proteome_cover/prediction
#BFVD=<FOLDSEEK_DB>/bfvd_logan/bfvd
BFVD=<FOLDSEEK_DB>/bfvd/bfvd
OUT=proteome_cover

foldseek easy-search $PROTEOME/*unrelaxed_rank_001*.pdb $BFVD $OUT/proteome_bfvd_logan.tsv $SCRATCH --format-mode 0 --format-output "query,target,evalue,pident,fident,bits,qcov,tcov,qtmscore,ttmscore" --remove-tmp-files false

mv $SCRATCH $OUT/tmp_logan

