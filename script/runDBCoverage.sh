#!/bin/bash


QDB=<FOLDSEEK_DB>/bfvd_logan/bfvd
AFDB=<FOLDSEEK_DB>/afdb/afdb50
PDB=<FOLDSEEK_DB>/pdb/pdb100
OUT=<BFVD_ANALYSIS>/db_cover

foldseek easy-search $QDB $PDB $OUT/bfvd_logan_pdb100_aln $SCRATCH/pdb100_tmp --remove-tmp-files false --format-output query,target,evalue,qtmscore,ttmscore,lddtfull --greedy-best-hits
foldseek easy-search $QDB $AFDB $OUT/bfvd_logan_afdb50_aln $SCRATCH/afdb50_tmp --remove-tmp-files false --format-output query,target,evalue,qtmscore,ttmscore,lddtfull --greedy-best-hits

mv $SCRATCH $OUT/tmp_bfvd_logan

