#!/bin/bash

#SBATCH -J DB_coverage
#SBATCH -p compute
#SBATCH -t 0
#SBATCH -o log/sbatch_bfvd_logan_db_coverage.out
#SBATCH -w super005
#SBATCH -c 128

#QDB=/fast/databases/foldseek/bfvd_logan/bfvd
QDB=/fast/databases/foldseek/bfvd/bfvd
AFDB=/fast/databases/foldseek/afdb/afdb50
PDB=/fast/databases/foldseek/pdb/pdb100
OUT=db_cover

foldseek easy-search $QDB $PDB $OUT/bfvd_logan_pdb100_aln $SCRATCH/pdb100_tmp --remove-tmp-files false --format-output query,target,evalue,qtmscore,ttmscore,lddtfull --greedy-best-hits
foldseek easy-search $QDB $AFDB $OUT/bfvd_logan_afdb50_aln $SCRATCH/afdb50_tmp --remove-tmp-files false --format-output query,target,evalue,qtmscore,ttmscore,lddtfull --greedy-best-hits

mv $SCRATCH $OUT/tmp_bfvd_logan

