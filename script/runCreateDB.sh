#!/bin/bash

#SBATCH -J BFVD_DB
#SBATCH -t 0
#SBATCH -p compute
#SBATCH -w super005
#SBATCH -c 64
#SBATCH -o log/sbatch_bfvd_logan_createdb.out

PDB=pdbs_bfvd_logan
PDBOLD=pdbs
#FSDB=/fast/databases/foldseek/bfvd_logan
FSDB=/fast/databases/foldseek/bfvd
FCDB=/fast/databases/foldcomp/bfvd_logan
foldseek createdb $PDB $FSDB/bfvd
foldseek createdb $PDBOLD $FSDB/bfvd_base
foldcomp compress $PDB $FCDB/bfvd -t 128 --skip-discontinuous --db
