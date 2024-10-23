#!/bin/bash

#SBATCH -J BFVD_DB
#SBATCH -t 0
#SBATCH -p compute
#SBATCH -w super005
#SBATCH -c 64
#SBATCH -o log/sbatch_bfvd_logan_createdb.out

#PDB=/home/seamustard52/bfvd-analysis/pdbs
#PATHFILE=/home/seamustard52/bfvd-analysis/metadata/entry_filepath_rep.tsv
#FSDB=/fast/databases/foldseek/bfvd
#FCDB=/fast/databases/foldcomp/bfvd

#foldseek createdb $PDB $FSDB/bfvd
#foldcomp compress $PDB $FCDB/bfvd -t 128 --skip-discontinuous --db

#awk -F"\t" 'NR==FNR {if ($3==0) {wounk[$1]=$0;next}} $1 in wounk  {print $2}' metadata/bfvd-model_length_unk_nmsa_plddt_ptm.tsv metadata/entry_filepath_rep.tsv > wounk.tsv
#foldseek createdb wounk.tsv  $FSDB/bfvd_wounk
#foldcomp compress -t 128 --file wounk.tsv --skip-discontinuous --db 1
#rm wounk.tsv

PDB=/home/seamustard52/bfvd-analysis/pdbs_bfvd_logan
FSDB=/fast/databases/foldseek/bfvd_logan
FCDB=/fast/databases/foldcomp/bfvd_logan
foldseek createdb $PDB $FSDB/bfvd
foldcomp compress $PDB $FCDB/bfvd -t 128 --skip-discontinuous --db
