#!/bin/bash


#PDB=<BFVD_ANALYSIS>/pdbs
#PATHFILE=<BFVD_ANALYSIS>/metadata/entry_filepath_rep.tsv
#FSDB=<FOLDSEEK_DB>/bfvd
#FCDB=<FOLDCOMP_DB>/bfvd

#foldseek createdb $PDB $FSDB/bfvd
#foldcomp compress $PDB $FCDB/bfvd -t 128 --skip-discontinuous --db

#awk -F"\t" 'NR==FNR {if ($3==0) {wounk[$1]=$0;next}} $1 in wounk  {print $2}' metadata/bfvd-model_length_unk_nmsa_plddt_ptm.tsv metadata/entry_filepath_rep.tsv > wounk.tsv
#foldseek createdb wounk.tsv  $FSDB/bfvd_wounk
#foldcomp compress -t 128 --file wounk.tsv --skip-discontinuous --db 1
#rm wounk.tsv

PDB=<BFVD_ANALYSIS>/pdbs_bfvd_logan
FSDB=<FOLDSEEK_DB>/bfvd_logan
FCDB=<FOLDCOMP_DB>/bfvd_logan
foldseek createdb $PDB $FSDB/bfvd
foldcomp compress $PDB $FCDB/bfvd -t 128 --skip-discontinuous --db
