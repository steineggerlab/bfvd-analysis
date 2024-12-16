#!/bin/bash


PDB=pdbs_bfvd_logan
PDBOLD=pdbs
#FSDB=<FOLDSEEK_DB>/bfvd_logan
FSDB=<FOLDSEEK_DB>/bfvd
FCDB=<FOLDCOMP_DB>/bfvd_logan
foldseek createdb $PDB $FSDB/bfvd
foldseek createdb $PDBOLD $FSDB/bfvd_base
foldcomp compress $PDB $FCDB/bfvd -t 128 --skip-discontinuous --db
