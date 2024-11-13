#!/bin/bash


PDB=pdbs_bfvd_logan
FSDB=<FOLDSEEK_DB>/bfvd_logan
FCDB=<FOLDCOMP_DB>/bfvd_logan
foldseek createdb $PDB $FSDB/bfvd
foldcomp compress $PDB $FCDB/bfvd -t 128 --skip-discontinuous --db
