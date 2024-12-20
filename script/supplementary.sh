#!/bin/bash

METADIR=metadata
ETC=etc

old=$METADIR/bfvd-model_length_unk_nmsa_plddt_ptm.tsv
new=$METADIR/logan-model_nmsa_plddt_ptm.tsv

awk -F"\t" 'BEGIN {OFS="\t"} NR==FNR {logan[$1]=$2;next} ($1 in logan) {print $1,$4,logan[$1]}' $new $old > $ETC/nmsa_change_logan-model_before_after.tsv
