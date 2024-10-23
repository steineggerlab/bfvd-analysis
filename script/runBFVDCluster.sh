#!/bin/bash


# conda env foldseek-latest
## BFVD without Logan
#DB=<BFVD_ANALYSIS>/pdbs
#CLUDIR=<BFVD_ANALYSIS>/bfvd_cluster
#C=(0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1)

## BFVD with Logan
DB=<BFVD_ANALYSIS>/pdbs_bfvd_logan
CLUDIR=<BFVD_ANALYSIS>/bfvd_cluster
C=(0.7 0.9 0.8 0.6 0.5 0.4 0.3 0.2 0.1)
for c in "${C[@]}"
do
    coverage=$(bc <<< "scale=1; $c * 100" | cut -d . -f 1)
    foldseek easy-cluster $DB $CLUDIR/bfvd_logan_clu_$coverage $SCRATCH/c$coverage -c $c --remove-tmp-files false
done

mv $SCRATCH $CLUDIR/tmp/tmp_bfvd_logan
