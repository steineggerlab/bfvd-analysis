#!/bin/bash

#SBATCH -J bfvd-cluster
#SBATCH -o log/sbatch_bfvdlogan_cluster.out
#SBATCH -p compute
#SBATCH -w super005
#SBATCH -c 96
#SBATCH -t 0

# conda env foldseek-latest

DB=pdbs_bfvd_logan
CLUDIR=bfvd_cluster
C=(0.7 0.9 0.8 0.6 0.5 0.4 0.3 0.2 0.1)
for c in "${C[@]}"
do
    coverage=$(bc <<< "scale=1; $c * 100" | cut -d . -f 1)
    foldseek easy-cluster $DB $CLUDIR/bfvd_logan_clu_$coverage $SCRATCH/c$coverage -c $c --remove-tmp-files false
done

mv $SCRATCH $CLUDIR/tmp/tmp_bfvd_logan
