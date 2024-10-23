#!/bin/bash

#SBATCH -J phage_annotation
#SBATCH -p compute
#SBATCH -t 0
#SBATCH -o log/sbatch_phage_annotation_logan.out
#SBATCH -w hulk
#SBATCH -c 64

PHAGE=/home/seamustard52/virusDB/application/Small_ACCs
#BFVD=/fast/databases/foldseek/bfvd/bfvd
#OUT=/home/seamustard52/bfvd-analysis/phage_annotation/base
BFVD=/fast/databases/foldseek/bfvd_logan/bfvd
OUT=/home/seamustard52/bfvd-analysis/phage_annotation/logan

for sample in $(find $PHAGE -maxdepth 1 -type d -name gac* );
do 
    for f in $(find $sample/INHERIT_phages/colabfold_predictions -maxdepth 1 -mindepth 1 -type d)
    do 
        foldseek easy-search $f/*unrelaxed_rank_001*.pdb $BFVD $OUT/${f##*/}.tab $SCRATCH --format-mode 0 --format-output "query,target,evalue,pident,fident,bits,qcov,tcov,alntmscore" -e 0.001; sort -u -k1,1 $OUT/${f##*/}.tab > $OUT/${f##*/}_besthits.tab
    done
done

mv $SCRATCH $OUT/tmp

### Reference
#nohup bash -c 'for f in colabfold_predictions/*; do \
#  foldseek easy-search $f/*unrelaxed_rank_1*.pdb /Volumes/data2/database/afdb ./foldseek_hits/${f##*/}.tab ./foldseek_hits/tmp --format-mode 0 --format-output "query,target,evalue,pident,fident,bits,qcov,tcov,alntmscore" -e 0.001; sort -u -k1,1 ./foldseek_hits/${f##*/}.tab > ./foldseek_hits/${f##*/}_besthits.tab ; \
#  done' &> ./foldseek_hits/foldseek.log &
