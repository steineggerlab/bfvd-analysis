#!/bin/bash

#SBATCH --job-name=cf_search
#SBATCH --nodelist=super005
#SBATCH --partition=compute
#SBATCH --time=0
#SBATCH --cpus-per-task=128
#SBATCH --output="log/cf_search_proteome.out"

SEQ_FILE=proteome_cover/concat_processed.fasta
MSA_DIR=proteome_cover/msa

colabfold_search --db1 uniref30_2202_db \
--db3 colabfold_envdb_202108_db \
--use-env 1 --use-templates 0 --filter 1 -s 7 --threads 64 --db-load-mode 1 \
$SEQ_FILE /storage/databases/colabfold_db_all/ $SCRATCH

mv $SCRATCH $MSA_DIR
