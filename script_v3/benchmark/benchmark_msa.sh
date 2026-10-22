#!/bin/bash
set -euo pipefail
if [ -z ${SLURM_SUBMIT_DIR+x} ]; then
    cd $(dirname $0)
    source ../common/env.sh
else 
    cd $SLURM_SUBMIT_DIR
    source ./common/env.sh
fi
source activate localcolabfold

OUTDIR=${BENCHMARK_DIR}/input

NCPUS=$(nproc)

if [ ! -d $OUTDIR ]; then
    mkdir -p $OUTDIR
else
   echo "Output directory $OUTDIR already exists. Results can be overwritten."
fi

if [ ! -d $SCRATCH ]; then
    mkdir -p $SCRATCH
fi

date +"%Y-%m-%d %T Starting job on %N"

fasta=${BENCHMARK_DIR}/candidate_filtered_rep_seq.fasta
colabfold_search \
    --db1 uniref30_2302_db --db2 bfvd_2023_02_v2 --db3 colabfold_envdb_202108_db \
    --use-env 1 --use-templates 1 --filter 1 --threads $NCPUS --db-load-mode 2 -s 7 \
    --mmseqs "$MMSEQS_COLABFOLD" \
    "$fasta" "$COLABFOLD_DB" "$SCRATCH"

find $SCRATCH -name *.a3m | xargs -I {} -P $NCPUS mv {} $OUTDIR
find $SCRATCH -name *.m8 | xargs -I {} -P $NCPUS mv {} $OUTDIR

date +"%Y-%m-%d %T Job completed on %N"
