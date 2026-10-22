#!/bin/bash


set -euo pipefail

cd $SLURM_SUBMIT_DIR
source ./common/env.sh
source activate localcolabfold

INDIR=${WKDIR}/fasta/batch500k
MSADIR=${WKDIR}/msa_base
TEMPLATEDIR=${WKDIR}/template_base

NCPUS=$(nproc)

if [ ! -d $INDIR ]; then
    echo "Input directory $INDIR does not exist. Exiting."
    exit 1
fi

if [ ! -d $MSADIR ]; then
    mkdir -p $MSADIR
#else
#    echo "Output directory $MSADIR already exists. Exiting to prevent overwriting."
#    exit 1
fi

if [ ! -d $TEMPLATEDIR ]; then
    mkdir -p $TEMPLATEDIR
#else
#    echo "Output directory $TEMPLATEDIR already exists. Exiting to prevent overwriting."
#    exit 1
fi

if [ ! -d $SCRATCH ]; then
    mkdir -p $SCRATCH
fi

date +"%Y-%m-%d %T Starting job on %N"
if false;then
for fasta in $INDIR/*.fa*; do
    file=$(basename $fasta)
    batch_id=${file%%.fa*}
    echo "Processing batch: $batch_id"

    colabfold_search \
        --db1 uniref30_2302_db --db2 bfvd_2023_02_v2 --db3 colabfold_envdb_202108_db \
        --use-env 1 --use-templates 1 --filter 1 --threads $NCPUS --db-load-mode 2 -s 7 \
        --mmseqs "$MMSEQS_COLABFOLD" \
        "$fasta" "$COLABFOLD_DB" "$SCRATCH"
    
    find $SCRATCH -name *.a3m | xargs -I {} -P $NCPUS mv {} $MSADIR
    find $SCRATCH -name *.m8 | xargs -I {} -P $NCPUS mv {} $TEMPLATEDIR
done
fi
MSADIR=${WKDIR}/msa_base_3rd
TEMPLATEDIR=${WKDIR}/template_base_3rd
file=${WKDIR}/fasta/uniprot_2025_03_virus_1line-ab3000_nosplit_3rd.fasta
batch_id=${file%%.fa*}
echo "Processing batch: $batch_id"

colabfold_search \
    --db1 uniref30_2302_db --db2 bfvd_2023_02_v2 --db3 colabfold_envdb_202108_db \
    --use-env 1 --use-templates 1 --filter 1 --threads $NCPUS --db-load-mode 2 -s 7 \
    --mmseqs "$MMSEQS_COLABFOLD" \
    "$file" "$COLABFOLD_DB" "$SCRATCH"

find $SCRATCH -name *.a3m | xargs -I {} -P $NCPUS mv {} $MSADIR
find $SCRATCH -name *.m8 | xargs -I {} -P $NCPUS mv {} $TEMPLATEDIR

date +"%Y-%m-%d %T Job completed on %N"
