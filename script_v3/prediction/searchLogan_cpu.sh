#!/bin/bash


set -euo pipefail

if [ -z ${SLURM_SUBMIT_DIR+x} ]; then
    cd $(dirname $0)
    source ../common/env.sh
else
    cd $SLURM_SUBMIT_DIR
    source ./common/env.sh
fi

INDIR=${WKDIR}/fasta/tmp #${WKDIR}/fasta/batch500k
OUTDIR=${WKDIR}/logan_search_3rd #logan_msa

date +"%Y-%m-%d %T Starting job on %N"
if [ ! -d $LOGANDB ]; then
    mkdir -p $LOGANDB
    # wget "#TODO" # Download Logan DB files
fi

if [ ! -d $OUTDIR/aln ]; then
    mkdir -p $OUTDIR
    mkdir -p $OUTDIR/aln
    mkdir -p $OUTDIR/db
fi

if [ ! -d $SCRATCH/result ]; then
    mkdir -p $SCRATCH
    mkdir -p $SCRATCH/db
    mkdir -p $SCRATCH/tmp
    mkdir -p $SCRATCH/result
fi

# ulimit -l unlimited
#$MMSEQS touchdb $LOGANDB/logan_db --touch-lock &
#TOUCH_PID=$!

for fasta in $INDIR/*.fa*; do
    file=$(basename $fasta)
    batch_id=${file%%.fa*}
    echo "Processing batch: $batch_id"

   $MMSEQS createdb $fasta $SCRATCH/db/${batch_id}_qdb
   MMSEQS_IGNORE_INDEX=1 $MMSEQS search ${SCRATCH}/db/${batch_id}_qdb ${LOGANDB}/logan_db ${SCRATCH}/${batch_id}_alndb ${SCRATCH}/tmp -e 0.1 --max-seqs 1000 --db-load-mode 2 --num-iterations 2 -a
   $MMSEQS mvdb ${SCRATCH}/${batch_id}_alndb ${OUTDIR}/aln/${batch_id}_alndb
   $MMSEQS mvdb ${SCRATCH}/db/${batch_id}_qdb  ${OUTDIR}/db/${batch_id}_qdb
   $MMSEQS mvdb ${SCRATCH}/db/${batch_id}_qdb_h  ${OUTDIR}/db/${batch_id}_qdb_h

    # $MMSEQS result2msa ${SCRATCH}/db/${batch_id}_qdb ${LOGANDB}/logan_db ${SCRATCH}/${batch_id}_alndb ${SCRATCH}/result/${batch_id}_msadb --msa-format-mode 6 --db-load-mode 2
#    cp ${SCRATCH}/db/${batch_id}_qdb.lookup ${SCRATCH}/result/${batch_id}_msadb.lookup
#    $MMSEQS mvdb ${SCRATCH}/result/${batch_id}_msadb ${OUTDIR}/${batch_id}_msadb
    # $MMSEQS unpackdb ${SCRATCH}/result/${batch_id}_msadb ${OUTDIR} --unpack-suffix _logan.a3m #TODO: adjust
done
#kill "$TOUCH_PID"
#wait
date +"%Y-%m-%d %T Job finished on %N"
