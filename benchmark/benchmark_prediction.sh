#!/bin/bash

set -euo pipefail
# set -x
if [ -z ${SLURM_SUBMIT_DIR+x} ]; then
    cd $(dirname $0)
    source ../common/env.sh
else 
    cd $SLURM_SUBMIT_DIR
    source ./common/env.sh
fi
source activate colabfold-dev

# Parse parameters: recycles, template, msa
while getopts r:m:t flag
do
    case "${flag}" in
        r) num_recycle=${OPTARG};;
        m) msa_mode=${OPTARG};;
        t) use_template=true;;
    esac
done

INDIR=${BENCHMARK_DIR}/input
PREDDIR=${BENCHMARK_DIR}/prediction

use_template=${use_template:-false}
num_recycle=${num_recycle:-3}
msa_mode=${msa_mode:-"base"}
num_models=${num_models:-5}

if [ "$msa_mode" = "base" ] || [ "$msa_mode" = "logan" ] ; then
    :
else
    echo "Invalid msa mode: $msa_mode"
    exit 1
fi

if [ "$use_template" = true ] ; then
    template_args="--templates --pdb-hit-suffix _bfvd_2023_02_v2.m8  --local-pdb-path ${WKDIR}/template/bfvd_cif"
    num_models=2
else
    template_args=""
fi
INDIR=${INDIR}/${msa_mode}
PREDDIR=${PREDDIR}/${msa_mode}_r${num_recycle}$( [ "$use_template" = true ] && echo "_template" || echo "" )

echo "Running benchmark with settings: msa_mode=$msa_mode, num_recycle=$num_recycle, use_template=$use_template"
echo "Input dir: $INDIR"
echo "Prediction dir: $PREDDIR"

colabfold_batch $INDIR $PREDDIR \
    --data ${COLABFOLD_DATA_DIR} \
    --num-recycle $num_recycle --num-models $num_models $template_args