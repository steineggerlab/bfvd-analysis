#!/bin/bash
# Site configuration for the pipeline scripts (prediction/, benchmark/, analysis/clusterV2V3.sh).
# Copy to common/env.sh and edit:
#   cp common/env.example.sh common/env.sh
# Scripts source ../common/env.sh (or ./common/env.sh when run from the repository root).

# --- working directories ---
export WKDIR=${WKDIR:-$(pwd)}                         # large working area (MSAs, predictions, DBs)
export METADIR=${METADIR:-${WKDIR}/metadata}          # tables written by benchmark_analysis.sh
export MSADIR=${MSADIR:-${WKDIR}/msa_base}            # base MSAs (colabfold_search)
export TEMPLATEDIR=${TEMPLATEDIR:-${WKDIR}/template_base}
export TEMPLATE_SUFFIX=${TEMPLATE_SUFFIX:-_bfvd_2023_02_v2.m8}
export LOGAN_DIR=${LOGAN_DIR:-${WKDIR}/logan_search}  # Logan50 search results / MSAs
export BENCHMARK_DIR=${BENCHMARK_DIR:-${WKDIR}/benchmark}
export DBDIR=${DBDIR:-${WKDIR}/db}                    # Foldseek/MMseqs2 databases built from the structures
export SCRATCH=${SCRATCH:-${WKDIR}/scratch}           # fast local disk for temporary files
export BFVD_BATCH_DIR=${BFVD_BATCH_DIR:-${WKDIR}/data}
export BFVD_SCRATCH_DIR=${BFVD_SCRATCH_DIR:-${SCRATCH}}
export BFVD_BATCH_URL=${BFVD_BATCH_URL:-}             # optional: URL to fetch prediction input batches from

# --- databases ---
export COLABFOLD_DB=${COLABFOLD_DB:-/path/to/colabfold_db}        # UniRef30 2302 + ColabFold envdb, see setup_databases.sh in ColabFold
export LOGANDB=${LOGANDB:-/path/to/logan50_db}                     # MMseqs2 database of Logan50
export BFVD_V2=${BFVD_V2:-/path/to/bfvd_2023_02_v2/bfvd}           # previous release (BFVD v2) Foldseek DB, used as template source
export BFVD=${BFVD:-${DBDIR}/bfvd_v3}                              # this release (BFVD v3) Foldseek DB
export COLABFOLD_DATA_DIR=${COLABFOLD_DATA_DIR:-/path/to/colabfold/params}  # AlphaFold2 weights for colabfold_batch

# --- binaries / environments ---
export MMSEQS=${MMSEQS:-mmseqs}                       # MMseqs2 18.8cc5c
export MMSEQS_COLABFOLD=${MMSEQS_COLABFOLD:-mmseqs}   # MMseqs2 build used by colabfold_search
export COLABFOLD_ENV_NAME="colabfold-bfvd"            # environment with colabfold_batch (ColabFold 1.6.0)
export MAMBA_ROOT_PREFIX="$HOME/.local/micromamba"
export MICROMAMBA_BIN="$HOME/.local/bin/micromamba"
export JAX_CACHE_DIR="$SCRATCH/.jax_cache"
