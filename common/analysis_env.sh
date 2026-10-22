#!/usr/bin/env bash
# Common variables shared by scripts under scripts/*/. Source this near the top
# of a run script, e.g.:
#   REPO_ROOT="."
#   source "$REPO_ROOT/common/analysis_env.sh"

REPO_ROOT="${REPO_ROOT:-.}"

DATA_DIR="${DATA_DIR:-$REPO_ROOT/data}"
WORK_ROOT="${WORK_ROOT:-$REPO_ROOT/work}"
ANALYSES_ROOT="${ANALYSES_ROOT:-$REPO_ROOT/analyses}"

THREADS="${THREADS:-$(nproc)}"
FOLDSEEK="${FOLDSEEK:-foldseek}"
MMSEQS="${MMSEQS:-mmseqs}"

STAMP=$(cat "$REPO_ROOT/.run_stamp" 2>/dev/null || echo "manual")
