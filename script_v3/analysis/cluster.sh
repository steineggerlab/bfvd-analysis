#!/usr/bin/env bash
set -euo pipefail

# === INPUTS ===
# work/db/bfvd_v3 (repo root) — BFVD v3 foldseek DB
# =============
#
# Clusters the BFVD v3 protein sequences with mmseqs and writes a cluster.tsv
# (representative, member), the same format foldseek easy-cluster produces.

cd "$(dirname "$0")/.."  # repository root
source ./common/analysis_env.sh

SEQDB="$WORK_ROOT/db/bfvd_v3"
if [ ! -f "${SEQDB}.dbtype" ]; then
    echo "sequence db not found at ${SEQDB} - build it first (foldseek createdb of the BFVD v3 structures)" >&2
    exit 1
fi

WORK="$WORK_ROOT/cluster"
CLUDB="$WORK/bfvd_v3_seqclu"
OUT_TSV="$WORK/bfvd_v3_seqclu.tsv"
TMP="$WORK/tmp_seqclu"

MIN_SEQ_ID="${MIN_SEQ_ID:-0.3}"
COV="${COV:-0.9}"

if [ -f "$OUT_TSV" ]; then
    echo "cluster tsv already exists at ${OUT_TSV}, skipping (delete it to force a re-run)"
    exit 0
fi

mkdir -p "$TMP"
echo "=== $(date) : clustering $SEQDB with mmseqs (min-seq-id=$MIN_SEQ_ID, -c=$COV) ==="
"$MMSEQS" cluster "$SEQDB" "$CLUDB" "$TMP" --min-seq-id "$MIN_SEQ_ID" -c "$COV" --threads "$THREADS"

echo "--- writing cluster tsv ---"
"$MMSEQS" createtsv "$SEQDB" "$SEQDB" "$CLUDB" "$OUT_TSV"

rm -rf "$TMP"
echo "=== $(date) : done: $OUT_TSV ==="
