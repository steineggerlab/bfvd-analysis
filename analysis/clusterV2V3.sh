#!/usr/bin/bash
set -euo pipefail

cd "$(dirname "$0")"
source ../common/env.sh

CONCATDB=${DBDIR}/bfvd_concat_v2_v3/bfvd_concat
if [ ! -f "${CONCATDB}.dbtype" ]; then
    echo "db not found at ${CONCATDB} - create concatentated db first" >&2
    exit 1
fi

WORK=${WKDIR}/bfvd_concat_cluster

MIN_SEQ_ID="${MIN_SEQ_ID:-0.3}"
SEQCOV="${SEQCOV:-0.9}"
STRCOV="${STRCOV:-0.9}"

OUT_TSV="$WORK/bfvd_concat_seqclu.tsv"
if [ ! -f "$OUT_TSV" ]; then
echo "=== $(date) : clustering $CONCATDB with mmseqs (min-seq-id=$MIN_SEQ_ID, -c=$SEQCOV) ==="
"$MMSEQS" cluster "$CONCATDB" "$SCRATCH/seqclu" "$SCRATCH/seqtmp" --min-seq-id "$MIN_SEQ_ID" -c "$SEQCOV"

echo "--- writing cluster tsv ---"
"$MMSEQS" createtsv "$CONCATDB" "$CONCATDB" "$SCRATCH/seqclu" "$OUT_TSV"

echo "=== $(date) : sequence clustering done: $OUT_TSV ==="
else
    echo "cluster tsv already exists at ${OUT_TSV}, skipping (delete it to force a re-run)"
fi

OUT_TSV="$WORK/bfvd_concat_strclu.tsv"
if [ ! -f "$OUT_TSV" ]; then
echo "=== $(date) : clustering $CONCATDB with foldseek (-c=$STRCOV) ==="
foldseek cluster "$CONCATDB" "$SCRATCH/strclu" "$SCRATCH/strtmp" -c "$STRCOV"

echo "--- writing cluster tsv ---"
foldseek createtsv "$CONCATDB" "$CONCATDB" "$SCRATCH/strclu" "$OUT_TSV"

echo "=== $(date) : structure clustering done: $OUT_TSV ==="
else
    echo "cluster tsv already exists at ${OUT_TSV}, skipping (delete it to force a re-run)"
fi
