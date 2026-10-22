#!/usr/bin/env bash
# Usage: ./compute_neff.sh MSA_DIR OUT_TSV [JOBS]
set -euo pipefail

IN="$1"
OUT="$2"
JOBS="${3:-$(nproc)}"
: "${SCRATCH:?SCRATCH not set}"

TMP="$(mktemp -p "$SCRATCH" neff.XXXXXX.tsv)"
trap 'rm -f "$TMP"' EXIT

find "$IN" -type f -name '*.a3m' -print0 \
| xargs -0 -n1 -P "$JOBS" -I{} bash -c '
    f="$1"
    id=$(basename "$f" .a3m)
    neff=$(hhmake -i "$f" -o stdout -v 0 2>/dev/null | awk "/^NEFF/ {print \$2; exit}")
    printf "%s\t%s\n" "$id" "${neff:-NA}"
  ' _ {} > "$TMP"

mkdir -p "$(dirname "$OUT")"
cp "$TMP" "$OUT"

