#!/usr/bin/bash

set -euo pipefail
set +u

if [ -z ${SLURM_SUBMIT_DIR+x} ]; then
  cd "$(dirname "$0")"
  source ../common/env.sh
else
  cd "$SLURM_SUBMIT_DIR"
  source ./common/env.sh
fi

NUM_GPUS=8
BATCH_START=""
BATCH_END=""
DRY_RUN=false
COMPRESS=false

print_usage() {
  cat <<'EOF'
Usage:
  ./prediction/run_prediction.sh [OPTIONS]

Options:
  --start    N     First batch index to process (required)
  --end      N     Last batch index to process, inclusive (required)
  --compress       Compress prediction output after each batch
  --dry-run        Print commands without executing them
EOF
}

# ── Parse flags ───────────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
  case "$1" in
    --start)    BATCH_START="$2"; shift 2 ;;
    --end)      BATCH_END="$2";   shift 2 ;;
    --compress) COMPRESS=true;    shift   ;;
    --dry-run)  DRY_RUN=true;     shift   ;;
    *)
      echo "Unknown argument: $1" >&2
      print_usage
      exit 1
      ;;
  esac
done

# ── Validate ──────────────────────────────────────────────────────────────────
if [[ -z "$BATCH_START" || -z "$BATCH_END" ]]; then
  echo "Error: --start and --end are required." >&2
  print_usage
  exit 1
fi
if ! [[ "$BATCH_START" =~ ^[0-9]+$ && "$BATCH_END" =~ ^[0-9]+$ ]]; then
  echo "Error: --start and --end must be non-negative integers." >&2
  exit 1
fi
if [[ "$BATCH_START" -gt "$BATCH_END" ]]; then
  echo "Error: --start ($BATCH_START) must be <= --end ($BATCH_END)." >&2
  exit 1
fi
# ─────────────────────────────────────────────────────────────────────────────

TOTAL=$(( BATCH_END - BATCH_START + 1 ))
DRY_RUN_FLAG=""
[[ "$DRY_RUN" == "true" ]] && DRY_RUN_FLAG="--dry-run"
COMPRESS_FLAG=""
[[ "$COMPRESS" == "true" ]] && COMPRESS_FLAG="--compress"

mkdir -p "${BFVD_BATCH_DIR}/logs"

echo "Batch range:  b$(printf '%04d' "$BATCH_START") → b$(printf '%04d' "$BATCH_END")"
echo "Total:        $TOTAL"
echo "GPUs:         $NUM_GPUS"
[[ "$DRY_RUN" == "true" ]] && echo "[dry-run] No commands will be executed."
echo

# ── GPU slot tracker ──────────────────────────────────────────────────────────
declare -a GPU_PIDS=()
for (( g=0; g<NUM_GPUS; g++ )); do
  GPU_PIDS[$g]=""
done

# ── Dispatch batches across GPUs ──────────────────────────────────────────────
NEXT_IDX=$BATCH_START

# Initial fill: assign one batch to each GPU
for (( g=0; g<NUM_GPUS && NEXT_IDX<=BATCH_END; g++ )); do
  BATCH_NAME="$(printf 'b%04d' "$NEXT_IDX")"
  echo "[GPU $g] Assigning batch '$BATCH_NAME' ($((NEXT_IDX - BATCH_START + 1)) / $TOTAL)..."
  ./prediction/predict.sh --batch "$BATCH_NAME" --gpu "$g" $COMPRESS_FLAG $DRY_RUN_FLAG \
    >> "${BFVD_BATCH_DIR}/logs/${BATCH_NAME}_gpu${g}.log" 2>&1 &
  GPU_PIDS[$g]=$!
  NEXT_IDX=$(( NEXT_IDX + 1 ))
done

# As each GPU finishes, assign the next unprocessed batch
while [[ $NEXT_IDX -le $BATCH_END ]]; do
  ASSIGNED=false
  for (( g=0; g<NUM_GPUS; g++ )); do
    if [[ -n "${GPU_PIDS[$g]}" ]] && ! kill -0 "${GPU_PIDS[$g]}" 2>/dev/null; then
      set +e; wait "${GPU_PIDS[$g]}"; EXIT_CODE=$?; set -e
      [[ $EXIT_CODE -ne 0 ]] && echo "Warning: GPU $g exited with code $EXIT_CODE. Check logs." >&2

      BATCH_NAME="$(printf 'b%04d' "$NEXT_IDX")"
      echo "[GPU $g] Assigning batch '$BATCH_NAME' ($((NEXT_IDX - BATCH_START + 1)) / $TOTAL)..."
      ./prediction/predict.sh --batch "$BATCH_NAME" --gpu "$g" $COMPRESS_FLAG $DRY_RUN_FLAG \
        >> "${BFVD_BATCH_DIR}/logs/${BATCH_NAME}_gpu${g}.log" 2>&1 &
      GPU_PIDS[$g]=$!
      NEXT_IDX=$(( NEXT_IDX + 1 ))
      ASSIGNED=true

      [[ $NEXT_IDX -gt $BATCH_END ]] && break
    fi
  done
  # Only sleep if no GPU was free this iteration
  [[ "$ASSIGNED" == "false" ]] && sleep 2
done

# ── Wait for all remaining jobs ───────────────────────────────────────────────
echo "All batches dispatched. Waiting for remaining jobs..."
for (( g=0; g<NUM_GPUS; g++ )); do
  if [[ -n "${GPU_PIDS[$g]}" ]]; then
    set +e; wait "${GPU_PIDS[$g]}"; EXIT_CODE=$?; set -e
    [[ $EXIT_CODE -ne 0 ]] && echo "Warning: GPU $g exited with code $EXIT_CODE. Check logs." >&2
    echo "[GPU $g] Done."
  fi
done

echo
echo "All batches complete."
# ─────────────────────────────────────────────────────────────────────────────
