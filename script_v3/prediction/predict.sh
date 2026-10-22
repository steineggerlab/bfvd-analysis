#!/usr/bin/env bash
set -euo pipefail

BATCH_NAME=""
GPU_ID=""
DRY_RUN=false
COMPRESS="${COMPRESS:-false}"

print_usage() {
  cat <<'EOF'
Usage:
  ./prediction/predict.sh [OPTIONS]

Options:
  --batch     NAME    Batch name to process (e.g. b0000, required)
  --gpu       N       GPU device index to use (required)
  --compress          Compress prediction output after completion
  --dry-run           Print commands without executing them
EOF
}

# ── Parse flags ───────────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
  case "$1" in
    --batch)     BATCH_NAME="$2"; shift 2 ;;
    --gpu)       GPU_ID="$2";     shift 2 ;;
    --compress)  COMPRESS=true;   shift   ;;
    --dry-run)   DRY_RUN=true;    shift   ;;
    *)
      echo "Unknown argument: $1" >&2
      print_usage
      exit 1
      ;;
  esac
done

# ── Validate ──────────────────────────────────────────────────────────────────
if [[ -z "$BATCH_NAME" || -z "$GPU_ID" ]]; then
  echo "Error: --batch and --gpu are required." >&2
  print_usage
  exit 1
fi
# ─────────────────────────────────────────────────────────────────────────────

have_cmd() {
  command -v "$1" >/dev/null 2>&1
}

download_file() {
  local url="$1"
  local output="$2"
  if have_cmd curl; then
    curl -fL "$url" -o "$output"
  elif have_cmd wget; then
    wget "$url" -O "$output"
  else
    echo "Error: curl or wget is required to download files." >&2
    exit 1
  fi
}

ensure_batch() {
  local batch_dir="${BFVD_BATCH_DIR}/input/${BATCH_NAME}"
  local download_retries=3
  local download_retry_delay=5

  if [[ -d "$batch_dir" ]]; then
    echo "[$BATCH_NAME] Input already present. Skipping download."
    return
  fi

  if [[ "$DRY_RUN" == "true" ]]; then
    echo "[dry-run] Would download ${BFVD_BATCH_URL}/${BATCH_NAME}.tar.gz → $batch_dir"
    return
  fi

  if [[ -z "$BFVD_BATCH_URL" ]]; then
    echo "Error: BFVD_BATCH_URL is not set and '$batch_dir' does not exist." >&2
    exit 1
  fi

  local url="${BFVD_BATCH_URL}/${BATCH_NAME}.tar.gz"
  local archive
  archive="$(mktemp --suffix=.tar.gz)"
  trap 'rm -f "$archive"' RETURN

  local attempt=1
  while [[ $attempt -le $download_retries ]]; do
    echo "[$BATCH_NAME] Download attempt $attempt / $download_retries from '$url'..."
    if download_file "$url" "$archive"; then
      echo "[$BATCH_NAME] Download succeeded. Decompressing..."
      mkdir -p "${BFVD_BATCH_DIR}/input"
      tar -xzf "$archive" -C "${BFVD_BATCH_DIR}/input"
      echo "[$BATCH_NAME] Decompression complete."
      return
    fi
    rm -f "$archive"
    if [[ $attempt -lt $download_retries ]]; then
      echo "[$BATCH_NAME] Retrying in ${download_retry_delay}s..."
      sleep "$download_retry_delay"
    fi
    attempt=$(( attempt + 1 ))
  done

  echo "Error: Failed to download '$BATCH_NAME' after $download_retries attempts." >&2
  exit 1
}

if have_cmd conda; then
  ENV_RUN_CMD=(conda run -n "$COLABFOLD_ENV_NAME")
else
  ENV_RUN_CMD=("$MICROMAMBA_BIN" run -n "$COLABFOLD_ENV_NAME")
  export MAMBA_ROOT_PREFIX
fi

format_duration() {
  local secs="$1"
  printf '%02d:%02d:%02d' $(( secs / 3600 )) $(( (secs % 3600) / 60 )) $(( secs % 60 ))
}

# ── Download and predict ──────────────────────────────────────────────────────
BATCH_INDIR="${BFVD_BATCH_DIR}/input/${BATCH_NAME}"
BATCH_PREDDIR="${BFVD_SCRATCH_DIR}/prediction/${BATCH_NAME}"

echo "Batch:      $BATCH_NAME"
echo "GPU:        $GPU_ID"
echo "Input dir:  $BATCH_INDIR"
echo "Output dir: $BATCH_PREDDIR"

ensure_batch

if [[ ! -d "$BATCH_INDIR" && "$DRY_RUN" == "false" ]]; then
  echo "Error: input dir '$BATCH_INDIR' not found after download." >&2
  exit 1
fi

if [[ "$DRY_RUN" == "true" ]]; then
  echo "[dry-run] Would run:"
  echo "  CUDA_VISIBLE_DEVICES=$GPU_ID ${ENV_RUN_CMD[*]} colabfold_batch $BATCH_INDIR $BATCH_PREDDIR \\"
  echo "    --num-recycle 3 \\"
  echo "    --num-models 1 --model-order 4"
  if [[ "$COMPRESS" == "true" ]]; then
    echo "[dry-run] Would compress: ${BATCH_PREDDIR}.tar.gz"
  fi
  exit 0
fi

mkdir -p "$BATCH_PREDDIR"

BATCH_START_TIME=$(date +%s)

JAX_COMPILATION_CACHE_DIR=$JAX_CACHE_DIR  CUDA_VISIBLE_DEVICES="$GPU_ID" "${ENV_RUN_CMD[@]}" colabfold_batch "$BATCH_INDIR" "$BATCH_PREDDIR" \
  --num-recycle 3 --skip-output msa,plots \
  --num-models 1 --model-order 4

# Compress output if requested
if [[ "$COMPRESS" == "true" ]]; then
  # Remove input
  rm -rf "$BATCH_INDIR"
  tar -czf "${BATCH_PREDDIR}.tar.gz" -C "${BFVD_BATCH_DIR}/prediction" "$BATCH_NAME"
fi

BATCH_END_TIME=$(date +%s)
BATCH_ELAPSED=$(( BATCH_END_TIME - BATCH_START_TIME ))
echo "[$BATCH_NAME] Complete. Time: $(format_duration $BATCH_ELAPSED)"
# ─────────────────────────────────────────────────────────────────────────────
