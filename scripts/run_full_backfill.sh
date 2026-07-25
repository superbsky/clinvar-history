#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${PROJECT_ROOT}"

START_DATE="${CLINVAR_START_DATE:-2020-01-01}"
END_DATE="${CLINVAR_END_DATE:-2026-07-24}"
DATA_ROOT="${CLINVAR_DATA_ROOT:-data/full}"
RESULTS_ROOT="${CLINVAR_RESULTS_ROOT:-results/full-v1}"
MANIFEST="${DATA_ROOT}/manifest/releases.parquet"
UV_BIN="${CLINVAR_UV_BIN:-uv}"
RETRY_DELAY_SECONDS="${CLINVAR_RETRY_DELAY_SECONDS:-60}"

mkdir -p "$(dirname "${MANIFEST}")" "${RESULTS_ROOT}"

stage() {
  printf '[%s] %s\n' "$(date -u '+%Y-%m-%dT%H:%M:%SZ')" "$1"
}

run_with_retry() {
  local label="$1"
  shift
  local attempt=1
  until "$@"; do
    local status=$?
    stage "${label} failed with status ${status}; retrying in ${RETRY_DELAY_SECONDS}s (attempt ${attempt})"
    attempt=$((attempt + 1))
    sleep "${RETRY_DELAY_SECONDS}"
  done
}

parquet_ready() {
  local output="$1"
  [[ -s "${output}" ]] || return 1
  "${UV_BIN}" run python -c \
    'import sys; import pyarrow.parquet as pq; pq.ParquetFile(sys.argv[1])' \
    "${output}" >/dev/null 2>&1
}

run_parquet_stage() {
  local label="$1"
  local output="$2"
  shift 2
  if parquet_ready "${output}"; then
    stage "Reusing existing ${label} ${output}"
  else
    "$@"
  fi
}

if [[ -f "${MANIFEST}" ]]; then
  stage "Reusing existing manifest ${MANIFEST}"
else
  stage "Building release manifest for ${START_DATE} through ${END_DATE}"
  run_with_retry \
    "Manifest discovery" \
    "${UV_BIN}" run clinvar-history build-manifest \
    --start "${START_DATE}" \
    --end "${END_DATE}" \
    --output "${MANIFEST}"
fi

stage "Ingesting weekly GRCh38 VCF releases"
run_with_retry \
  "VCF ingestion" \
  "${UV_BIN}" run clinvar-history ingest-vcf \
  --manifest "${MANIFEST}" \
  --root "${DATA_ROOT}"

stage "Ingesting monthly submission summaries"
run_with_retry \
  "Submission ingestion" \
  "${UV_BIN}" run clinvar-history ingest-submissions \
  --manifest "${MANIFEST}" \
  --root "${DATA_ROOT}"

stage "Attributing aggregate events"
run_parquet_stage \
  "attributions" \
  "${RESULTS_ROOT}/attributions.parquet" \
  "${UV_BIN}" run clinvar-history attribute-events \
    --root "${DATA_ROOT}" \
    --output "${RESULTS_ROOT}/attributions.parquet" \
    --max-lag-days 35

stage "Running offline analyses"
run_parquet_stage \
  "Renovo analysis" \
  "${RESULTS_ROOT}/renovo_quarterly.parquet" \
  "${UV_BIN}" run python analysis/01_reproduce_renovo.py \
    --root "${DATA_ROOT}" \
    --start 2020-01-06 \
    --end "${END_DATE}" \
    --output "${RESULTS_ROOT}/renovo_quarterly.parquet"
run_parquet_stage \
  "sampling-bias analysis" \
  "${RESULTS_ROOT}/sampling_bias.parquet" \
  "${UV_BIN}" run python analysis/02_sampling_bias.py \
    --root "${DATA_ROOT}" \
    --output "${RESULTS_ROOT}/sampling_bias.parquet"
run_parquet_stage \
  "transition-matrix analysis" \
  "${RESULTS_ROOT}/transition_matrix.parquet" \
  "${UV_BIN}" run python analysis/03_transition_matrix.py \
    --root "${DATA_ROOT}" \
    --output "${RESULTS_ROOT}/transition_matrix.parquet"
run_parquet_stage \
  "VCV time-to-resolution analysis" \
  "${RESULTS_ROOT}/vus_time_to_resolution.parquet" \
  "${UV_BIN}" run python analysis/04_time_to_resolution.py \
    --root "${DATA_ROOT}" \
    --study-end "${END_DATE}" \
    --output "${RESULTS_ROOT}/vus_time_to_resolution.parquet"
run_parquet_stage \
  "event-mechanisms analysis" \
  "${RESULTS_ROOT}/event_mechanisms.parquet" \
  "${UV_BIN}" run python analysis/05_event_mechanisms.py \
    --root "${DATA_ROOT}" \
    --attributions "${RESULTS_ROOT}/attributions.parquet" \
    --output "${RESULTS_ROOT}/event_mechanisms.parquet"
run_parquet_stage \
  "SCV time-to-resolution analysis" \
  "${RESULTS_ROOT}/scv_vus_time_to_resolution.parquet" \
  "${UV_BIN}" run python analysis/06_scv_time_to_resolution.py \
    --root "${DATA_ROOT}" \
    --study-end "${END_DATE}" \
    --output "${RESULTS_ROOT}/scv_vus_time_to_resolution.parquet"

stage "Full backfill completed"
