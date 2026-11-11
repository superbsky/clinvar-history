#!/usr/bin/env bash
set -euo pipefail
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${PROJECT_ROOT}"
exec "${CLINVAR_UV_BIN:-uv}" run clinvar-history run-study \
  --root "${CLINVAR_DATA_ROOT:-data/full}" \
  --results "${CLINVAR_RESULTS_ROOT:-results/full-v1}" \
  --start "${CLINVAR_START_DATE:-2020-01-01}" \
  --end "${CLINVAR_END_DATE:-2026-07-24}" \
  --max-lag-days "${CLINVAR_MAX_LAG_DAYS:-35}" \
  --attempts "${CLINVAR_RETRY_ATTEMPTS:-3}" \
  --retry-delay "${CLINVAR_RETRY_DELAY_SECONDS:-10}" "$@"
