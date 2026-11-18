#!/usr/bin/env bash
# DATABASE=单库 或 DATABASES_FILE=逐行绝对路径清单；多库必须显式SQLITE_INDEX。
set -euo pipefail
DATABASE="${DATABASE:-}"
DATABASES_FILE="${DATABASES_FILE:-}"
SQLITE_INDEX="${SQLITE_INDEX:-${DATABASE:+${DATABASE}.sqlite}}"
MIN_DEPTH="${MIN_DEPTH:-0}"
MIN_SAMPLES="${MIN_SAMPLES:-1}"
VERIFY_PAYLOADS="${VERIFY_PAYLOADS:-1}"
SKIP_INVALID_SAMPLES="${SKIP_INVALID_SAMPLES:-0}"
SKIPPED_SAMPLES_TSV="${SKIPPED_SAMPLES_TSV:-}"
PROJECT="${PROJECT:-/home/lianglifeng/01.WorkDir/01.Project/03.CoverageMap}"
PYTHON="${PYTHON:-/home/lianglifeng/01.WorkDir/env/conda-envs/coverage-map/bin/python}"
INDEXER="$PROJECT/source/bin/reference_index.py"
[[ -n "$SQLITE_INDEX" && "$SQLITE_INDEX" = /* ]] || { echo '请设置绝对路径SQLITE_INDEX' >&2; exit 2; }
[[ "$VERIFY_PAYLOADS" = 0 || "$VERIFY_PAYLOADS" = 1 ]] || { echo 'VERIFY_PAYLOADS必须为0或1' >&2; exit 2; }
[[ "$SKIP_INVALID_SAMPLES" = 0 || "$SKIP_INVALID_SAMPLES" = 1 ]] || { echo 'SKIP_INVALID_SAMPLES必须为0或1' >&2; exit 2; }
options=()
if [[ -n "$DATABASES_FILE" ]]; then
  [[ -z "$DATABASE" ]] || { echo 'DATABASE与DATABASES_FILE只能设置一个' >&2; exit 2; }
  while IFS= read -r database || [[ -n "$database" ]]; do
    [[ -n "$database" ]] || continue
    [[ "$database" = /* ]] || { echo '库清单中必须是绝对路径' >&2; exit 2; }
    options+=(--database "$database")
  done < "$DATABASES_FILE"
else
  [[ "$DATABASE" = /* ]] || { echo '请设置绝对路径DATABASE或DATABASES_FILE' >&2; exit 2; }
  options+=(--database "$DATABASE")
fi
(( ${#options[@]} > 0 )) || { echo '库清单为空' >&2; exit 2; }
if [[ "$VERIFY_PAYLOADS" = 1 ]]; then options+=(--verify-payloads); fi
if [[ "$SKIP_INVALID_SAMPLES" = 1 ]]; then options+=(--skip-invalid-samples); fi
if [[ -n "$SKIPPED_SAMPLES_TSV" ]]; then options+=(--skipped-samples-tsv "$SKIPPED_SAMPLES_TSV"); fi
"$PYTHON" -B "$INDEXER" build "${options[@]}" --out "$SQLITE_INDEX" --min-depth "$MIN_DEPTH" --min-samples "$MIN_SAMPLES"
"$PYTHON" -B "$INDEXER" inspect --index "$SQLITE_INDEX"
"$PYTHON" -B "$INDEXER" report --index "$SQLITE_INDEX" --outdir "${SQLITE_INDEX}.report"
