#!/usr/bin/env bash
# 原始参考库：不预建抽样深度。先填输出目录；默认仅预览。
set -euo pipefail

# ==================== 用户填写区 ====================
CPUS="${CPUS:-4}"
COLLECT_JOBS="${COLLECT_JOBS:-4}" # 最终校验并发上限
MEMORY="${MEMORY:-}"             # 留空：不显式申请内存
DATABASE="${DATABASE:-}"         # 必填：尚不存在的绝对输出目录
WINDOW_SIZE="${WINDOW_SIZE:-1000}" # bp，每条contig独立划窗，保留末尾短窗口
JOBS="${JOBS:-20}"
DRY_RUN="${DRY_RUN:-1}"           # 1仅生成计划/脚本；0执行并构建SQLite
PARTITION="${PARTITION:-cpu}"
PROJECT="${PROJECT:-/home/lianglifeng/01.WorkDir/01.Project/03.CoverageMap}"
PYTHON="${PYTHON:-/home/lianglifeng/01.WorkDir/env/conda-envs/coverage-map/bin/python}"
RECORDS="${RECORDS:-$PROJECT/01.demo_data/DemoCramRecord_add_GCF.csv}"
REFERENCES="${REFERENCES:-$PROJECT/01.demo_data/fai_path.csv}"
SQLITE_INDEX="${SQLITE_INDEX:-${DATABASE}.sqlite}"
# ==================== 配置结束 ====================
: "${DATABASE:?请填写DATABASE}"
[[ "$DATABASE" = /* && "$SQLITE_INDEX" = /* ]] || { echo '输出必须为绝对路径' >&2; exit 2; }
[[ ! -e "$DATABASE" && ! -L "$DATABASE" && ! -e "$SQLITE_INDEX" && ! -L "$SQLITE_INDEX" ]] || { echo '输出已存在，请使用新目录和新索引路径' >&2; exit 2; }
for value in "$CPUS" "$JOBS" "$WINDOW_SIZE" "$COLLECT_JOBS"; do
  [[ "$value" =~ ^[1-9][0-9]*$ ]] || { echo 'CPU/JOBS/WINDOW_SIZE必须是正整数' >&2; exit 2; }
done
[[ "$DRY_RUN" = 0 || "$DRY_RUN" = 1 ]] || { echo 'DRY_RUN必须为0或1' >&2; exit 2; }
OPTIONS=()
[[ -z "$MEMORY" ]] || OPTIONS+=(--memory "$MEMORY")
[[ "$DRY_RUN" = 0 ]] || OPTIONS+=(--dry-run)
"$PYTHON" -B "$PROJECT/source/bin/ref_depth_database.py" run --mode slurm \
  --records "$RECORDS" --references "$REFERENCES" --outdir "$DATABASE" \
  --cpus "$CPUS" --jobs "$JOBS" --collect-jobs "$COLLECT_JOBS" --window-size "$WINDOW_SIZE" \
  --partition "$PARTITION" --geometry-format packed "${OPTIONS[@]}"
if [[ "$DRY_RUN" = 1 ]]; then
  printf '仅预览，未执行/投递：%s/dry_run.json\n' "$DATABASE"
  exit 0
fi
"$PYTHON" -B "$PROJECT/source/bin/reference_index.py" build \
  --database "$DATABASE" --out "$SQLITE_INDEX" --verify-payloads
printf '下游读取：--reference-sqlite %s\n资源：%s/RESOURCE_INDEX.tsv\n' "$SQLITE_INDEX" "$DATABASE"
# 同一计划续跑：execute --database "$DATABASE" --mode local --jobs "$JOBS"
# Slurm预览后的正式运行：submit --database "$DATABASE" --jobs "$JOBS"
# 然后 status --database "$DATABASE" --watch --scheduler --refill --jobs "$JOBS"
# 续跑完成后 collect，再执行上方reference_index.py build命令。
