# Sample–GCF 查询导出

`bin/reference_export.py` 仅使用 Python 标准库，只读查询当前 schema-3 SQLite。
不访问原始 CRAM/BAM、geometry 或外部 skipped_samples.tsv，不重建或重新筛选索引。

参数：

- `--index FILE`：已有 SQLite，必填。
- `--out FILE`：新 TSV 路径，必填，父目录需存在；拒绝覆盖已有文件或符号链接。
- `--sample ID`、`--gcf ID`：区分大小写、完整 ID 精确匹配，各可重复；同类 OR，不同条件 AND。省略表示不限。
- `--min-depth N`、`--max-depth N`：原始完整参考基因组平均深度 `full_mean_depth`，含两端；可单独指定。必须有限、非负，且下限不高于上限。
- `--selected-only`：只取当前索引 `selected=1`。默认导出全部有效登记记录，包括筛后未保留者；不是重新应用深度/最少标本数条件，也不能找回建索引时已跳过的错误样本。

输出按 Sample、GCFID 的二进制字典序稳定排序，一行一个 Sample–GCF：
`Sample GCFID Depth CramFilePath DatabaseRoot Selected EligibleDepth ExclusionReason`（制表符分隔）。
Selected 为当前索引保留标记；EligibleDepth 仅为当前索引建库筛选深度达标标记，不能替代 Selected，也不随本次查询区间变化。
DatabaseRoot 为此记录主来源；重复副本不会重复输出。CramFilePath 为 builder 的 sources[0].path 登记值，CRAM/BAM 以外或缺失时留空，终端 JSON 的 missing_alignment_paths 汇总缺失数，不承诺路径当前可读。
零匹配成功生成仅表头 TSV；终端 JSON 包含记录数、去重 Sample/GCF 数、查询条件、原索引筛选参数及 plan SHA。
TSV 流式写临时文件，成功后排他原子发布；错误不发布半份 TSV。连接 mode=ro/query_only，单一读事务保证元数据快照。未知/旧 schema 显式拒绝，不猜测筛选状态。

NH 示例（在 NH shell 执行，变量赋值一次）：

```bash
PYTHON=/home/lianglifeng/01.WorkDir/env/conda-envs/coverage-map/bin/python
EXPORT=/home/lianglifeng/01.WorkDir/01.Project/03.CoverageMap/source/bin/reference_export.py
INDEX=/home/lianglifeng/share/lianglifeng/CoverageMap/ref_depth_database_build/RefPool/RefPool_IH_MAG_vOTU_V20261022_all.sqlite
"$PYTHON" "$EXPORT" --index "$INDEX" --min-depth 1 --max-depth 5 --out "$PWD/Sample_GCF_1to5X.tsv"
"$PYTHON" "$EXPORT" --index "$INDEX" --sample XNYY01-A035 --out "$PWD/XNYY01-A035.tsv"
"$PYTHON" "$EXPORT" --index "$INDEX" --gcf GCF_002082765.1 --out "$PWD/GCF_002082765.1.tsv"
"$PYTHON" "$EXPORT" --index "$INDEX" --sample XNYY01-A035 --gcf GCF_002082765.1 --min-depth 1 --max-depth 5 --out "$PWD/combined.tsv"
"$PYTHON" "$EXPORT" --index "$INDEX" --sample XNYY01-A035 --sample XNYY01-A036 --out "$PWD/two_samples.tsv"
```

只查已有 5×/3 筛选库的保留记录：

```bash
"$PYTHON" "$EXPORT" --index /home/lianglifeng/share/lianglifeng/CoverageMap/ref_depth_database_build/RefPool/RefPool_IH_MAG_vOTU_V20261022_d3s3.sqlite --selected-only --out "$PWD/retained.tsv"
```

注意文件名 d3s3 的已核实参数实际为深度≥5×、达标样本≥3；以程序输出 index_selection 为准。
在 all.sqlite 上加 selected-only 读取的是 all.sqlite 自身标记，不会自动切换到过滤库。
