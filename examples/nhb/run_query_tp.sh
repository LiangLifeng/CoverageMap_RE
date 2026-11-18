#!/usr/bin/env bash
# User-owned formal execution. Append -bg for native background mode, -resume to recover.
set -euo pipefail
export PATH="/home/lianglifeng/01.WorkDir/env/conda-envs/coverage-map/bin:$PATH"
export NXF_SYNTAX_PARSER=v2
unset NXF_BIN NXF_DIST NXF_LAUNCHER
export NXF_HOME=/home/lianglifeng/01.WorkDir/env/app/nextflow-official-26.04.6/runtime
export NXF_VER=26.04.6
CODE=/home/lianglifeng/01.WorkDir/01.Project/03.CoverageMap/runs/coveragemap-re-vnext-20261118/code
mkdir -p /home/lianglifeng/01.WorkDir/01.Project/03.CoverageMap/runs/coveragemap-re-vnext-20261118/launch_tp
cd /home/lianglifeng/01.WorkDir/01.Project/03.CoverageMap/runs/coveragemap-re-vnext-20261118/launch_tp
/home/lianglifeng/01.WorkDir/env/app/nextflow-official-26.04.6/bin/nextflow \
    -C "$CODE/raw_query/nextflow.config,$CODE/examples/nhb/raw_query.config" \
    run "$CODE/raw_query/main.nf" \
    -profile slurm,nhb -qs 100 \
    --index /home/lianglifeng/share/lianglifeng/CoverageMap/Data/step.1.ref_depth_database/0.Raw_Cram_SQL/reference_cram.sqlite \
    --manifest /home/share/lianglifeng/CoverageMap/Data/TestPool/query_manifests_viral_G1000_S10_tp_fp/TP_ClinicalVerified_FeverPatients_nonreference.tsv \
    --outdir /home/lianglifeng/01.WorkDir/01.Project/03.CoverageMap/runs/coveragemap-re-vnext-20261118/query_tp \
    -work-dir /home/lianglifeng/01.WorkDir/01.Project/03.CoverageMap/runs/coveragemap-re-vnext-20261118/work_tp \
    --geometry_cache_dir /home/lianglifeng/01.WorkDir/01.Project/03.CoverageMap/runs/coveragemap-re-vnext-20261118/shared_geometry \
    --sampling_unit read --sampling_tolerance_bp 150 --reference_draws 100 --bin_count 60 \
    --alpha 0.05 --coverage_unit percent --coverage_alignment_policy primary_only \
    --min_samples 9 --depth_threshold 5 --seed 401 --query_cpus 4 \
    --jaccard_methods minimum --bootstrap_replicates 1000 --bootstrap_confidence 0.95 \
    --jaccard_coverage_threshold 0 --diagnostics true "$@"
