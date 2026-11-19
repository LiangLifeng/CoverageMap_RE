process RAW_QUERY_PLAN {
    tag 'query-plan'
    label 'raw_query_control'

    input:
    path manifest, stageAs: 'queries.tsv'
    path index, stageAs: 'reference.sqlite'
    path code, stageAs: 'query_code/bin/*'
    path fai_map, stageAs: 'fai_map/*'

    output:
    path 'plan/tasks.tsv', emit: tasks
    path 'plan/roster.json', emit: roster
    path 'versions.yml', emit: versions

    script:
    def fai = params.fai_map ? "--fai-map '${params.fai_map}'" : (params.fai_dir ? "--fai-dir '${params.fai_dir}'" : '')
    def coverage = params.coverage_unit ? "--coverage-unit ${params.coverage_unit}" : ''
    def grid = params.jaccard_coverage_thresholds ? "--jaccard-coverage-thresholds ${params.jaccard_coverage_thresholds}" : ''
    def windows = params.bin_count != null ? "--bin-count ${params.bin_count}" : "--window-size ${params.window_size}"
    def geometryCache = new File(params.geometry_cache_dir ?: "${params.outdir}/geometry_cache").absolutePath
    """
    export PATH="\$PWD/query_code/bin:\$PATH"
    raw_cram_query.py plan --manifest "${manifest}" --index "${index}" --output plan \
        ${fai} ${windows} --alpha ${params.alpha} --reference-draws ${params.reference_draws} \
        --seed ${params.seed} --depth-threshold ${params.depth_threshold} --min-reference-depth ${params.min_reference_depth} \
        --sampling-tolerance-bp ${params.sampling_tolerance_bp} \
        --insufficient-reference-policy ${params.insufficient_reference_policy} \
        --geometry-cache-dir "${geometryCache}" \
        --sampling-unit ${params.sampling_unit} --min-mapq ${params.min_mapq} \
        --min-samples ${params.min_samples} --threads ${params.query_cpus} \
        --jaccard-coverage-threshold ${params.jaccard_coverage_threshold} ${grid} ${coverage} \
        --coverage-alignment-policy ${params.coverage_alignment_policy} \
        --jaccard-methods "${params.jaccard_methods}" \
        --bootstrap-replicates ${params.bootstrap_replicates} --bootstrap-confidence ${params.bootstrap_confidence}
    cat <<-VERSIONS > versions.yml
    "${task.process}":
        raw_cram_query.py: \$(raw_cram_query.py --version)
    VERSIONS
    """
    stub:
    def fai = params.fai_map ? "--fai-map '${params.fai_map}'" : (params.fai_dir ? "--fai-dir '${params.fai_dir}'" : '')
    def coverage = params.coverage_unit ? "--coverage-unit ${params.coverage_unit}" : ''
    def grid = params.jaccard_coverage_thresholds ? "--jaccard-coverage-thresholds ${params.jaccard_coverage_thresholds}" : ''
    def windows = params.bin_count != null ? "--bin-count ${params.bin_count}" : "--window-size ${params.window_size}"
    def geometryCache = new File(params.geometry_cache_dir ?: "${params.outdir}/geometry_cache").absolutePath
    """
    export PATH="\$PWD/query_code/bin:\$PATH"
    raw_cram_query.py plan --manifest "${manifest}" --index "${index}" --output plan \
        ${fai} ${windows} --alpha ${params.alpha} --reference-draws ${params.reference_draws} \
        --seed ${params.seed} --depth-threshold ${params.depth_threshold} --min-reference-depth ${params.min_reference_depth} \
        --sampling-tolerance-bp ${params.sampling_tolerance_bp} \
        --insufficient-reference-policy ${params.insufficient_reference_policy} \
        --geometry-cache-dir "${geometryCache}" \
        --sampling-unit ${params.sampling_unit} --min-mapq ${params.min_mapq} \
        --min-samples ${params.min_samples} --threads ${params.query_cpus} \
        --jaccard-coverage-threshold ${params.jaccard_coverage_threshold} ${grid} ${coverage} \
        --coverage-alignment-policy ${params.coverage_alignment_policy} \
        --jaccard-methods "${params.jaccard_methods}" \
        --bootstrap-replicates ${params.bootstrap_replicates} --bootstrap-confidence ${params.bootstrap_confidence}
    cat <<-VERSIONS > versions.yml
    "${task.process}":
        raw_cram_query.py: "stub"
    VERSIONS
    """
}

