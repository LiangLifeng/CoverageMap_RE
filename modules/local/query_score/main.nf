process QUERY_SCORE {
    tag "${meta.id}"
    label 'query_worker'

    input:
    val meta
    path manifest, stageAs: 'input_manifest.tsv'
    path reference_index, stageAs: 'reference_index.sqlite'
    path code, stageAs: 'query_code/bin/*'
    path tools_file, stageAs: 'query_code/Tools.py'

    output:
    tuple val(meta.id), path("${meta.id}.json"), emit: receipts
    tuple val(meta), path('query_result'), optional: true, emit: details
    path 'versions.yml', emit: versions

    script:
    def windows = params.window_size == null ? "--bin-count ${params.bin_count.toInteger()}" : "--window-size ${params.window_size.toInteger()}"
    def sampling = params.sampling_unit ? "--sampling-unit ${params.sampling_unit}" : ''
    def alignments = params.coverage_alignment_policy ? "--coverage-alignment-policy ${params.coverage_alignment_policy}" : ''
    def grid = params.jaccard_coverage_thresholds ? "--jaccard-coverage-thresholds ${params.jaccard_coverage_thresholds}" : ''
    def coverage = params.coverage_unit ? "--coverage-unit ${params.coverage_unit}" : ''
    """
    export PATH="\$PWD/query_code/bin:\$PATH"
    coverage_query_flow.py score \
        --manifest "${manifest}" --index "${reference_index}" \
        --query-id "${meta.id}" --input-stamp "${meta.input_stamp}" \
        --jaccard-coverage-threshold ${params.jaccard_coverage_threshold.toDouble()} ${grid} \
        --jaccard-methods "${params.jaccard_methods}" --bootstrap-replicates ${params.bootstrap_replicates.toInteger()} \
        --bootstrap-confidence ${params.bootstrap_confidence.toDouble()} --reference-draws ${params.reference_draws.toInteger()} \
        --threads ${task.cpus} --alpha ${params.alpha.toDouble()} ${windows} ${coverage} ${sampling} ${alignments} \
        --depth-threshold ${params.depth_threshold.toDouble()} --min-samples ${params.min_samples.toInteger()} \
        --seed ${params.seed.toInteger()} --cache-mib ${params.cache_mib.toInteger()} --work-mib ${params.work_mib.toInteger()} \
        --query-outdir query_result --receipt "${meta.id}.json"
    cat <<-VERSIONS > versions.yml
    "${task.process}":
        coverage_query_flow.py: \$(coverage_query_flow.py --version)
        query_manifest_v2.py: \$(query_manifest_v2.py --version)
    VERSIONS
    """

    stub:
    def windows = params.window_size == null ? "--bin-count ${params.bin_count.toInteger()}" : "--window-size ${params.window_size.toInteger()}"
    def sampling = params.sampling_unit ? "--sampling-unit ${params.sampling_unit}" : ''
    def alignments = params.coverage_alignment_policy ? "--coverage-alignment-policy ${params.coverage_alignment_policy}" : ''
    def grid = params.jaccard_coverage_thresholds ? "--jaccard-coverage-thresholds ${params.jaccard_coverage_thresholds}" : ''
    def coverage = params.coverage_unit ? "--coverage-unit ${params.coverage_unit}" : ''
    """
    export PATH="\$PWD/query_code/bin:\$PATH"
    coverage_query_flow.py score \
        --manifest "${manifest}" --index "${reference_index}" \
        --query-id "${meta.id}" --input-stamp "${meta.input_stamp}" \
        --jaccard-coverage-threshold ${params.jaccard_coverage_threshold.toDouble()} ${grid} \
        --jaccard-methods "${params.jaccard_methods}" --bootstrap-replicates ${params.bootstrap_replicates.toInteger()} \
        --bootstrap-confidence ${params.bootstrap_confidence.toDouble()} --reference-draws ${params.reference_draws.toInteger()} \
        --threads ${task.cpus} --alpha ${params.alpha.toDouble()} ${windows} ${coverage} ${sampling} ${alignments} \
        --depth-threshold ${params.depth_threshold.toDouble()} --min-samples ${params.min_samples.toInteger()} \
        --seed ${params.seed.toInteger()} --cache-mib ${params.cache_mib.toInteger()} --work-mib ${params.work_mib.toInteger()} \
        --query-outdir query_result --receipt "${meta.id}.json" --stub
    cat <<-VERSIONS > versions.yml
    "${task.process}":
        coverage_query_flow.py: "stub"
    VERSIONS
    """
}
