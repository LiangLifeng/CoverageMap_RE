process RAW_QUERY_COLLECT {
    tag 'query-summary'
    label 'raw_query_control'

    input:
    path roster, stageAs: 'roster.json'
    path results, stageAs: 'results/*'
    path code, stageAs: 'query_code/bin/*'

    output:
    path 'collected', emit: summary
    path 'versions.yml', emit: versions

    script:
    def grid = params.jaccard_coverage_thresholds ? "--jaccard-coverage-thresholds ${params.jaccard_coverage_thresholds}" : ''
    """
    export PATH="\$PWD/query_code/bin:\$PATH"
    raw_cram_query.py collect --roster "${roster}" --results-root results --output collected \
        --jaccard-methods "${params.jaccard_methods}" \
        --jaccard-coverage-threshold ${params.jaccard_coverage_threshold} ${grid}
    cat <<-VERSIONS > versions.yml
    "${task.process}":
        raw_cram_query.py: \$(raw_cram_query.py --version)
    VERSIONS
    """
    stub:
    def grid = params.jaccard_coverage_thresholds ? "--jaccard-coverage-thresholds ${params.jaccard_coverage_thresholds}" : ''
    """
    export PATH="\$PWD/query_code/bin:\$PATH"
    raw_cram_query.py collect --roster "${roster}" --results-root results --output collected \
        --jaccard-methods "${params.jaccard_methods}" \
        --jaccard-coverage-threshold ${params.jaccard_coverage_threshold} ${grid}
    cat <<-VERSIONS > versions.yml
    "${task.process}":
        raw_cram_query.py: "stub"
    VERSIONS
    """
}

