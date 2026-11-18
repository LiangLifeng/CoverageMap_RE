process QUERY_COLLECT {
    tag 'all queries'
    label 'query_control'

    input:
    path roster
    path receipts, stageAs: 'receipts/*'
    path code, stageAs: 'query_code/bin/*'
    path tools_file, stageAs: 'query_code/Tools.py'

    output:
    path 'collected', emit: results
    path 'versions.yml', emit: versions

    script:
    def grid = params.jaccard_coverage_thresholds ? "--jaccard-coverage-thresholds ${params.jaccard_coverage_thresholds}" : ''
    """
    export PATH="\$PWD/query_code/bin:\$PATH"
    coverage_query_flow.py collect --roster "${roster}" --receipts receipts --outdir collected --jaccard-methods "${params.jaccard_methods}" --jaccard-coverage-threshold ${params.jaccard_coverage_threshold.toDouble()} ${grid}
    cat <<-VERSIONS > versions.yml
    "${task.process}":
        coverage_query_flow.py: \$(coverage_query_flow.py --version)
    VERSIONS
    """

    stub:
    def grid = params.jaccard_coverage_thresholds ? "--jaccard-coverage-thresholds ${params.jaccard_coverage_thresholds}" : ''
    """
    export PATH="\$PWD/query_code/bin:\$PATH"
    coverage_query_flow.py collect --roster "${roster}" --receipts receipts --outdir collected --jaccard-methods "${params.jaccard_methods}" --jaccard-coverage-threshold ${params.jaccard_coverage_threshold.toDouble()} ${grid}
    cat <<-VERSIONS > versions.yml
    "${task.process}":
        coverage_query_flow.py: "stub"
    VERSIONS
    """
}
