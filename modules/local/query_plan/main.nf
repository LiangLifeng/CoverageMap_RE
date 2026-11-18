process QUERY_PLAN {
    tag 'manifest'
    label 'query_control'

    input:
    path manifest, stageAs: 'input_manifest.tsv'
    path index, stageAs: 'reference.sqlite'
    path code, stageAs: 'query_code/bin/*'
    path tools_file, stageAs: 'query_code/Tools.py'

    output:
    path 'query_tasks.tsv', emit: tasks
    path 'query_roster.json', emit: roster
    path 'versions.yml', emit: versions

    script:
    def coverage = params.coverage_unit ? "--coverage-unit ${params.coverage_unit}" : ''
    """
    export PATH="\$PWD/query_code/bin:\$PATH"
    coverage_query_flow.py plan --manifest "${manifest}" --index "${index}" \
        --tasks query_tasks.tsv --roster query_roster.json ${coverage}
    cat <<-VERSIONS > versions.yml
    "${task.process}":
        coverage_query_flow.py: \$(coverage_query_flow.py --version)
    VERSIONS
    """

    stub:
    def coverage = params.coverage_unit ? "--coverage-unit ${params.coverage_unit}" : ''
    """
    export PATH="\$PWD/query_code/bin:\$PATH"
    coverage_query_flow.py plan --manifest "${manifest}" --index "${index}" \
        --tasks query_tasks.tsv --roster query_roster.json ${coverage}
    cat <<-VERSIONS > versions.yml
    "${task.process}":
        coverage_query_flow.py: "stub"
    VERSIONS
    """
}
