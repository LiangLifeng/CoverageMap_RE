process RAW_QUERY_SCORE {
    tag "${meta.id}"
    label 'raw_query_score'

    input:
    tuple val(meta), val(payload), path(dependencies, stageAs: 'tracked????/*')
    path code, stageAs: 'query_code/bin/*'

    output:
    tuple val(meta), path("${meta.id}"), emit: results
    path 'versions.yml', emit: versions

    script:
    """
    export PATH="\$PWD/query_code/bin:\$PATH"
    printf '%s' '${payload}' | base64 -d > spec.json
    raw_cram_query.py score --spec spec.json --output "${meta.id}" --failure-audit-dir "${task.workDir}/failure_audit"
    cat <<-VERSIONS > versions.yml
    "${task.process}":
        raw_cram_query.py: \$(raw_cram_query.py --version)
    VERSIONS
    """

    stub:
    """
    export PATH="\$PWD/query_code/bin:\$PATH"
    printf '%s' '${payload}' | base64 -d > spec.json
    raw_cram_query.py score --spec spec.json --output "${meta.id}" --stub
    cat <<-VERSIONS > versions.yml
    "${task.process}":
        raw_cram_query.py: "stub"
    VERSIONS
    """
}

