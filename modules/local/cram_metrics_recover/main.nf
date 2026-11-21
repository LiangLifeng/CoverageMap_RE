process CRAM_METRICS_RECOVER {
    tag "${meta.id}"
    label 'cram_metrics_recover'
    input:
    tuple val(meta), val(rows), path(dependencies, stageAs: 'dependencies????/*')
    path code, stageAs: 'metric_code/bin/*'
    output:
    tuple val(meta.id), path("${meta.id}.sqlite"), emit: metrics
    path 'versions.yml', emit: versions
    script:
    def payload = groovy.json.JsonOutput.toJson(rows)
    """
    export PATH="\$PWD/metric_code/bin:\$PATH"
    cat > rows.json <<'ROWS_JSON'
    ${payload}
    ROWS_JSON
    cram_index_metrics.py compute --strict-mates --chunk rows.json --output "${meta.id}.sqlite" --jobs ${task.cpus} --min-mapq ${params.index_min_mapq} --m5-policy ${params.index_m5_policy}
    cat <<-VERSIONS > versions.yml
    "${task.process}":
        cram_index_metrics.py: \$(cram_index_metrics.py --version)
    VERSIONS
    """
    stub:
    """
    touch "${meta.id}.sqlite"
    echo 'cram_index_metrics.py: stub' > versions.yml
    """
}
