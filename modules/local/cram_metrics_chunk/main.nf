process CRAM_METRICS_CHUNK {
    tag "${meta.id}"
    label 'cram_metrics_chunk'
    input:
    tuple val(meta), path(chunk), path(dependencies, stageAs: 'dependencies????/*')
    path code, stageAs: 'metric_code/bin/*'
    output:
    path "${meta.id}.sqlite", emit: metrics
    path 'versions.yml', emit: versions
    script:
    """
    export PATH="\$PWD/metric_code/bin:\$PATH"
    cram_index_metrics.py compute --chunk "${chunk}" --output "${meta.id}.sqlite" --jobs ${task.cpus} --min-mapq ${params.index_min_mapq} --m5-policy ${params.index_m5_policy}
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
