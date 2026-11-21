process CRAM_METRICS_MERGE {
    tag 'cram-index-merge'
    label 'cram_metrics_merge'
    input:
    path records, stageAs: 'build_input.csv'
    path shards, stageAs: 'shards/*'
    path fna, stageAs: 'fna.csv'
    path fai, stageAs: 'fai_path.csv'
    path code, stageAs: 'metric_code/bin/*'
    output:
    path 'reference_cram.sqlite', emit: index
    path 'versions.yml', emit: versions
    script:
    """
    export PATH="\$PWD/metric_code/bin:\$PATH"
    cram_index_metrics.py merge --records "${records}" --chunk-dir shards --fna "${fna}" --fai "${fai}" --output reference_cram.sqlite --min-mapq ${params.index_min_mapq} --m5-policy ${params.index_m5_policy}
    cat <<-VERSIONS > versions.yml
    "${task.process}":
        cram_index_metrics.py: \$(cram_index_metrics.py --version)
    VERSIONS
    """
    stub:
    """
    touch reference_cram.sqlite
    echo 'cram_index_metrics.py: stub' > versions.yml
    """
}
