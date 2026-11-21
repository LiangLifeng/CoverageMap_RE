process CRAM_METRICS_PLAN {
    tag 'cram-index-plan'
    label 'cram_metrics_control'
    input:
    path records, stageAs: 'build_input.csv'
    path fna, stageAs: 'fna.csv'
    path fai, stageAs: 'fai_path.csv'
    path code, stageAs: 'metric_code/bin/*'
    output:
    path 'plan', emit: bundle
    path 'plan/plan.json', emit: plan
    path 'versions.yml', emit: versions
    script:
    """
    export PATH="\$PWD/metric_code/bin:\$PATH"
    cram_index_metrics.py plan --records "${records}" --fna "${fna}" --fai "${fai}" --chunk-size ${params.index_chunk_size} --outdir plan
    cat <<-VERSIONS > versions.yml
    "${task.process}":
        cram_index_metrics.py: \$(cram_index_metrics.py --version)
    VERSIONS
    """
    stub:
    """
    mkdir -p plan/chunks
    touch plan/chunks/chunk-000001.jsonl plan/plan.json
    printf 'id\tchunk\tdependencies\nchunk-000001\tchunks/chunk-000001.jsonl\t%s\n' "\$(readlink -f '${fna}')" > plan/tasks.tsv
    echo 'cram_index_metrics.py: stub' > versions.yml
    """
}
