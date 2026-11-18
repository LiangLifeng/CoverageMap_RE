process REFERENCE_PLAN {
    tag 'global-plan'
    label 'reference_control'

    input:
    path records, stageAs: 'records.csv'
    path references, stageAs: 'references.csv'
    path code, stageAs: 'reference_code/bin/*'
    path tools_file, stageAs: 'reference_code/Tools.py'

    output:
    path 'plan', emit: plan
    path 'plan/chunks/chunk-*.json', emit: chunks
    path 'versions.yml', emit: versions

    script:
    def bins = params.reference_bin_count == null ? "--window-size ${params.reference_window_size.toInteger()}" : "--bin-count ${params.reference_bin_count.toInteger()}"
    """
    export PATH="\$PWD/reference_code/bin:\$PATH"
    ref_build_pipeline.py plan --records "${records}" --references "${references}" \
        --outdir plan --chunk-size ${params.chunk_size.toInteger()} --cpus ${task.cpus} \
        --smoke-records ${params.smoke_records.toInteger()} ${bins}
    ref_build_pipeline.py --version > versions.yml
    """

    stub:
    """
    export PATH="\$PWD/reference_code/bin:\$PATH"
    mkdir -p plan/chunks
    touch plan/chunks/chunk-000000.json
    ref_build_pipeline.py --version > versions.yml
    """
}
