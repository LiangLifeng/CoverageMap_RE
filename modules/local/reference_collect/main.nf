process REFERENCE_COLLECT {
    tag 'single-writer-collector'
    label 'reference_collect'

    input:
    path plan
    path receipts, stageAs: 'receipts/*'
    val output_root
    path code, stageAs: 'reference_code/bin/*'
    path tools_file, stageAs: 'reference_code/Tools.py'

    output:
    path 'collected.json', emit: receipt
    path 'versions.yml', emit: versions

    script:
    """
    export PATH="\$PWD/reference_code/bin:\$PATH"
    ref_build_pipeline.py collect --plan "${plan}" --receipts receipts \
        --database "${output_root}/database" --cpus ${task.cpus} --receipt collected.json
    ref_build_pipeline.py --version > versions.yml
    """

    stub:
    """
    export PATH="\$PWD/reference_code/bin:\$PATH"
    touch collected.json
    ref_build_pipeline.py --version > versions.yml
    """
}
