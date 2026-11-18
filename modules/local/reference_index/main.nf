process REFERENCE_INDEX {
    tag 'single-writer-sqlite'
    label 'reference_collect'

    input:
    path collection_receipt
    val output_root
    path code, stageAs: 'reference_code/bin/*'
    path tools_file, stageAs: 'reference_code/Tools.py'

    output:
    path 'reference_build.json', emit: receipt
    path 'versions.yml', emit: versions

    script:
    """
    export PATH="\$PWD/reference_code/bin:\$PATH"
    ref_build_pipeline.py index --database "${output_root}/database" \
        --cpus ${task.cpus} --receipt reference_build.json
    ref_build_pipeline.py --version > versions.yml
    """

    stub:
    """
    export PATH="\$PWD/reference_code/bin:\$PATH"
    touch reference_build.json
    ref_build_pipeline.py --version > versions.yml
    """
}
