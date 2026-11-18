process REFERENCE_CHUNK {
    tag "${meta.id}"
    label 'reference_chunk'

    input:
    tuple val(meta), path(chunk)
    val output_root
    path code, stageAs: 'reference_code/bin/*'
    path tools_file, stageAs: 'reference_code/Tools.py'

    output:
    tuple val(meta), path("${meta.id}.json"), emit: receipts
    path 'versions.yml', emit: versions

    script:
    """
    export PATH="\$PWD/reference_code/bin:\$PATH"
    ref_build_pipeline.py chunk --chunk "${chunk}" --cpus ${task.cpus} \
        --checkpoints "${output_root}/batches/${meta.id}/checkpoints" \
        --database "${output_root}/batches/${meta.id}/database" --receipt "${meta.id}.json"
    ref_build_pipeline.py --version > versions.yml
    """

    stub:
    """
    export PATH="\$PWD/reference_code/bin:\$PATH"
    touch "${meta.id}.json"
    ref_build_pipeline.py --version > versions.yml
    """
}
