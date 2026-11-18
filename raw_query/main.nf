#!/usr/bin/env nextflow
include { RAW_CRAM_QUERIES } from '../workflows/raw_cram_queries'

workflow {
    if (!params.manifest || !params.index || params.alpha == null)
        error 'Specify --manifest, --index and --alpha'
    if (params.fai_map != null && params.fai_dir != null)
        error 'Specify at most one of --fai_map or --fai_dir; SQL FNA/FAI catalog needs neither'
    if ((params.bin_count == null) == (params.window_size == null))
        error 'Specify exactly one of --bin_count or --window_size'
    if (!(params.diagnostics.toString().toLowerCase() in ['true', 'false']))
        error '--diagnostics must be true or false'
    if (params.diagnostics_cpus.toInteger() < 1 || params.query_cpus.toInteger() < 1)
        error 'CPU parameters must be positive'
    RAW_CRAM_QUERIES(
        channel.value(file(params.manifest, checkIfExists: true)),
        channel.value(file(params.index, checkIfExists: true)),
        channel.value(files("${projectDir}/../bin/*.py", checkIfExists: true)),
        channel.value(params.fai_map ? [file(params.fai_map, checkIfExists: true)] : []),
        channel.value(file("${projectDir}/../bin/diagnostics/query_similarity_plot.py", checkIfExists: true))
    )
}

