#!/usr/bin/env nextflow
include { CRAM_INDEX_BUILD } from './workflows/cram_index_build'

workflow {
    if (!params.records)
        error 'Specify --records with the BuildInput CSV'
    if (!params.fna || !params.fai)
        error 'Specify --fna and --fai with the reference mapping CSVs'
    CRAM_INDEX_BUILD(
        channel.value(file(params.records, checkIfExists: true)),
        channel.value(file(params.fna, checkIfExists: true)),
        channel.value(file(params.fai, checkIfExists: true)),
        channel.value(file("${projectDir}/bin/build_cram_index.py", checkIfExists: true))
    )
}

