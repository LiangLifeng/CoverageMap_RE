#!/usr/bin/env nextflow
include { CRAM_INDEX_BUILD } from './workflows/cram_index_build'

workflow {
    if (!params.records || !params.fna || !params.fai)
        error 'Specify --records, --fna and --fai'
    if (params.mate_audit)
        error 'Computed indexing scans mate labels itself; omit --mate_audit'
    CRAM_INDEX_BUILD(
        channel.value(file(params.records, checkIfExists: true)),
        channel.value(file(params.fna, checkIfExists: true)),
        channel.value(file(params.fai, checkIfExists: true)),
        channel.value(['cram_index_metrics.py', 'build_cram_index.py', 'raw_cram_pool.py', 'coverage_core.py', '__init__.py'].collect { name -> file("${projectDir}/bin/${name}", checkIfExists: true) })
    )
}
