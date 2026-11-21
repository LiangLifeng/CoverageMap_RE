#!/usr/bin/env nextflow
include { CRAM_INDEX_BUILD } from './workflows/cram_index_mates_build'

workflow {
    if (!params.records || !params.fna || !params.fai)
        error 'Specify --records, --fna and --fai'
    if (params.mate_audit)
        error 'Computed indexing scans mate labels itself; omit --mate_audit'
    if ((params.index_chunk_size as int) < 1 || (params.index_chunk_buckets as int) < 1)
        error 'index_chunk_size and index_chunk_buckets must be positive'
    if (!params.cached_index_code)
        error 'Specify --cached_index_code pointing to the original computed index checkout'
    def scripts = ['cram_index_metrics.py', 'build_cram_index.py', 'raw_cram_pool.py', 'coverage_core.py', '__init__.py']
    CRAM_INDEX_BUILD(
        channel.value(file(params.records, checkIfExists: true)),
        channel.value(file(params.fna, checkIfExists: true)),
        channel.value(file(params.fai, checkIfExists: true)),
        channel.value(scripts.collect { name -> file("${params.cached_index_code}/bin/${name}", checkIfExists: true) }),
        channel.value(scripts.collect { name -> file("${projectDir}/bin/${name}", checkIfExists: true) })
    )
}
