#!/usr/bin/env nextflow
include { REFERENCE_BUILD } from './workflows/reference_build'

workflow {
    if (!params.records || !params.references)
        error 'Specify complete --records and --references CSV files'
    if ([params.chunk_size, params.cpus_per_chunk, params.collect_cpus, params.plan_cpus].any { value -> value.toInteger() < 1 })
        error 'Chunk size, concurrency and CPU budgets must be positive'
    def output_root = file(params.outdir).toAbsolutePath().normalize().toString()
    if (output_root == projectDir.toString() || output_root.startsWith(projectDir.toString() + '/'))
        error 'Reference output must be outside the source checkout'
    if (params.max_chunks != null)
        log.warn '--max_chunks is retired and ignored; use executor.queueSize or -qs for this run'
    log.info "Reference build: CPUs per chunk=${params.cpus_per_chunk}; submissions use native Nextflow queue settings"
    log.info "Preprocessing full input into whole-GCF byte-balanced batches; target records/batch=${params.chunk_size}"
    REFERENCE_BUILD(
        channel.value(file(params.records, checkIfExists: true)),
        channel.value(file(params.references, checkIfExists: true)),
        channel.value(output_root),
        channel.value(file("${projectDir}/bin/*.py", checkIfExists: true)),
        channel.value(file("${projectDir}/Tools.py", checkIfExists: true))
    )
}
