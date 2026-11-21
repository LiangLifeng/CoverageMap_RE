include { CRAM_METRICS_PLAN } from '../modules/local/cram_metrics_plan/main'
include { CRAM_METRICS_CHUNK } from '../modules/local/cram_metrics_chunk/main'
include { CRAM_METRICS_MERGE } from '../modules/local/cram_metrics_merge/main'

workflow CRAM_INDEX_BUILD {
    take:
    ch_records
    ch_fna
    ch_fai
    ch_code

    main:
    CRAM_METRICS_PLAN(ch_records, ch_fna, ch_fai, ch_code)
    chunks = CRAM_METRICS_PLAN.out.bundle.flatMap { bundle ->
        file("${bundle}/tasks.tsv").readLines().drop(1).collect { line ->
            def fields = line.split('\t', -1)
            tuple([id: fields[0]], file("${bundle}/${fields[1]}", checkIfExists: true), fields[2].split('\\|').collect { source -> file(source, checkIfExists: true) })
        }
    }
    CRAM_METRICS_CHUNK(chunks, ch_code)
    // Only small scalar SQLite shards are collected, never CRAM or coverage arrays.
    CRAM_METRICS_MERGE(CRAM_METRICS_PLAN.out.plan, CRAM_METRICS_CHUNK.out.metrics.collect(), ch_fna, ch_fai, ch_code)

    emit:
    index = CRAM_METRICS_MERGE.out.index
    versions = CRAM_METRICS_PLAN.out.versions.mix(CRAM_METRICS_CHUNK.out.versions, CRAM_METRICS_MERGE.out.versions)
}
