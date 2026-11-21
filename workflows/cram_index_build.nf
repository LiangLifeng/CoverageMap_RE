include { CRAM_METRICS_CHUNK } from '../modules/local/cram_metrics_chunk/main'
include { CRAM_METRICS_MERGE } from '../modules/local/cram_metrics_merge/main'

workflow CRAM_INDEX_BUILD {
    take:
    ch_records
    ch_fna
    ch_fai
    ch_code

    main:
    // CSV rows enter the channel before grouping. Numeric source fields are
    // excluded from computation/cache keys; all original cells reach merge.
    rows = ch_records.splitCsv(header: true).map { row -> tuple(row.GCFID, [row.Sample, row.CramFilePath, row.GCFID]) }
    fnas = ch_fna.splitCsv(header: true).map { row -> tuple(row.GCF_ID, file(row.GCF_path, checkIfExists: false)) }.unique()
    fais = ch_fai.splitCsv(header: true).map { row -> tuple(row.GCF_ID, file(row.FAI_path, checkIfExists: false)) }.unique()
    catalog = fnas.combine(fais, by: 0)
    chunks = rows.combine(catalog, by: 0).map { gcf, row, fna, fai ->
        def bucket = Math.floorMod(row.join('\u001f').hashCode(), params.index_chunk_buckets as int)
        tuple(bucket, [row[0], row[1], row[2], fna.toString(), fai.toString()])
    }.groupTuple().flatMap { bucket, grouped ->
        def ordered = grouped.unique().sort { left, right -> left.join('\u001f') <=> right.join('\u001f') }
        ordered.collate(params.index_chunk_size as int).withIndex().collect { batch ->
            def dependencies = batch[0].collectMany { row -> [row[1], row[3], row[4]] }.unique().sort().collect { source -> file(source, checkIfExists: true) }
            tuple([id: "bucket-${bucket}-part-${batch[1]}"], batch[0], dependencies)
        }
    }
    CRAM_METRICS_CHUNK(chunks, ch_code)
    // Only compact scalar SQLite shards are collected; no CRAM/coverage arrays.
    CRAM_METRICS_MERGE(ch_records, CRAM_METRICS_CHUNK.out.metrics.collect(), ch_fna, ch_fai, ch_code)

    emit:
    index = CRAM_METRICS_MERGE.out.index
    versions = CRAM_METRICS_CHUNK.out.versions.mix(CRAM_METRICS_MERGE.out.versions)
}
