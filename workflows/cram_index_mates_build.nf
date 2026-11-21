include { CRAM_METRICS_CHUNK } from '../modules/local/cram_metrics_chunk/main'
include { CRAM_METRICS_RECOVER } from '../modules/local/cram_metrics_recover/main'
include { CRAM_METRICS_MATE_MERGE } from '../modules/local/cram_metrics_mate_merge/main'

workflow CRAM_INDEX_BUILD {
    take:
    ch_records
    ch_fna
    ch_fai
    ch_legacy_code
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
    // Keep the original process name, tuple and staged code for native -resume.
    def known = params.mate_recovery_buckets ? params.mate_recovery_buckets.split(',').toList() : []
    candidates = chunks.filter { meta, chunk_rows, deps -> !known.contains(meta.id.toString()) }
    known_failed = chunks.filter { meta, chunk_rows, deps -> known.contains(meta.id.toString()) }
    CRAM_METRICS_CHUNK(candidates, ch_legacy_code)
    successful = CRAM_METRICS_CHUNK.out.metrics.map { shard -> tuple(shard.baseName, shard) }
    failed = candidates.map { meta, chunk_rows, deps -> tuple(meta.id.toString(), meta, chunk_rows, deps) }
        .join(successful, remainder: true)
        .filter { row -> row.size() == 4 || row[4] == null }
        .map { row -> tuple(row[1], row[2], row[3]) }
    CRAM_METRICS_RECOVER(known_failed.mix(failed), ch_code)
    shards = CRAM_METRICS_CHUNK.out.metrics.mix(CRAM_METRICS_RECOVER.out.metrics.map { id, shard -> shard })
    CRAM_METRICS_MATE_MERGE(ch_records, shards.collect().ifEmpty { error 'No completed CRAM metric shards' }, ch_fna, ch_fai, ch_code)

    emit:
    index = CRAM_METRICS_MATE_MERGE.out.index
    exclusions = CRAM_METRICS_MATE_MERGE.out.exclusions
    versions = CRAM_METRICS_CHUNK.out.versions.mix(CRAM_METRICS_RECOVER.out.versions, CRAM_METRICS_MATE_MERGE.out.versions)
}
