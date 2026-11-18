include { RAW_QUERY_PLAN } from '../modules/local/raw_query_plan/main'
include { RAW_QUERY_SCORE } from '../modules/local/raw_query_score/main'
include { RAW_QUERY_COLLECT } from '../modules/local/raw_query_collect/main'
include { QUERY_DIAGNOSTICS as RAW_QUERY_DIAGNOSTICS } from '../modules/local/query_diagnostics/main'

workflow RAW_CRAM_QUERIES {
    take:
    ch_manifest  // value: original TSV
    ch_index     // value: raw-CRAM metadata SQLite
    ch_code      // value: tracked Python entry points and every local helper
    ch_fai_map   // value: optional FAI mapping file list
    ch_plotter   // value: legacy diagnostics renderer

    main:
    RAW_QUERY_PLAN(ch_manifest, ch_index, ch_code, ch_fai_map)
    ch_queries = RAW_QUERY_PLAN.out.tasks
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            tuple([id: row.query_id], row.payload,
                  row.inputs ? row.inputs.split('\\|').collect { name -> file(name) } : [])
        }
    RAW_QUERY_SCORE(ch_queries, ch_code)
    ch_results = RAW_QUERY_SCORE.out.results
        .toSortedList { left, right -> left[0].id <=> right[0].id }
        .map { rows -> rows.collect { row -> row[1] } }
    RAW_QUERY_COLLECT(RAW_QUERY_PLAN.out.roster, ch_results, ch_code)
    if (params.diagnostics.toString().toBoolean()) {
        ch_nulls = RAW_QUERY_SCORE.out.results
            .map { meta, directory -> tuple(meta.id, directory) }
            .filter { id, directory -> directory.resolve('null.tsv.gz').exists() }
            .toSortedList { left, right -> left[0] <=> right[0] }
            .map { rows -> tuple(rows.collect { row -> row[0] }, rows.collect { row -> row[1] }) }
        RAW_QUERY_DIAGNOSTICS(RAW_QUERY_COLLECT.out.summary, ch_nulls, ch_plotter)
    }

    emit:
    results = RAW_QUERY_SCORE.out.results
    summary = RAW_QUERY_COLLECT.out.summary
    versions = RAW_QUERY_PLAN.out.versions.mix(RAW_QUERY_SCORE.out.versions, RAW_QUERY_COLLECT.out.versions)
}

