include { QUERY_PLAN } from '../modules/local/query_plan/main'
include { QUERY_SCORE } from '../modules/local/query_score/main'
include { QUERY_COLLECT } from '../modules/local/query_collect/main'
include { QUERY_DIAGNOSTICS } from '../modules/local/query_diagnostics/main'

workflow COVERAGE_QUERY {
    take:
    ch_manifest  // value: original query TSV
    ch_index     // value: native reference SQLite
    ch_code      // value: Python entry points and local dependencies
    ch_tools     // value: Tools.py
    ch_plotter   // value: independently staged rendering code

    main:
    QUERY_PLAN(ch_manifest, ch_index, ch_code, ch_tools)
    ch_query_meta = QUERY_PLAN.out.tasks
        .splitCsv(header: true, sep: '\t')
        .map { row -> [id: row.query_id, input_stamp: row.input_stamp] }

    QUERY_SCORE(ch_query_meta, ch_manifest, ch_index, ch_code, ch_tools)
    // Canonicalize by query ID before hashing the collector's file list.
    // Keep the empty list: an all-failed cohort still needs its error ledger.
    ch_receipts = QUERY_SCORE.out.receipts
        .toSortedList { left, right -> left[0] <=> right[0] }
        .map { rows -> rows.collect { row -> row[1] } }
    QUERY_COLLECT(QUERY_PLAN.out.roster, ch_receipts, ch_code, ch_tools)
    if (params.diagnostics.toString().toBoolean()) {
        ch_nulls = QUERY_SCORE.out.details
            .map { meta, result -> tuple(meta.id, result.resolve("queries/${meta.id}")) }
            .filter { id, query_dir -> query_dir.resolve("null.tsv.gz").exists() }
            .toSortedList { left, right -> left[0] <=> right[0] }
            .map { rows -> tuple(rows.collect { row -> row[0] }, rows.collect { row -> row[1] }) }
        QUERY_DIAGNOSTICS(QUERY_COLLECT.out.results, ch_nulls, ch_plotter)
    }

    emit:
    results = QUERY_COLLECT.out.results
    details = QUERY_SCORE.out.details
    versions = QUERY_PLAN.out.versions.mix(QUERY_SCORE.out.versions, QUERY_COLLECT.out.versions)
}
