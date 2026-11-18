include { REFERENCE_PLAN } from '../modules/local/reference_plan/main'
include { REFERENCE_CHUNK } from '../modules/local/reference_chunk/main'
include { REFERENCE_COLLECT } from '../modules/local/reference_collect/main'
include { REFERENCE_INDEX } from '../modules/local/reference_index/main'

workflow REFERENCE_BUILD {
    take:
    ch_records
    ch_references
    ch_output_root
    ch_code
    ch_tools

    main:
    REFERENCE_PLAN(ch_records, ch_references, ch_code, ch_tools)
    ch_chunks = REFERENCE_PLAN.out.chunks
        .flatten()
        .map { chunk -> tuple([id: chunk.baseName], chunk) }
    REFERENCE_CHUNK(ch_chunks, ch_output_root, ch_code, ch_tools)
    ch_receipts = REFERENCE_CHUNK.out.receipts
        .toSortedList { left, right -> left[0].id <=> right[0].id }
        .map { rows -> rows.collect { row -> row[1] } }
    REFERENCE_COLLECT(REFERENCE_PLAN.out.plan, ch_receipts, ch_output_root, ch_code, ch_tools)
    REFERENCE_INDEX(REFERENCE_COLLECT.out.receipt, ch_output_root, ch_code, ch_tools)

    emit:
    receipt = REFERENCE_INDEX.out.receipt
    chunks = REFERENCE_CHUNK.out.receipts
    versions = REFERENCE_PLAN.out.versions.mix(REFERENCE_CHUNK.out.versions, REFERENCE_COLLECT.out.versions, REFERENCE_INDEX.out.versions)
}
