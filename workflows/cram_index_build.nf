include { BUILD_CRAM_INDEX } from '../modules/local/build_cram_index/main'

workflow CRAM_INDEX_BUILD {
    take:
    ch_records  // value: BuildInput CSV
    ch_fna  // value: full FNA catalog CSV
    ch_fai  // value: full FAI catalog CSV
    ch_mate_audit // value: optional complete mate-label inventory
    ch_builder  // value: tracked standalone Python entry point

    main:
    BUILD_CRAM_INDEX(ch_records, ch_fna, ch_fai, ch_mate_audit, ch_builder)

    emit:
    index = BUILD_CRAM_INDEX.out.index
    versions = BUILD_CRAM_INDEX.out.versions
}

