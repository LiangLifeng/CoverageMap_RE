process BUILD_CRAM_INDEX {
    tag 'raw-cram-index'
    label 'cram_index'

    input:
    path records, stageAs: 'build_input.csv'
    path fna, stageAs: 'fna.csv'
    path fai, stageAs: 'fai_path.csv'
    path mate_audit, stageAs: 'mate_inventory/*'
    path builder, stageAs: 'index_code/build_cram_index.py'

    output:
    path 'reference_cram.sqlite', emit: index
    path 'versions.yml', emit: versions

    script:
    def auditFile = mate_audit instanceof List ? (mate_audit ? mate_audit.first() : null) : mate_audit
    def mateAudit = auditFile ? "--mate-audit '${auditFile}'" : ''
    """
    export PATH="\$PWD/index_code:\$PATH"
    build_cram_index.py --records "${records}" --fna "${fna}" --fai "${fai}" ${mateAudit} --output reference_cram.sqlite
    cat <<-VERSIONS > versions.yml
    "${task.process}":
        build_cram_index.py: \$(build_cram_index.py --version)
    VERSIONS
    """

    stub:
    """
    touch reference_cram.sqlite
    cat <<-VERSIONS > versions.yml
    "${task.process}":
        build_cram_index.py: "stub"
    VERSIONS
    """
}

