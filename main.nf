#!/usr/bin/env nextflow
include { COVERAGE_QUERY } from './workflows/coverage_query'

workflow {
    if (!params.manifest || !params.index || params.alpha == null)
        error 'Specify --manifest, --index and --alpha'
    if ((params.window_size == null) == (params.bin_count == null))
        error 'Specify exactly one of --window_size or --bin_count'
    if (!(params.alpha.toDouble() > 0 && params.alpha.toDouble() < 1))
        error '--alpha must be between 0 and 1'
    def methods = params.jaccard_methods.toString().split(',').collect { method -> method.trim() }
    if (!methods || methods.unique(false).size() != methods.size() || methods.any { method -> !(method in ['normal', 'empirical', 'minimum', 'bootstrap']) })
        error '--jaccard_methods must be distinct names from normal,empirical,minimum,bootstrap'
    if (!(params.reference_draws.toString() ==~ /[0-9]+/) || params.reference_draws.toLong() < 2)
        error '--reference_draws must be an integer at least 2'
    if (params.bootstrap_replicates.toInteger() < 10 || params.bootstrap_replicates.toInteger() > 10000)
        error '--bootstrap_replicates must be within 10..10000'
    if (!(params.bootstrap_confidence.toDouble() > 0 && params.bootstrap_confidence.toDouble() < 1))
        error '--bootstrap_confidence must be between 0 and 1'
    if (!(params.diagnostics.toString().toLowerCase() in ['true', 'false']))
        error '--diagnostics must be true or false'
    if (!(params.jaccard_coverage_threshold.toDouble() >= 0 && params.jaccard_coverage_threshold.toDouble() <= 1))
        error '--jaccard_coverage_threshold must be within [0,1]'
    if (params.jaccard_coverage_thresholds && !(params.jaccard_coverage_thresholds == 'default' || params.jaccard_coverage_thresholds ==~ /[0-9.eE+,-]+/))
        error '--jaccard_coverage_thresholds must be default or comma-separated fractions'
    if (params.query_cpus.toInteger() < 2)
        error '--query_cpus must be at least 2'
    if (params.coverage_unit && !(params.coverage_unit in ['fraction', 'percent']))
        error '--coverage_unit must be fraction or percent'

    if (params.sampling_unit && !(params.sampling_unit in ['read', 'fragment']))
        error '--sampling_unit must be read or fragment'
    if (params.coverage_alignment_policy && !(params.coverage_alignment_policy in ['primary_only', 'primary_plus_supplementary']))
        error 'Unsupported --coverage_alignment_policy'

    // CramPath is the legacy column name; both BAM and CRAM are decoded by content.
    // Unindexed BAMs from sample_cram_modes.py are streamed sequentially.
    log.info 'Query inputs: CramPath accepts BAM or CRAM; no BAI/CRAI is required'
    log.info 'V2 windows: GCF-concat; verified legacy native geometry is rebinned read-only'
    log.info "Query source: ${projectDir}; reference_draws=${params.reference_draws}; bootstrap_replicates=${params.bootstrap_replicates}"

    COVERAGE_QUERY(
        channel.value(file(params.manifest, checkIfExists: true)),
        channel.value(file(params.index, checkIfExists: true)),
        channel.value(file("${projectDir}/bin/*.py", checkIfExists: true)),
        channel.value(file("${projectDir}/Tools.py", checkIfExists: true)),
        channel.value(file("${projectDir}/bin/diagnostics/query_similarity_plot.py", checkIfExists: true))
    )
}
