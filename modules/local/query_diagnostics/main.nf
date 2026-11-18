process QUERY_DIAGNOSTICS {
    tag 'similarity distributions'
    label 'query_diagnostics'

    input:
    path collected
    // Query directories retain null scores, coverage arrays, window lengths and draw identities.
    tuple val(ids), path(nulls, stageAs: 'nulls/input??????/*')
    path plotter, stageAs: 'plot_bin/*'

    output:
    path 'diagnostics', emit: plots
    path 'versions.yml', emit: versions

    script:
    def runtime = params.diagnostics_bin_dir ? "${params.diagnostics_bin_dir}:" : ''
    """
    export PATH="\$PWD/plot_bin:${runtime}\$PATH"
    export MPLCONFIGDIR="\$PWD/matplotlib-config" MPLBACKEND=Agg
    printf '%s\\n' '${ids.join('\n')}' > query_ids.txt
    mkdir diagnostics
    rc=0
    query_similarity_plot.py --summary "${collected}/summary.tsv" --decision-methods "${collected}/decision_methods.tsv" --methods "${params.jaccard_methods}" --outdir diagnostics \\
        --staged-nulls nulls --staged-ids query_ids.txt --jobs ${task.cpus} \\
        2>diagnostics/plot_stderr.log || rc=\$?
    printf '%s\\n' "\$rc" > diagnostics/plot_exit_code.txt
    version=\$(query_similarity_plot.py --version)
    printf '"%s":\\n    query_similarity_plot.py: "%s"\\n' "${task.process}" "\$version" > versions.yml
    exit \$rc
    """

    stub:
    def runtime = params.diagnostics_bin_dir ? "${params.diagnostics_bin_dir}:" : ''
    """
    export PATH="\$PWD/plot_bin:${runtime}\$PATH"
    printf '%s\\n' '${ids.join('\n')}' > query_ids.txt
    mkdir diagnostics
    query_similarity_plot.py --summary "${collected}/summary.tsv" --decision-methods "${collected}/decision_methods.tsv" --methods "${params.jaccard_methods}" --outdir diagnostics \\
        --staged-nulls nulls --staged-ids query_ids.txt --stub
    printf '0\\n' > diagnostics/plot_exit_code.txt
    printf '"%s":\\n    query_similarity_plot.py: "stub"\\n' "${task.process}" > versions.yml
    """
}
