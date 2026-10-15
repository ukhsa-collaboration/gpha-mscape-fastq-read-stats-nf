#!/usr/bin/env nextflow

process ONYX_WRITE {
    /*
       This process adds pneumokity results to an onyx analysis table.

        Inputs:
            - meta: Sample metadata
            - result_json: Path to json containing analysis fields for onyx
            - server: Server working on - one of mscape or synthscape

        Outputs:
            - analysis_id: Tuple of sample metadata and path to file with onyx analysis ID
            - logs: Path to log file for process

    */
    container 'ghcr.io/ukhsa-collaboration/gpha-mscape-onyx-analysis-helper:0.6.5'
    cpus 1
    memory '1GB'
    tag "${meta.id}"
    publishDir "${params.outdir}/${meta.id}/onyx", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(result_json)
    val server

    output:
    tuple val(meta), path("${meta.id}.onyx_analysis.write.analysis_id.txt"), emit: analysis_id
    path "${meta.id}.onyx_analysis.write.log.txt", emit: logs

    script:
    """
    onyx_analysis write \\
    -i ${meta.id} \\
    -s $server \\
    -o . \\
    -j $result_json \\
    --prod
    """
    }
