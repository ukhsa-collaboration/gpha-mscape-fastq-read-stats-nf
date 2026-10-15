#!/usr/bin/env nextflow

process ONYX_PUBLISH {
    /*
       This process publishes a complete onyx analysis table.

        Inputs:
            - server: Server working on - one of mscape or synthscape
            - meta: Sample metadata
            - analysis_id: File containing analysis ID

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
    val server
    tuple val(meta), path(analysis_id)

    output:
    tuple val(meta), path("${meta.id}.onyx_analysis.publish.analysis_id.txt"), emit: analysis_id
    path "${meta.id}.onyx_analysis.publish.log.txt", emit: logs

    script:
    """
    onyx_analysis publish \\
    -i ${meta.id} \\
    -s $server \\
    -o . \\
    -a $analysis_id \\
    --prod
    """
    }
