#!/usr/bin/env nextflow

process ONYX_UPDATE {
    /*
       This process updates an unpublished onyx analysis table with s3 file paths.

        Inputs:
            - server: Server working on - one of mscape or synthscape
            - meta: Sample metadata
            - analysis_id: File containing analysis ID
            - update_json: Onyx analysis json with s3 locations of result files

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
    tuple val(meta), path(analysis_id), path(update_json)
    val server

    output:
    tuple val(meta), path("${meta.id}.onyx_analysis.update.analysis_id.txt"), emit: analysis_id
    path "${meta.id}.onyx_analysis.update.log.txt", emit: logs

    script:
    """
    onyx_analysis update \\
    -i ${meta.id} \\
    -s $server \\
    -o . \\
    -a $analysis_id \\
    -j $update_json \\
    --prod
    """
    }
