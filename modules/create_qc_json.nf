#!/usr/bin/env nextflow

process CREATE_READQC_ONYX_JSON {
    /*
       This process creates a json file used to populate an onyx analysis table.

        Inputs:
            - meta: Sample metadata
            - readQC_status: Status of read QC run
            - csv: Paths to unclassified read QC output table
            - server: Server working on - one of mscape or synthscape
            - context: str of orange box version

        Outputs:
            - stats_summary - tuple of sample metadata and serotyping results file in onyx analysis table format
            - logs: Path to log file for the process

    */
    container 'ghcr.io/ukhsa-collaboration/gpha-mscape-onyx-analysis-helper:0.6.1'
    cpus 2
    memory '2GB'
    tag "${meta}"
    publishDir "${params.outdir}/${meta}/onyx", mode: params.publish_dir_mode

    input:
    // update to unclassified QC 
    tuple val(meta), val(readQC_status), path(csv)
    val server
    val context

    output:
    tuple val(meta), path("${meta}.unclassified_readQC.analysis_fields.json"), emit: stats_summary
    path "${meta}.unclassified.readQC.analysis_fields.log.txt", emit: logs

    script:
    """
    create_onyx_json.py \\
    --store-onyx \\
    -i ${meta} \\
    -s $server \\
    -v 1.0 \\
    -c $context \\
    -o . \\
    -p "${workflow.manifest.name}, ${workflow.manifest.version}, ${workflow.manifest.homePage}" \\
    -f $csv \\
    -r $readQC_status \\
    """
    }
