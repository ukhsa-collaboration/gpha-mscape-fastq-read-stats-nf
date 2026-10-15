#!/usr/bin/env nextflow

process S3_UPLOAD {
    /*
       This process uploads the files produced by pneumokity to s3.

        Inputs:
            - meta: Sample metadata
            - *_csv: Paths to pneumokity files for upload
            - server: Server working on - one of mscape or synthscape
            - bucket: Bucket to upload data to
            - analysis_id: Analysis ID to use in s3 key

        Outputs:
            - s3_locations: Tuple of sample metadata and path to s3 locations file
            - logs: Path to log file for process
    */
    container 'ghcr.io/ukhsa-collaboration/gpha-mscape-onyx-analysis-helper:0.6.5'
    cpus 1
    memory '1GB'
    tag "${meta.id}"
    stageInMode 'copy'
    publishDir "${params.outdir}/${meta.id}/onyx", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(alldata_csv), path(quality_csv), path(result_csv), path(analysis_id)
    val server
    val bucket

    output:
    tuple val(meta), path("${meta.id}.onyx_analysis.s3_upload.analysis_fields.json"), emit: s3_locations
    path "${meta.id}.onyx_analysis.s3_upload.log.txt", emit: logs

    script:
    """
    onyx_analysis s3_upload \\
    -i ${meta.id} \\
    -s $server \\
    -b $bucket \\
    -o . \\
    -a $analysis_id \\
    -f "${alldata_csv},${quality_csv},${result_csv}" \\
    --prod
    """
    }
