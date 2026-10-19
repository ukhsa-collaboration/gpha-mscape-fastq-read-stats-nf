#!/usr/bin/env nextflow

include { CREATE_READQC_ONYX_JSON } from '../modules/create_qc_json'
include { ONYX_WRITE } from '../modules/onyx_write'
include { S3_UPLOAD } from '../modules/s3_upload'
include { ONYX_UPDATE } from '../modules/onyx_update'
include { ONYX_PUBLISH } from '../modules/onyx_publish'

/*******************************************************************************
 * COLLATE_QC_RESULTS
 *
 * Expects ch_readQC to be a channel of tuples: [climb_id, readQC_complete (file), result_csv (file or list)]
 * Branching now uses readQC_complete.text safely.
 *******************************************************************************/


workflow COLLATE_QC_RESULTS {
    take:
    ch_readQC   // channel of tuples: (climb_id, readQC_complete_file, stats_file(s))
    server      // val: name of server running pipeline on
    bucket      // val: name of bucket to upload results to
    context     // val: string context from the sample sheet

    main:

    // Normalize: ensure result_csv is a single file (if *.stats matched multiple, pick first)
    normalized_readQC = ch_readQC
        .map { climb_id, readQC_complete, result_csv ->
            def stats_file = (result_csv instanceof List) ? result_csv[0] : result_csv
            tuple(climb_id, readQC_complete, stats_file)
        }

    // Branch into pass/fail based on the contents of readQC_complete.txt
    normalized_readQC
        .branch { climb_id, readQC_complete, result_csv ->
            pass: readQC_complete.text.contains("True")
            fail: readQC_complete.text.contains("False")
        }
        .set { readQC_status }
    normalized_readQC
        .map { meta, readQC_complete, result_csv -> tuple(meta, readQC_complete.text, result_csv) }
        .set { normalized_readQC_complete }
    // Create onyx JSON summary from the readQC tuples
    CREATE_READQC_ONYX_JSON(normalized_readQC_complete, server, context)

    // Write to onyx (example usage)
    ONYX_WRITE(CREATE_READQC_ONYX_JSON.out.stats_summary, server)

    // Prepare S3 upload inputs: join the readQC tuples with the analysis id from ONYX_WRITE
    ch_s3_inputs = normalized_readQC
        .join(ONYX_WRITE.out.analysis_id)
         //Filter to only successful runs before uploading
        .filter { climb_id, readQC_complete, result_csv, analysis_id -> readQC_complete.text.contains("True") }
        .map { climb_id, readQC_complete, result_csv, analysis_id ->
            tuple(climb_id, result_csv, analysis_id)
        }

    // Upload results to S3
    S3_UPLOAD(ch_s3_inputs, server, bucket)

    // Update onyx with S3 locations
    ch_update_inputs = ONYX_WRITE.out.analysis_id
        .join(S3_UPLOAD.out.s3_locations)
    ONYX_UPDATE(ch_update_inputs, server)

    // Publish (combine write/update outputs for publish)
    ch_publish = channel.of()
        .concat(ONYX_WRITE.out.analysis_id, ONYX_UPDATE.out.analysis_id.ifEmpty([[], []]))
        .collect()
        .map { meta_write, id_write, meta_update, id_update ->
            tuple(meta_write, id_write)
        }
    ONYX_PUBLISH(server, ch_publish)


    emit:
    readQC_summary = CREATE_READQC_ONYX_JSON.out.stats_summary   // channel [val(meta), path(readQC_json)]
    readQC_analysis_id = ONYX_PUBLISH.out.analysis_id
}