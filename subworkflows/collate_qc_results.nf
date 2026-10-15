#!/usr/bin/env nextflow

include { CREATE_READQC_ONYX_JSON } from '../modules/create_qc_json'
include { ONYX_WRITE } from '../modules/onyx_write'
include { S3_UPLOAD } from '../modules/s3_upload'
include { ONYX_UPDATE } from '../modules/onyx_update'
include { ONYX_PUBLISH } from '../modules/onyx_publish'


workflow COLLATE_QC_RESULTS {
    take:
    ch_readQC // channel [val(meta), path(txt), path(csv)
    server // val: name of server running pipeline on
    bucket // val: name of bucket to upload results to
    context // val: string of context from the samplsheet (orange box version)

    main:
    // this section needs to be amended
    ch_readQC
        .branch { meta, readQC_complete, result_csv ->
            pass: readQC_complete.text.contains("True")
            fail: readQC_complete.text.contains("False")
    }
    .set { readQC_status }
    ch_readQC
        .map { meta, readQC_complete, result_csv -> tuple(meta, readQC_complete, result_csv) }
        .set { ch_readQC_complete }

    CREATE_READQC_ONYX_JSON(ch_readQC,  server, context)



    ONYX_WRITE(CREATE_READQC_ONYX_JSON.out.stats_summary, server)
    ch_s3_inputs = ch_readQC
            .join(ONYX_WRITE.out.analysis_id)
            .filter { meta, readQC_complete, result_csv, analysis_id -> readQC_complete.contains("True") }
            .map { meta, readQC_complete, result_csv, analysis_id
                -> tuple(meta, result_csv, analysis_id)
            }

    S3_UPLOAD(ch_s3_inputs, server, bucket)
    ch_update_inputs = ONYX_WRITE.out.analysis_id
        .join(S3_UPLOAD.out.s3_locations)
    ONYX_UPDATE(ch_update_inputs, server)
    ch_publish = channel.of()
        .concat(ONYX_WRITE.out.analysis_id, ONYX_UPDATE.out.analysis_id.ifEmpty([[], []]))
        .collect()
        .map { meta_write, id_write, meta_update, id_update ->
            tuple(meta_write, id_write)
         }
    ONYX_PUBLISH(server, ch_publish)

    emit:
    
    // this section needs to be amended
    readQC_summary = CREATE_READQC_ONYX_JSON.out.stats_summary // channel: [val(meta), path(readQC_json)]
    readQC_analysis_id = ONYX_PUBLISH.out.analysis_id

}
