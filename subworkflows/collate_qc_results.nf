#!/usr/bin/env nextflow

include { CREATE_READQC_ONYX_JSON } from '../modules/create_qc_json'

/*******************************************************************************
 * COLLATE_QC_RESULTS
 *
 * Expects ch_readQC to be a channel of tuples: [climb_id, readQC_complete (file), result_csv (file or list)]
 * Branching now uses readQC_complete.text safely.
 *******************************************************************************/


workflow COLLATE_QC_RESULTS {
    take:
    collate_qc_results_ch   // channel of tuples: (climb_id, readQC_complete_file, stats_file(s))
    server      // val: name of server running pipeline on
    context     // val: string context from the sample sheet

    main:

    // Unwrap the nested tuple shape:
    // (climb_id, [readQC_complete, [stats_file_1, stats_file_2]])
    normalized_qc = collate_qc_results_ch
        .map { climb_id, qc_bundle ->
            def readQC_complete = qc_bundle[0]
            def stats_files = qc_bundle[1]

            // ensure stats_files is a list
            if (!(stats_files instanceof List)) {
                stats_files = [stats_files]
            }

            tuple(climb_id, readQC_complete, stats_files)
        }

    normalized_qc
        .branch { climb_id, readQC_complete, stats_files ->
            pass: readQC_complete.toString().contains("True")
            fail: readQC_complete.toString().contains("False")
        }
        .set { readQC_status }

    normalized_qc
        .map { climb_id, readQC_complete, stats_files ->
            def qc_text = readQC_complete instanceof Path || readQC_complete instanceof File
                ? readQC_complete.text
                : readQC_complete.toString()

            tuple(climb_id, qc_text, stats_files)
        }
        .set { normalized_readQC_complete }

    CREATE_READQC_ONYX_JSON(normalized_readQC_complete, server, context)
    
    emit:
    readQC_summary = CREATE_READQC_ONYX_JSON.out.stats_summary   // path(readQC_json)]
}
