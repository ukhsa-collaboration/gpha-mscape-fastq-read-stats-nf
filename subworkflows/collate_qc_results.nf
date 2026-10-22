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

    collate_qc_results_ch
        .map { climb_id, readqc_file, stats_files ->
            def qc_text = readqc_file.text
            tuple(climb_id, qc_text, stats_files)
        }
        .set { normalized_readQC_complete }

    CREATE_READQC_ONYX_JSON(normalized_readQC_complete, server, context)
    
    emit:
    readQC_summary = CREATE_READQC_ONYX_JSON.out.stats_summary   // path(readQC_json)]
}
