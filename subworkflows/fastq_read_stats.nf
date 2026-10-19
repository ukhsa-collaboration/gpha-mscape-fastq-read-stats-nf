#!/usr/bin/env nextflow

// module includes (adjust relative paths as needed)
include { CREATE_READQC_ONYX_JSON } from '../modules/create_qc_json'
include { ONYX_WRITE } from '../modules/onyx_write'
include { S3_UPLOAD } from '../modules/s3_upload'
include { ONYX_UPDATE } from '../modules/onyx_update'
include { ONYX_PUBLISH } from '../modules/onyx_publish'


/*******************************************************************************
 * Processes: generate_fastq_stats (single & paired)
 *
 * These now emit a single tuple channel "readQC" containing:
 *  - val(climb_id)
 *  - path('readQC_complete.txt')   (file containing "True" or "False")
 *  - path('*.stats')               (one or more stats files)
 *******************************************************************************/

process generate_fastq_stats {
    label "process_tiny"
    container "community.wave.seqera.io/library/pyfastx_pip_pyfunctional_pytest:10dc317fc97efcc2"
    publishDir "fastq_read_stats_results/", mode: 'copy'

    input:
    tuple val(climb_id), path(fastqgz)

    // Emit a single tuple channel so downstream receives (id, completeFile, statsFile(s))
    output:
    tuple val(climb_id), path("readQC_complete.txt"), path("*.stats"), emit: readQC

    script:
    """
    set +e
    fastq_read_stats.py ${fastqgz} > ${fastqgz}.stats 2> fastq_stdout.txt
    exit_status=\$?
    if [ \$exit_status -eq 0 ] ; then
        echo "True" > readQC_complete.txt
        exit 0
    elif grep -q "ERROR:" fastq_stdout.txt ; then
        echo "False" > readQC_complete.txt
        exit 0
    else
        exit \$exit_status
    fi
    """
}

process generate_fastq_stats_paired {
    label "process_tiny"
    container "community.wave.seqera.io/library/pyfastx_pip_pyfunctional_pytest:10dc317fc97efcc2"
    publishDir "fastq_read_stats_results/", mode: 'copy'

    input:
    tuple val(climb_id), path(fastqgz1), path(fastqgz2)

    output:
    tuple val(climb_id), path("readQC_complete.txt"), path("*.stats"), emit: readQC

    script:
    """
    set +e
    fastq_read_stats.py ${fastqgz1} > ${fastqgz1}.stats &>> fastq_stdout.txt
    fastq_read_stats.py ${fastqgz2} > ${fastqgz2}.stats &>> fastq_stdout.txt
    exit_status=\$?
    if [ \$exit_status -eq 0 ] ; then
        echo "True" > readQC_complete.txt
        exit 0
    elif grep -q "ERROR:" fastq_stdout.txt ; then
        echo "False" > readQC_complete.txt
        exit 0
    else
        exit \$exit_status
    fi
    """
}


/*******************************************************************************
 * Workflows that wrap the processes and expose a single readQC channel
 *******************************************************************************/

workflow generate_stats {
    take:
    single_end_ch

    main:
    // stats is a workflow result object; .readQC is the emitted channel
    stats = single_end_ch | generate_fastq_stats

    emit:
    readQC = stats.readQC
}

workflow generate_stats_paired {
    take:
    paired_end_ch

    main:
    stats = paired_end_ch | generate_fastq_stats_paired

    emit:
    readQC = stats.readQC
}
