process generate_fastq_stats {
    label "process_tiny"
    
    container "community.wave.seqera.io/library/pyfastx_pip_pyfunctional_pytest:10dc317fc97efcc2"

    publishDir "fastq_read_stats_results/", mode: 'copy'

    input:
    tuple val(climb_id), path(fastqgz)

    // include climb_id with the emitted file so downstream can join on it
    output:
    path("*.stats"), emit: readQC_csv
    path("readQC_complete.txt"), optional: true, emit: readQC_complete
    val(climb_id), emit: climb_id
    
    script:
    """
    set +e
    fastq_read_stats.py ${fastqgz} > ${fastqgz}.stats
    exit_status=\$?
    if [ \$exit_status == 0 ] ; then
        result="True"
        echo \$result > readQC_complete.txt
        exit \$exit_status
    elif grep "ERROR:" fastq_stdout.txt ; then
        result="False"
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

    // include climb_id so both stats files can be matched to the sample
    output:
    path("*.stats"), emit: readQC_csv
    path("readQC_complete.txt"), optional: true, emit: readQC_complete
    val(climb_id), emit: climb_id
    
    script:
    """
    set +e
    fastq_read_stats.py ${fastqgz1} > ${fastqgz1}.stats &>> fastq_stdout.txt
    fastq_read_stats.py ${fastqgz2} > ${fastqgz2}.stats &>> fastq_stdout.txt
    exit_status=\$?
    if [ \$exit_status == 0 ] ; then
        result="True"
        echo \$result > readQC_complete.txt
        exit \$exit_status
    elif grep "ERROR:" fastq_stdout.txt ; then
        result="False"
    else
        exit \$exit_status
    fi
    """
}

workflow generate_stats{
    take:
    single_end_ch
    
    main:
    stats = single_end_ch | generate_fastq_stats
    
    emit:
    readQC_complete = stats.readQC_complete
    readQC_csv      = stats.readQC_csv
    climb_id = stats.climb_id
}

workflow generate_stats_paired{
    take:
    paired_end_ch
    
    main:
    stats = paired_end_ch | generate_fastq_stats_paired
    
    emit:
    readQC_complete = stats.readQC_complete
    readQC_csv      = stats.readQC_csv
    climb_id = stats.climb_id
} 
