process generate_fastq_stats {
    label "process_tiny"
    
    container "community.wave.seqera.io/library/pyfastx_pip_pyfunctional_pytest:10dc317fc97efcc2"

    publishDir "fastq_read_stats_results/", mode: 'copy'

    input:
    tuple val(climb_id), path(fastqgz)

    output:
    path "*.stats"
    
    script:
    """
    fastq_read_stats.py ${fastqgz} > ${fastqgz}.stats
    """
}

process generate_fastq_stats_paired {
    label "process_tiny"
    
    container "community.wave.seqera.io/library/pyfastx_pip_pyfunctional_pytest:10dc317fc97efcc2"

    publishDir "fastq_read_stats_results/", mode: 'copy'

    input:
    tuple val(climb_id), path(fastqgz1), path(fastqgz2)

    output:
    path "*.stats"
    
    script:
    """
    fastq_read_stats.py ${fastqgz1} > ${fastqgz1}.stats
    fastq_read_stats.py ${fastqgz2} > ${fastqgz2}.stats
    """
}

workflow generate_stats{
    take:
    single_end_ch
    
    main:
    single_end_ch | generate_fastq_stats
}

workflow generate_stats_paired{
    take:
    paired_end_ch
    
    main:
    paired_end_ch | generate_fastq_stats_paired
} 