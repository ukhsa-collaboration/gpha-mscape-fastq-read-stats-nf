#!/usr/bin/env nextflow

include { GENERATE_SAMPLESHEET } from './modules/local/samplesheet'
include { generate_stats      } from './subworkflows/fastq_read_stats'
include { generate_stats_paired } from './subworkflows/fastq_read_stats'
// include { ONYX_UPLOAD          } from "./modules/local/onyx_upload"


workflow {
    // Handle either samplesheet or climb id
    if (params.samplesheet && params.unique_id){
        exit(1, "Please specify one of --unique_id or --samplesheet. Not both.")
    }
    else if(params.samplesheet){
        log.info "Samplesheet input: ${params.samplesheet}"
        samplesheet_ch = channel.fromPath(params.samplesheet)
    } else if (params.unique_id) {
        log.info "Unique ID input: ${params.unique_id}"
        sample_ch = Channel.of(tuple (params.unique_id, params.samplesheet_columns))
        samplesheet_ch = GENERATE_SAMPLESHEET(sample_ch).samplesheet
    }
    else{
        exit(1, "Please specify either --unique_id or --samplesheet")
    }
    // Split csv into channels
    samples = samplesheet_ch.splitCsv(header: true, quote: '\"')
        .map { row ->
            def climb_id = row.climb_id
            // def taxon_report_dir = row.taxon_reports
            def fastq1 = row.unclassified_reads_1
            def fastq2 = row.containsKey('unclassified_reads_2') ? row.unclassified_reads_2 : null
            return fastq2 ? tuple(climb_id, fastq1, fastq2) : tuple(climb_id, fastq1)
        }
        // Split channel by the number of reads
        .branch{ v -> 
            paired_end: v.size() == 3
            single_end: v.size() == 2
        }
        .set { ch_fastqs }
            
    generate_stats(ch_fastqs.single_end)
    
    generate_stats_paired(ch_fastqs.paired_end)
}