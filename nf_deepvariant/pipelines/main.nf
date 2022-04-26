nextflow.preview.dsl=2

include { sort_cram } from '../modules/sort_cram.nf'
include { remap_cram } from '../modules/remap_cram.nf'
include { markDuplicates } from '../modules/markDuplicates.nf'
include { coord_sort_cram } from '../modules/coord_sort_cram.nf'
include { bam_to_cram } from '../modules/bam_to_cram.nf'
include { deepvariant } from '../modules/deepvariant.nf'
include { gatk_haplotypecaller } from '../modules/gatk_haplotypecaller.nf'

workflow {
    main:
    
    log.info "params: ${params}"
    log.info "params.tsv_file: $params.tsv_file"
    
    channel_inputs_bams = Channel
        .fromPath(params.tsv_file)
        .splitCsv(header: true, sep: '\t')
        .map{row->tuple(row.sample, row.object, row.object_index)}
        .take(params.samples_to_process)
    
    // channel_inputs_bams.view()
    
    if (params.run_mode == "remap_inputs") {
        remap_cram(channel_inputs_bams) 
        markDuplicates(remap_cram.out.remapped_sample_bam)
        coord_sort_cram(markDuplicates.out.markdup_sample_cram)
        bam_to_cram(coord_sort_cram.out.markdup_sample_cram_crai)
        deepvariant(coord_sort_cram.out.markdup_sample_cram_crai)
        gatk_haplotypecaller(coord_sort_cram.out.markdup_sample_cram_crai)
    }
    else if (params.run_mode == "sort_inputs") {
        sort_cram(channel_inputs_bams)
        markDuplicates(sort_cram.out.sorted_sample_cram)
        coord_sort_cram(markDuplicates.out.markdup_sample_cram)
        bam_to_cram(coord_sort_cram.out.markdup_sample_cram_crai)
        deepvariant(coord_sort_cram.out.markdup_sample_cram_crai)
        gatk_haplotypecaller(coord_sort_cram.out.markdup_sample_cram_crai)
    }
    else if (params.run_mode == "sort_inputs") {
        deepvariant(channel_inputs_bams)
        gatk_haplotypecaller(channel_inputs_bams)
    }
    else {
	log.info "params.run_mode must be either \"remap_inputs\" \"sort_inputs\" or \"no_remap_no_sort\""	
    }
    
    emit:
    deepvariant_out = deepvariant.out
}

workflow.onError {
    log.info "Pipeline execution stopped with the following message: ${workflow.errorMessage}" }
workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Command line: $workflow.commandLine"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
