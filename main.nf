nextflow.enable.dsl=2

params.runtag = 'deepvariant'
//params.index_crams = false
//params.run_deepvariant = true
//params.outdir = "${baseDir}/../../outputs"
//params.ref_dir = "/lustre/scratch118/humgen/resources/ref/Homo_sapiens/HS38DH/"
//params.bed_dir = "/lustre/scratch118/humgen/resources/exome/Homo_sapiens/"

include { imeta_study } from './nf_deepvariant/modules/imeta_study.nf'
include { imeta_study_lane } from './nf_deepvariant/modules/imeta_study_lane.nf'
include { iget_study_cram } from './nf_deepvariant/modules/iget_study_cram.nf'
include { sort_cram } from './nf_deepvariant/modules/sort_cram.nf' //params(run: true, outdir: params.outdir)
include { remap_cram } from './nf_deepvariant/modules/remap_cram.nf' //params(run: true, outdir: params.outdir)
include { markDuplicates } from './nf_deepvariant/modules/markDuplicates.nf' //params(run: true, outdir: params.outdir)
include { coord_sort_cram } from './nf_deepvariant/modules/coord_sort_cram.nf' //params(run: true, outdir: params.outdir)
include { bam_to_cram } from './nf_deepvariant/modules/bam_to_cram.nf' //params(run: true, outdir: params.outdir)
include { deepvariant } from './nf_deepvariant/modules/deepvariant.nf' //params(run: true, outdir: params.outdir)
include { gatk_haplotypecaller } from './nf_deepvariant/modules/gatk_haplotypecaller' //params(run: true, outdir: params.outdir)

workflow {

    log.info "${params}"

    if (params.run_deepvariant) {
	
//    if (params.tsv_file) {
//	ch_cram_file = Channel
//		.fromPath(params.tsv_file)
//		.splitCsv(header: true, sep: '\t')
//                .map{row->tuple(row.study_id, row.sample, row.object)}
//                .filter { it[2] =~ /.cram$/ } // Need to check for bam too?
//		.take(params.samples_to_process)
		//.view()
//    }

     main:
        
        log.info "${params.ref_dir}"

    if (params.run_mode == "study_id") {
        if (params.study_id_mode.input_study_lanes) {
            imeta_study_lane( [params.study_id_mode.input_studies, params.study_id_mode.input_study_lanes] )
            samples_irods_tsv = imeta_study_lane.out.irods_samples_tsv
            work_dir_to_remove = imeta_study_lane.out.work_dir_to_remove
        }
        else if (params.study_id_mode.exclude_samples_file) {
            imeta_study( [params.study_id_mode.input_studies, params.study_id_mode.exclude_samples_file])
            samples_irods_tsv = imeta_study.out.irods_samples_tsv
            work_dir_to_remove = imeta_study.out.work_dir_to_remove
        }
	else{
            imeta_study(Channel.from(params.study_id_mode.input_studies))
	    samples_irods_tsv = imeta_study.out.irods_samples_tsv
	    work_dir_to_remove = imeta_study.out.work_dir_to_remove
        } 
        iget_study_cram(
            samples_irods_tsv
                .map{study_id, samples_tsv -> samples_tsv}
                .splitCsv(header: true, sep: '\t')
                .map{row->tuple(row.study_id, row.sample, row.object)}
                .filter { it[2] =~ /.cram$/ } // Need to check for bam too?
                .take(params.samples_to_process)
                .dump()
                .unique())
    }
    else if (params.tsv_file) {
        cram_file = Channel
                .fromPath(params.tsv_file)
                .splitCsv(header: true, sep: '\t')
                .map{row->tuple(row.study_id, row.sample, row.object, row.object_index)}
                .filter { it[2] =~ /.cram$/ } // Need to check for bam too?
                .take(params.samples_to_process)
                .dump()
                //.view()
    }

        if (params.remap) {
           remap_cram(cram_file) 
           markDuplicates(remap_cram.out.remapped_sample_bam)
           coord_sort_cram(markDuplicates.out.markdup_sample_cram)
           bam_to_cram(coord_sort_cram.out.markdup_sample_cram_crai)
           deepvariant(coord_sort_cram.out.markdup_sample_cram_crai)
           gatk_haplotypecaller(coord_sort_cram.out.markdup_sample_cram_crai)
        }
        else if (params.run_sort_cram){
                 if (params.tsv_file) {
			sort_cram(cram_file)
                 }
                 else {
			sort_cram(iget_study_cram.out.study_sample_cram_crai)
                 }
		markDuplicates(sort_cram.out.sorted_sample_cram)
		coord_sort_cram(markDuplicates.out.markdup_sample_cram)
		bam_to_cram(coord_sort_cram.out.markdup_sample_cram_crai)
		deepvariant(coord_sort_cram.out.markdup_sample_cram_crai)
		gatk_haplotypecaller(coord_sort_cram.out.markdup_sample_cram_crai)
                 }
              else {
                 if (params.tsv_file) {
                     deepvariant(cram_file)
                     gatk_haplotypecaller(cram_file)}
                 else {
                     deepvariant(iget_study_cram.out.study_sample_cram_crai)
                     gatk_haplotypecaller(iget_study_cram.out.study_sample_cram_crai)}
              }  
     emit:
        my_data = deepvariant.out
        
   }
}

workflow.onError {
    log.info "Pipeline execution stopped with the following message: ${workflow.errorMessage}" }
workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Command line: $workflow.commandLine"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"

    if (params.on_complete_uncache_irods_search) {
        log.info "You have selected \"on_complete_uncache_irods_search = true\"; will therefore attempt to remove Irods work dirs to forcefully uncache them even if successful."
        if (! file("${params.outdir}/irods_work_dirs_to_remove.csv").isEmpty()) {
            log.info "file ${params.outdir}/irods_work_dirs_to_remove.csv exists and not empty ..."
            file("${params.outdir}/irods_work_dirs_to_remove.csv")
                .eachLine {  work_dir ->
                if (file(work_dir).isDirectory()) {
                    log.info "removing work dir $work_dir ..."
                    file(work_dir).deleteDir()
                } } } }

    if (params.on_complete_remove_workdir_failed_tasks) {
        log.info "You have selected \"on_complete_remove_workdir_failed_tasks = true\"; will therefore remove work dirs of all tasks that failed (.exitcode file not 0)."
        // work dir and other paths are hardcoded here ... :
        def proc = "bash ./nextflow_ci/bin/del_work_dirs_failed.sh ${workDir}".execute()
        def b = new StringBuffer()
        proc.consumeProcessErrorStream(b)
        log.info proc.text
        log.info b.toString() }
}
