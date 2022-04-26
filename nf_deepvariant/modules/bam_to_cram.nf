//params.ref_dir = "/lustre/scratch118/humgen/resources/ref/Homo_sapiens/HS38DH/"

process bam_to_cram {
    tag "$sample"
    publishDir path: "${params.outdir}/bam_to_cram/",
	mode: "${params.copy_mode}",
	overwrite: "true"
    publishDir "${params.cram_output_dir}", mode: "${params.copy_mode}"

    when:
    params.run_coord_sort_cram
     
    input:
    tuple val(sample), path(cram_file_sorted_dups_coord), path(cram_file_sorted_dups_coord_index)

    output:
    tuple val(sample), path("${sample}.sorted.dups.coord.cram"), path("${sample}.sorted.dups.coord.cram.crai"), emit: processed_sample_cram_crai

    script:
""" 
samtools view ${cram_file_sorted_dups_coord} \\
  -o ${sample}.sorted.dups.coord.cram -O CRAM \\
  -T ${params.ref_dir}/${params.ref_file}

samtools index ${sample}.sorted.dups.coord.cram
"""
}

