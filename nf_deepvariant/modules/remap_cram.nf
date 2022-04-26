//params.ref_dir = "/lustre/scratch118/humgen/resources/ref/Homo_sapiens/HS38DH/"

process remap_cram {
    tag "$sample"
    publishDir path: "${params.outdir}/remap_cram/",
	mode: "${params.copy_mode}",
	overwrite: "true"

    input:
    tuple val(sample), path(cram_file), path(cram_file_index)

    output:
    tuple val(sample), path("${sample}.merged.bam"), emit: remapped_sample_bam

    script:
"""
singularity run -B ${params.ref_dir}:${params.ref_dir} -B \$PWD:/home --pwd "/home" --no-home /software/hgi/containers/oqfe_remap.sif -1 /home/${cram_file} --sample ${sample} --cram-reference-fasta ${params.ref_dir}/${params.ref_file} -c -j 4
"""
}
