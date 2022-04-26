//params.ref_dir = "/lustre/scratch118/humgen/resources/ref/Homo_sapiens/HS38DH/"
//params.bed_dir = "/lustre/scratch118/humgen/resources/exome/Homo_sapiens/"

process gatk_haplotypecaller {
    tag "$sample"
    publishDir path: "${params.outdir}/gatk_haplotypecaller/",
	mode: "${params.copy_mode}",
	overwrite: "true",
	pattern: "*gatk.g.vcf.gz*"
    
    when:
    params.run_haplotypecaller
     
    input:
    tuple val(sample), path(cram_file_sorted_dups_coord), path(cram_file_sorted_dups_coord_index)

    output:
    tuple path("${cram_file_sorted_dups_coord}.gatk.g.vcf.gz"), path("${cram_file_sorted_dups_coord}.gatk.g.vcf.gz.tbi")

    script:
""" 
/gatk/gatk --java-options \"-Xms6g -Xmx6g -XX:+UseSerialGC\" \\
  HaplotypeCaller \\
  -I ${cram_file_sorted_dups_coord} \\
  -O ${cram_file_sorted_dups_coord}.gatk.g.vcf.gz \\
  -R ${params.ref_dir}/${params.ref_file} \\
  -L ${params.bed_dir}/${params.bed_file_deepvariant} \\
  -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \\
  -ERC GVCF \\
  -G StandardAnnotation \\
  -G StandardHCAnnotation \\
  -G AS_StandardAnnotation
"""
}

