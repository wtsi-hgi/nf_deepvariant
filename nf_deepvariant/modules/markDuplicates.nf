//params.ref_dir = "/lustre/scratch118/humgen/resources/ref/Homo_sapiens/HS38DH/"

process markDuplicates {
    tag "$sample"
    publishDir path: "${params.outdir}/markDuplicates/",
	mode: "${params.copy_mode}",
	overwrite: "true"

    when:
    params.run_markDuplicates
     
    input:
    tuple val(sample), path(cram_file_sorted)

    output:
    tuple val(sample), path("${cram_file_sorted}.dups"), emit: markdup_sample_cram
    path "${cram_file_sorted}.dups"

    script:
""" 
/gatk/gatk --java-options \"-Xms4g -Xmx4g  -XX:+UseSerialGC\" 
  MarkDuplicates \\
  -I ${cram_file_sorted} \\
  -O ${cram_file_sorted}.dups \\
  --METRICS_FILE ${cram_file_sorted}.dups.metrics \\
  -R /ref/hs38DH.fa \\
  --VALIDATION_STRINGENCY SILENT \\
  --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \\
  --ASSUME_SORT_ORDER \"queryname\" \\
  --CLEAR_DT \"true\" \\
  --ADD_PG_TAG_TO_READS false \\
  --TMP_DIR ./tmp
"""
}

