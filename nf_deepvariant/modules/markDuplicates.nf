//params.run = true

//params.ref_dir = "/lustre/scratch118/humgen/resources/ref/Homo_sapiens/HS38DH/"

process markDuplicates {
    memory '8G'
    tag "$sample"
    //cpus 1
    disk '20 GB'
    //time '100m'
    //queue 'normal'
    container  = 'file:///software/hgi/containers/gatk_4.2.4.1.sif'
    containerOptions = "--bind /lustre --bind ${params.ref_dir}:/ref --bind /tmp:/tmp"
    // errorStrategy 'terminate'
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    //publishDir "${params.outdir}/cram_index/", mode: 'symlink', overwrite: true, pattern: "${cram_file}.crai"
    maxRetries 3

    when:
    params.run_markDuplicates
     
    input:
    tuple val(study_id), val(sample), path(cram_file_sorted)
    //path cram_file_sorted

    output:
    tuple val(study_id), val(sample), path("${cram_file_sorted}.dups"), emit: markdup_sample_cram
    path "${cram_file_sorted}.dups"
    //tuple file("${cram_file}.sorted"), emit: indexes

    script:
""" 
/gatk/gatk --java-options "-Xms4g -Xmx4g  -XX:+UseSerialGC" MarkDuplicates -I ${cram_file_sorted} -O ${cram_file_sorted}.dups --METRICS_FILE ${cram_file_sorted}.dups.metrics -R /ref/hs38DH.fa --VALIDATION_STRINGENCY SILENT --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 --ASSUME_SORT_ORDER "queryname" --CLEAR_DT "true" --ADD_PG_TAG_TO_READS false --TMP_DIR ./tmp
"""
}

