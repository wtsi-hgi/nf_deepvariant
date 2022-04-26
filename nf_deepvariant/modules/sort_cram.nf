//params.ref_dir = "/lustre/scratch118/humgen/resources/ref/Homo_sapiens/HS38DH/"

process sort_cram {
    tag "$sample"
    publishDir path: "${params.outdir}/sort_cram/",
	mode: "${params.copy_mode}",
	overwrite: "true"

    input:
    tuple val(sample), path(cram_file), path(cram_file_index)

    output:
    tuple val(sample), path("${cram_file}.sorted"), emit: sorted_sample_cram

    script:
""" 
/opt/samtools/bin/samtools view -b ${cram_file} \\
  | sambamba sort -p -m 7GB -n --tmpdir ./tmp /dev/stdin -o ${cram_file}.sorted

/opt/samtools/bin/samtools quickcheck ${cram_file}.sorted
"""
}

