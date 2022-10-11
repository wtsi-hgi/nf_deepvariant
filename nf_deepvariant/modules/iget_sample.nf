params.copy_mode = "symlink"

process 'iget_sample' {
    tag "$sample"
    publishDir "${params.outdir}/iget/no_study_id/${sample}/", mode: "${params.copy_mode}"

    input:
    val(sample)
    
  output:
    tuple val(sample), path("*.cram"), path("*.crai"), emit: sample_cram_crai optional true

  script:
    """
imeta qu -z seq -d sample = ${sample} and target = \"library\" and manual_qc = 1 | grep collection | awk -F ' ' '{print \$2}' > collection.txt || true
imeta qu -z seq -d sample = ${sample} and target = \"library\" and manual_qc = 1 | grep dataObj | awk -F ' ' '{print \$2}' > dataObj.txt || true
paste -d '/' collection.txt dataObj.txt | grep '.cram\$' > to_iget.txt || true
if [[ -f "to_iget.txt" && -s "to_iget.txt" ]]; then 
    echo "target library exist and not empty"
else 
    echo "target library not exist or empty" 
    echo try with target 1 instead
    imeta qu -z seq -d sample = ${sample} and target = 1 and manual_qc = 1 | grep collection | awk -F ' ' '{print \$2}' > collection.txt || true
    imeta qu -z seq -d sample = ${sample} and target = 1 and manual_qc = 1 | grep dataObj | awk -F ' ' '{print \$2}' > dataObj.txt || true
    paste -d '/' collection.txt dataObj.txt | grep '.cram\$' > to_iget.txt || true
fi
if [[ -f "to_iget.txt" && -s "to_iget.txt" ]]; then 
    echo "target 1 exist and not empty"
else 
    echo "target 1 not exist or empty" 
    exit 1
fi
cat to_iget.txt | while read line
do
    filename=\$(echo \${line} | sed s'/^.*\\///'g)
    echo filename is \$filename
    iget -K -f -v \${line} ${sample}__\$filename
    # get index file if exists:
    iget -K -f -v \${line}.crai ${sample}__\${filename}.crai || true
done
   """
}
