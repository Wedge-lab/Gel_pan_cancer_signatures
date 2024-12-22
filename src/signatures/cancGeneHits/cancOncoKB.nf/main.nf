#!/usr/bin/env nextflow

// help message
def helpMessage() {
    log.info"""
    Mandatory arguments:
      --filename_oncokb_list  List of sample OncoKB files to extract
      --filename_genes        Gene loci to extract - tsv file with chromosome, start, end, ID
      --dir_output            Directory to save results

    Optional arguments:
      --projectName	     LSF project name that jobs are submitted with (e.g. re_gecip_cancer_colorectal for CRC GeCIP). Default: re_gecip_cancer_pan
      --max_signatures   Max number of signatures to run in the method - use a small number when testing. Default: 1000
      --n_samples        Number of samples to process (put in a small number if doing a trial run)
      --chunk_size       Number of samples per nextflow process
    """.stripIndent()
}
if (params.help){
    helpMessage()
    exit 0
}


// parameters
file_samples = file(params.filename_oncokb_list)

process split {
    tag "Split files"
    label 'standard'//'short_job'

    input:
    file(samples) from file_samples

    output:
    path 'sample_chunks_*' into file_chunks

    """
    HEADER=\$(head -1 ${samples})

    tail -n +2 ${samples} | split -l ${params.chunk_size} - sample_chunks_
    for i in sample_chunks_*; do
        sed -i -e "1i\$HEADER" "\$i"
    done
    """
}

process mutations {
  tag "Extract mutations"
  label 'standard'

  input:
  file(oncokb_list) from file_chunks.flatten()

  output:
  file("OncoKB_somatic_hits.tsv") into somatic_hits_ch

  module = 'lang/Anaconda3/2021.11'

  """
  python ${baseDir}/scripts/collectOncoKB.py ./ ${oncokb_list} ${params.filename_genes}
  """
}
somatic_hits_ch.collectFile(name: "OncoKB_somatic_cancer-gene_hits.tsv", keepHeader: true, skip: 1, storeDir: "${params.dir_output}")
