#!/usr/bin/env nextflow
// treatmentInducedMutations.nf pipeline
// written by Andy Everall
//

// help message
def helpMessage() {
    log.info"""
    Mandatory arguments:
      --filename_battenberg_list   List of battenberg files - tsv file with tumour_sample_platekey, filename_batt
      --filename_sample_list  List of samples to extract
      --filename_genes        Gene loci to extract - tsv file with chrom, start, end, ID
      --dir_output            Directory to save results

    Optional arguments:
      --projectName	     LSF project name that jobs are submitted with (e.g. re_gecip_cancer_colorectal for CRC GeCIP). Default: re_gecip_cancer_pan
      --run_name         Identifying name given to the run (usually would explain why it is different from previous runs). Default: testrun
      --run_meta         Description of the run and what makes it different to previous runs. Default: ''
    """.stripIndent()
}
if (params.help){
    helpMessage()
    exit 0
}

// Samples channel
file_samples=file(params.filename_sample_list)

// run time 1-10 hours with 9500 samples
// Channel.fromPath(params.filename_sample_list).splitCsv(by: params.chunk_size).randomSample(2).set{ chunk_ch }
process split_samples {
  tag "Split samples"
  label 'standard'

  input:
  file(sample_list) from file_samples

  output:
  path 'sample_chunk_*' into sample_chunks_ch

  """
  split -l ${params.chunk_size} --numeric-suffixes ${sample_list} sample_chunk_
  """

}
sample_chunks_ch2 = sample_chunks_ch.flatMap()//.randomSample(2)

process loh {
  tag "Extract LOH mutations"
  label 'standard' // 'short_job'

  input:
  file(sample_list) from sample_chunks_ch2

  output:
  // file("loh_matrix.tsv") into loh_hits_ch
  file("loh_fraction.tsv") into loh_frac_ch
  file("n_overlapping_cnv.tsv") into n_cnv_ch
  file("loh_number_size.tsv") into number_size_ch

  module 'bio/BCFtools:singularity/3.2.1:lang/Anaconda3/2021.11'

  """
  # Run LOH python script - searches from Battenberg output

  python ${baseDir}/scripts/loh.py \
      ${sample_list} \
      ${params.filename_battenberg_list} \
      ${params.filename_genes}
  """
}
// loh_hits_ch.collectFile(name: "lohclonal_gene_sample_matrix.tsv", keepHeader: true, skip: 1, storeDir: "${params.dir_output}")
loh_frac_ch.collectFile(name: "lohfrac_gene_sample_matrix.tsv", keepHeader: true, skip: 1, storeDir: "${params.dir_output}")
n_cnv_ch.collectFile(name: "ncnv_gene_sample_matrix.tsv", keepHeader: true, skip: 1, storeDir: "${params.dir_output}")
number_size_ch.collectFile(name: "loh_number_size.tsv", keepHeader: true, skip: 1, storeDir: "${params.dir_output}")
