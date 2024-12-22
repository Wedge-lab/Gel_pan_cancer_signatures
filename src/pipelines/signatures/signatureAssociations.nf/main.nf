#!/usr/bin/env nextflow
// signatureAssociations.nf pipeline
// written by Andy Everall
//

// help message
def helpMessage() {
    log.info"""
    Mandatory arguments:
      --filename_samples     TSV file, first col 'sample_id', second col 'hypermutation', third col 'category' and covariates for all samples
      --filename_signatures  TSV file, first col 'sample_id' and signature activities
      --filename_targets     TSV file, first col 'sample_id' and variables to be fit
      --filename_tests       TSV file, columns "signature", "group", "target", "type"
      --filename_tests_binary       TSV file, columns "signature", "group", "target"
      --output_dir           Directory name into which results will be saved

    Optional arguments:
      --projectName	     LSF project name that jobs are submitted with (e.g. re_gecip_cancer_colorectal for CRC GeCIP). Default: re_gecip_cancer_pan
      --max_tests        Max number of signatures to run in the method - use a small number when testing. Default: 1000
    """.stripIndent()
}
if (params.help){
    helpMessage()
    exit 0
}

// Get covariates file channel
Channel.from(file(params.filename_tests)).into{ tests_ch; test_hits_ch }
tests_ch.splitCsv(by: 1, skip: 1, sep: "\t", limit: params.max_tests ) //header: ['group', 'signature', 'target'],
        .groupTuple(by: [0,1]).into{ test_ch; test_lognmax_ch; test_binary_ch; test_nb_ch }

process zinb {
    tag "Run fits using glm in R"
    label 'short_job'//'standard'//

    input:
    tuple val(signature), val(group), val(targets), val(types) from test_ch

    output:
    file("*_glm_results.csv") into glm_results_ch

    """
    echo "${group}, ${signature}"

    # Get rows of test file for the given group signature combination
    echo -e "signature\\tgroup\\ttarget\ttype" > tests_subset.tsv
    grep -P "${signature}\\t${group}" ${params.filename_tests} >> tests_subset.tsv

    Rscript --vanilla ${projectDir}/scripts/zinb_regress.r \
                ./ \
                ${projectDir}/scripts \
                tests_subset.tsv \
                ${params.filename_samples} \
                ${params.filename_signatures} \
                ${params.filename_targets} \
                ${params.resampling_method} \
                ${params.power_analysis} \
                true
    mv glm_results.csv "${group}_${signature}_glm_results.csv"
    """
}
// glm_summary_ch = glm_results_ch.collectFile(name: "signature_target_assoc_zinb.csv", keepHeader: true, skip: 1, storeDir: "${params.output_dir}")

process negbin {
    tag "Run fits using glm in R"
    label 'short_job'//'standard'//

    input:
    tuple val(signature), val(group), val(targets), val(types) from test_nb_ch

    output:
    file("*_nb_results.csv") into nb_results_ch

    """
    echo "${group}, ${signature}"

    # Get rows of test file for the given group signature combination
    echo -e "signature\\tgroup\\ttarget\ttype" > tests_subset.tsv
    grep -P "${signature}\\t${group}" ${params.filename_tests} >> tests_subset.tsv

    Rscript --vanilla ${projectDir}/scripts/zinb_regress.r \
                ./ \
                ${projectDir}/scripts \
                tests_subset.tsv \
                ${params.filename_samples} \
                ${params.filename_signatures} \
                ${params.filename_targets} \
                ${params.resampling_method} \
                ${params.power_analysis} \
                false
    mv glm_results.csv "${group}_${signature}_nb_results.csv"
    """
}
// glm_summary_ch = glm_results_ch.collectFile(name: "signature_target_assoc_zinb.csv", keepHeader: true, skip: 1, storeDir: "${params.output_dir}")


// // Get covariates file channel
// Channel.from(file(params.filename_tests_binary)).set{ tests_binary_ch }
// tests_binary_ch.splitCsv(by: 1, skip: 1, sep: "\t", limit: 100000 ) //header: ['group', 'signature', 'target'],
//         .groupTuple(by: [0,1]).set{ test_binary_ch }

process logistic {
    tag "Run fits using logistic regression in R"
    label 'short_job'//'standard'//

    input:
    tuple val(signature), val(group), val(targets), val(types) from test_binary_ch

    output:
    file("*_logistic_results.csv") into binary_results_ch

    """
    echo "${group}, ${signature}"


    # Get rows of test file for the given group signature combination
    echo -e "signature\\tgroup\\ttarget\ttype" > tests_subset.tsv
    grep -P "${signature}\\t${group}" ${params.filename_tests} >> tests_subset.tsv

    Rscript --vanilla ${projectDir}/scripts/logistic_regress.r \
                      ./ \
                      ${projectDir}/scripts \
                      tests_subset.tsv \
                      ${params.filename_samples} \
                      ${params.filename_signatures} \
                      ${params.filename_targets} \
                      ZERO50 \
                      NULL
    mv logistic_results.csv "${group}_${signature}_logistic_results.csv"
    """
}



// Get covariates file channel
test_hits_ch.splitCsv(by: 1, skip: 1, sep: "\t", limit: params.max_tests ) //header: ['group', 'signature', 'target'],
             .groupTuple(by: [2,3])
             .set{ test_hit_ch }

process hithit {
    tag "Regressing hits against one another in R"
    label 'short_job'//'standard'//

    input:
    tuple val(signatures), val(group), val(target), val(type) from test_hit_ch

    output:
    file("*_target_results.csv") into target_results_ch

    """
    echo "${group}, ${target}, ${type}"

    # Get rows of test file for the given group signature combination
    echo -e "group\\ttarget" > tests_subset.tsv
    cat ${params.filename_tests} | tail -n +2 | cut -d \$'\\t' -f2,3 | sort | uniq >> tests_subset.tsv
    Rscript --vanilla ${projectDir}/scripts/target_regress.r \
                ${projectDir}/scripts \
                tests_subset.tsv \
                ${params.filename_samples} \
                ${params.filename_targets} \
                ${target} ${type}
    mv target_results.csv "${target}_target_results.csv"
    """
}

// Concatenate all results together into single files
glm_allresults_ch = glm_results_ch.collect()
nb_allresults_ch = nb_results_ch.collect()
binary_allresults_ch = binary_results_ch.collect()
// lognmax_allresults_ch = lognmax_results_ch.collect()
target_allresults_ch = target_results_ch.collect()
process concatenate {
    tag "Run fits using glm in R"
    label 'short_job'//'standard'//

    publishDir "${params.output_dir}", mode: "copy"

    input:
    file(glm_results) from glm_allresults_ch
    file(nb_results) from nb_allresults_ch
    file(logistic_results) from binary_allresults_ch
    // file(lognmax_results) from lognmax_allresults_ch
    file(target_results) from target_allresults_ch

    output:
    tuple file("signature_target_assoc_zinb.csv"), file("signature_target_assoc_logistic.csv"), file("signature_target_assoc_nb.csv"), file("target_target_assoc.csv")  into summary_ch
    // file("signature_target_assoc_lognmax.csv"),

    """
    Rscript --vanilla ${projectDir}/scripts/concatenate.r "." "*_glm_results.csv" "signature_target_assoc_zinb.csv"
    Rscript --vanilla ${projectDir}/scripts/concatenate.r "." "*_nb_results.csv" "signature_target_assoc_nb.csv"
    Rscript --vanilla ${projectDir}/scripts/concatenate.r "." "*_logistic_results.csv" "signature_target_assoc_logistic.csv"
    # Rscript --vanilla ${projectDir}/scripts/concatenate.r "." "*_lognmax_results.csv" "signature_target_assoc_lognmax.csv"
    Rscript --vanilla ${projectDir}/scripts/concatenate.r "." "*_target_results.csv" "target_target_assoc.csv"
    """
}

process analysis {
    tag "Run fits using glm in R"
    label 'short_job'//'standard'//

    publishDir "${params.output_dir}/figures", pattern: "*.png", mode: "copy"

    errorStrategy "terminate"

    input:
    tuple file(zinb_file), file(logistic_file), file(nb_file), file(target_file) from summary_ch
    // file(lognmax_file),

    output:
    file("*.png") into figures_ch

    """
    module load lang/Anaconda3/2021.11

    mkdir -p "${params.output_dir}/figures"
    python3 ${projectDir}/scripts/analysis.py ${zinb_file} ${logistic_file}
    # \${lognmax_file}
    """
}
