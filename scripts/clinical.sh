#!/bin/bash
parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd $parent_path
source ../.env

# directories and filenames
run_name=clinical_tumour-group_InDelQC_nb
run_meta="Clinical variable analysis against signature activities by tumour group. 06/04/2022."
resampling_method='dCRT'

# Where the input files for the nextflow pipeline are stored
dir_input=../data/clinicalSigs/${run_name}/input
# Where nextflow will save files
dir_output=../data/clinicalSigs/${run_name}/output
dir_output_cohort=../data/clinicalSigs/${run_name}/output_cohort

filename_samples=${dir_input}/samples.tsv
filename_signatures=${dir_input}/signatures.tsv
filename_targets=${dir_input}/targets.tsv
filename_tests=${dir_input}/tests.tsv
filename_tests_binary=${dir_input}/tests_binary.tsv

filename_samples_cohort=${dir_input}/samples_cohort.tsv
filename_targets_cohort=${dir_input}/targets_cohort.tsv
filename_tests_cohort=${dir_input}/tests_cohort.tsv
filename_tests_cohort_binary=${dir_input}/tests_cohort_binary.tsv

# create output directory
mkdir -p $dir_output
mkdir -p $dir_output_cohort
mkdir -p $dir_input

# Direct a copy of the output to a log file
LOGFILE="${dir_output}/clinical.log"
echo -e $run_meta | tee -a $LOGFILE

# Put covariates variables together and subset samples. Set tumour_sample_platekey as index
# Should also include category column and hypermutation column
python3 ../src/signatures/associations/clinical/samples.py ${filename_samples} | tee -a $LOGFILE
# Put all signatures together. Set tumour_sample_platekey as index.
python3 ../src/signatures/associations/clinical/signatures.py ${filename_samples} ${filename_signatures} | tee -a $LOGFILE

python3 ../src/signatures/associations/clinical/targets_stage.py ${filename_samples} ${filename_targets} | tee -a $LOGFILE

python3 ../src/signatures/associations/clinical/histology.py ${filename_samples} ${filename_samples_cohort} ${filename_targets_cohort} | tee -a $LOGFILE

python3 ../src/signatures/associations/clinical/tests.py ${filename_samples} ${filename_samples_cohort} ${filename_signatures} ${filename_targets} ${filename_targets_cohort} \
                    ${filename_tests} ${filename_tests_cohort} ${filename_tests_binary} ${filename_tests_cohort_binary} | tee -a $LOGFILE


# create working directory and set as current
dir_wd=${PROJECT_DIR}
dir_wd_old=$( pwd )
mkdir -p $dir_wd
cd $dir_wd

nextflow run $parent_path/../src/signatures/signatureAssociations.nf -with-trace -dsl1 \
        --filename_tests ${filename_tests} \
        --filename_tests_binary ${filename_tests_binary} \
        --filename_signatures ${filename_signatures} \
        --filename_targets ${filename_targets} \
        --filename_samples ${filename_samples} \
        --max_tests 100000 \
        --output_dir ${dir_output} \
        --resampling_method $resampling_method \
        --executor 'local' \
        --power_analysis false

nextflow run $parent_path/../src/signatures/signatureAssociations.nf -with-trace -dsl1 \
        --filename_tests ${filename_tests_cohort} \
        --filename_tests_binary ${filename_tests_cohort_binary} \
        --filename_signatures ${filename_signatures} \
        --filename_targets ${filename_targets_cohort} \
        --filename_samples ${filename_samples_cohort} \
        --max_tests 100000 \
        --output_dir ${dir_output_cohort} \
        --resampling_method $resampling_method \
        --executor 'local' \
        --power_analysis false

# move back to original working directory
cd $dir_wd_old
