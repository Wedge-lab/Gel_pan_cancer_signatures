#!/bin/bash
#!/bin/bash
parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd $parent_path

# parameters
PROJECT_DIR=../workdir

# directories and filenames
resampling_method=dCRT
target_style=multibinom2
run_name=twohit_lognloh

# Where the input files for the nextflow pipeline are stored
dir_input=../data/cancGeneHits/${run_name}/input
# Where nextflow will save files
dir_output=../data/cancGeneHits/${run_name}/output

# Files to be created
filename_samples=${dir_input}/samples.tsv
filename_signatures=${dir_input}/signatures.tsv
filename_targets=${dir_input}/targets.tsv
filename_tests=${dir_input}/tests.tsv
filename_tests_binary=${dir_input}/tests_logistic.tsv

# create output directory
mkdir -p $dir_output
mkdir -p $dir_input

# Direct a copy of the output to a log file
LOGFILE="${dir_output}/twohit_assoc.log"
echo -e "Combined signature activities vs Two hit cancer gene hits with dCRT - multiple binomials for germline, somatic, loh." | tee -a $LOGFILE
echo -e "Tumour-group. SBS,DBS,ID,CNV,SV." | tee -a $LOGFILE

# Put covariates variables together and subset samples. Set tumour_sample_platekey as index
n=2 # Two hit additive model
# Should also include category column and hypermutation column
python3 ../src/signatures/associations/cancGeneHit/samples.py ${filename_samples} twohit | tee -a $LOGFILE
# Put all signatures together. Set tumour_sample_platekey as index.
python3 ../src/signatures/associations/cancGeneHit/signatures.py ${filename_samples} ${filename_signatures} | tee -a $LOGFILE
# Put all target variables together. Set tumour_sample_platekey as index
python3 ../src/signatures/associations/cancGeneHit/targets.py ${filename_samples} ${filename_targets} $n | tee -a $LOGFILE
# Generate list of signature-group-target tests to run
python3 ../src/signatures/associations/cancGeneHit/tests.py ${filename_samples} ${filename_signatures} ${filename_targets} ${filename_tests} ${filename_tests_binary} $n | tee -a $LOGFILE

# create working directory and set as current
dir_wd=${PROJECT_DIR} # ${dir_output}/nextflow
dir_wd_old=$( pwd )
mkdir -p $dir_wd
cd $dir_wd

nextflow run $parent_path/../src/signatures/signatures/signatureAssociations.nf -with-trace -dsl1 -resume \
        --filename_tests ${filename_tests} \
        --filename_tests_binary ${filename_tests_binary} \
        --filename_signatures ${filename_signatures} \
        --filename_targets ${filename_targets} \
        --filename_samples ${filename_samples} \
        --max_tests 1000000 \
        --output_dir ${dir_output} \
        --resampling_method $resampling_method \
        --power_analysis false

# move back to original working directory
cd $dir_wd_old
