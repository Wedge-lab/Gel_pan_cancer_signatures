# run combineSignatures for ID83
# Iteratively add cohort extracted signatures to pan-cancer COSMIC list

#!/bin/bash
parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd $parent_path
source ../../.env

min_stability=1.0
# Signature type
sig_type=ID83
max_sigs=25
# Directory containing signatures
sig_dir=${SIG_DIR}/${sig_type}
# COSMIC signatures file - column headered Type with each mutation type
# Input signatures should be in additional columns
reference=${REF_SIGNATURES_DIR}/COSMIC_v3.3_ID_GRCh37.txt
sample_file=${SAMPLE_LIST}

dir_output=${COMBINED_SIGS_DIR}/combinedSignatures_${sig_type}
cohort_file=$dir_output/cohort_list.tsv

# create input and output data directories
mkdir -p $dir_output

log_file=$dir_output/readme.txt
echo "Combined Signatures ${sig_type}" > $log_file
echo "min_stability ${min_stability}" >> $log_file
echo "max_sigs ${max_sigs}" >> $log_file
echo "Min similarity selection" >> $log_file
echo "AIC solution selection" >> $log_file
echo "COSMIC reference ${cosmic}" >> $log_file
echo "Cohort signatures ${sig_dir}" >> $log_file

# Make cohort list
echo -e "cohort\tcohort_dir" > $cohort_file
for cohort in $(ls $sig_dir)
do
    if [ true ]
    then
        echo $cohort
        cohort_dir=${sig_dir}/${cohort}/${sig_type}_random_max${max_sigs}sigs_500nmf_reps
        echo -e "${cohort//CANCERTYPE_/}\t${cohort_dir}" >> $cohort_file
    fi
done

# Aggregate file paths for SigProfiler-type output
echo "getFilesSigProfiler..."
python ../../src/signatures/combineSignatures/_01_getFilesSigProfiler.py $sig_dir $cohort_file $dir_output ${sig_type} ${max_sigs}

# # Get best solution for each cohort
echo "bestSolution..."
python ../../src/signatures/combineSignatures/_02_bestSolutionAIC.py ${dir_output}/all_cohort_stats.tsv ${max_sigs} $min_stability $dir_output ${sig_type}

# Combine signatures between cohorts
echo "combineSignatures..."
python ../../src/signatures/combineSignatures/_03_combineSignatures.py ${dir_output}/best_solutions.tsv $reference ${dir_output} ${sig_type} tumour_sample_platekey

# Plot results of signature combiner
echo "plotResults..."
python ../../src/signatures/combineSignatures/_04_plotResults.py ${dir_output} ${sig_type} ${sample_file}

chmod 755 -R $dir_output
