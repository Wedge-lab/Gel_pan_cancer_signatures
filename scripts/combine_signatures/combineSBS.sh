# run combineSignatures for SBS288
# Iteratively add cohort extracted signatures to pan-cancer COSMIC list

#!/bin/bash

min_stability=1.0
# Signature type
sig_type=SBS288
max_sigs=30
# Directory containing signatures
sig_dir=/re_gecip/cancer_pan/fprefect/botl/results/SIGmats/v2_draft/${sig_type}
# COSMIC signatures file - column headered Type with each mutation type
# Input signatures should be in additional columns
cosmic=${DATA_DIR}/COSMIC_v3.3.1_SBS_GRCh38.txt
sample_file=${DATA_DIR}/sample_lists_incl_SEGs/sample_list_2021_06_29_incl_SEGs.tsv

dir_output=../../data/combinedSignatures_${sig_type}
cohort_file=$dir_output/cohort_list.tsv

# create input and output data directories
mkdir -p $dir_output

log_file=$dir_output/readme.txt
echo "Combined Signatures ${sig_type}" > $log_file
echo "min_stability ${min_stability}" >> $log_file
echo "max_sigs ${max_sigs}" >> $log_file
echo "Min similarity selection" >> $log_file
echo "COSMIC reference ${cosmic}" >> $log_file
echo "Cohort signatures ${sig_dir}" >> $log_file

# Make cohort list
echo -e "cohort\tcohort_dir" > $cohort_file
for cohort in $(ls $sig_dir)
do
    echo $cohort
    cohort_dir=${sig_dir}/${cohort}/${sig_type}_random_max${max_sigs}sigs_500nmf_reps
    echo -e "${cohort//CANCERTYPE_/}\t${cohort_dir}" >> $cohort_file
done

# Aggregate file paths for SigProfiler-type output
echo "getFilesSigProfiler..."
python ../../src/signatures/combineSignatures/_01_getFilesSigProfiler.py $sig_dir $cohort_file $dir_output ${sig_type} ${max_sigs}

# Get best solution for each cohort
echo "bestSolution..."
python ../../src/signatures/combineSignatures/_02_bestSolution.py ${dir_output}/all_cohort_stats.tsv ${max_sigs} $min_stability $dir_output ${sig_type}

# Combine signatures between cohorts
echo "combineSignatures..."
python ../../src/signatures/combineSignatures/_03_combineSignatures.py ${dir_output}/best_solutions.tsv $cosmic ${dir_output} ${sig_type} pid_tumour_germline

# Plot results of signature combiner
echo "plotResults..."
python ../../src/signatures/combineSignatures/_04_plotResults.py ${dir_output} ${sig_type} ${sample_file}

chmod 755 -R $dir_output