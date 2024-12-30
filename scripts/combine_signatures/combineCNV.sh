# run combineSignatures.py
# Iteratively add cohort extracted signatures to pan-cancer COSMIC list

#!/bin/bash
parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd $parent_path

min_stability=1.0
# Signature type
sig_type=CNV48
# Directory containing signatures
sig_dir=${DATA_DIR}/copy_number_signatures/SigProfilerCNV48
# COSMIC signatures file - column headered Type with each mutation type
# Input signatures should be in additional columns
# If there aren't any COSMIC reference signatures or extracting signatures deNovo, only include the Type column
reference=${DATA_DIR}/COSMIC_v3.3_CN_GRCh37.txt
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
    if [ $cohort != "eofiles" ] && [ $cohort != "Summary" ] && [ $cohort != "Pan" ]
    then
        echo $cohort
        cohort_dir=${sig_dir}/${cohort}/sigs_1_to_15_bestrun_inclFailed/CNV48
        echo -e "${cohort}\t${cohort_dir}" >> $cohort_file
    fi
done

# # Aggregate file paths for SigProfiler-type output
echo "getFilesSigProfiler..."
python ../../src/signatures/combineSignatures/_01_getFilesSigProfiler.py $sig_dir $cohort_file $dir_output SBSCNV 15

# # Get best solution for each cohort
echo "bestSolution..."
python ../../src/signatures/combineSignatures/_02_bestSolution.py ${dir_output}/all_cohort_stats.tsv 15 $min_stability $dir_output CNV48

# # Combine signatures between cohorts
echo "combineSignatures..."
python ../../src/signatures/combineSignatures/_03_combineSignatures.py ${dir_output}/best_solutions.tsv $reference ${dir_output} CNV48 pid_germline_tumour

# Plot results of signature combiner
echo "plotResults..."
python ../../src/signatures/combineSignatures/_04_plotResults.py ${dir_output} ${sig_type} ${sample_file}


chmod 755 -R $dir_output
