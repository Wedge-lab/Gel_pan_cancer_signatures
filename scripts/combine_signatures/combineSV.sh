# run combineSignatures for SV32
# Iteratively add cohort extracted signatures to pan-cancer Reference list

#!/bin/bash

min_stability=1.0
# Signature type
sig_type=SV32
max_sigs=15
# Directory containing signatures
sig_dir=${SIG_DIR}/${sig_type}
# Reference signatures file - column headered Type with each mutation type
# Input signatures should be in additional columns
reference=${REF_SIGNATURES_DIR}/Breast560_rearrangement.signatures.tsv
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
echo "Reference ${reference}" >> $log_file
echo "Cohort signatures ${sig_dir}" >> $log_file

# Make cohort list
echo -e "cohort\tcohort_dir" > $cohort_file
for cohort_dir in $(ls -d $sig_dir/*_S1_S15_R500_PoisNoisRes_T/)
do
    cohort=$(echo $cohort_dir | rev | cut -d'/' -f2 | rev)
    if [ $cohort != "AllSamples_S1_S15_R500_PoisNoisRes_T" ]
    then
        echo $cohort
        cohort_dir=${sig_dir}/${cohort}/${sig_type}
        echo -e "${cohort//_S1_S15_R500_PoisNoisRes_T/}\t${cohort_dir}" >> $cohort_file
    fi
done

# Construct reference file with index only
reference=$dir_output/empty_reference.tsv
cohort_dir=$(ls -d $sig_dir/*S1_S15_R500_PoisNoisRes_T/ | tail -n 1)
cut -d$'\t' -f1 ${cohort_dir}/SV32/All_Solutions/SBSSV_1_Signatures/Signatures/SBSSV_S1_Signatures.txt | sed 's/MutationType/Type/' > $reference

# Aggregate file paths for SigProfiler-type output
echo "getFilesSigProfiler..."
python ../../src/signatures/combineSignatures/_01_getFilesSigProfiler.py $sig_dir $cohort_file $dir_output SBSSV ${max_sigs}

# Get best solution for each cohort
echo "bestSolution..."
python ../../src/signatures/combineSignatures/_02_bestSolutionAIC.py ${dir_output}/all_cohort_stats.tsv ${max_sigs} $min_stability $dir_output ${sig_type}

# Combine signatures between cohorts
echo "combineSignatures..."
python ../../src/signatures/combineSignatures/_03_combineSignatures.py ${dir_output}/best_solutions.tsv $reference ${dir_output} ${sig_type} pid_tumour_germline

# Plot results of signature combiner
echo "plotResults..."
python ../../src/signatures/combineSignatures/_04_plotResults.py ${dir_output} ${sig_type} ${sample_file}

chmod 755 -R $dir_output
