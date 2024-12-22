
# Run OncoKB
# Runs OncoKB for every single individual in the cohort. We later perform additional processing
#


username=aeverall
project_dir=/re_scratch/re_gecip/cancer_pan/${username}/projects/cancgene_sigs
output_dir=/re_gecip/cancer_pan/${username}/data/

module purge
export PATH=/usr/share/lsf/10.1/linux3.10-glibc2.17-x86_64/etc:/usr/share/lsf/10.1/linux3.10-glibc2.17-x86_64/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/home/rculliford/.local/bin:/home/rculliford/bin

cd ${project_dir}
# /re_gecip/cancer_renal_cell/analysisResults/6.driverAnalysis/Updated_OncoKB_Annotation/VEP/nextflow_work_dirs

module load singularity/3.2.1

/re_gecip/shared_allGeCIPs/bkinnersley/nextflow/nextflow run /re_gecip/cancer_renal_cell/rculliford/driver_analysis/Scripts/Initiation_Scripts/run_vep_CADD_UTRannotate_pancan_UsingRAWCADD.nf \
--input_samples /re_gecip/cancer_renal_cell/analysisResults/6.driverAnalysis/Updated_OncoKB_Annotation/SampleLists/RCC_input_VEP.tsv \
--output_dir /re_gecip/cancer_renal_cell/analysisResults/6.driverAnalysis/VEP_Annotated_VCF/ \
-c /re_gecip/shared_allGeCIPs/bkinnersley/intogen/rerun_VEP/run_vep.conf -resume -with-report -with-trace


#module load lang/R/4.1.0-foss-2019b

mkdir -p /re_gecip/cancer_renal_cell/analysisResults/6.driverAnalysis/VEP_Annotated_VCF/
mkdir -p /re_gecip/cancer_renal_cell/analysisResults/6.driverAnalysis/Updated_OncoKB_Annotation/VEP/nextflow_work_dirs

module purge
export PATH=/usr/share/lsf/10.1/linux3.10-glibc2.17-x86_64/etc:/usr/share/lsf/10.1/linux3.10-glibc2.17-x86_64/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/home/rculliford/.local/bin:/home/rculliford/bin

cd /re_gecip/cancer_renal_cell/analysisResults/6.driverAnalysis/Updated_OncoKB_Annotation/VEP/nextflow_work_dirs

module load singularity/3.2.1

/re_gecip/shared_allGeCIPs/bkinnersley/nextflow/nextflow run /re_gecip/cancer_renal_cell/rculliford/driver_analysis/Scripts/Initiation_Scripts/run_vep_CADD_UTRannotate_pancan_UsingRAWCADD.nf \
--input_samples /re_gecip/cancer_renal_cell/analysisResults/6.driverAnalysis/Updated_OncoKB_Annotation/SampleLists/RCC_input_VEP.tsv \
--output_dir /re_gecip/cancer_renal_cell/analysisResults/6.driverAnalysis/VEP_Annotated_VCF/ \
-c /re_gecip/shared_allGeCIPs/bkinnersley/intogen/rerun_VEP/run_vep.conf -resume -with-report -with-trace

## --input_samples /re_gecip/cancer_colorectal/analysisResults/13.driverAnalysis/OncoKB_annotation/CRC_v8_ALL_mut_input_for_nextflow.tsv

mkdir -p /re_gecip/cancer_renal_cell/analysisResults/6.driverAnalysis/Updated_OncoKB_Annotation/output
mkdir -p /re_gecip/cancer_renal_cell/analysisResults/6.driverAnalysis/Updated_OncoKB_Annotation/nextflow_work_dirs

cd /re_gecip/cancer_renal_cell/analysisResults/6.driverAnalysis/Updated_OncoKB_Annotation/nextflow_work_dirs

module load singularity/3.2.1

/re_gecip/shared_allGeCIPs/bkinnersley/nextflow/nextflow run \
/re_gecip/shared_allGeCIPs/bkinnersley/intogen/OncoKB/run_OncoKB_mutations.nf -resume -with-report -with-trace \
--input_samples /re_gecip/cancer_renal_cell/analysisResults/6.driverAnalysis/Updated_OncoKB_Annotation/SampleLists/RCC_input_oncokb.tsv \
--output_dir /re_gecip/cancer_renal_cell/analysisResults/6.driverAnalysis/Updated_OncoKB_Annotation/output \
-c /re_gecip/shared_allGeCIPs/bkinnersley/intogen/OncoKB/run_OncoKB.conf
