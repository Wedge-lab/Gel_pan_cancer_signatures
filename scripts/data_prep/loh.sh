#!/bin/bash
# run cancGeneSignatures.nf

parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd $parent_path
source ../../.env

# Run parameters
chunk_size=100 # Number of samples to run on a single process

# directories and filenames
dir_analysis=../../data/cancGeneHits/somatic
filename_battenberg_list=${dir_analysis}/input/battenberg_list.tsv #sample_data.tsv
filename_sample_list=${dir_analysis}/input/loh_sample_list.tsv #sample_data.tsv
filename_genes=${dir_analysis}/input/gene_loci.tsv # Download from https://cancer.sanger.ac.uk/cosmic/download

dir_output=${dir_analysis}/output

# create output directory
mkdir -p $dir_analysis/input
mkdir -p $dir_output

# Battenberg file list
awk 'NR==1 {for (i=1; i<=NF; i++) {f[$i] = i}}{ print $(f["tumour_sample_platekey"])"\t"$(f["filename_cna"]) }' ${SAMPLE_LIST} | sed 's/filename_cna/filename_batt/' > ${filename_battenberg_list}

# Gene regions from
echo -e "chrom\tstart\tend\tgene_id" > ${filename_genes}
sed -e 's/\t/_/g' ../../data/cancer_gene_census.csv | awk -F'","' '{print $4"\t"$1}' | sed -e 's/"//g' -e 's/:/\t/' -e 's/-/\t/' | tail -n +2 | sed -e 's/^/chr/'  >> ${filename_genes}

# Make sample list
python make_somatic_input.py $dir_analysis/input $filename_sample_list

#### code below this point should not need to be changed ####

# create working directory and set as current
dir_wd=${PROJECT_DIR} # ${dir_output}/nextflow
dir_wd_old=$( pwd )
mkdir -p $dir_wd
cd $dir_wd

# create and set temporary directory
dir_tmp=${PROJECT_DIR}/tmp # /re_scratch/${username}/tmp_run_gwas_${project_code}
mkdir -p $dir_tmp
export NXF_TEMP=$dir_tmp

echo "Run nextflow"

# Run pipeline to get function mutations in cancer genes
nextflow run ../../src/pipelines/cancGeneHits/cancLOH.nf -with-trace \
  --run_name "test"\
  --filename_battenberg_list ${filename_battenberg_list}\
  --filename_sample_list ${filename_sample_list}\
  --filename_genes ${filename_genes}\
  --chunk_size ${chunk_size}\
  --dir_output ${dir_output}\
  --projectName ${project_code}

# move back to original working directory
cd $dir_wd_old
