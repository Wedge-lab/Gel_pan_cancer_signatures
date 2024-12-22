#!/bin/bash
# Combine results from OncoKB

parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd $parent_path

source ../.env

# parameters
PROJECT_DIR=../../workdir

# directories and filenames
dir_analysis=../../data/cancGeneHits/somatic
filename_oncokb_list=${dir_analysis}/input/oncokb_list.tsv
filename_genes=${dir_analysis}/input/gene_loci.tsv # Download from https://cancer.sanger.ac.uk/cosmic/download

dir_output=${dir_analysis}/output

chunk_size=100
n_samples=200

# create output directory
mkdir -p $dir_analysis/input
mkdir -p $dir_output

# OncoKB file list
echo "OncoKB file list"
oncokb_dir=ONCOKB_DIR
echo -e "tumour_sample_platekey\tfile" > ${filename_oncokb_list}
ls ${oncokb_dir} | awk -v dir=$oncokb_dir '{ print $1"\t"dir"/"$1"/"$1"_combined_mutations_OncoKB_query_flatfile.tsv" }' >> ${filename_oncokb_list}

# Gene regions
echo "Gene regions"
echo -e "chrom\tstart\tend\tgene_id" > ${filename_genes}
sed -e 's/\t/_/g' ../../data/cancer_gene_census.csv | awk -F'","' '{print $4"\t"$1}' | sed -e 's/"//g' -e 's/:/\t/' -e 's/-/\t/' | tail -n +2 | sed -e 's/^/chr/'  >> ${filename_genes}

#### code below this point should not need to be changed ####

BASEDIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
echo $BASEDIR

# create working directory and set as current
dir_wd=${project_dir} # ${dir_output}/nextflow
dir_wd_old=$( pwd )
mkdir -p $dir_wd
cd $dir_wd

# create and set temporary directory
dir_tmp=${project_dir}/tmp
mkdir -p $dir_tmp
export NXF_TEMP=$dir_tmp

# Run pipeline to get function mutations in cancer genes
nextflow run ../../src/pipelines/cancGeneHits/cancOncoKB.nf -with-trace \
  --projectName ${project_code} \
  --filename_oncokb_list ${filename_oncokb_list} \
  --filename_genes ${filename_genes} \
  --chunk_size ${chunk_size} \
  --n_samples ${n_samples} \
  --dir_output ${dir_output}

# move back to original working directory
cd $dir_wd_old
