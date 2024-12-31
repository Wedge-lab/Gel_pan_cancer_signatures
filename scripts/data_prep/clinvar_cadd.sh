#!/bin/bash
# run cancGeneHits.nf

parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd $parent_path
source ../../.env

# directories and filenames
dir_analysis=CLINVAR_DIR
cadd_cmd=../../data/cancGeneHits/clinvar
filename_genes=../../data/cancer_gene_census.csv # Download from https://cancer.sanger.ac.uk/cosmic/download
filename_clinvar_vcf=../../data/clinvar/clinvar.vcf.gz

dir_output=${dir_analysis}/output

#### code below this point should not need to be changed ####

BASEDIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
echo $BASEDIR

# create input and output data directories
mkdir -p $dir_output
mkdir -p $dir_analysis/input

# Gene regions
sed -e 's/\t/_/g' ${filename_genes} | awk -F'","' '{print $4"\t"$1}' | sed -e 's/"//g' -e 's/:/\t/' -e 's/-/\t/' | tail -n +2 | sed -e 's/^/chr/' \
| awk 'BEGIN { FS = OFS = "\t" } { for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = 0 }; 1'  > ${dir_analysis}/input/gene_loci.tsv

# create working directory and set as current
dir_wd=${project_dir} # ${dir_output}/nextflow
dir_wd_old=$( pwd )
mkdir -p $dir_wd
cd $dir_wd

# create and set temporary directory
dir_tmp=${project_dir}/tmp # /re_scratch/${username}/tmp_run_gwas_${project_code}
mkdir -p $dir_tmp
export NXF_TEMP=$dir_tmp

# Run pipeline to get function mutations in cancer genes
module load bio/nextflow
nextflow run ../../src/pipelines/cancGeneHits/clinvarCadd.nf -with-trace -resume \
    --filename_genes ${dir_analysis}/input/gene_loci.tsv\
    --filename_clinvar_vcf ${filename_clinvar_vcf}\
    --dir_output ${dir_output}\
    --projectName ${project_code}

# move back to original working directory
cd $dir_wd_old
