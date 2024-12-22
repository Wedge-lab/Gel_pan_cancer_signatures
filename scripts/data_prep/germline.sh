#!/bin/bash
# run cancGeneHits.nf

parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd $parent_path

source ../../.env

# parameters
PROJECT_DIR=../workdir

# Run parameters
n_files=10000 # Max number of files to run (input a small number if just testing)
PHRED_threshold=20

# directories and filenames
dir_analysis=../../data/cancGeneHits/germline_${PHRED_threshold} #test_data
germline_dir=GERMLINE_DIR
cadd_cmd=CLINVAR_CADD_CMD
filename_germline_vcfs=${dir_analysis}/input/germline_vcfs.tsv
filename_sample_list=${dir_analysis}/input/germline_sample_list.tsv
filename_genes=../../data/cancer_gene_census.csv # Download from https://cancer.sanger.ac.uk/cosmic/download
filename_clinvar_vcf=../../data/clinvar/clinvar.vcf.gz

dir_output=${dir_analysis}/output

#### code below this point should not need to be changed ####

BASEDIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
echo $BASEDIR

# create input and output data directories
mkdir -p $dir_output
mkdir -p $dir_analysis/input

if true; then
    # Germline vcf file list
    ls ${germline_dir} > ${dir_analysis}/germline_file_list.txt
    cat ${dir_analysis}/germline_file_list.txt | rev | cut -d'/' -f1 | rev | cut -d'_' -f4 > ${dir_analysis}/chr_list.txt
    echo -e "chromosome\tfilename_vcf" > ${filename_germline_vcfs}
    paste -d$'\t' ${dir_analysis}/chr_list.txt ${dir_analysis}/germline_file_list.txt >> ${filename_germline_vcfs}
    rm ${dir_analysis}/germline_file_list.txt ${dir_analysis}/chr_list.txt

    # Gene regions
    sed -e 's/\t/_/g' ${filename_genes} | awk -F'","' '{print $4"\t"$1}' | sed -e 's/"//g' -e 's/:/\t/' -e 's/-/\t/' | tail -n +2 | sed -e 's/^/chr/' \
    | awk 'BEGIN { FS = OFS = "\t" } { for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = 0 }; 1'  > ${dir_analysis}/input/gene_loci.tsv

    # Make sample list
    python make_germline_input.py $dir_analysis/input $filename_sample_list
fi

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
module load bio/nextflow/21.04.3
nextflow run ../../src/pipelines/cancGeneHits/clinvarCadd.nf -with-trace -resume \
    --filename_germline_vcfs ${filename_germline_vcfs}\
    --filename_sample_list ${filename_sample_list}\
    --filename_genes ${dir_analysis}/input/gene_loci.tsv\
    --filename_clinvar_vcf ${filename_clinvar_vcf}\
    --dir_output ${dir_output}\
    --projectName ${project_code}\
    --n_files ${n_files}\
    --PHRED_threshold ${PHRED_threshold}

# move back to original working directory
cd $dir_wd_old
