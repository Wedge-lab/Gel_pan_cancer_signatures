#!/usr/bin/env nextflow
// cancGermlineHits.nf pipeline
// written by Andy Everall using code from Ben Kinnersley
//

// help message
def helpMessage() {
    log.info"""
    Mandatory arguments:
      --filename_germline_vcfs tsv file with chromosome, vcf file path columns
      --filename_sample_list  List of samples to extract
      --filename_genes        Gene loci to extract - tsv file with chromosome, start, end, ID
      --filename_clinvar_vcf  vcf file of clinvar mutations
      --dir_output            Directory to save results

    Optional arguments:
      --projectName	     LSF project name that jobs are submitted with (e.g. re_gecip_cancer_colorectal for CRC GeCIP). Default: re_gecip_cancer_pan
      --max_signatures   Max number of signatures to run in the method - use a small number when testing. Default: 1000
      --n_files          Number of vcf files to process (put in a small number if doing a trial run)
    """.stripIndent()
}
if (params.help){
    helpMessage()
    exit 0
}


// Input channels
chrom_ch = Channel.from(1..22).concat(Channel.of('X'))
file_chunks = file(params.filename_germline_vcfs)
file_genes = file(params.filename_genes)

// run time 1-10 hours with 9500 samples
Channel.from(file_chunks)
       .splitCsv(header: true, sep: '\t', limit: params.n_files )
       .map { row -> [row.chromosome, row.filename_vcf] }
       .set { chunk_ch }

process clinvar {
    tag "Filter ClinVar"
    label 'short_job'

    // errorStrategy { task.exitStatus==123  ? 'ignore':'terminate'}

    input:
    val(chromosome) from chrom_ch

    output:
    set env(chr_label), file("clinvar_filtered_*.vcf") into clinvar_filtered_ch

    module 'bio/BCFtools'

    """
    chr_label='chr${chromosome}'
    echo -e '${chromosome} chr${chromosome}\n' > rename_chrs.txt

    # Select gene regions on chromosome
    grep -P '^chr${chromosome}\t' ${params.filename_genes} | sed 's/^chr//g' > chr${chromosome}_gene_regions.txt
    # Filter by pathogenic or likely pathogenic or at least one study says pathogenic/likely pathogenic and none say benign
    bcftools view -i 'CLNSIG="Pathogenic/i"|CLNSIG="Likely_pathogenic"|(CLNSIGCONF~"Pathogenic/i"&CLNSIGCONF!~"Benign/i")' \
                  -R chr${chromosome}_gene_regions.txt ${params.filename_clinvar_vcf} | \
                  bcftools annotate --rename-chrs rename_chrs.txt > clinvar_filtered_chr${chromosome}.vcf
    """
}

process cadd {
    tag "Run CADD"
    label 'short_job'

    // errorStrategy { task.exitStatus==123  ? 'ignore':'terminate'}

    input:
    set chromosome, vcf_fname from chunk_ch
    file(gene_regions) from file_genes

    output:
    set val(chromosome), file("*_regions.vcf"), file("*.cadd.txt.gz") into mutations_ch //

    module 'bio/BCFtools:tools/snakemake/5.7.1-foss-2019a-Python-3.7.2:lang/Anaconda3/2021.11'

    """
    TODO: use conda environment from .env

    # Get unique ID of file chunk
    uid=\$(echo ${vcf_fname} | awk -F'aggV2_' '{print \$2}' | cut -d'.' -f1)

    CADD_OUT="${params.dir_output}/cadd_output/\${uid}.cadd.txt.gz"

    vcf_prefix=\$(echo ${vcf_fname} | cut -d"." -f1)

    # For renameing chromosomes (CADD only uses numbers without chr)
    chrom=${chromosome}
    echo -e "${chromosome} \${chrom#chr}" > rename_chrs.txt

    echo "Filter samples and genes with bcftools"
    # if [ -f "\$CADD_OUT" ]; then
    if false; then
        ln -s \$CADD_OUT
    else
        # Find all the variants that fall in the gene regions, keep only pass variants and only alt variants (vcf file contains pathogenic variants but alt genotype may not actually be present)
        echo "bcftools to filter regions"
        bcftools view -i 'FILTER="PASS" && GT="alt"' -R ${gene_regions} -S ${params.filename_sample_list} ${vcf_fname} > \${uid}_regions.vcf
        echo "bcftools annotate to change chromosome label format from chrX to X"
        bcftools annotate --rename-chrs rename_chrs.txt \${uid}_regions.vcf > \${uid}_reformat.vcf

        # bug in cadd - cannot have "chr" in chromosome name, needs to tab separated
        # TODO set core number equal to number of cpus
        echo "Run CADD"
        bash ${params.cadd_cmd} -c 30 -g GRCh38 -o \${uid}.cadd.txt.gz \${uid}_reformat.vcf
    fi
    """
}

// Group genotypes by chromosome
mutations_ch.groupTuple().join(clinvar_filtered_ch)
                                 .into { germline_grouped_ch; view_ch }

process germline_hits {
    tag "Extract mutations"
    label 'short_job'

    // publishDir "${params.dir_output}/mutations", mode: "copy", overwrite: true

    input:
    tuple val(chromosome), file(list_of_cadds), file(list_of_vcfs), file(clinvar_vcf) from germline_grouped_ch

    output:
    file("aggv2-germline_cancer-gene_hits_*.tsv") into germline_hits_ch
    file("CADD_ClinVar_filter_df_*.tsv") into filter_ch

    module 'bio/BCFtools:singularity/3.2.1:lang/Anaconda3/2021.11'

    """
    echo "${chromosome}, ${list_of_vcfs}"

    # Activate python virtual environment
    # TODO: use conda environment from .env

    # Make file of vcfs
    ls *_regions.vcf > vcf_list.txt
    ls *.cadd.txt.gz > cadd_list.txt

    echo "Get ClinVar data table from vcf"
    bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/CLNSIG\t%INFO/CLNSIGCONF\n' ${clinvar_vcf} > clinvar.tsv

    # Import to Hail format
    python3 ${baseDir}/scripts/germline_hail.py \
        ${chromosome} \
        'vcf_list.txt' \
        ${params.filename_sample_list} \
        'cadd_list.txt' \
        ${params.filename_genes} \
        clinvar.tsv \
        1 \
        'GRCh38' \
        ${params.PHRED_threshold}

    """
}
germline_hits_ch.collectFile(name: "aggv2-germline_cancer-gene_hits.tsv", keepHeader: true, skip: 1, storeDir: "${params.dir_output}")
filter_ch.collectFile(name: "CADD_ClinVar_filter_df.tsv", keepHeader: true, skip: 1, storeDir: "${params.dir_output}")
