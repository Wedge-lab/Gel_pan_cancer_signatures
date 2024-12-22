#!/usr/bin/env nextflow
// cancGermlineHits.nf pipeline
// written by Andy Everall using code from Ben Kinnersley
//

// help message
def helpMessage() {
    log.info"""
    Mandatory arguments:
      --filename_genes        Gene loci to extract - tsv file with chromosome, start, end, ID
      --filename_clinvar_vcf  vcf file of clinvar mutations
      --dir_output            Directory to save results

    Optional arguments:
      --projectName	     LSF project name that jobs are submitted with (e.g. re_gecip_cancer_colorectal for CRC GeCIP). Default: re_gecip_cancer_pan
      --n_files          Number of vcf files to process (put in a small number if doing a trial run)
    """.stripIndent()
}
if (params.help){
    helpMessage()
    exit 0
}

// Input channels
chrom_ch = Channel.from(1..22).concat(Channel.of('X'))
file_genes = file(params.filename_genes)

process clinvar {
    tag "Filter ClinVar"
    label 'standard'

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
    bcftools view -R chr${chromosome}_gene_regions.txt ${params.filename_clinvar_vcf} | \
                  bcftools annotate --rename-chrs rename_chrs.txt > clinvar_filtered_chr${chromosome}.vcf
    """
    // -i 'CLNSIG="Pathogenic/i"|CLNSIG="Likely_pathogenic"|(CLNSIGCONF~"Pathogenic/i"&CLNSIGCONF!~"Benign/i")'
}

process cadd {
    tag "Run CADD"
    label 'standard'

    input:
    set chromosome, vcf_file from clinvar_filtered_ch

     output:
     set chromosome, vcf_file, file("*.cadd.txt.gz") into caddout_ch

    module 'bio/BCFtools:tools/snakemake/5.7.1-foss-2019a-Python-3.7.2:lang/Anaconda3/2021.11'

    """
    CADD_OUT="${params.dir_output}/cadd_output/${chromosome}.cadd.txt.gz"

    # For renameing chromosomes (CADD only uses numbers without chr)
    chrom=${chromosome}
    echo -e "${chromosome} \${chrom#chr}" > rename_chrs.txt

    echo "bcftools annotate to change chromosome label format from chrX to X"
    bcftools annotate --rename-chrs rename_chrs.txt ${vcf_file} > clinvar_${chromosome}_reformat.vcf

    # TODO set core number equal to number of cpus
    echo "Run CADD"
    bash ${params.cadd_cmd} -c 30 -g GRCh38 -o clinvar_${chromosome}.cadd.txt.gz clinvar_${chromosome}_reformat.vcf
    """
}

// Group genotypes by chromosome
// caddout_ch.collect().into { cadd_grouped_ch; view_ch }

process germline_hits {
    tag "Extract mutations"
    label 'standard'

    input:
    tuple val(chromosome), file(clinvar_vcf), file(clinvar_cadd) from caddout_ch

    output:
    file("clinvar.tsv") into clinvar_tsv_ch
    file("CADD_ClinVar.tsv") into cadd_gene_ch

    module 'bio/BCFtools:singularity/3.2.1:lang/Anaconda3/2021.11'

    """
    echo "${chromosome}"
    # Activate python virtual environment
    TODO: use conda environment from .env

    echo "Get ClinVar data table from vcf"
    echo -e 'CHROM\tPOS\tID\tREF\tALT\tAF_ESP\tAF_EXAC\tAF_TGP\tALLELEID\tCLNDN\tCLNDNINCL\tCLNDISDB\t\
CLNDISDBINCL\tCLNHGVS\tCLNREVSTAT\tCLNSIG\tCLNSIGCONF\tCLNSIGINCL\t\
CLNVC\tCLNVCSO\tCLNVI\tDBVARID\tGENEINFO\tMC\tORIGIN\tRS\tSSR' > clinvar.tsv
    bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AF_ESP\t%INFO/AF_EXAC\t%INFO/AF_TGP\t%INFO/ALLELEID\t%INFO/CLNDN\t%INFO/CLNDNINCL\t%INFO/CLNDISDB\t%INFO/\
CLNDISDBINCL\t%INFO/CLNHGVS\t%INFO/CLNREVSTAT\t%INFO/CLNSIG\t%INFO/CLNSIGCONF\t%INFO/CLNSIGINCL\t%INFO/\
CLNVC\t%INFO/CLNVCSO\t%INFO/CLNVI\t%INFO/DBVARID\t%INFO/GENEINFO\t%INFO/MC\t%INFO/ORIGIN\t%INFO/RS\t%INFO/SSR\n' ${clinvar_vcf} >> clinvar.tsv

    # Import to Hail format
    # export PYSPARK_SUBMIT_ARGS='--driver-memory ${task.cpus}0G pyspark-shell --conf spark.port.maxRetries=50'
    # singularity exec -B ${params.bind} ${baseDir}/singularity/hail.img /usr/bin/
    python3 ${baseDir}/scripts/clinvar_hail.py \
        ${chromosome} \
        ${clinvar_cadd} \
        ${params.filename_genes} \
        1 \
        'GRCh38'
    """
}
cadd_gene_ch.collectFile(name: "CADD_ClinVar_all.tsv", keepHeader: true, skip: 1, storeDir: "${params.dir_output}")
clinvar_tsv_ch.collectFile(name: "clinvar_all.tsv", keepHeader: true, skip: 1, storeDir: "${params.dir_output}")
