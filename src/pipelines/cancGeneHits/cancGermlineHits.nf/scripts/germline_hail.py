#!/usr/bin/env python

## Implement the ALFRED analysis in GEL using the Battenberg calls to identify LOH
'''
 Input:

 Output:

'''
import hail as hl
import pandas as pd, numpy as np
import sys
import re

# Input arguments
chromosome, vcf_list_file, sample_file, cadd_list_file, gene_file, clinvar_file, n_cpu, reference, PHRED_threshold = sys.argv[1:]
PHRED_threshold = int(PHRED_threshold)
print(sys.argv[1:])

# Load list of VCF files
vcf_files = pd.read_csv(vcf_list_file, sep="\t", header=None)[0].to_list()
print(vcf_files)

# Load list of CADD files
cadd_files = pd.read_csv(cadd_list_file, sep="\t", header=None)[0].to_list()
cadd_df = pd.concat([pd.read_csv(cadd_file,sep="\t",compression='gzip',
                                comment='#', names=['chrom','pos',
                                                    'ref','alt',
                                                    'raw','phred']) for cadd_file in cadd_files])
cadd_df['chrom'] = "chr"+cadd_df.chrom.astype(str)

# Load clinvar
clinvar_df = pd.read_csv(clinvar_file, sep="\t", names=['chrom','pos','id','ref','alt','CLNSIG','CLNSIGCONF'])

# Read in samples
samples = pd.read_csv(sample_file, sep="\t", names=["germline_sample_platekey"])

# Read in genes - chrom,start,end,gene_id
genes = pd.read_csv(gene_file, sep="\t", names=['chrom','start','end','gene_id']).fillna(0)
genes['start'] = genes['start'].astype(int)
genes['end'] = genes['end'].astype(int)
genes = genes[genes.chrom==chromosome]
print(genes[genes.start==0])
genes = genes[genes.start!=0]
print(f"genes: {len(genes)}")


# CADD and ClinVar filters
filter_df = pd.merge(cadd_df, clinvar_df, on=('chrom','pos','ref','alt'), how='left')
# Apply filters
# Negative controls
benign_match = np.vectorize(lambda x:bool(re.match('.*benign.*', x, re.IGNORECASE)))
filter_df = filter_df[(filter_df.phred>PHRED_threshold)&\
                      (~benign_match(filter_df.CLNSIG.astype(str)))]
filter_df['alleles'] = filter_df[['ref','alt']].values.tolist()

# initialise Hail - 101 retries to cope with parallel processing``
hl.init(master="local[" + n_cpu + "]", default_reference=reference, spark_conf={"spark.port.maxRetries":"101"})

# Construct matrix_table in Hail
aggv2_mt = hl.import_vcf(vcf_files, reference_genome=reference)#, force_bgz=True)
print("All variants: ", aggv2_mt.count())

# Filter on CADD and ClinVar classification
filter_mt = hl.Table.from_pandas(filter_df[['chrom','pos','alleles']])
filter_mt = filter_mt.annotate(locus=hl.locus(filter_mt.chrom, filter_mt.pos)).key_by('locus','alleles')
aggv2_mt = aggv2_mt.semi_join_rows(filter_mt)
print("CADD & ClinVar filter: ", aggv2_mt.count())

germline_hits_arr = np.zeros((len(genes), len(samples)), dtype=np.int64)
for idx, row in genes.reset_index(drop=True).iterrows():

    # Get genotype for gene
    print(row.gene_id, f"{row['chrom']}:{row['start']}-{row['end']}")
    gene_mt = hl.filter_intervals(aggv2_mt, [hl.parse_locus_interval(f"{row['chrom']}:{row['start']}-{row['end']}",
                                                                     reference_genome=reference)])
    if gene_mt.count()[0]==0:
        print(f"{row['gene_id']} - empty")
        continue

    print(f"{row['gene_id']} - not empty")

    # Get sum of genotype for each sample within gene
    agg_sum = gene_mt.annotate_cols(n_variants_sum=hl.agg.sum(gene_mt.GT.n_alt_alleles()))
    mutated = agg_sum.cols().to_pandas()

    # Save non-zero indices to hdf5
    germline_hits_arr[idx] = pd.merge(samples, mutated, how='left',
                                      left_on='germline_sample_platekey', right_on='s').n_variants_sum

germline_hits_df = pd.DataFrame(germline_hits_arr, index=genes.gene_id, columns=samples.germline_sample_platekey)
germline_hits_df.to_csv(f"aggv2-germline_cancer-gene_hits_{chromosome}.tsv", sep="\t", header=True, index=True)

filter_df.to_csv(f"CADD_ClinVar_filter_df_{chromosome}.tsv", sep="\t", header=True, index=False)