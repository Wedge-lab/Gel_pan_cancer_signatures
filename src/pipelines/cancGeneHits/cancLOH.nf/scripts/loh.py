#!/usr/bin/env python

## Implement the ALFRED analysis in GEL using the Battenberg calls to identify LOH
'''
 Input:
    sample_list
    gene_list - four columns: chrom, start, end, gene_id
    battenberg_list - two columns: tumour_sample_platekey, filename_batt
 Output:
    File where each row is gene, column for participant_id, with 1=LoH in this gene, 0=no LoH
'''
import pandas as pd, numpy as np
import sys

#python 1_loh_calc.py $battenberg_id $updated_germline $updated_somatic ${params.ref_gene}
# participant_id, battenberg_dir, cadd_region_file, gene_db, gene_list = sys.argv[1:]
sample_list_file, battenberg_list_file, gene_list_file = sys.argv[1:]

# Read in samples
samples = pd.read_csv(sample_list_file, sep="\t", names=['tumour_sample_platekey'])
battenberg_list = pd.read_csv(battenberg_list_file, sep="\t")
battenberg_list = battenberg_list[battenberg_list.filename_batt.astype(str)!='nan']
samples = pd.merge(samples, battenberg_list, how='inner', on='tumour_sample_platekey')
print(f"Samples: {len(samples)}")

# Read in genes
genes = pd.read_csv(gene_list_file, sep="\t", usecols=['chrom', 'start', 'end', 'gene_id'])
genes = genes[genes.start.astype(str)!='nan']
genes = genes.astype({'start':int,'end':int})
genes.chrom = genes.chrom.str[3:]
genes.reset_index(drop=True)

# Results table
# loh_results=np.zeros((len(genes), len(samples)), dtype=np.int16)#, columns=genes.gene_id, index=samples.tumour_sample_platekey) # pd.DataFrame(index=genes.gid)
loh_fraction_arr=np.zeros((len(genes), len(samples)), dtype=np.float64)
n_cnv=np.zeros((len(genes), len(samples)), dtype=np.float64)
n_loh = np.zeros(len(samples), dtype=int)
n_loh50 = np.zeros(len(samples), dtype=int)
n_loh100 = np.zeros(len(samples), dtype=int)
size_loh = np.zeros(len(samples), dtype=int)

# Columns to load from battenberg
battenberg_cols = ['chr','startpos','endpos']\
                 +[f'nMaj{i}_A' for i in range(1,3)]\
                 +[f'nMin{i}_A' for i in range(1,3)]\
                 +[f'frac{i}_A' for i in range(1,3)]

for idx, row in samples.iterrows():

    # Import battenberg data
    battenberg_df = pd.read_csv(row['filename_batt'], sep="\t", usecols=battenberg_cols).fillna(0.)

    # Get intersections between battenberg mutations and genes
    batt_chrchr, gene_chrchr = np.meshgrid(battenberg_df.chr, genes.chrom)
    batt_startstart, gene_startstart = np.meshgrid(battenberg_df.startpos, genes.start)
    batt_endend, gene_endend = np.meshgrid(battenberg_df.endpos, genes.end)
    batt_gene_overlap = (batt_chrchr==gene_chrchr)&(batt_startstart<gene_endend)&(batt_endend>gene_startstart)

    # Number of mutations overlapping gene
    n_cnv[:,idx]=np.sum(batt_gene_overlap.astype(int), axis=1)

    # Fraction of cells with LOH
    loh_fraction = np.nansum( np.array([battenberg_df[f"frac{i}_A"]*((battenberg_df[f"nMaj{i}_A"]==0)|(battenberg_df[f"nMin{i}_A"]==0)).astype(float) \
                                    for i in range(1,3)]), axis=0)

    loh_fraction_arr[:,idx]=np.nanmax( batt_gene_overlap.astype(float) * loh_fraction[None,:], axis=1 )

    # Number and total size of LoH events
    n_loh[idx] = len(battenberg_df)
    n_loh50[idx] = np.sum(loh_fraction>0.5)
    n_loh100[idx] = np.sum(loh_fraction==1)
    size_loh[idx] = np.sum(battenberg_df['endpos']-battenberg_df['startpos'])

# loh_results = pd.DataFrame(loh_results.T, columns=genes.gene_id, index=samples.tumour_sample_platekey)
# loh_results.to_csv(f"loh_matrix.tsv", sep="\t", header=True, index=True)

loh_number_size = pd.DataFrame({"n_loh":n_loh, "n_loh50":n_loh50, "n_loh100":n_loh100, "size_loh":size_loh},
                                index=samples.tumour_sample_platekey)
loh_number_size.to_csv(f"loh_number_size.tsv", sep="\t", header=True, index=True)

loh_fraction_arr = pd.DataFrame(loh_fraction_arr.T, columns=genes.gene_id, index=samples.tumour_sample_platekey)
loh_fraction_arr.to_csv(f"loh_fraction.tsv", sep="\t", header=True, index=True)

n_cnv = pd.DataFrame(n_cnv.T, columns=genes.gene_id, index=samples.tumour_sample_platekey)
n_cnv.to_csv(f"n_overlapping_cnv.tsv", sep="\t", header=True, index=True)
