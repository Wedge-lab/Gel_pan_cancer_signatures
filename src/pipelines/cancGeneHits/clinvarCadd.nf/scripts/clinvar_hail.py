#!/usr/bin/env python

## Implement the ALFRED analysis in GEL using the Battenberg calls to identify LOH
import pandas as pd, numpy as np
import sys

# Input arguments
chromosome, cadd_file, gene_file, n_cpu, reference = sys.argv[1:]
print(sys.argv[1:])

# Load list of CADD files
cadd_df = pd.read_csv(cadd_file,sep="\t",compression='gzip',
                      comment='#', names=['chrom','pos','ref','alt','raw','phred'])
cadd_df['chrom'] = "chr"+cadd_df.chrom.astype(str)

# Read in genes - chrom,start,end,gene_id
genes = pd.read_csv(gene_file, sep="\t", names=['chrom','start','end','gene_id']).fillna(0)
genes['start'] = genes['start'].astype(int)
genes['end'] = genes['end'].astype(int)
genes = genes[genes.chrom==chromosome]

# Get gene labels for CADD tsv file
cadd_df_gene = pd.DataFrame(index=cadd_df.index)
for idx,row in genes.iterrows():
    print(row.gene_id)
    cadd_df_gene[row.gene_id] = np.where((cadd_df['chrom']==row.chrom)&(cadd_df['pos']>row.start)&(cadd_df['pos']<row.end),
                                         row.gene_id, "")
cadd_df['genes'] = [",".join(filter(lambda a: a!="", x)) for x in cadd_df_gene.values.tolist()]

# Save CADD results
cadd_df.to_csv(f"CADD_ClinVar.tsv", sep="\t", header=True, index=False)
