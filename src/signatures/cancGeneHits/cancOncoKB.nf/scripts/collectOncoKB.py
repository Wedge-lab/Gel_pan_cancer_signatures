#
# Collect OncoKB results from pan cancer analysis
#

import sys
import numpy as np, pandas as pd, re


if __name__=='__main__':

    print(sys.argv[1:])
    dir_output, oncokb_file, gene_list_file = sys.argv[1:]

    # Sample list
    sample_list = pd.read_csv(oncokb_file, sep="\t")
    print(sample_list.head())

    # Gene list
    genes = np.array(pd.read_csv(gene_list_file, sep="\t", usecols=['gene_id']).gene_id)

    # Iterate through samples and get gene mutations
    mutation_df = pd.DataFrame(np.zeros((len(genes), len(sample_list)), dtype=int), columns=sample_list.tumour_sample_platekey, index=genes)
    unique_entries = np.array([])
    unique_oncogenic = np.array([])
    oncogenic = np.vectorize(lambda x:bool(re.match('.*oncogenic.*', x, re.IGNORECASE)))
    for i, row in sample_list.iterrows():

        tumour_sample_platekey = row['tumour_sample_platekey']
        annotations = pd.read_csv(row['file'], sep="\t",
                                  usecols=['hugoSymbol','oncogenic','knownEffect'])\
                        .rename({'hugoSymbol':'gene'}, axis=1)
        if len(annotations)==0:
            continue
        unique_entries = np.unique(np.hstack((unique_entries, annotations.oncogenic)))
        # Filter to only rows which are oncogenic
        annotations = annotations[oncogenic(annotations.oncogenic)]
        unique_oncogenic = np.unique(np.hstack((unique_oncogenic, annotations.oncogenic)))

        # mutation_df[tumour_sample_platekey] = np.zeros(len(genes), dtype=int)
        unique_genes, unique_counts = np.unique(annotations, return_counts=True)
        match_genes, unique_indices, _ = np.intersect1d(unique_genes, genes, return_indices=True)
        mutation_df[tumour_sample_platekey].loc[match_genes] = unique_counts[unique_indices]

        print(f"{int(100*i/len(sample_list))}%", end="\r")
    print("")

    print("Unique entries: ", unique_entries)
    print("Oncogenic entries: ", unique_oncogenic)

    mutation_df.T.to_csv(f"{dir_output}/OncoKB_somatic_hits.tsv", header=True, index=True, sep="\t")
