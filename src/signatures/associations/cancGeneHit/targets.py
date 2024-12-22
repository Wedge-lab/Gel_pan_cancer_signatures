# Import standard modules
print("Importing packages")
import sys, os, numpy as np, pandas as pd
import scipy.stats
from dotenv import load_dotenv

load_dotenv()
DATA_DIR = os.getenv('DATA_DIR')

if __name__=="__main__":

    samples_file, targets_file, n = sys.argv[1:]
    n=int(n)

    # Sample list
    sample_df = pd.read_csv(samples_file, sep='\t')
    ids = sample_df.sample_id.map(lambda x: [x.split("_")[0], "_".join(x.split("_")[1:3]), "_".join(x.split("_")[3:])])
    sample_df = pd.concat((sample_df,
                           pd.DataFrame(ids.tolist(), index= ids.index,
                                        columns=["participant_id", "tumour_sample_platekey", "germline_sample_platekey"])),
                          axis=1)

    # Germline
    germline_df = pd.read_csv(f"{DATA_DIR}/aggv2-germline_cancer-gene_hits.tsv",#"germline_30/output/aggv2-germline_cancer-gene_hits.tsv",
                             sep="\t").set_index("gene_id").T

    # OncoKB
    somatic_df = pd.read_csv(f"{DATA_DIR}/OncoKB_somatic_cancer-gene_hits.tsv",
                               sep="\t", index_col=0)

    # LOH
    loh_df = pd.read_csv(f"{DATA_DIR}/lohfrac_gene_sample_matrix.tsv",
                             sep="\t").set_index("tumour_sample_platekey")
    # Add in platekeys which failed battenberg and assume no LoH
    failed_battenberg_platekeys = np.setxor1d(np.intersect1d(loh_df.index, somatic_df.index), somatic_df.index)
    loh_df = pd.concat((loh_df,
                        pd.DataFrame(np.zeros((len(failed_battenberg_platekeys), len(loh_df.keys()))),
                                     index=failed_battenberg_platekeys,
                                     columns=loh_df.keys())))

    # Only tumour suppressor genes
    cancer_genes = pd.read_csv(f"{DATA_DIR}/cancer_gene_census.csv")
    genes = cancer_genes['Gene Symbol'][cancer_genes['Role in Cancer'].map(lambda x: "TSG" in str(x))]

    # Get binary two hit model data
    tumour_sample_platekeys = list(sample_df.tumour_sample_platekey)
    germline_sample_platekeys = list(sample_df.germline_sample_platekey)
    germline_df = germline_df.loc[germline_sample_platekeys][genes].set_index(sample_df.sample_id)
    loh_df = loh_df.loc[tumour_sample_platekeys][genes].set_index(sample_df.sample_id)
    somatic_df = somatic_df.loc[tumour_sample_platekeys][genes].set_index(sample_df.sample_id)
    twohit_df = pd.DataFrame(np.array(germline_df > 0).astype(int) + \
                             np.array(loh_df>0.5).astype(int) + \
                             np.array(somatic_df),
                             index=sample_df.sample_id,
                             columns=genes)
    print(np.unique(twohit_df, return_counts=True))
    twohit_df.index.name = 'sample_id'
    if n==1:
        # Value=1 when both alleles have been knocked out
        twohit_df = (twohit_df>1).astype(int)
    elif n==2:
        # Value=2 when two of germline/somatic/loh have been hit
        twohit_df[twohit_df>2] = 2
    else:
        raise ValueError(f"n={n} but n must be 1 or 2")

    # Add 0,1,2 columns for germline, somatic, loh
    germline_df.rename({gene:gene+"_0" for gene in genes}, axis=1, inplace=True)
    germline_df[germline_df>1] = 1
    somatic_df.rename({gene:gene+"_1" for gene in genes}, axis=1, inplace=True)
    # somatic_df[somatic_df>2] = 2
    loh_df.rename({gene:gene+"_2" for gene in genes}, axis=1, inplace=True)
    loh_df = (loh_df==1).astype(int)
    twohit_df = pd.concat([twohit_df, germline_df, loh_df, somatic_df], axis=1)

    # DNA types from database
    DNA_repair = pd.merge(pd.read_csv(f"{DATA_DIR}/human-dna-repair-genes.tsv", sep="\t"),
                              pd.DataFrame(index=np.array(genes)),
                              left_on="Gene", right_index=True, how='inner')
    DNA_repair_types = {'Chromatin Structure and Modification':'CSM',
                 'DNA polymerases ':'POL',
                 'Direct reversal of damage':'DRD',
                 'Editing and processing nucleases':'EPN',
                 'Fanconi anemia':'FA',
                 'Genes defective in diseases associated with sensitivity to DNA damaging agents':'DDA',
                 'Homologous recombination':'HR',
                 'Other conserved DNA damage response genes':'OTHER'}
    DNA_repair.replace({'Type':DNA_repair_types}, inplace=True)
    for type in np.unique(DNA_repair['Type']):
        for i in range(3):
            type_hit = np.zeros(len(twohit_df), dtype=int)
            type_hit[np.sum(twohit_df[[f"{gene}_{i}" for gene in DNA_repair.Gene[DNA_repair.Type==type]]]>0, axis=1)>0] = 1
            type_hit[np.sum(twohit_df[[f"{gene}_{i}" for gene in DNA_repair.Gene[DNA_repair.Type==type]]]>1, axis=1)>0] = 2
            twohit_df[f"{type}_{i}"] = type_hit
        total_hits = twohit_df[f'{type}_0'] + twohit_df[f'{type}_1'] + twohit_df[f'{type}_2']
        twohit_df[f"{type}"] = np.where(total_hits>n, 2, total_hits)

    # Generate mock targets - for 1%,5%,20%,50% hit rate
    for hit_rate in [5,15,50]:
        twohit_df[f'MOCK{hit_rate:02d}_0'] = scipy.stats.binom.rvs(1,hit_rate/100,size=len(twohit_df))
        twohit_df[f'MOCK{hit_rate:02d}_1'] = scipy.stats.poisson.rvs(hit_rate/100,size=len(twohit_df))
        twohit_df[f'MOCK{hit_rate:02d}_2'] = scipy.stats.binom.rvs(1,hit_rate/100,size=len(twohit_df))
        total_hits = twohit_df[f'MOCK{hit_rate:02d}_0'] + twohit_df[f'MOCK{hit_rate:02d}_1'] + twohit_df[f'MOCK{hit_rate:02d}_2']
        twohit_df[f'MOCK{hit_rate:02d}'] = np.where(total_hits>n, 2, total_hits)

    # Generate mock targets - perturbation of genes within cohorts
    for gene in ["NTHL1","BRCA2","MGMT","POL","FA","HR","MMR"]:
        indices = np.arange(len(twohit_df))
        MOCKgene = np.zeros((len(twohit_df),4), dtype=int)
        for group in np.unique(sample_df.group):
            sample = np.random.choice(indices[sample_df.group==group], size=np.sum(sample_df.group==group), replace=False)
            MOCKgene[sample_df.group==group,0] = twohit_df[gene][sample]
            MOCKgene[sample_df.group==group,1] = twohit_df[gene+"_0"][sample]
            MOCKgene[sample_df.group==group,2] = twohit_df[gene+"_1"][sample]
            MOCKgene[sample_df.group==group,3] = twohit_df[gene+"_2"][sample]
        twohit_df[f'MOCK{gene}'] = MOCKgene[:,0]
        twohit_df[f'MOCK{gene}_0'] = MOCKgene[:,1]
        twohit_df[f'MOCK{gene}_1'] = MOCKgene[:,2]
        twohit_df[f'MOCK{gene}_2'] = MOCKgene[:,3]

    # Save targets file
    twohit_df.to_csv(targets_file, sep="\t")
