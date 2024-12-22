# Import standard modules
print("Importing packages")
import sys, pandas as pd

# signature       group   target  type

if __name__=="__main__":

    samples_file, samples_cohort_file, targets_cohort_file = sys.argv[1:]

    # Load samples file
    sample_df = pd.read_csv(samples_file, sep="\t")

    # Get group tissues
    sample_df['tumour_tissue'] = sample_df.group.map(lambda x: x.split("-")[0])
    group_sizes = sample_df.groupby(['tumour_tissue', 'group']).size()
    # Generate new samples file with groups as tumour tissues
    sample_df.drop('group', axis=1).rename({'tumour_tissue':'group'}, axis=1).to_csv(samples_cohort_file, sep="\t", index=False)

    # Number of groups per tissue
    tissue_group_counts = group_sizes.reset_index().groupby('tumour_tissue').size()
    bigroup = tissue_group_counts[tissue_group_counts==2].index
    polygroup = tissue_group_counts[tissue_group_counts>2].index

    # Test groups
    test_groups = list(group_sizes.reset_index().set_index('tumour_tissue').loc[polygroup].group)
    for i in range(len(bigroup)):
        test_groups += [group_sizes.reset_index().set_index('tumour_tissue').loc[bigroup[i]].group.iloc[0],]

    # Generate variables
    for group in test_groups:
        sample_df[group] = (sample_df.group==group).astype(int)

    # Save sample dataframe
    sample_df[['sample_id']+test_groups].to_csv(targets_cohort_file, sep="\t", index=False)
