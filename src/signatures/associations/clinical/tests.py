# Import standard modules
import sys, pandas as pd

if __name__=="__main__":

    samples_file, samples_cohort_file, signatures_file, targets_file, targets_cohort_file, \
    tests_file, tests_cohort_file, tests_file_binary, tests_cohort_file_binary = sys.argv[1:]
    min_sig_count=10
    min_tgt_count=5

    # Load in signatures, targets and groups
    signatures = pd.read_csv(signatures_file, sep="\t").set_index('sample_id')
    # Convert to binary
    signatures = (signatures>0).astype(int).T

    for i in range(2):
        targets = pd.read_csv([targets_file, targets_cohort_file][i], sep="\t").set_index('sample_id')
        groups = pd.read_csv([samples_file, samples_cohort_file][i], usecols=['sample_id', 'group'], sep="\t")

        targets_0 = (targets==0).astype(int).T
        targets_1 = (targets==1).astype(int).T

        # Number of non-zero activities for each signature-group pair
        signature_groups = pd.merge(signatures.T, groups.set_index('sample_id'),
                                    left_index=True, right_index=True, how='inner')\
                            .groupby("group")

        signature_counts = pd.merge(signature_groups.sum().reset_index()\
                                     .melt(id_vars=['group'], value_name='count_sig', var_name='signature'),
                                  signature_groups.count().reset_index()\
                                     .melt(id_vars=['group'], value_name='total', var_name='signature'))

        # Number of mutated samples in each gene-group pair
        target_group_0 = pd.merge(targets_0.T, groups.set_index('sample_id'),
                                    left_index=True, right_index=True, how='inner')\
                            .groupby("group").sum().reset_index()\
                                         .melt(id_vars=['group'], value_name='count_tgt_0', var_name='target')
        target_group_1 = pd.merge(targets_1.T, groups.set_index('sample_id'),
                                    left_index=True, right_index=True, how='inner')\
                            .groupby("group").sum().reset_index()\
                                         .melt(id_vars=['group'], value_name='count_tgt_1', var_name='target')

        # Only perform tests with at least min_sig_count non-zero signatures and min_tgt_count mutated and non-mutated samples
        set_counts = pd.merge(signature_counts,
                     pd.merge(target_group_0, target_group_1,
                              on=('group','target'), how='inner'),
                              on='group', how='inner')
        set_counts = set_counts[(set_counts['count_sig']>min_sig_count)&\
                                (set_counts['count_tgt_0']>min_tgt_count)&\
                                (set_counts['count_tgt_1']>min_tgt_count)]

        # Target is binomial with n total possible hits
        set_counts['type'] = "binomial"

        # Save test file
        print(f"{len(set_counts)} tests")
        header = ['signature', 'group', 'target', 'type']
        set_counts[header].sort_values(header).to_csv([tests_file, tests_cohort_file][i], sep="\t", index=False)

        # Test file for non-zero logistic regression
        set_counts = set_counts[(set_counts['total']-set_counts['count_sig']>min_sig_count)]
        print(f"{len(set_counts)} binary tests")
        header = ['signature', 'group', 'target']
        set_counts[header].sort_values(header).to_csv([tests_file_binary, tests_cohort_file_binary][i], sep="\t", index=False)
