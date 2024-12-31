import sys

import pandas as pd

if __name__ == "__main__":
    sample_file, groups_file = sys.argv[1:]

    samples_df = pd.read_csv(sample_file, usecols=["group", "vital_status"], sep="\t")

    # Count number of deaths
    group_df = samples_df.groupby("group").sum().copy()

    # Keep groups with at least 5 deaths recorded
    group_df = group_df[group_df.vital_status > 5].reset_index().copy()

    # Save group list
    group_df[["group"]].to_csv(groups_file, sep="\t", index=False)
