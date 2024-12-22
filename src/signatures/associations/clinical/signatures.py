# Import standard modules
print("Importing packages")
import sys, numpy as np, pandas as pd
from signatures.plotting.combinedSignatures import sig_dirs

if __name__=="__main__":

    samples_file, activities_file = sys.argv[1:]

    # Load samples file
    samples_df = pd.read_csv(samples_file, sep="\t")
    print("Samples: ", len(samples_df))

    print("Load signature activity files")
    activities_df = pd.DataFrame(index=samples_df['sample_id'])
    signatures_files = {key: f"{sig_dir}/Combined_Solution_Activities.tsv" for key,sig_dir in sig_dirs.items()}
    for sig_set, sig_file in signatures_files.items():
        sig_df = pd.read_csv(sig_file, sep="\t").set_index("Samples")
        sig_df[f"T{sig_set}"] = np.sum(sig_df, axis=1)
        activities_df = pd.merge(activities_df, sig_df,
                                 how="left", left_index=True, right_index=True).fillna(0).astype(int)
    activities_df.index.rename("sample_id", inplace=True)
    print(activities_df.head())
    print("Activities: ", activities_df.shape)

    # Get activities for sample set
    activities_df = pd.merge(samples_df[["sample_id"]], activities_df,
                            left_on="sample_id", right_index=True,
                            how="inner").set_index("sample_id")
    print("Combined: ", activities_df.shape)

    # Save activities file
    activities_df.to_csv(activities_file, sep="\t")
