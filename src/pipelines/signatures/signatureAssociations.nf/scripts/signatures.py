# Import standard modules
print("Importing packages")
import sys, os, numpy as np, pandas as pd
from multiprocessing import Pool
import re

verbose=True
def print_verbose(message):
    if verbose: print(message)

new_path = '/'.join(os.path.abspath(__file__).split('/')[:-2] + ['python'])
print_verbose(new_path)
sys.path.append(new_path)


if __name__=='__main__':

    # runname to go in filename
    activities_file = sys.argv[1]

    # Import signature activities
    sig_df = pd.read_csv(activities_file, sep="\t")

    sig_df['Samples'] = sig_df['Samples'].map(lambda x: "tumo"+x.split("_")[1]+"_"+x.split("_")[2]+"_norm"+x.split("_")[3]+"_"+x.split("_")[4])
    sig_df.set_index('Samples', inplace=True)

    for key in sig_df.keys():
        sig_df[key] = sig_df[key].astype(float)
        sig_df[key][np.isnan(sig_df[key])] = 0

        sig_df[key] = np.round(sig_df[key]).astype(int)

    sig_df['TSMC'] = np.sum(np.array(sig_df), axis=1)

    try:
        # Add new signatures
        sig_df['SBSapo'] = sig_df.SBS2 + sig_df.SBS13

        # Add together signatures which have been sub-split
        signatures = sig_df.keys()
        split_signatures = np.unique([
            re.search("(SBS[0-9]+?)[a-z]+", sig).group(1) for sig in signatures if re.search("SBS([0-9]+)[a-z]", sig) is not None
        ])
        for split_sig in split_signatures:
            sig_df[split_sig] = np.zeros(len(sig_df))
        for sig in signatures:
            if re.search("SBS([0-9]+)[a-z]", sig) is not None:
                split_sig = re.search("(SBS[0-9]+?)[a-z]+", sig).group(1)
                sig_df[split_sig] += sig_df[sig]
    except AttributeError:
        sig_df['TSMC'] = np.sum(np.array(sig_df), axis=1)
    print(sig_df.head())


    sig_df.to_csv("activities_extended.tsv", index_label=False, header=True, sep="\t")
