"Combine signatures between cohorts"
import sys, os
import pandas as pd
import SigProfilerAssignment.Analyzer as Analyze

def decompose_and_analyze(args):

    output, samples_file, signatures_file, sig_type, reference = args

    if not os.path.exists(output): os.makedirs(output)

    print(output.split("/")[-2:])

    Analyze.decompose_fit(samples_file,
                                   output,
                                   signatures=signatures_file,
                                   signature_database=reference,
                                   genome_build="GRCh38",
                                   verbose=False,
                                   new_signature_thresh_hold=0.8,
                                   make_plots=False,
                                   collapse_to_SBS96=sig_type=='SBS288'
                                   )

if __name__=='__main__':

    output_dir, stats_file, sig_type, max_sigs, reference = sys.argv[1:]
    max_sigs = int(max_sigs)

    # Load file locations
    stats_df = pd.read_csv(stats_file, sep="\t")
    stats_df = stats_df[stats_df.signatures<=stats_df.max_sigs]

    # Set output directories
    stats_df['output'] = output_dir+"/"+stats_df.cohort+"/"+stats_df.signatures.astype(str)
    stats_df['reference'] = reference
    stats_df['sig_type'] = sig_type

    stats_list = list(stats_df[['output', 'samples_file', 'signatures_file', 'sig_type', 'reference']].itertuples(index=False, name=None))

    for i, args in enumerate(stats_list):
        decompose_and_analyze(args)