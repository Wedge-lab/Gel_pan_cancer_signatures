# Import standard modules
import sys, numpy as np, pandas as pd
from signatures.associations.sampleCuration.pan_cancer_sample import getSamples

if __name__=='__main__':

    # filename to save covariates into
    filename_output = sys.argv[1]

    # GET DATA - include hypermutated samples
    print("Load samples")
    sample_df = getSamples(PC=True, signatures=True, primary=True, ethnicity=False, hyper=True, ER=False, MSI=False)

    # Remove samples where sex doesn't match organ
    sample_df = sample_df[ ( ~(sample_df.tumour_tissue.map(lambda x: x in ["Breast","Ovary","Uterus"])&(sample_df.is_female==-1)) ) &\
                           ( ~(sample_df.tumour_tissue.map(lambda x: x in ["Testis", "Prostate"])&(sample_df.is_female==1)) ) ]
    print("Removed male samples for Breast,Ovary,Uterus and female for Testis,Prostate")
    print(f"Number of samples: {len(sample_df)}")

    # Confounding variables
    variables = ['log_age', 'is_female'] + [f"pc{i}" for i in range(1,4)] # logit_purity',

    # Construct X matrix
    print("Construct X")
    X = pd.DataFrame({variable:sample_df[variable] for variable in variables+['tumour_tissue', 'tumour_group',
                                                                              'participant_id', 'tumour_sample_platekey', 'germline_sample_platekey']})

    # Get non-nan values
    [print(key, np.sum(~np.isnan(X[key]))) for key in X.columns if X[key].dtype in [float, int]]
    cuts = np.prod(np.array([~np.isnan(X[key]) for key in X.columns if X[key].dtype in [float, int]]),
                   axis=0).astype(bool)
    print(f"Non-nan: {np.sum(cuts)}")

    # Apply cuts
    X = X[cuts]

    # Save
    X['sample_id'] = X.participant_id.astype(str)+"_"+X.tumour_sample_platekey+"_"+X.germline_sample_platekey
    X.rename(columns={"tumour_tissue":"group"}, inplace=True)
    print(X.columns)
    X[["sample_id", "group"] + variables].to_csv(filename_output, index=False, sep="\t")
