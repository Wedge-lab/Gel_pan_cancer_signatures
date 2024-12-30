import numpy as np
import re
import pandas as pd

def orderSignatures(signatures):

    """orderSignatures - sort signatures into a sensible ordering
    1) mutation type - SBS, DBS, ID
    2) signature number - 1,2,3,...
    3) novel signature cohort - Bladder, Breast...
    4) signature letter - 7a,7b,7c
    """

    mut_type  = np.array([re.search("([A-Za-z]+)([0-9]+)([a-zA-Z]*)", sig).group(1) for sig in signatures])
    mut_index = np.array([re.search("([A-Za-z]+)([0-9]+)([a-zA-Z]*)", sig).group(2) for sig in signatures]).astype(int)
    mut_exten = np.array([re.search("([A-Za-z]+)([0-9]+)([a-zA-Z]*)", sig).group(3) for sig in signatures])
    mut_tissue= np.array([re.search("([A-Za-z]+)([0-9]+)([a-zA-Z]*)-([a-zA-Z_]+)", sig).group(4) \
                          if '-' in sig else 'Z' for sig in signatures])

    order_exten = np.zeros(len(mut_exten))
    order_exten[np.argsort(mut_exten)] = np.arange(len(mut_exten))

    order_tissue = np.zeros(len(mut_tissue))
    order_tissue = np.arange(len(np.unique(mut_tissue)))[
        np.unique(mut_tissue, return_inverse=True)[1]
    ]

    print(mut_type)

    order = np.argsort(pd.Series(mut_type).replace(['SBS','DBS','ID','CN','CNV','RefSigR','RS','SV'],[1,2,3,4,5,6,7,8]).astype(int) \
                       + mut_index/1e3 \
                       + order_tissue/1e5 \
                       + order_exten/1e7)
    return order