import os
import sys

import numpy as np
import pandas as pd

if __name__ == "__main__":
    dir_output, sample_list_file = sys.argv[1:]

    # Get all samples from previous projects
    cancer_analysis_tables = [
        f"/re_gecip/shared_allGeCIPs/labkey_tables/{version}/cancer_analysis.tsv"
        for version in ["V8", "V11/V11_reheadered", "v14/v14_reheadered"]
    ]
    tumour_platekey_list = np.array([])
    for table in cancer_analysis_tables:
        tumour_platekey_list = np.union1d(
            tumour_platekey_list,
            np.array(pd.read_csv(table, sep="\t", usecols=["tumour_sample_platekey"])),
        )

    print(f"{len(tumour_platekey_list)} tumour_samples")
    pd.DataFrame(tumour_platekey_list, columns=["tumour_sample_platekey"]).to_csv(
        f"{sample_list_file}", header=False, index=False
    )
