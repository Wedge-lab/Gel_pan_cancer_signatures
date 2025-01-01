import os
import sys

import numpy as np
import pandas as pd

from signatures.config import load_environment

load_environment()
AGGV2_SAMPLE_LIST = os.getenv("AGGV2_SAMPLE_LIST")

if __name__ == "__main__":
    dir_output, sample_list_file = sys.argv[1:]

    # Get all samples from previous projects
    cancer_analysis_tables = [
        f"/re_gecip/shared_allGeCIPs/labkey_tables/{version}/cancer_analysis.tsv"
        for version in ["V8", "V11/V11_reheadered", "v14/v14_reheadered"]
    ]
    germline_platekey_list = np.array([])
    for table in cancer_analysis_tables:
        germline_platekey_list = np.union1d(
            germline_platekey_list,
            np.array(
                pd.read_csv(table, sep="\t", usecols=["germline_sample_platekey"])
            ),
        )

    # Crossmatch with aggV2
    aggv2_sample_list = np.array(
        pd.read_csv(AGGV2_SAMPLE_LIST, sep="\t", header=None)[0]
    )

    germline_platekey_list = np.intersect1d(germline_platekey_list, aggv2_sample_list)
    print(f"{len(germline_platekey_list)} germline samples")
    pd.DataFrame(germline_platekey_list, columns=["germline_sample_platekey"]).to_csv(
        f"{sample_list_file}", header=False, index=False
    )
