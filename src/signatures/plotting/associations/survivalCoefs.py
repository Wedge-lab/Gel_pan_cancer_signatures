import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
import scipy.stats
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test

from signatures.plotting.combinedSignatures import loadSignatures, signatureRenamer

mpl.rcParams["mathtext.fontset"] = "stix"
mpl.rcParams["font.family"] = "STIXGeneral"
plt.rc("axes", labelsize=16)
plt.rc("xtick", labelsize=16)
plt.rc("ytick", labelsize=16)
plt.rc("legend", fontsize=16)


def survival(sample_df, target, covariates, transformation="lognorm"):
    target_arr = np.array(sample_df[target].copy())

    median = np.median(target_arr, axis=0)
    target_arr[target_arr <= median] = 0
    target_arr[target_arr > median] = 1
    subset_split = sample_df[target] < np.mean(sample_df[target])

    logrank = logrank_test(
        np.array(sample_df[target_arr == 0].survival_time),
        np.array(sample_df[target_arr == 1].survival_time),
        event_observed_A=np.array(sample_df[target_arr == 0].vital_status),
        event_observed_B=np.array(sample_df[target_arr == 1].vital_status),
    )

    ### KM curve split by target
    plt.figure()
    kmf = KaplanMeierFitter()
    kmf.fit(
        durations=np.array(sample_df[subset_split].survival_time),
        event_observed=np.array(sample_df[subset_split].vital_status),
    )
    kmf.plot_survival_function(label=f"No {target}")
    kmf.fit(
        durations=np.array(sample_df[~subset_split].survival_time),
        event_observed=np.array(sample_df[~subset_split].vital_status),
    )
    kmf.plot_survival_function(label=f"{target}")
    plt.xlabel("Time since sampling")
    plt.ylim(0, 1)
    plt.xlim(0, np.max(sample_df.survival_time))

    try:
        cph_df = sample_df.copy()
        if transformation == "lognorm":
            cph_df[target] = np.log10(cph_df[sig] + 1)
            cph_df[target] = (cph_df[target] - np.mean(cph_df[target])) / np.std(
                cph_df[target]
            )

        # Cox proportional hazard
        cph_alt = CoxPHFitter()
        alt_df = cph_df[
            ["survival_time", "vital_status"]
            + covariates
            + [
                target,
            ]
        ]
        cph_alt.fit(alt_df, duration_col="survival_time", event_col="vital_status")
    except:
        cph_alt = None
        print(f"CPH failed for {target}")

    return cph_alt, logrank


def result_string(sig, cph_alt):
    """Get string of results to print out"""
    mu = cph_alt.params_.loc[sig]
    z = mu / cph_alt.standard_errors_.loc[sig]
    p_wald = scipy.stats.chi2.sf(z**2, df=1)
    return f"Wald={p_wald:.2e}, z={z:.2e}, mu={mu:.2e} CI=[{cph_alt.confidence_intervals_['95% lower-bound'].loc[sig]:.2e}, {cph_alt.confidence_intervals_['95% upper-bound'].loc[sig]:.2e}]"


if __name__ == "__main__":
    group = "Breast-DuctalCA"
    print(f"Using group {group}")

    # Get renaming dictionary
    combined_acts, combined_sigs, relabel_map = loadSignatures()
    rename_dict = signatureRenamer(relabel_map)

    # Import sample data
    sample_df = pd.read_csv(
        "/re_gecip/shared_allGeCIPs/pancancer_signatures/results/survival/survival_log_grade-4_unique/input/samples.tsv",
        sep="\t",
        index_col=0,
    )
    sample_df.drop_duplicates(inplace=True)
    sample_df = sample_df[~np.isnan(sample_df.grade)]
    covariates = [
        key
        for key in sample_df.keys()
        if not key in ["group", "survival_time", "vital_status"]
    ]

    # Add tumour type
    tumour_type = pd.read_csv(
        "/re_gecip/shared_allGeCIPs/pancancer_signatures/results/sample_lists_incl_SEGs/sample_list_2021_06_29.tsv",
        sep="\t",
        usecols=[
            "participant_id",
            "tumour_sample_platekey",
            "germline_sample_platekey",
            "tumour_type",
        ],
    )
    tumour_type["sample_id"] = (
        tumour_type.participant_id.astype(str)
        + "_"
        + tumour_type["tumour_sample_platekey"]
        + "_"
        + tumour_type["germline_sample_platekey"]
    )
    sample_df = pd.merge(sample_df, tumour_type, on="sample_id", how="inner").set_index(
        "sample_id"
    )

    # Import tgtnature data
    target_df = pd.read_csv(
        "/re_gecip/shared_allGeCIPs/pancancer_signatures/results/survival/survival_log_grade-4_unique/input/signatures.tsv",
        sep="\t",
        index_col=0,
    )
    target_df.rename(rename_dict, axis=1, inplace=True)
    targets = target_df.keys()

    # Merge targets with samples
    subset = np.array(sample_df.group == group)
    subset_df = pd.merge(
        sample_df[subset], target_df, how="inner", left_index=True, right_on="sample_id"
    ).sort_values("survival_time")

    # Add clinical data
    clinical_data = pd.read_csv(
        "/re_gecip/shared_allGeCIPs/pancancer_signatures/results/associations/clinicalSigs/clinical_tumour-group_InDelQC_nb/input/targets.tsv",
        sep="\t",
        index_col=0,
    )
    clinical_data = clinical_data.loc[subset_df.index]
    subset_df = pd.merge(
        subset_df,
        clinical_data[["pr_status", "er_status", "her2_status"]],
        left_index=True,
        right_index=True,
        how="inner",
    )

    subset_df["HRD"] = (
        (
            (subset_df["SBS3"] > 0).astype(int)
            + (subset_df["ID6"] > 0).astype(int)
            + (subset_df["CN17"] > 0).astype(int)
        )
        > 1
    ).astype(int)

    print("Association between signatures and survival")
    for sig in ["HRD", "SBS3", "DBS6", "ID6", "CN17", "SV3"]:
        cph_alt, logrank = survival(subset_df, sig, covariates)
        print(f"{sig}: {result_string(sig, cph_alt)}")

    print("Controlling for ER status")
    subset = ~pd.isnull(clinical_data.er_status)
    for sig in ["HRD", "SBS3", "DBS6", "ID6", "CN17", "SV3"]:
        cph_alt, logrank = survival(subset_df[subset], sig, covariates + ["er_status"])
        print(f"{sig}: {result_string(sig, cph_alt)}")

    print("Survival vs clinical status")
    for tgt in ["er_status", "pr_status", "her2_status"]:
        subset = ~pd.isnull(clinical_data[tgt])
        cph_alt, logrank = survival(
            subset_df[subset], tgt, covariates, transformation="none"
        )
        print(f"{tgt}: {result_string(tgt, cph_alt)}")

    group = "CNS-GBM-IDHwt"
    print(f"Using group {group}")

    # Merge targets with samples
    subset = np.array(sample_df.group == group)
    subset_df = pd.merge(
        sample_df[subset], target_df, how="inner", left_index=True, right_on="sample_id"
    ).sort_values("survival_time")

    print("Association between signatures and survival")
    for sig in ["ID8"]:
        cph_alt, logrank = survival(subset_df, sig, covariates)
        print(f"{sig}: {result_string(sig, cph_alt)}")

    print("Only primary tumours")
    subset = subset_df.tumour_type == "PRIMARY"
    print(subset_df.tumour_type.value_counts())
    for sig in ["ID8"]:
        cph_alt, logrank = survival(subset_df[subset], sig, covariates)
        print(f"{sig}: {result_string(sig, cph_alt)}")
