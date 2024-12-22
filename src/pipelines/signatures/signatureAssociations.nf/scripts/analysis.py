import sys
import scipy, scipy.stats, re
import pandas as pd, numpy as np

import matplotlib as mpl
mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.family'] = 'STIXGeneral'
import matplotlib.pyplot as plt

plt.rc('axes', labelsize=16)
plt.rc('xtick',labelsize=16)
plt.rc('ytick',labelsize=16)
plt.rc('legend',fontsize=16)




if __name__=='__main__':

    assoc_zinb, assoc_log0 = sys.argv[1:] # , assoc_log

    run_name="twohit_dCRT_binary"
    results_zinb = pd.read_csv(assoc_zinb)
    # results_log = pd.read_csv(assoc_log)
    results_log0 = pd.read_csv(assoc_log0)

    for i, results_table in enumerate([results_zinb]):#, results_log]):


        # QQ plot
        plt.figure(figsize=(8,8))

        subset = (results_table.resample_pvalue>0)&(results_table.model_alt.isin(['negbin', 'glm.binom', 'binomial']))&(np.abs(results_table.target_means_alt)<1000)
        n_test = np.sum(subset)
        k = np.arange(1,n_test+1)

        print(results_table.head())
        print("N_test", n_test)
        print(np.unique(results_table.resample_pvalue))

        pv_expected = np.linspace(0,1-1/n_test,n_test)+1/(2*n_test)
        pv_expected = np.arange(1/n_test,1+1e-10,1/n_test)
        percentiles = scipy.stats.beta.ppf(np.array([0.05,0.5,0.95])[:,None], k[None,:],n_test+1-k[None,:])

        plt.scatter(-np.log10(percentiles[1][::-1]),
                    np.sort(-np.log10(results_table[subset]['wilks_pvalue'])),
                    s=5, label="Wilks", vmax=15, marker='D',
                    c=np.sort(-np.log10(results_table[subset]['wilks_pvalue'])),
                    cmap='viridis_r')
        plt.scatter(-np.log10(percentiles[1][::-1]),
                    np.sort(-np.log10(results_table[subset]['resample_pvalue'])),
                    s=20, label="Resampled", marker='+',
                    c=-np.log10(np.array(results_table[subset]['wilks_pvalue']))[
                        np.argsort(-np.log10(results_table[subset]['resample_pvalue']))
                    ],
                    cmap='viridis_r', vmax=15)

        plt.plot(-np.log10(percentiles[1][::-1]),-np.log10(percentiles[1][::-1]))
        plt.fill_between(-np.log10(percentiles[1][::-1]),
                         -np.log10(percentiles[0][::-1]),
                         -np.log10(percentiles[2][::-1]), alpha=0.6)

        #plt.ylim(0, np.max(-np.log10(results_table[subset]['resample_pvalue']))+1)
        plt.ylim(0, min( np.max(-np.log10(results_table[subset]['resample_pvalue']))*2,
                         np.max(-np.log10(results_table[subset]['wilks_pvalue'])))*1.1)
        plt.gca().set_xlim(left=0)

        plt.legend(loc="upper left")
        plt.xlabel(r"$-\log_{10}(\mathrm{pvalue}_\mathrm{expected})$")
        plt.ylabel(r"$-\log_{10}(\mathrm{pvalue}_\mathrm{observed})$")

        plt.savefig(f"sig_germline_qqplot_{['zinb','log', 'log0'][i]}.png", dpi=200, bbox_inches='tight', facecolor='w', transparent=False)


        # Z-score - Q plot
        plt.figure(figsize=(8,5))

        subset = (results_table.resample_pvalue>0)&\
                 (results_table.model_alt.isin(['negbin', 'glm.binom', 'binomial']))&\
                 (np.abs(results_table.target_means_alt)<1000)
        zscore = results_table.target_means_alt/np.sqrt(results_table.target_covs_alt)
        plt.scatter(zscore[subset],
                     -np.log10(results_table.resample_pvalue[subset]),
                     c=-np.log10(results_table.wilks_pvalue[subset]),
                     s=4, cmap='viridis_r', vmax=15)

        plt.xlim(-np.max(np.abs(zscore[subset])), np.max(np.abs(zscore[subset])))
        plt.gca().set_ylim(bottom=0)

        cbar = plt.colorbar()
        cbar.set_label(r"$-\log_{10}(\mathrm{pvalue}_\mathrm{Wilks})$")

        plt.xlabel(r"$z-\mathrm{score}$")
        plt.ylabel(r"$-\log_{10}(\mathrm{pvalue}_\mathrm{resampled})$")

        plt.savefig(f"sig_germline_zqplot_{['zinb','log', 'log0'][i]}.png", dpi=200, bbox_inches='tight', facecolor='w', transparent=False)


    # Combine logistic and ZINB results
    keys = ['target', 'signature', 'group', 'wilks_pvalue', 'model_alt', 'target_means_alt', 'target_covs_alt']
    results_combined = pd.merge(results_zinb[keys+['resample_pvalue',]],
                                results_log0[keys],
                     on=['target', 'signature', 'group'],
                     how='inner', suffixes=("_zinb", "_log"))

    plt.figure(figsize=(8,5))

    subset = (results_combined.resample_pvalue>0)&\
             (results_combined.wilks_pvalue_log>0)&\
             (results_combined.model_alt_zinb=='negbin')&\
             (np.abs(results_combined.target_means_alt_log)<1000)

    plt.scatter(-np.log10(results_combined[subset].wilks_pvalue_log),
                -np.log10(results_combined[subset].resample_pvalue),
                 c=-np.log10(results_combined[subset].wilks_pvalue_zinb),
                 s=4, cmap='viridis_r', vmax=15)
    plt.gca().set_xlim(left=0)
    plt.gca().set_ylim(bottom=0)

    cbar = plt.colorbar()
    cbar.set_label(r"$-\log_{10}(\mathrm{pvalue}_\mathrm{ZINB Wilks})$")

    plt.xlabel(r"$-\log_{10}(\mathrm{pvalue}_\mathrm{LOG Wilks})$")
    plt.ylabel(r"$-\log_{10}(\mathrm{pvalue}_\mathrm{ZINB resampled})$")

    plt.savefig(f"sig_germline_ppplot_zinb-log.png", dpi=200, bbox_inches='tight', facecolor='w', transparent=False)


    # Negative controls
    vmatch = np.vectorize(lambda x:bool(re.match('.*mock.*', x, re.IGNORECASE)))

    default_colours = plt.rcParams['axes.prop_cycle'].by_key()['color']
    markers = ['.', '+', '^', 'x', '>', '<']

    for j, results_table in enumerate([results_zinb, results_log0]):

        if j==0: pvalue_key = 'resample_pvalue'
        else: pvalue_key = 'wilks_pvalue'

        plt.figure(figsize=(8,8))

        targets = np.unique(results_table.target)
        mock_targets = targets[vmatch(targets)]
        signatures = np.unique(results_table.signature)
        mock_signatures = signatures[vmatch(signatures)]
        mock_set = np.hstack((mock_targets, mock_signatures))
        max_n_test = 0

        for i,mock in enumerate(mock_set):

            subset = (results_table[pvalue_key]>0)&\
                     (results_table.model_alt.isin(['negbin', 'glm.binom', 'binomial']))&\
                     ((results_table['target']==mock)|(results_table['signature']==mock))

            n_test = np.sum(subset)
            max_n_test = max(n_test, max_n_test)
            k = np.arange(1,n_test+1)

            pv_expected = np.linspace(0,1-1/n_test,n_test)+1/(2*n_test)
            pv_expected = np.arange(1/n_test,1+1e-10,1/n_test)
            percentiles = scipy.stats.beta.ppf(np.array([0.05,0.5,0.95])[:,None], k[None,:],n_test+1-k[None,:])

            plt.scatter(-np.log10(percentiles[1][::-1]),
                        np.sort(-np.log10(results_table[subset]['wilks_pvalue'])),
                        s=10, label=f"{mock} Wilks",
                        c=default_colours[3], marker=markers[i%len(markers)])

            if j==0:
                plt.scatter(-np.log10(percentiles[1][::-1]),
                            np.sort(-np.log10(results_table[subset]['resample_pvalue'])),
                            s=10, label=f"{mock} dCRT",
                            c=default_colours[0], marker=markers[i%len(markers)])

        k = np.arange(1,max_n_test+1)
        percentiles = scipy.stats.beta.ppf(np.array([0.05,0.5,0.95])[:,None],
                                           k[None,:],
                                           max_n_test+1-k[None,:])
        plt.plot(-np.log10(percentiles[1][::-1]),-np.log10(percentiles[1][::-1]), c='k')
        plt.fill_between(-np.log10(percentiles[1][::-1]),
                         -np.log10(percentiles[0][::-1]),
                         -np.log10(percentiles[2][::-1]), alpha=0.3, color='k')

        plt.gca().set_xlim(left=0)
        plt.gca().set_ylim(bottom=0)

        plt.legend(loc="upper left")
        plt.xlabel(r"$-\log_{10}(\mathrm{pvalue}_\mathrm{expected})$")
        plt.ylabel(r"$-\log_{10}(\mathrm{pvalue}_\mathrm{observed})$")

        plt.savefig(f"sig_germline_mock-qqplot_{['zinb', 'log0'][j]}.png", dpi=200, bbox_inches='tight', facecolor='w', transparent=False)
