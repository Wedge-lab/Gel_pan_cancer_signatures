"Find best solution for each cohort"
import sys, os
import pandas as pd, numpy as np
import inspect
import matplotlib as mpl
import matplotlib.pyplot as plt
from SigProfilerExtractor import subroutines as sub

mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.family'] = 'STIXGeneral'
plt.rc('axes', labelsize=16)
plt.rc('xtick',labelsize=16)
plt.rc('ytick',labelsize=16)
plt.rc('legend',fontsize=16)


print(f"SigProfilerAssignment file location: {inspect.getfile(sub.calculate_similarities)}")

def estimate_solution(stats_df,
                        genomes="Samples.txt",
                        output="results",
                        title="Selection_Plot",
                        stability=0.8,
                        min_stability=0.2,
                        combined_stability=1.0,
                        statistics=True,
                        select=None,
                        sequence="genome",
                        max_signatures=50,
                         mtype="SBS96"):

    """estimate_solution
    SigProfilerExtractor.estimate_best_solution.estimate_solution hacked to allow a max signature cutoff
    """

    # Create output directory if it doesn't exist
    try:
        if not os.path.exists(output): os.makedirs(output)
    except:
        print (f"The {output} folder could not be created")

    # Solution statistics file
    stats_df.sort_values("signatures", inplace=True)
    stats_df = stats_df[:max_signatures]

    # Mutation count file
    genomes=pd.read_csv(genomes, sep="\t", index_col=0)

    # Get signature and activity matrices
    all_similarities_list=[]
    for i, row in stats_df.iterrows():
        w=np.array(pd.read_csv(row.signatures_file, sep="\t", index_col=0))
        h=np.array(pd.read_csv(row.activities_file, sep="\t", index_col=0)).T

        all_similarities, _ = sub.calculate_similarities(np.array(genomes), np.dot(w,h), genomes.columns)
        all_similarities_list.append(all_similarities)

    # Rename columns for input to stabVsRError
    columns=['signatures', 'min_stability', 'matrix_frobenius', 'avg_stability']
    stats_df = stats_df[columns].reset_index(drop=True)\
                .rename(dict(zip(columns,
                                 ["Total Signatures",  "Stability",  "Matrix Frobenius%",  "avgStability"])),
                       axis=1)

    print("Run stabVsRError")
    solution, all_stats= sub.stabVsRError(stats_df, output, title,
                                            all_similarities_list,
                                            input_type="dataframe",
                                            stability=stability, min_stability=min_stability, combined_stability=combined_stability,
                                            mtype=mtype,
                                            statistics=statistics,
                                            select=select,
                                            sequence=sequence)

    all_stats.insert(1, 'Stability (Avg Silhouette)', stats_df["avgStability"])
    all_stats=all_stats.set_index(["Signatures"])
    all_stats.to_csv(output+"/All_solutions_stat.csv", sep = ",")

    return solution



if __name__=='__main__':

    cohort_stats_file, max_sigs, min_stability, output_dir, sig_type = sys.argv[1:]
    min_stability=float(min_stability)
    print(output_dir, min_stability)

    if not os.path.exists(f"{output_dir}/figs"): os.mkdir(f"{output_dir}/figs")

    # Load Cohorts
    all_cohort_stats = pd.read_csv(cohort_stats_file, sep="\t")
    cohorts = np.unique(all_cohort_stats.cohort)

    # Initialise dataframes
    best_solutions_df = pd.DataFrame()
    all_cohort_stats_aic = pd.DataFrame()

    # Iterate through cohorts to find best solutions
    for cohort in cohorts:
        print(f"Best solutions for {cohort}")
        max_sigs = all_cohort_stats[all_cohort_stats.cohort==cohort].max_sigs.iloc[0]
        cohort_stats = all_cohort_stats[all_cohort_stats.cohort==cohort].sort_values('signatures')[:max_sigs]

        # Estimate best solution for signatures
        best_solution = estimate_solution(stats_df=cohort_stats, #f"{cohort_dir}/All_solutions_stat.csv",
                                            genomes=cohort_stats.samples_file.iloc[0], #f"{cohort_dir}/Samples.txt",
                                            output=f"{output_dir}/{cohort}",
                                            title="Selection_Plot",
                                            stability=0.8,
                                            min_stability=0.2,
                                            combined_stability=min_stability,
                                            statistics=True,
                                            select=None,
                                            sequence="genome",
                                            mtype=sig_type,
                                            max_signatures=max_sigs)

        # Get Sample mutations
        mut_df=pd.read_csv(cohort_stats['samples_file'].iloc[0], sep="\t", index_col=0).T
        if len(mut_df.T)==288: mut_df = mut_df.T.groupby(mut_df.T.index.str[2:9]).sum().T
        m = len(mut_df.T)

        # Estimate BIC and AIC for each k signatures tried
        bic = np.zeros(max_sigs)
        aic = np.zeros(max_sigs)
        decomposed_sigs_file_exists = True
        for i in range(1,max_sigs+1):
            try:
                sig_df = pd.read_csv(cohort_stats.set_index('signatures')['signatures_decomposed_file'].loc[i],
                                        sep="\t", index_col=0)
                act_df = pd.read_csv(cohort_stats.set_index('signatures')['activities_decomposed_file'].loc[i],
                                        sep="\t", index_col=0)
            except:
                decomposed_sigs_file_exists = False
                sig_df = pd.read_csv(cohort_stats.set_index('signatures')['signatures_file'].loc[i],
                                        sep="\t", index_col=0)
                act_df = pd.read_csv(cohort_stats.set_index('signatures')['activities_file'].loc[i],
                                        sep="\t", index_col=0)
            mut_mean = (act_df @ sig_df.T).loc[mut_df.index][mut_df.keys()]
            mut_mean[mut_mean<1e-5] = 1e-5

            bic[i-1] = (m+len(mut_df))*i*np.log(m*len(mut_df)) - 2*np.sum(np.array(mut_df*np.log(mut_mean) - mut_mean))
            aic[i-1] = 2*(m+len(mut_df))*i - 2*np.sum(np.array(mut_df*np.log(mut_mean) - mut_mean))
            if np.isnan(aic[i-1]):
                print(m, len(mut_df), np.sum(np.log(np.array(mut_mean))), np.sum(np.array(mut_df*np.log(mut_mean))))
                print(mut_df.head())
                print(mut_mean.head())
        cohort_stats['aic'] = aic
        cohort_stats['bic'] = bic
        if not decomposed_sigs_file_exists:
            print("Decomposed sigs file does not exist! Running AIC on deNovo activities")
            print(cohort_stats.set_index('signatures')['signatures_decomposed_file'].loc[i])

        # Plot Stability and Best solution
        fig,ax=plt.subplots(1,1,figsize=(16,10))
        default_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

        # Plot stability
        ln1=plt.plot(np.arange(1,max_sigs+1), cohort_stats["avg_stability"],
                     "s-.", c=default_colors[0], label="Stability")
        plt.plot([0,max_sigs+1], [0.8,0.8], '--k', alpha=0.5)
        plt.ylabel("Stability")

        # Plot mean cosine similarity
        ln_cos=plt.plot(np.arange(1,max_sigs+1), cohort_stats["Mean Cosine Distance"],
                    "o:", c=default_colors[3], label="Mean Cosine Dist")
        plt.fill_betweenx(ax.get_ylim(), [best_solution-0.25,]*2, [best_solution+0.25,]*2, color=default_colors[3], alpha=0.25)
        plt.ylim(0,1)
        ax.tick_params(axis='y', labelcolor=default_colors[0])

        # Generate a new Axes instance, on the twin-X axes (same position)
        stable_sol = cohort_stats["avg_stability"]>0.8
        best_solution_aic_stable = np.arange(1,max_sigs+1)[stable_sol][np.argmin(aic[stable_sol])]
        best_solution_aic = np.argmin(aic)+1
        ax2 = ax.twinx()
        ln2=plt.plot(np.arange(1,max_sigs+1), aic,
                     "*--", label="AIC", c=default_colors[9])
        plt.plot([best_solution_aic_stable,]*2, ax2.get_ylim(), color=default_colors[9], alpha=0.25)
        plt.fill_betweenx(ax2.get_ylim(), [best_solution_aic-0.2,]*2, [best_solution_aic+0.2,]*2,
                          color=default_colors[9], alpha=0.25)
        plt.ylabel("AIC")
        ax2.tick_params(axis='y', labelcolor=default_colors[9])
        IQR = np.nanpercentile(aic, 84)-np.nanpercentile(aic, 16)
        plt.ylim(np.nanmin(aic)-IQR/5, np.nanmedian(aic)+IQR/2)
        # Add legend for lines
        lns = ln1+ln2+ln_cos
        labs = [l.get_label() for l in lns]
        ax.legend(lns, labs, loc='lower left', bbox_to_anchor=(0.,1.))

        plt.xlim(0,max_sigs+1)
        plt.title(cohort, fontsize=24)
        # Save figure
        plt.savefig(f"{output_dir}/figs/{cohort}_aic.png",
               bbox_inches='tight', dpi=200, facecolor='w', transparent=False)

        print(cohort, best_solution, best_solution_aic)

        solution_df=pd.DataFrame(cohort_stats.set_index('signatures').loc[[best_solution_aic]])
        solution_df['cohort']=cohort
        solution_df['samples_file']=cohort_stats.samples_file.iloc[0]
        best_solutions_df = pd.concat((best_solutions_df, solution_df))[['cohort','samples_file','signatures_file','signatures_SE_file','activities_file']]

        cohort_stats[cohort] = cohort
        all_cohort_stats_aic = pd.concat((all_cohort_stats_aic, cohort_stats))

    # Save best solutions and signatures for each cohort
    best_solutions_df.to_csv(f"{output_dir}/best_solutions.tsv",sep="\t",index=True)
    all_cohort_stats_aic.to_csv(f"{output_dir}/all_cohort_stats_aic.tsv",sep="\t",index=True)
