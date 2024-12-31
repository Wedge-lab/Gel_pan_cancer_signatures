"Find best solution for each cohort"
import sys, os
import pandas as pd, numpy as np
import inspect
from SigProfilerExtractor import subroutines as sub

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

        all_similarities, cosine_similarities = sub.calculate_similarities(np.array(genomes), np.dot(w,h),
                                                                           genomes.columns)
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
    max_sigs=int(max_sigs)
    print(output_dir, min_stability)

    # Load Cohorts
    cohort_stats = pd.read_csv(cohort_stats_file, sep="\t")
    cohorts = np.unique(cohort_stats.cohort)

    # Initialise dataframes
    best_solutions_df = pd.DataFrame()

    # Iterate through cohorts to find best solutions
    for cohort in cohorts:
        # cohort_dir = f"{sig_dir}/CANCERTYPE_{cohort}/{sig_set}_random_max{max_sigs}sigs_500nmf_reps"
        print(f"Best solutions for {cohort}")

        # Estimate best solution for signatures
        best_solution = estimate_solution(stats_df=cohort_stats[cohort_stats.cohort==cohort], #f"{cohort_dir}/All_solutions_stat.csv",
                                            genomes=cohort_stats[cohort_stats.cohort==cohort].samples_file.iloc[0], #f"{cohort_dir}/Samples.txt",
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

        print(cohort, best_solution)

        solution_df=pd.DataFrame(cohort_stats[cohort_stats.cohort==cohort].set_index('signatures').loc[[best_solution]])
        solution_df['cohort']=cohort
        solution_df['samples_file']=cohort_stats[cohort_stats.cohort==cohort].samples_file.iloc[0]
        best_solutions_df = pd.concat((best_solutions_df, solution_df))[['cohort','samples_file','signatures_file','signatures_SE_file','activities_file']]

    # Save best solutions and signatures for each cohort
    best_solutions_df.to_csv(f"{output_dir}/best_solutions.tsv",sep="\t",index=True)
