suppressMessages(library(tidyverse))
suppressMessages(library(ggplot2, viridis, glmnet))

# require(pscl)
require(MASS)
require(boot)
library(pscl)


args = commandArgs(trailingOnly=TRUE)
base_dir = args[1]
test_file = args[2]
covariate_file = args[3]
target_file = args[4]
target_variable = args[5]
target_type = args[6]

# source(paste(base_dir, "logistic_functions.r", sep="/"))
source(paste(base_dir, "logistic_regress.r", sep="/"))
source(paste(base_dir, "regression_functions.r", sep="/"))

# Set of tests to run
test_df <- read_delim(test_file, delim="\t")
print(test_df[1:10,])

# Load in samples
sample_df <- read_delim(covariate_file, delim="\t")
covariates <- colnames(sample_df)[! colnames(sample_df) %in% c("sample_id", "group")]
# covariates <- c("log_age", "is_female", "pc1", "pc2", "pc3")
print(c("Covariates: ", covariates))


# Load in target variables
# Trait file should be csv with "germline_platekey" and a column for each trait
target_df <- read_delim(target_file, delim="\t")

# Join tables together
sample_df <- inner_join(sample_df, target_df, by=c("sample_id"), suffix=c("", "_tgt")) %>%mutate_if(is.numeric, ~replace_na(., NA) %>%replace(., is.infinite(.), NA))


# Iterate through tumour groups for the given signature
groups <- unique(test_df$group)
for (grp in groups) {

    # Get all variable targets for the given group and signature
    targets <- unique(test_df[test_df$group==grp,]$target)

    print(c("Group: ", grp))
    print(c("Target variable: ", target_variable))
    subsample_df <- sample_df %>% filter(group==grp)

    # Filter out instances where only one target value
    nonsingular <- c()
    for (target in targets) { nonsingular<-c(nonsingular, nrow(unique(subsample_df[,target]))>1) }
    targets <- targets[nonsingular]
    print(c("Targets: ", targets))
    if (length(targets)==0) { stop("No non-unique targets") }

    print("Subsample info: ")
    print(length(subsample_df))
    subsample_df <- subsample_df %>% mutate(count=.data[[target_variable]]) %>% dplyr::select(count, targets, covariates) %>% filter(!is.na(count))
    print(dim(subsample_df))
    print(head(subsample_df))

    if (dim(subsample_df)[1]>2){
        results_df <- model_fit(data.frame(subsample_df), covariates, gsub("-", ".", targets), family=target_type, power.analysis=F)
        results_df["group"] <- grp
        results_df["dependent"] <- target_variable

        full_results_df <- tryCatch( bind_rows(full_results_df, results_df),
                                    error=function (e) {results_df} )
    }

}

if (nrow(full_results_df)==0) { stop("Failed all runs!!!") }

# Statistical tests
full_results_df["wilks_pvalue"] <- pchisq( -2 * (full_results_df$loglike_null - full_results_df$loglike_alt), df=1, lower.tail=FALSE )

print(head(full_results_df))

# Save results
write_delim(full_results_df %>% dplyr::select(group, dependent, target, wilks_pvalue, everything()), "target_results.csv", delim=",")
