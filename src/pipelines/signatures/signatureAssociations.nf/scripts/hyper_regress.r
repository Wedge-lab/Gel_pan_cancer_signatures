suppressMessages(library(tidyverse))
suppressMessages(library(ggplot2, viridis, glmnet))
library(rhdf5)

# require(pscl)
require(MASS)
require(boot)
library(pscl)



logistic_fit <- function(data.df, covariates) {

  # Run glm.nb - Initialise theta=1 if fit with auto calculating theta fails
  results <- glm(k ~ . , data=data.df, family='binomial')

  count_means <- c(results$coefficients[c("(Intercept)", covariates)])
  count_covs <- c(diag(vcov(results))[c("(Intercept)", covariates)])

  loglike <- length(results$coefficients)-results$aic/2
  # pred = results$coefficients[1] + results$coefficients[2:length(results$coefficients)] %*% t(data.matrix(variables.df[,2:dim(variables.df)[2]]))
  # pred <- 1/(1+exp(-pred))
  # loglike <- sum(variables.df$k * log(pred) + (1-variables.df$k) * log(1-pred))

  fit_results = list(count_means=count_means, count_covs=count_covs,
                     converged=results$converged, loglike=loglike, model='glm.nb',
                     dof=1)

  return(fit_results)

}


model_fit <- function(variables.df, covariates, test_variables) {


  # Renormalise covariates
  print("Renormalising")
  means = c()
  stds = c()
  for ( key in c(covariates, test_variables) ) {
    means = c(means, mean(variables.df[, key]))
    stds = c(stds, sd(variables.df[, key]))
    variables.df[key] = ( variables.df[, key]-mean(variables.df[, key]) )/sd(variables.df[, key])
  }
  # print(variables.df %>%top_n(10))
  # Select covariates with non-zero variance
  covariates = covariates[stds[1:length(covariates)-1]>0]

  # print(select_(data.frame(variables.df), .dots = c('k', covariates, 'prs')))
  # print("Null fit")
  fit_results_null <- logistic_fit(select_(data.frame(variables.df), .dots = c('k', covariates)), covariates)

  for (test_variable in test_variables) {

    # print(paste("Alternative fit", test_variable, sep=" "))
    fit_results_alt   <- logistic_fit(select_(data.frame(variables.df), .dots = c('k', covariates, test_variable)), c(covariates, test_variable))

    # print("Collecting results")
    results.df <- tibble(trait_fit=fit_results_null["converged"][1], null_fit=fit_results_alt["converged"][1],
                         loglike_null=fit_results_null["loglike"][[1]], loglike_alt=fit_results_alt["loglike"][[1]],
                         model_null=fit_results_null["model"][[1]], model_alt=fit_results_alt["model"][[1]],
                         dof=fit_results_null["dof"][[1]],
                         sample_size=nrow(variables.df), nz_size=sum(variables.df$k>0),
                         trait=test_variable)

    # Count paramters
    # print("Processing results")
    names <- c("c", covariates, "prs", "theta")
    for (i in 1:length(names)){
        results.df[paste(names[i], "means", "alt", sep="_")] <- fit_results_alt$count_means[i]
        results.df[paste(names[i], "covs", "alt", sep="_")] <- fit_results_alt$count_covs[i]
    }
    names <- c("c", covariates, "theta")
    for (i in 1:length(names)){
        results.df[paste(names[i], "means", "null", sep="_")] <- fit_results_null$count_means[i]
        results.df[paste(names[i], "covs", "null", sep="_")] <- fit_results_null$count_covs[i]
    }

    # Renormalisation values
    # print("Processing renormalisation values")
    for (i in 1:length(covariates)) {
        results.df[paste("Xmean", covariates[i], sep="_")] <- means[i]
        results.df[paste("Xsd", covariates[i], sep="_")] <- stds[i]
    }

    stacked_results.df <- tryCatch( bind_rows(stacked_results.df, results.df),
                                  error=function (e) {results.df} )
  }

  return(stacked_results.df)
}


args = commandArgs(trailingOnly=TRUE)
trait_file = args[1]
covariate_file = args[2]
category_col = args[3]

# Load in data
covariates.df <- data.frame(read_csv(covariate_file))
categories <- unique(covariates.df$category)
covariates <- colnames(covariates.df)[-c(1:3)]
print(c("Covar", nrow(covariates.df)))

# Trait file should be csv with "germline_platekey" and a column for each trait
trait.df <- data.frame(read_csv(trait_file))
traits <- colnames(trait.df)[-c(1)]

data.df <- inner_join(covariates.df, trait.df, by=c("tumour_sample_platekey"="tumour_sample_platekey")) %>%mutate_if(is.numeric, ~replace_na(., NA) %>%replace(., is.infinite(.), NA))
data.df <- mutate(data.df, is_female=as.double(is_female), category=category+1)
print(c("Data", nrow(data.df)))

print(categories)
print(covariates)
print(traits)


data.df$k = as.integer(data.df$hypermutation>0)
k_list = list(data.df$hypermutation>0,
              data.df$hypermutation>1)
dependent = c("hyper1", "hyper2")

for ( i_k in c( 1:length(k_list) ) ) {
  print(dependent[i_k])
  data.df$k <- as.integer(k_list[[i_k]])

  results.df <- model_fit(data.frame(data.df), covariates, traits)
  results.df["category"] <- "_all_"
  results.df["dependent_variable"] <- dependent[i_k]
  full_results.df <- tryCatch(bind_rows(full_results.df, results.df),
                                    error=function (e) { results.df } )

  for (i_category in 1:length(categories)) {

    category_data.df <- data.frame(data.df) %>% filter(category==i_category)

    print(categories[i_category])
    print(sum(category_data.df$k>0))

    if ( sum(category_data.df$k>0)<5 ) {
      next
    } else {
      results.df <- model_fit(data.frame(category_data.df), covariates, traits)
    }
    results.df["category"] <- categories[i_category]
    results.df["dependent_variable"] <- dependent[i_k]

    full_results.df <- bind_rows(full_results.df, results.df)
  }
}


print(full_results.df)
# if ( !exists(full_results.df) ) { quit("no", 120) }

# Statistical tests
full_results.df["wilks"] <- pchisq( -2 * (full_results.df$loglike_null - full_results.df$loglike_alt), df=full_results.df$dof, lower.tail=FALSE )

# Save results
write_delim(full_results.df %>%relocate(where(is.character)), paste("hyper_results.csv", sep="_"), delim=",")
