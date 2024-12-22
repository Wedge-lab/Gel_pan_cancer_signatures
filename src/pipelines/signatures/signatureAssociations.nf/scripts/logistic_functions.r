suppressMessages(library(tidyverse))
suppressMessages(library(ggplot2, viridis, glmnet))

# require(pscl)
require(MASS)
require(boot)
library(pscl)



logistic_fit <- function(data.df, covariates) {

  # Run glm.nb - Initialise theta=1 if fit with auto calculating theta fails
  formula <- as.formula( paste(c("k ~ 0", covariates), collapse=" + ") )
  # print(c("Formula: ", formula))
  suppressWarnings(results <- glm(formula, data=data.df, family='binomial'))

  count_means <- c(results$coefficients[covariates])
  count_covs <- c(diag(vcov(results))[covariates])

  loglike <- length(results$coefficients)-results$aic/2
  # pred = results$coefficients[1] + results$coefficients[2:length(results$coefficients)] %*% t(data.matrix(variables.df[,2:dim(variables.df)[2]]))
  # pred <- 1/(1+exp(-pred))
  # loglike <- sum(variables.df$k * log(pred) + (1-variables.df$k) * log(1-pred))

  fit_results = list(count_means=count_means, count_covs=count_covs,
                     converged=results$converged, loglike=loglike, model='glm.binom',
                     dof=1)

  return(fit_results)

}

resample_logistic_fit <- function(data.df, covariates, null_prob, nulliter=100) {

    # Alternative fit
    fit_results_alt   <- logistic_fit(data.df, covariates)
    z_star <- fit_results_alt$count_means[length(covariates)]/sqrt(fit_results_alt$count_covs[length(covariates)])
    print(data.df$k)
    print(sum(data.df$k))

    # Resample n times to get null fits
    z_null <- c()
    for (i in c(1:nulliter)) {
        data.df$k <- rbinom(length(null_prob), 1, null_prob)

        # print(paste("Alternative fit", test_variable, sep=" "))
        fit_resample_null <- logistic_fit(data.df, covariates)
        z_null <- c(z_null, fit_resample_null$count_means[length(covariates)]/sqrt(fit_resample_null$count_covs[length(covariates)]) )

    }

    # Get p-value target parameter
    pv_score <- nullSamplePV(z_null, z_star)

    return(append(fit_results_alt, list(z_star=z_star, z_null=z_null, pv_score=pv_score)))

}


zero_model_fit <- function(variables.df, covariates, test_variables) {

  # Convert signature to binary 1/0
  variables.df$k <- variables.df$k>0

  # Renormalise covariates
  print("Renormalising")
  means = c()
  stds = c()
  for ( key in c(covariates) ) {
    means = c(means, mean(variables.df[, key]))
    stds = c(stds, sd(variables.df[, key]))
    variables.df[key] = ( variables.df[, key]-mean(variables.df[, key]) )/sd(variables.df[, key])
  }
  # print(variables.df %>%top_n(10))
  # Select covariates with non-zero variance
  variables.df$intercept <- 1.0
  covariates = c("intercept", covariates[stds[1:length(covariates)-1]>0])

  # print(select_(data.frame(variables.df), .dots = c('k', covariates, 'target')))
  # print("Null fit")
  fit_results_null <- logistic_fit(variables.df, covariates)

  for (test_variable in test_variables) {

    # print(paste("Alternative fit", test_variable, sep=" "))
    fit_results_alt   <- logistic_fit(variables.df, c(covariates, test_variable))

    # print("Collecting results")
    results.df <- tibble(target=test_variable,
                         alt_fit=fit_results_null["converged"][[1]], null_fit=fit_results_alt["converged"][[1]],
                         loglike_null=fit_results_null["loglike"][[1]], loglike_alt=fit_results_alt["loglike"][[1]],
                         model_null=fit_results_null["model"][[1]], model_alt=fit_results_alt["model"][[1]],
                         dof=fit_results_null["dof"][[1]],
                         sample_size=nrow(variables.df), nz_size=sum(variables.df$k>0))

    # Count paramters
    # print("Processing results")
    names <- c(covariates, "target", "theta")
    for (i in 1:length(names)){
        results.df[paste(names[i], "means", "alt", sep="_")] <- fit_results_alt$count_means[i]
        results.df[paste(names[i], "covs", "alt", sep="_")] <- fit_results_alt$count_covs[i]
    }
    names <- c(covariates, "theta")
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


nmax_model_fit <- function(variables.df, covariates, test_variables, n_resample=100) {

  # Renormalise covariates
  print("Renormalising")
  means = c()
  stds = c()
  for ( key in c(covariates) ) {
    means = c(means, mean(variables.df[, key]))
    stds = c(stds, sd(variables.df[, key]))
    variables.df[key] = ( variables.df[, key]-mean(variables.df[, key]) )/sd(variables.df[, key])
  }
  # print(variables.df %>%top_n(10))
  # Select covariates with non-zero variance
  variables.df$intercept <- 1.0
  covariates = c("intercept", covariates[stds[1:length(covariates)-1]>0])

  for (test_variable in test_variables) {

    # Generate new variable from max n samples in target
    n_hit = min( sum(variables.df[, test_variable]), sum(variables.df$k>0) )
    hit = order(- (variables.df$k + runif(length(variables.df$k))*0.1) )[1:n_hit]
    new_k = rep(0,length(variables.df$k))
    new_k[hit] <- 1
    new_variables.df <- variables.df %>% mutate(k = new_k)

    # print(select_(data.frame(variables.df), .dots = c('k', covariates, 'target')))
    # print("Null fit")
    fit_results_null <- logistic_fit(new_variables.df, covariates)
    nullprob <- 1/( 1+exp( -(as.matrix(new_variables.df[covariates]) %*% as.vector(fit_results_null$count_means))[,1] ) )

    # print(paste("Alternative fit", test_variable, sep=" "))
    fit_results_alt <- resample_logistic_fit(new_variables.df, c(covariates, test_variable), nullprob)
    fit_results_alt <- resampleGLM(new_variables.df %>% dplyr::select(-count), new_variables.df$count, test_variable,
                                    covariates=covariates, nulliter=500,
                                    family=family, target_family=target_families[target_types[i]][[1]],
                                    zeroinfl=F, notarget_result=fit_notarget, resampling=resampling_method)

    results.df <- data.frame(target=test_variable,
                             resample_pvalue=fit_results_alt$pv_score,
                             null_fit=fit_results_null$converged, alt_fit=fit_results_alt$converged,
                             loglike_null=fit_results_null$loglike, loglike_alt=fit_results_alt$loglike,
                             model_null=fit_results_null$model, model_alt=fit_results_alt$model,
                             dof=fit_results_null$dof,
                             sample_size=nrow(new_variables.df), nz_size=sum(new_variables.df$k>0), nz_target=sum(new_variables.df[, test_variable]))

    # Count paramters
    # print("Processing results")
    names <- c(covariates, "target", "theta")
    for (i in 1:length(names)){
        results.df[paste(names[i], "means", "alt", sep="_")] <- fit_results_alt$count_means[i]
        results.df[paste(names[i], "covs", "alt", sep="_")] <- fit_results_alt$count_covs[i]
    }
    names <- c(covariates, "theta")
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


#' Estimate p-value from null samples and fitted coefficient.
#'
#' Estimate p-value from null samples and fitted coefficient.
#' @param null_samples=vector. Vector of z-scores of null fits.
#' @param fitted_value=float. Value of z-score for true run.
#' @param model = str. Null model to be used. ('chisquare' (D), 'skewt'). Only chisquare currently implemented
#' @return float. P-value of fitted value against null samples under null model.
#' @examples
#' null_samples <- rnorm(100)
#' fitted_value <- 5
#' pvalue <- nullSamplePV(null_samples, fitted_value, model='chisquare')
nullSamplePV <- function(null_samples, fitted_value, model='chisquare') {

    if (model=='chisquare') {
          normal_var = sd(null_samples)^2
          normal_mean = mean(null_samples)
          pv = pchisq((fitted_value-normal_mean)^2/normal_var, df=1, lower.tail=F)
    } else {
          stop('Unknown model for null sample passed to nullSamplePV')
    }

    return(pv)

}