suppressMessages(library(tidyverse))
suppressMessages(library(ggplot2, viridis, glmnet))

# require(pscl)
require(MASS)
require(boot)
library(pscl)

model_fit <- function(variables.df, covariates, targets, family = "binomial", power.analysis = T) {
  # Convert signature to binary 1/0
  variables.df$count <- variables.df$count > 0

  # Renormalise covariates
  print("Renormalising")
  means <- c()
  stds <- c()
  for (key in c(covariates)) {
    means <- c(means, mean(variables.df[, key]))
    stds <- c(stds, sd(variables.df[, key]))
    variables.df[key] <- (variables.df[, key] - mean(variables.df[, key])) / sd(variables.df[, key])
  }
  # print(variables.df %>%top_n(10))
  # Select covariates with non-zero variance
  variables.df$intercept <- 1.0
  print(c("Covariates: ", covariates))
  print(c("SD of covariates: ", stds))
  covariates <- list(count = c("intercept", covariates[stds[1:length(covariates)] > 0]))
  print(covariates)

  # print(select_(data.frame(variables.df), .dots = c('k', covariates, 'target')))
  # print("Null fit")
  # fit_null <- logistic_fit(variables.df, covariates)

  fit_notarget <- glmFit(variables.df, covariates, family = family, zeroinfl = F)

  for (target in targets) {
    print(target)

    print("Removing nan target values - to make null fit fair (glmFit does this automatically when doing alt)")
    # print(head(variables.df[,0:5]))
    # print(dim(variables.df)[1])
    variables_nona.df <- variables.df %>% filter(!is.na(get(target)))
    # print(head(variables_nona.df[,0:5]))
    # print(dim(variables_nona.df)[1])

    if (dim(variables_nona.df)[1] != dim(variables.df)[1]) {
      fit_notarget_nona <- glmFit(variables_nona.df, covariates, family = family, zeroinfl = F)
    } else {
      fit_notarget_nona <- fit_notarget
    }

    # print(paste("Alternative fit", target, sep=" "))
    fit_alt <- glmFit(variables_nona.df, list(count = c(covariates$count, target)), family = family, zeroinfl = F)

    # print("Collecting results")
    results.df <- data.frame(
      target = target,
      loglike_null = fit_notarget_nona$loglike, loglike_alt = fit_alt$loglike,
      model_null = fit_notarget_nona$model, model_alt = fit_alt$model,
      sample_size = nrow(variables_nona.df), nz_size = sum(variables_nona.df$count > 0), nz_target = sum(variables_nona.df[, target])
    )

    # Count paramters
    # print("Processing results")
    names <- c(covariates$count, "target")
    for (j in 1:length(names)) {
      results.df[paste(names[j], "means", "alt", sep = "_")] <- fit_alt$betas[j]
      results.df[paste(names[j], "covs", "alt", sep = "_")] <- diag(fit_alt$covar)[j]
    }
    names <- c(covariates$count)
    for (j in 1:length(names)) {
      results.df[paste(names[j], "means", "null", sep = "_")] <- fit_notarget_nona$betas[j]
      results.df[paste(names[j], "covs", "null", sep = "_")] <- diag(fit_notarget_nona$covar)[j]
    }

    # Renormalisation values
    # print("Processing renormalisation values")
    for (j in 1:length(covariates$count)) {
      results.df[paste("Xmean", covariates$count[j], sep = "_")] <- means[j]
      results.df[paste("Xsd", covariates$count[j], sep = "_")] <- stds[j]
    }

    ### Power Analysis
    power.percentiles <- c(20, 50, 80)
    if (power.analysis) {
      # Set parameters to known ones with target effect size of log(2)
      # Resample from Zero-Inflated Negative Binomial N times from new parameters
      # Apply resampleGLM method to simulated signature
      # Fit beta distribution to set of p-values
      effect.size <- log(2)
      power.ntest <- 10

      # New count probability
      mu <- expit(as.matrix(variables.df[c(covariates$count, target)]) %*% c(fit_alt$betas[1:length(covariates$count)], effect.size))

      pvalue_list <- c()
      for (j in c(1:power.ntest)) {
        counts <- rbinom(dim(mu)[1], size = 1, prob = mu)
        fit_notarget_power <- glmFit(variables.df %>% mutate(count = counts), covariates, family = family, zeroinfl = F)
        fit_alt_power <- glmFit(variables.df %>% mutate(count = counts), list(count = c(covariates$count, target)), family = family, zeroinfl = F)

        pvalue_power <- pchisq(-2 * (fit_notarget_power$loglike - fit_alt_power$loglike), df = 1, lower.tail = FALSE)
        pvalue_list <- c(pvalue_list, pvalue_power)
        print(c("Logistic fit", fit_notarget_power$loglike, fit_alt_power$loglike))
        print(c("Pass: ", pvalue_power))
      }

      for (percentile in power.percentiles) {
        results.df[paste0("power", percentile)] <- quantile(na.omit(pvalue_list), percentile / 100)[[1]]
      }
      results.df[paste0("powerMean")] <- mean(na.omit(pvalue_list))
      results.df[paste0("powerVar")] <- sd(na.omit(pvalue_list))^2
      results.df[paste0("powerN")] <- length(na.omit(pvalue_list))
    } else {
      for (percentile in power.percentiles) {
        results.df[paste0("power", percentile)] <- -99
      }
      results.df[paste0("powerMean")] <- -99
      results.df[paste0("powerVar")] <- -99
      results.df[paste0("powerN")] <- -99
    }

    stacked_results.df <- tryCatch(bind_rows(stacked_results.df, results.df),
      error = function(e) {
        results.df
      }
    )
  }

  return(stacked_results.df)
}


nmax_model_fit <- function(variables.df, covariates, targets, target_types, resampling_method) {
  # Set up model and backup for when model fails
  family <- "binomial"

  # Renormalise covariates
  print("Renormalising")
  means <- c()
  stds <- c()
  for (key in c(covariates)) {
    means <- c(means, mean(variables.df[, key]))
    stds <- c(stds, sd(variables.df[, key]))
    variables.df[key] <- (variables.df[, key] - mean(variables.df[, key])) / sd(variables.df[, key])
  }
  # print(variables.df %>%top_n(10))
  # Select covariates with non-zero variance
  variables.df$intercept <- 1.0
  # covariates = c("intercept", covariates[stds[1:length(covariates)-1]>0])
  covariates <- list(
    count = c("intercept", covariates[stds[1:length(covariates) - 1] > 0]),
    target = c("intercept", covariates[stds[1:length(covariates) - 1] > 0])
  )

  for (i in c(1:length(targets))) {
    # Generate new variable from max n samples in target
    n_hit <- min(sum(variables.df[, targets[i]]), sum(variables.df$count > 0))
    hit <- order(-(variables.df$count + runif(length(variables.df$count)) * 0.1))[1:n_hit]
    new_count <- rep(0, length(variables.df$count))
    new_count[hit] <- 1
    new_variables.df <- variables.df %>% mutate(count = new_count)

    # print(select_(data.frame(variables.df), .dots = c('k', covariates, 'target')))
    # print("Null fit")
    # fit_null <- logistic_fit(new_variables.df, covariates)
    # nullprob <- 1/( 1+exp( -(as.matrix(new_variables.df[covariates]) %*% as.vector(fit_null$count_means))[,1] ) )
    fit_notarget <- glmFit(new_variables.df, covariates, family = family, zeroinfl = F)

    # print(paste("Alternative fit", target, sep=" "))
    fit_alt <- resampleGLM(new_variables.df %>% dplyr::select(-count), new_variables.df$count, targets[i],
      covariates = covariates, nulliter = 500,
      family = family,
      target_family = target_types[i],
      zeroinfl = F, notarget_result = fit_notarget, resampling = resampling_method
    )

    # print("Collecting results")
    results.df <- data.frame(
      target = targets[i], resample_pvalue = fit_alt$pv_score,
      loglike_null = fit_notarget$loglike, loglike_alt = fit_alt$loglike,
      model_null = fit_notarget$model, model_alt = fit_alt$model,
      sample_size = nrow(new_variables.df), nz_size = sum(new_variables.df$count > 0), nz_target = sum(new_variables.df[, targets[i]])
    )


    # Count paramters
    # print("Processing results")
    names <- c(covariates$count, "target")
    for (j in 1:length(names)) {
      results.df[paste(names[j], "means", "alt", sep = "_")] <- fit_alt$betas[j]
      results.df[paste(names[j], "covs", "alt", sep = "_")] <- diag(fit_alt$covar)[j]
    }
    names <- c(covariates$count)
    for (j in 1:length(names)) {
      results.df[paste(names[j], "means", "null", sep = "_")] <- fit_notarget$betas[j]
      results.df[paste(names[j], "covs", "null", sep = "_")] <- diag(fit_notarget$covar)[j]
    }

    # Renormalisation values
    # print("Processing renormalisation values")
    for (j in 1:length(covariates$count)) {
      results.df[paste("Xmean", covariates$count[j], sep = "_")] <- means[j]
      results.df[paste("Xsd", covariates$count[j], sep = "_")] <- stds[j]
    }

    stacked_results.df <- tryCatch(bind_rows(stacked_results.df, results.df),
      error = function(e) {
        results.df
      }
    )
  }

  return(stacked_results.df)
}


if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  output_dir <- args[1]
  base_dir <- args[2]
  test_file <- args[3]
  covariate_file <- args[4]
  signature_file <- args[5]
  target_file <- args[6]
  logistic_method <- args[7]
  resample_method <- args[8]

  # source(paste(base_dir, "logistic_functions.r", sep="/"))
  source(paste(base_dir, "regression_functions.r", sep = "/"))

  # covariates = c("log_age", "logit_purity", "is_female", "pc1", "pc2", "pc3")

  # Set of tests to run
  test_df <- read_delim(test_file, delim = "\t")
  print(test_df[1:10, ])
  signatures <- unique(test_df$signature)

  # Load in samples
  sample_df <- read_delim(covariate_file, delim = "\t")
  covariates <- colnames(sample_df)[!colnames(sample_df) %in% c("sample_id", "group")]
  print(covariates)
  # covariates <- c("log_age", "is_female", "pc1", "pc2", "pc3")

  # Load in signatures
  print(c("sample_id", signatures))
  signature_df <- read_delim(signature_file, delim = "\t") %>% dplyr::select(c("sample_id", signatures))
  signature_df[is.na(signature_df)] <- 0

  # Load in target variables
  # Trait file should be csv with "germline_platekey" and a column for each trait
  target_df <- read_delim(target_file, delim = "\t")

  # Join tables together
  sample_df <- inner_join(sample_df, target_df, by = c("sample_id"), suffix = c("", "_tgt")) %>% mutate_if(is.numeric, ~ replace_na(., NA) %>% replace(., is.infinite(.), NA))


  # Iterate through all signatures
  for (sig in signatures) {
    # Iterate through tumour groups for the given signature
    groups <- unique(test_df[test_df$signature == sig, ]$group)
    for (grp in groups) {
      # Get all variable targets for the given group and signature
      targets <- test_df[test_df$signature == sig & test_df$group == grp, ]$target
      target_types <- test_df[test_df$signature == sig & test_df$group == grp, ]$type

      print(sig)
      print(grp)
      print(targets)
      subsample_df <- sample_df %>% filter(group == grp)

      # Filter out instances where only one target value
      nonsingular <- c()
      for (target in targets) {
        nonsingular <- c(nonsingular, nrow(unique(subsample_df[, target])) > 1)
      }
      targets <- targets[nonsingular]
      print(targets)
      if (length(targets) == 0) {
        stop("No non-unique targets")
      }

      subsample_df <- inner_join(subsample_df,
        signature_df %>% dplyr::select(sample_id, sig) %>% rename(count = sig),
        by = c("sample_id")
      ) %>% dplyr::select(count, targets, covariates)

      if (logistic_method == "ZERO") {
        results_df <- model_fit(data.frame(subsample_df), covariates, gsub("-", ".", targets))
      } else if (logistic_method == "ZERO50") {
        # Take top 50% of activities if more than 50% non-zero
        n_hit <- min(floor(nrow(subsample_df) * 0.5), sum(subsample_df$count > 0))
        hit <- order(-(subsample_df$count + runif(length(subsample_df$count)) * 0.1))[1:n_hit]
        new_count <- rep(0, length(subsample_df$count))
        new_count[hit] <- 1
        subsample_df <- subsample_df %>% mutate(count = new_count)
        results_df <- model_fit(data.frame(subsample_df), covariates, gsub("-", ".", targets))
      } else if (logistic_method == "NMAX") {
        results_df <- nmax_model_fit(data.frame(subsample_df), covariates, gsub("-", ".", targets), target_types, resample_method)
      } else {
        stop(c("logistic_method ", logistic_method, " not known"))
      }
      results_df["signature"] <- sig
      results_df["group"] <- grp

      full_results_df <- tryCatch(bind_rows(full_results_df, results_df),
        error = function(e) {
          results_df
        }
      )
    }
  }

  # Statistical tests
  full_results_df["wilks_pvalue"] <- pchisq(-2 * (full_results_df$loglike_null - full_results_df$loglike_alt), df = 1, lower.tail = FALSE)

  # print(head(full_results_df))

  # Save results
  write_delim(full_results_df %>% dplyr::select(group, signature, target, wilks_pvalue, everything()), "logistic_results.csv", delim = ",")
}
