suppressMessages(library(tidyverse, rhdf5))


model_fit <- function(sample.df, covariate_list, targets, target_types, resampling_method, power.analysis = true, zeroinfl = true) {
    # Set up model and backup for when model fails
    family <- "negbin"
    backup <- "poisson"
    zeroinfl_model <- (sum(sample.df$count == 0) > 0) & (zeroinfl)
    # target_families <- list(linear="gaussian", logistic="binomial")

    # Renormalise covariate_list
    # print("Renormalising")
    means <- c()
    stds <- c()
    for (key in c(covariate_list, targets[target_types == "gaussian"])) {
        means <- c(means, mean(sample.df[, key]))
        stds <- c(stds, sd(sample.df[, key]))
        sample.df[key] <- (sample.df[, key] - mean(sample.df[, key])) / sd(sample.df[, key])
    }

    # Select covariate_list with non-zero variance
    sample.df$intercept <- 1.0
    covariates <- list(
        count = c("intercept", covariate_list[stds[1:length(covariate_list)] > 0]),
        zero = c("intercept", covariate_list[stds[1:length(covariate_list)] > 0]),
        target = c("intercept", covariate_list[stds[1:length(covariate_list)] > 0])
    )

    # Run GLM regression without target variable
    print("No target")
    print(covariates)
    print(covariate_list)
    print(stds[1:length(covariate_list) - 1])
    print(stds)
    fit_notarget <- tryCatch(glmFit(sample.df, covariates, family = family, zeroinfl = zeroinfl_model),
        error = function(e) {
            glmFit(sample.df, covariates, family = backup, zeroinfl = zeroinfl_model)
        }
    )
    family <- fit_notarget$model
    print("Done")

    # For each target, run GLM regression with target variable and null resampling
    for (i in c(1:length(targets))) {
        print(c(i, length(targets), " target: ", targets[i]))

        print("Removing nan target values - to make null fit fair (glmFit does this automatically when doing alt)")
        print(head(sample.df))
        print(dim(sample.df)[1])
        sample_nona.df <- sample.df %>% filter(!is.na(get(targets[i])))
        print(head(sample_nona.df))
        print(dim(sample_nona.df)[1])
        if (dim(sample_nona.df)[1] != dim(sample.df)[1]) {
            print("Fit nona")
            fit_notarget_nona <- tryCatch(glmFit(sample_nona.df, covariates, family = family, zeroinfl = zeroinfl_model),
                error = function(e) {
                    list(
                        loglike = -99, model = "NA", theta = -99,
                        betas = rep(-99, length(covariates$count) + 1), betas_zero = rep(-99, length(covariates$zero)),
                        covar = matrix(-99, ncol = length(covariates$count) + 1, nrow = length(covariates$count) + 1)
                    )
                }
            )
            family_nona <- fit_notarget$model
        } else {
            fit_notarget_nona <- fit_notarget
            family_nona <- family
        }

        print(c("Resampling: ", resampling_method))
        fit_alt <- tryCatch(
            resampleGLM(sample_nona.df %>% dplyr::select(-count), sample_nona.df$count, targets[i],
                covariates = covariates, nulliter = 500,
                family = family_nona, target_family = target_types[i],
                zeroinfl = zeroinfl_model, notarget_result = fit_notarget_nona, resampling = resampling_method
            ),
            error = function(e) {
                list(
                    pv_score = -99, loglike = -99, model = "NA",
                    betas = rep(-99, length(covariates$count) + 1), betas_zero = rep(-99, length(covariates$zero)),
                    covar = matrix(-99, ncol = length(covariates$count) + 1, nrow = length(covariates$count) + 1),
                    z_star = -99, z_null = c(-99, -99)
                )
            }
        )

        print(c("Z ", fit_alt$z_star, mean(na.omit(fit_alt$z_null)), sd(na.omit(fit_alt$z_null))))

        # print("Collecting results")
        # print(c('ll', fit_notarget_nona['loglike']))
        # print(c('ll', fit_notarget_nona$loglike))
        wilks <- pchisq(-2 * (fit_notarget_nona$loglike - fit_alt$loglike), df = 1, lower.tail = FALSE)
        print(list(
            target = targets[i], resample_pvalue = fit_alt$pv_score,
            theta_null = fit_notarget_nona$theta, theta_alt = fit_alt$theta,
            loglike_null = fit_notarget_nona$loglike, loglike_alt = fit_alt$loglike, wilks = wilks,
            model_null = fit_notarget_nona$model, model_alt = fit_alt$model,
            sample_size = nrow(sample_nona.df), nz_size = sum(sample_nona.df$count > 0),
            alt_params = fit_alt$betas
        ))
        results.df <- data.frame(
            target = targets[i], resample_pvalue = fit_alt$pv_score,
            loglike_null = fit_notarget_nona$loglike, loglike_alt = fit_alt$loglike,
            model_null = fit_notarget_nona$model, model_alt = fit_alt$model,
            sample_size = nrow(sample_nona.df), nz_size = sum(sample_nona.df$count > 0),
            target_z_means_alt = fit_alt$z_star, target_z_means_null = mean(na.omit(fit_alt$z_null)), target_z_covs_null = sd(na.omit(fit_alt$z_null))^2
        )
        results.df$theta_null <- fit_notarget_nona$theta
        results.df$theta_alt <- fit_alt$theta


        # Count paramters
        names <- c(covariates$count, "target")
        for (j in c(1:length(names))) {
            results.df[paste(names[j], "means", "alt", sep = "_")] <- fit_alt$betas[j]
            results.df[paste(names[j], "covs", "alt", sep = "_")] <- diag(fit_alt$covar)[j]
        }
        for (j in c(1:length(covariates$count))) {
            results.df[paste(covariates$count[j], "means", "null", sep = "_")] <- fit_notarget_nona$betas[j]
            results.df[paste(covariates$count[j], "covs", "null", sep = "_")] <- diag(fit_notarget_nona$covar)[j]
        }
        for (j in c(1:length(covariate_list))) {
            results.df[paste("Xmean", covariate_list[j], sep = "_")] <- means[j]
            results.df[paste("Xsd", covariate_list[j], sep = "_")] <- stds[j]
        }
        for (j in c(1:length(covariates$zero))) {
            results.df[paste(covariates$zero[j], "zero", "means", "alt", sep = "_")] <- fit_alt$betas_zero[j]
            # results.df[paste(covariates$zero[j], "zero", "covs", "alt", sep="_")] <- diag(fit_alt$covar)[j]
            results.df[paste(covariates$zero[j], "zero", "means", "null", sep = "_")] <- fit_notarget_nona$betas_zero[j]
            # results.df[paste(covariates$zero[j], "zero", "covs", "null", sep="_")] <- diag(fit_notarget_nona$covar)[j]
        }


        ### Power Analysis
        power.percentiles <- c(20, 50, 80)
        if ((power.analysis) & !(is.na(fit_alt$pv_score)) & (fit_alt$pv_score != -99)) {
            # Set parameters to known ones with target effect size of log(2)
            # Resample from Zero-Inflated Negative Binomial N times from new parameters
            # Apply resampleGLM method to simulated signature
            # Fit beta distribution to set of p-values
            effect.size <- log(2)
            power.ntest <- 10

            fit_target <- targetFit(sample.df %>% dplyr::select(-count), list(count = covariates$target), targets[i],
                family = target_types[i], zeroinfl = F, maxit = 100
            )


            mu <- exp(as.matrix(sample.df[c(covariates$count, targets[i])]) %*% c(fit_alt$betas[1:length(covariates$count)], effect.size))
            print(c("Input betas ", c(fit_alt$betas[1:length(covariates$count)], effect.size)))
            if (zeroinfl_model) {
                # mu zero is the expected probability of sample having activity forced to zero
                muzero <- 1 - expit(as.matrix(sample.df[covariates$zero]) %*% fit_alt$betas_zero)
            }
            pvalue_list <- c()
            rpv_list <- c()
            for (j in c(1:power.ntest)) {
                # Reorder true counts by drawn counts
                # counts <- sort(sample.df$count)[order(order( drawFromModel(mu, muzero, theta=fit_alt$theta, zeroinfl=zeroinfl_model, family=family) ))]
                # counts <- drawFromModel(mu, muzero, theta=fit_alt$theta, zeroinfl=zeroinfl_model, family=family)
                counts <- drawFromModel(mu, muzero, zeroinfl = zeroinfl_model, family = "poisson")
                fit_power <- tryCatch(
                    resampleGLM(sample.df %>% dplyr::select(-count), counts, targets[i],
                        covariates = covariates, nulliter = 500,
                        family = family_nona, target_family = target_types[i],
                        zeroinfl = zeroinfl_model, resampling = resampling_method,
                        fittarget_result = fit_target
                    ),
                    error = function(e) {
                        "failed"
                    }
                )
                # notarget_result=fit_notarget_nona, fittarget_result=fit_target,
                if (fit_power == "failed") {
                    next
                }
                # print(fit_power$z_star)
                # print(fit_power$z_null)
                print(mean(fit_power$z_null))
                print(sd(fit_power$z_null))
                wilks_pv <- pchisq(-2 * (fit_power$null_loglike - fit_power$loglike), df = 1, lower.tail = FALSE)
                pvalue_list <- c(pvalue_list, wilks_pv)
                rpv_list <- c(rpv_list, fit_power$pv_score)
                print(c("GLMft: ", wilks_pv))
                print(c("Pass:  ", fit_power$pv_score))
            }

            for (percentile in power.percentiles) {
                results.df[paste0("power", percentile)] <- quantile(na.omit(pvalue_list), percentile / 100)[[1]]
            }
            for (percentile in power.percentiles) {
                results.df[paste0("resample_power", percentile)] <- quantile(na.omit(rpv_list), percentile / 100)[[1]]
            }
            results.df[paste0("powerMean")] <- mean(na.omit(pvalue_list))
            results.df[paste0("powerVar")] <- sd(na.omit(pvalue_list))^2
            results.df[paste0("powerN")] <- length(na.omit(pvalue_list))
        } else {
            for (percentile in power.percentiles) {
                results.df[paste0("power", percentile)] <- -99
            }
            for (percentile in power.percentiles) {
                results.df[paste0("resample_power", percentile)] <- -99
            }
            results.df[paste0("powerMean")] <- -99
            results.df[paste0("powerVar")] <- -99
            results.df[paste0("powerN")] <- -99
        }

        ### Return results
        stacked_results.df <- tryCatch(bind_rows(stacked_results.df, results.df),
            error = function(e) {
                results.df
            }
        )
    }

    return(stacked_results.df)
}

drawFromModel <- function(mu, muzero, theta = NULL, zeroinfl = F, family = "X") {
    # Draw from count process
    if (family == "negbin") {
        counts <- rnbinom(dim(mu)[1], mu = mu, size = theta)
    } else if (family == "poisson") {
        counts <- rpois(dim(mu)[1], lambda = mu)
    }

    # Draw zero inflation from binomial distribution
    if (zeroinfl) {
        nonzero <- rbinom(dim(mu)[1], size = 1, prob = muzero)
    } else {
        nonzero <- 1
    }

    return(nonzero * counts)
}

if (sys.nframe() == 0) {
    args <- commandArgs(trailingOnly = TRUE)
    output_dir <- args[1]
    base_dir <- args[2]
    test_file <- args[3]
    covariate_file <- args[4]
    signature_file <- args[5]
    target_file <- args[6]
    resampling_method <- args[7] # Must be dCRT or bootstrap
    power.analysis <- args[8] # true or not
    zeroinfl <- args[9] == "true"
    print("Zeroinfl? ")
    print(zeroinfl)

    print("Power analysis?")
    power.analysis <- power.analysis == "true"
    print(power.analysis)

    source(paste(base_dir, "regression_functions.r", sep = "/"))

    # covariates = c("log_age", "logit_purity", "is_female", "pc1", "pc2", "pc3")

    # Set of tests to run
    test_df <- read_delim(test_file, delim = "\t")
    print(test_df[1:10, ])
    signatures <- unique(test_df$signature)
    print(test_df)
    print(signatures)

    # Load in samples
    sample_df <- read_delim(covariate_file, delim = "\t")
    covariates <- colnames(sample_df)[!colnames(sample_df) %in% c("sample_id", "group")]
    # covariates <- c("log_age", "is_female", "pc1", "pc2", "pc3")

    # Load in signatures
    print(c("sample_id", signatures))
    signature_df <- read_delim(signature_file, delim = "\t") %>% dplyr::select(c("sample_id", signatures))
    signature_df[is.na(signature_df)] <- 0

    # Load in target variables
    # Trait file should be csv with "germline_platekey" and a column for each trait
    target_df <- read_delim(target_file, delim = "\t")

    # Join tables together
    # sample_df <- inner_join(sample_df, signature_df, by=c("sample_id"), suffix=c("", "_sig"))
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
            target_types <- target_types[nonsingular]
            print(targets)
            if (length(targets) == 0) {
                stop("No non-unique targets")
            }

            # Get all cols with are either target or target_i for multicomponent target resampling
            target_cols <- c()
            for (target in targets) {
                target_cols <- c(
                    target_cols,
                    colnames(subsample_df)[str_detect(colnames(subsample_df), regex(paste0(target, "[_0-9]*")))]
                )
            }

            subsample_df <- inner_join(subsample_df,
                signature_df %>% dplyr::select(sample_id, sig) %>% rename(count = sig),
                by = c("sample_id")
            ) %>% dplyr::select(count, target_cols, covariates)
            # subsample_df <- sample_df %>% filter(group==grp) %>% dplyr::select(sig, targets, covariates) %>% rename(count=sig)

            results_df <- model_fit(data.frame(subsample_df), covariates, gsub("-", ".", targets), target_types, resampling_method, power.analysis = power.analysis, zeroinfl = zeroinfl)
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

    # print(full_results_df)

    # Save results
    write_delim(full_results_df %>% dplyr::select(group, signature, target, wilks_pvalue, everything()), paste0(output_dir, "glm_results.csv"), delim = ",")
}
