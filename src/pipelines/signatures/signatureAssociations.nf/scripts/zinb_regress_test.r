suppressMessages(library(tidyverse, rhdf5))

library(tidyverse, MASS)
# , stats)

#' Fit negative binomial with MASS::glm.nb
#'
#' Fit negative binomial with MASS::glm.nb
#' @param data.df = tibble. DataFrame including count column and column for each covarites.
#' @param covariates = vector of strings. Columns in data.df to be used as covariates.
#' @param start = vector (of length covariates) or NULL. Starting values for covariate betas.
#' @param maxit = int. Maximum number of iterations when fitting glm.
#' @param offset = NULL,vector of int. Sample offsets to be used (for fixing parameter values - offset=exp(X beta)).
#' @return list. Output from glm. Includes betas, residuals, beta covariances, NB size parameter (theta) and offsets
#' @examples
#' data.df <- tibble(intercept = rep(c(1), 10), x1 = c(1:10), count = c(21:30))
#' covariates <- c("intercept", "x1")
#' results <- negBinGLM(data.df, covariates)
negBinGLM <- function(data.df, covariates, start = NULL, maxit = 100, offset = NULL, sample_weights = NULL) {
    # If starting parameters provided
    if (!is.null(start)) {
        start.theta <- start$theta
        start.count <- start$count
    } else {
        start.theta <- 1
        start.count <- NULL
    }

    # If there is an offset column, add to formula
    if (!is.null(offset)) {
        data.df <- data.df %>%
            mutate(offset = offset) %>%
            dplyr::select(count, covariates$count, offset)
        formula <- "count ~ 0 + offset(offset) + "
    } else {
        formula <- "count ~ 0 + "
    }

    formula <- as.formula(paste0(formula, paste(covariates$count, collapse = " + ")))
    # print("NB")
    # print(formula)
    # print(start.count)
    # print(start.theta)
    # print(maxit)
    # print(sample_weights[50:60])
    # print(data.df[50:60,])

    fit_nb <- MASS::glm.nb(formula,
        data = data.df,
        start = start.count, init.theta = start.theta,
        maxit = maxit, weights = sample_weights
    )
    # print(fit_nb)

    return(list(
        betas = fit_nb$coefficients,
        residual = data.df$count - as.numeric(fit_nb$fitted.values),
        covar = vcov(fit_nb),
        theta = fit_nb$theta,
        offset = fit_nb$fitted.values,
        loglike = fit_nb$twologlik / 2,
        model = "negbin"
    ))
}

#' Fit GLM with stats::glm
#'
#' Fit GLM with stats::glm
#' @param data.df = tibble. DataFrame including count column and column for each covarites.
#' @param covariates = vector of strings. Columns in data.df to be used as covariates.
#' @param family='negbin' = str. GLM family to be used.
#' @param start = vector (of length covariates) or NULL. Starting values for covariate betas.
#' @param maxit = int. Maximum number of iterations when fitting glm.
#' @param offset = NULL,vector of int. Sample offsets to be used (for fixing parameter values - offset=exp(X beta)).
#' @return list. Output from glm. Includes betas, residuals, beta covariances, NB size parameter (theta) and offsets
#' @examples
#' data.df <- tibble(intercept = rep(c(1), 10), x1 = c(1:10), count = c(21:30))
#' covariates <- c("intercept", "x1")
#' results <- negBinGLM(data.df, covariates)
genLinMod <- function(data.df, covariates, family = "X", start = NULL, maxit = 100, offset = NULL, sample_weights = NULL) {
    # If there is an offset column, add to formula
    if (!is.null(offset)) {
        data.df <- data.df %>% mutate(offset = offset)
        formula <- "count ~ 0 + offset(offset)"
    } else {
        formula <- "count ~ 0"
    }

    formula <- as.formula(paste(c(formula, covariates$count), collapse = " + "))

    fit_glm <- stats::glm(formula,
        data = data.df, family = family,
        start = start$count,
        maxit = maxit, weights = sample_weights
    )

    residual <- data.df$count - as.numeric(fit_glm$fitted.values)
    # Estimate standard error from Gaussian residuals
    # Only meaningful in Gaussian models?
    sigma <- sqrt(sum(residual^2) / (length(residual) - length(covariates$count)))

    return(list(
        betas = fit_glm$coefficients,
        residual = residual,
        sigma = sigma,
        covar = vcov(fit_glm),
        offset = fit_glm$fitted.values,
        loglike = logLik(fit_glm),
        model = family
    ))
}

#' Fit zero inflated GLM with pscl::zeroinfl
#'
#' Fit zero inflated GLM with pscl::zeroinfl
#' @param data.df = tibble. DataFrame including count column and column for each covarites.
#' @param covariates = list of vectors with count and zero entries. Columns in data.df to be used as covariates.
#' @param family='negbin' = str. GLM family to be used.
#' @param start = list or NULL. Starting values for parameters with count, zero and theta entries.
#' @param maxit = int. Maximum number of iterations when fitting glm.
#' @param offset = NULL or list of count:vector, zero:vector.
#' Sample offsets to be used (for fixing parameter values - offset=exp(X beta)).
#' @return list. Output from glm. Includes betas, residuals, beta covariances, NB size parameter (theta) and offsets
#' @examples
#' data.df <- tibble(intercept = rep(c(1), 10), x1 = c(1:10), count = c(21:30))
#' covariates <- c("intercept", "x1")
#' results <- negBinGLM(data.df, covariates, head(covariates, 1))
zeroInflGLM <- function(data.df, covariates, family = "X",
                        start = NULL, maxit = 100, offset = NULL, sample_weights = NULL) {
    if (length(covariates$zero) == 0) {
        if (family == "negbin") {
            return(negBinGLM(data.df, covariates,
                start = start,
                maxit = maxit, offset = offset$count, sample_weights = sample_weights
            ))
        } else {
            return(genLinMod(data.df, covariates,
                start = start, family = family,
                maxit = maxit, offset = offset$count, sample_weights = sample_weights
            ))
        }
    }

    # If there is an offset column, add to formula
    if (!is.null(offset)) {
        data.df <- data.df %>% mutate(offset_count = offset$count, offset_zero = offset$zero)
        formula_count <- "count ~ 0 + offset(offset_count)"
        formula_zero <- "0 + offset(offset_zero)"
    } else {
        formula_count <- "count ~ 0"
        formula_zero <- "0"
    }

    formula <- as.formula(paste0(
        paste(c(formula_count, covariates$count), collapse = " + "), " | ",
        paste(c(formula_zero, covariates$zero), collapse = " + ")
    ))
    print(formula)

    fit_nb <- pscl::zeroinfl(formula,
        data = data.df,
        dist = family, link = "logit", start = start,
        maxit = maxit, reltol = 1e-12, weights = sample_weights
    )

    return(list(
        betas = fit_nb$coefficients$count,
        betas_zero = fit_nb$coefficients$zero,
        residual = data.df$count - as.numeric(fit_nb$fitted.values),
        covar = vcov(fit_nb),
        theta = fit_nb$theta,
        offset = fit_nb$fitted.values,
        loglike = fit_nb$loglik,
        model = family
    ))
}

#' Fit GLM with stats::glm
#'
#' Fit GLM with stats::glm
#' @param data.df = tibble. DataFrame including count column and column for each covarites.
#' @param covariates = vector of strings. Columns in data.df to be used as covariates.
#' @param family='poisson' = str. GLM family to be used.
#' @param zeroinfl=Bool. Use zero inflated GLM.
#' @param start = vector (of length covariates) or NULL. Starting values for covariate betas.
#' @param maxit = int. Maximum number of iterations when fitting glm.
#' @param offset = NULL,vector of int. Sample offsets to be used (for fixing parameter values - offset=exp(X beta)).
#' @return list. Output from glm. Includes betas, offsets and beta covariances.
#' @examples
#' data.df <- tibble(intercept = rep(c(1), 10), x1 = c(1:10), count = c(21:30))
#' covariates <- c("intercept", "x1")
#' results <- glmFit(data.df, covariates, family = "poisson")
glmFit <- function(data.df, covariates, family = "X", zeroinfl = F,
                   start = NULL, maxit = 100, offset = NULL, sample_weights = NULL) {
    if (family == "negbin" & !zeroinfl) {
        results <- negBinGLM(data.df, covariates, start = start, offset = offset, sample_weights = sample_weights)
    } else if (!zeroinfl) {
        results <- genLinMod(data.df, covariates, family = family, start = start, offset = offset, sample_weights = sample_weights)
    } else {
        results <- zeroInflGLM(data.df, covariates, family = family, start = start, offset = offset, sample_weights = sample_weights)
    }

    return(results)
}

#' Fit GLM
#'
#' Fit GLM to counts and covariates
#' @param X = tibble. DataFrame of covariates.
#' @param y = vector of int >= 0. Target variable (count data)
#' @param covariates = list of vector of strings with entries count, zero (if zeroinfl model).
#' Columns in data.df to be used as covariates.
#' @param family='X' = str. GLM family to be used.
#' @param zeroinfl=Bool. Use zero inflated GLM.
#' @param start = list of vector (of length covariates) (count=, zero=) or NULL. Starting values for covariate betas.
#' @param maxit = int. Maximum number of iterations when fitting glm.
#' @param offset = NULL,vector of int. Sample offsets to be used (for fixing parameter values - offset=exp(X beta)).
#' @return list. Output from glm. Includes betas, offsets, beta covariances and other parameters.
#' @examples
#' X <- tibble(intercept = rep(c(1), 10), x1 = c(1:10))
#' y <- c(21:30)
#' results <- glmFit(X, y)
runGlmFit <- function(X, y, covariates = NULL, family = "X", zeroinfl = F, start = NULL,
                      maxit = 100, offset = NULL, sample_weights = NULL) {
    if (is.null(coariates)) {
        covariates <- list(count = colnames(X), zero = colnames(X))
    }
    data.df <- X %>% mutate(count = y)

    results <- glmFit(data.df, covariates,
        family = family, zeroinfl = zeroinfl,
        start = start, maxit = maxit, offset = offset, sample_weights = sample_weights
    )

    return(results)
}

#' Fit GLM with variable resampling
#'
#' Fit GLM to counts and covariates with resampling of the target variable
#' Resampling performed using resample with replacement
#' @param X = tibble. DataFrame of covariates.
#' @param y = vector of int >= 0. Target variable (count data)
#' @param target = str. Column name of variable which is being analysed.
#' @param covariates=NULL
#' @param nulliter=100 = int. Number of resamples to perform for null distribution.
#' @param family='negbin' = str. GLM family to be used.
#' @param zeroinfl=Bool. Use zero inflated GLM.
#' @param maxit = int. Maximum number of iterations when fitting glm.
#' @param use_offsets = Bool. Whether to fix non-target parameters when resampling --> faster but less accurate?
#' @param notarget_results=NULL
#' @param resampling='dCRT', str. 'dCRT' (distilled CRT), 'bootstrap' or 'permutation'
#' --> CRT defaults to permuation if no target covatiates.
#' @return list. Output from glm. Includes betas, offsets, beta covariances and other parameters.
#' @examples
#' X <- tibble(intercept = rep(c(1), 10), x1 = c(1:10))
#' y <- c(21:30)
#' results <- resample_glmFit(X, y, "x1", nulliter = 100, family = "negbin")
resampleGLM <- function(X, y, target, covariates = NULL, nulliter = 100,
                        family = "X", target_family = "X",
                        zeroinfl = F, maxit = 100, use_offsets = F, sample_weights = NULL,
                        notarget_result = NULL,
                        fittarget_result = NULL,
                        resampling = "dCRT") {
    # Covariates excluding the target parameter
    if (is.null(covariates)) {
        covariates <- list(count = c(colnames(X)[colnames(X) != target]))
        if (zeroinfl) {
            covariates$zero <- covariates$count
        }
    }
    data.df <- X %>% mutate(count = y)

    # Fit resampling model
    if ((is.null(fittarget_result)) & (!is.null(covariates$target))) {
        fittarget_result <- glmFit(X %>% rename(count = target), list(count = covariates$target),
            family = target_family, zeroinfl = F, maxit = 100
        )
    }
    # print(c("Fit to target: ", fittarget_result$betas))

    # Fit model without target variable
    if (is.null(notarget_result)) {
        notarget_result <- glmFit(data.df, covariates, family = family, zeroinfl = zeroinfl, maxit = maxit)
    }
    # print(c("No target: ", notarget_result$betas))

    # Get fixed starting parameters from model fit without target
    if (zeroinfl) {
        start <- list(count = c(notarget_result$betas, rnorm(1, 0, 0)), zero = notarget_result$betas_zero)
        offset <- list(
            count = as.numeric(as.vector(notarget_result$betas) %*% t(data.matrix(data.df[covariates$count]))),
            zero = as.numeric(as.vector(notarget_result$betas_zero) %*% t(data.matrix(data.df[covariates$zero])))
        )
        sample_weights <- 1 - expit(offset$zero)
    } else {
        start <- list(count = c(notarget_result$betas, rnorm(1, 0, 0)))
        offset <- list(count = as.numeric(as.vector(notarget_result$betas) %*% t(data.matrix(data.df[covariates$count])))) # log(notarget_result$offset)
        sample_weights <- rep(1, length(offset$count))
    }
    if (family == "negbin") {
        start$theta <- notarget_result$theta
    }

    # Fit model with target variable and hot start
    covariates$count <- c(covariates$count, target)
    alt_result <- glmFit(data.df, covariates, family = family, zeroinfl = zeroinfl, maxit = maxit, start = start)
    # print(c("Alt result: ", alt_result$betas))
    print("Star")
    z_star <- scoreTestZ(data.df[target], data.df$count, exp(offset$count),
        notarget_result$theta,
        sample_weights = sample_weights, family = family
    )

    # Generate resample-with-replacement of target parameter
    target_resampled <- paste0(target, "_resampled")
    covariates$count <- c(covariates$count[covariates$count != target], target_resampled)
    # Run resample and fit lots of times
    z_null <- c()
    # print(head(data.df))
    nulliter <- 2
    print("NT")
    print(c(notarget_result$loglike, notarget_result$betas))
    for (i in c(1:nulliter)) {
        if (resampling == "bootstrap") {
            # Bootstrap
            resample_idx <- sample(c(1:nrow(data.df)), size = nrow(data.df), replace = T)
            # Get null z value from score test method
            z_null <- c(z_null, scoreTestZ(data.df[resample_idx, ][target],
                data.df[resample_idx, ]$count,
                exp(offset$count)[resample_idx],
                notarget_result$theta,
                sample_weights = sample_weights[resample_idx], family = family
            ))
            next
        } else if (resampling == "dCRT" & is.null(covariates$target)) {
            # Permutation
            data.df[target_resampled] <- sample(X[, target], size = length(y), replace = FALSE)
        } else if (resampling == "dCRT") {
            # Distilled CRT
            data.df[target_resampled] <- drawFromGLM(X, covariates$target, fittarget_result, family = target_family)
        } else {
            stop(paste0("resampling method, ", resampling, ", unknown"))
        }

        covariates$count <- c(head(covariates$count, -1), target_resampled)
        # print(notarget_result$offset[1:50])
        # print(offset$count[1:50])
        null_result <- glmFit(data.df, list(count = c(target_resampled)),
            family = family, zeroinfl = zeroinfl,
            maxit = maxit,
            start = list(count = NULL, theta = start$theta),
            offset = offset, sample_weights = sample_weights
        )
        print("Null")
        print(c(null_result$loglike, null_result$betas))
        tgt_result <- glmFit(data.df, list(count = c(target)),
            family = family, zeroinfl = zeroinfl,
            maxit = maxit,
            start = list(count = NULL, theta = start$theta),
            offset = offset, sample_weights = sample_weights
        )
        print("TGT")
        print(c(tgt_result$loglike, tgt_result$betas))

        # Get null z value from score test method
        print("Null")
        z_null <- c(z_null, scoreTestZ(data.df[target_resampled], data.df$count, exp(offset$count),
            notarget_result$theta,
            sample_weights = sample_weights, family = family
        ))
    }
    stop("Interrupted ")

    # Get p-value target parameter
    if (resampling == "bootstrap") {
        z_null <- z_null - mean(z_null)
    }
    pv_score <- nullSamplePV(z_null, z_star)

    return(append(
        alt_result,
        list(
            z_star = z_star, z_null = z_null, pv_score = pv_score,
            null_loglike = notarget_result$loglike
        )
    ))
}

#' Draw sample from GLM
#'
#' Draw sample from GLM with given parameters
#' @param X=data.frame. Dataframe of covariates for model.
#' @param covariates=vector of str. Column headers of covariates used for model.
#' @param parameters=list(betas=vector, sigma?). Parameters of GLM fit.
#' @param family = str. GLM distribution family used. Currently binomial and gaussian implemented.
#' @return vector of same length as X. Samples from the model.
drawFromGLM <- function(X, covariates, parameters, family = "X") {
    K <- as.matrix(X[covariates]) %*% as.vector(parameters$betas)
    if (family == "binomial") {
        mu <- expit(K)
        sample <- rbinom(dim(X)[1], 1, mu)
    } else if (family == "gaussian") {
        sample <- rnorm(dim(X)[1], K, parameters$sigma)
    } else {
        print(paste0("No family ", family, " for scoreTestZ"))
    }

    return(sample)
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
#' pvalue <- nullSamplePV(null_samples, fitted_value, model = "chisquare")
nullSamplePV <- function(null_samples, fitted_value, model = "chisquare") {
    if (model == "chisquare") {
        normal_var <- sd(null_samples)^2
        normal_mean <- mean(null_samples)
        pv <- pchisq((fitted_value - normal_mean)^2 / normal_var, df = 1, lower.tail = F)
    } else {
        stop("Unknown model for null sample passed to nullSamplePV")
    }

    return(pv)
}

#' Calculate expit (aka sigmoid) of x
#'
#' Calculate expit (aka sigmoid) of x
#' @param x = numeric type. Value or array of values.
#' @return p. expit(x).
#' @examples
#' x <- c(-10, -0.1, 0., 0.1, 10)
#' p <- expit(x)
expit <- function(x) {
    p <- 1 / (1 + exp(-x))
    return(p)
}

#' Estimate z-score used in score test
#'
#' Estimate z-score used in score test for negative binomial or poisson distributions
#' @param count: y predicted value from negative binomial draw
#' @param offset: predicted values from NB fit to covariates
#' @param theta: theta parameter from NB fit
#' @param sample_weights: weights for samples in likelihood function
#' @param family="X":
#' @return z score estimate
#' @examples
#' x <- c(-10, -0.1, 0., 0.1, 10)
#' p <- expit(x)
scoreTestZ <- function(X, count, offset, theta = 1, sample_weights = 1, family = "X") {
    if (family == "negbin") {
        dlogL <- sum(sample_weights * X * (count - offset * (count + theta) / (offset + theta)))
        fisherI <- sum(sample_weights * X^2 * offset / (offset + theta))
        print(dlogL)
        print(fisherI)
        print(dlogL / sqrt(fisherI))
    } else if (family == "poisson") {
        dlogL <- sum(sample_weights * X * (count - offset))
        fisherI <- sum(sample_weights * X^2 * offset)
    } else {
        print(paste0("No family ", family, " for scoreTestZ"))
    }

    return(dlogL / sqrt(fisherI))
}


model_fit <- function(sample.df, covariates, targets, target_types, resampling_method) {
    # Set up model and backup for when model fails
    family <- "negbin"
    backup <- "poisson"
    zeroinfl_model <- sum(sample.df$count == 0) > 0
    target_families <- list(linear = "gaussian", logistic = "binomial")

    # Renormalise covariates
    # print("Renormalising")
    means <- c()
    stds <- c()
    for (key in c(covariates, targets[target_types == "linear"])) {
        means <- c(means, mean(sample.df[, key]))
        stds <- c(stds, sd(sample.df[, key]))
        sample.df[key] <- (sample.df[, key] - mean(sample.df[, key])) / sd(sample.df[, key])
    }

    # Select covariates with non-zero variance
    sample.df$intercept <- 1.0
    covariates <- list(
        count = c("intercept", covariates[stds[1:length(covariates) - 1] > 0]),
        zero = c("intercept", covariates[stds[1:length(covariates) - 1] > 0]),
        target = c("intercept", covariates[stds[1:length(covariates) - 1] > 0])
    )

    # Run GLM regression without target variable
    print("No target")
    # fit_notarget <- glmFit(sample.df, covariates, family=family, zeroinfl=zeroinfl_model, maxit=maxit)
    fit_notarget <- tryCatch(glmFit(sample.df, covariates, family = family, zeroinfl = zeroinfl_model, maxit = maxit),
        error = function(e) {
            glmFit(sample.df, covariates, family = backup, zeroinfl = zeroinfl_model, maxit = maxit)
        }
    )
    family <- fit_notarget$model
    print(fit_notarget)

    # For each target, run GLM regression with target variable and null resampling
    for (i in c(1:length(targets))) {
        print(c("target", targets[i]))

        print(c("Resampling: ", resampling_method))
        fit_alt <- resampleGLM(sample.df %>% dplyr::select(-count), sample.df$count, targets[i],
            covariates = covariates, nulliter = 10,
            family = family, target_family = target_families[target_types[i]][[1]],
            zeroinfl = zeroinfl_model, notarget_result = fit_notarget, resampling = resampling_method
        )
        print(
            list(
                target = targets[i], resample_pvalue = fit_alt$pv_score,
                loglike_null = fit_notarget$loglike, loglike_alt = fit_alt$loglike,
                model_null = fit_notarget$model, model_alt = fit_alt$model,
                sample_size = nrow(sample.df), nz_size = sum(sample.df$count > 0)
            ),
            z_star = fit_alt$z_star, z_null = fit_alt$z_null
        )
        results.df <- data.frame(
            target = targets[i], resample_pvalue = fit_alt$pv_score,
            loglike_null = fit_notarget$loglike, loglike_alt = fit_alt$loglike,
            model_null = fit_notarget$model, model_alt = fit_alt$model,
            sample_size = nrow(sample.df), nz_size = sum(sample.df$count > 0)
        )
        results.df$theta_null <- fit_notarget$theta
        results.df$theta_alt <- fit_alt$theta

        # Count paramters
        names <- c(covariates$count, "target")
        for (i in c(1:length(names))) {
            results.df[paste(names[i], "means", "alt", sep = "_")] <- fit_alt$betas[i]
            results.df[paste(names[i], "covs", "alt", sep = "_")] <- diag(fit_alt$covar)[i]
        }
        for (i in c(1:length(covariates$count))) {
            results.df[paste(covariates$count[i], "means", "null", sep = "_")] <- fit_notarget$betas[i]
            results.df[paste(covariates$count[i], "covs", "null", sep = "_")] <- diag(fit_notarget$covar)[i]

            results.df[paste("Xmean", covariates$count[i], sep = "_")] <- means[i]
            results.df[paste("Xsd", covariates$count[i], sep = "_")] <- stds[i]
        }
        for (i in c(1:length(covariates$zero))) {
            results.df[paste(covariates$zero[i], "zero", "means", "alt", sep = "_")] <- fit_alt$betas_zero[i]
            # results.df[paste(covariates$zero[i], "zero", "covs", "alt", sep="_")] <- diag(fit_alt$covar)[i]

            results.df[paste(covariates$zero[i], "zero", "means", "null", sep = "_")] <- fit_notarget$betas_zero[i]
            # results.df[paste(covariates$zero[i], "zero", "covs", "null", sep="_")] <- diag(fit_notarget$covar)[i]
        }

        stacked_results.df <- tryCatch(bind_rows(stacked_results.df, results.df),
            error = function(e) {
                results.df
            }
        )
    }

    return(stacked_results.df)
}


args <- commandArgs(trailingOnly = TRUE)
output_dir <- args[1]
base_dir <- args[2]
test_file <- args[3]
covariate_file <- args[4]
signature_file <- args[5]
target_file <- args[6]
resampling_method <- args[7] # Must be dCRT or bootstrap

# source(paste(base_dir, "regression_functions.r", sep="/"))

# covariates = c("log_age", "logit_purity", "is_female", "pc1", "pc2", "pc3")

# Set of tests to run
test_df <- read_delim(test_file, delim = "\t")
print(test_df[1:10, ])
signatures <- unique(test_df$signature)

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
        if (grp %in% c("Breast", "Ovary", "Uterus")) {
            subsample_df <- sample_df %>% filter(group == grp, is_female == 1)
        } else if (grp %in% list("Testis", "Prostate")) {
            subsample_df <- sample_df %>% filter(group == grp, is_female == 0)
        } else {
            subsample_df <- sample_df %>% filter(group == grp)
        }

        # Filter out instances where only one target value
        nonsingular <- c()
        for (target in targets) {
            nonsingular <- c(nonsingular, nrow(unique(subsample_df[, target])) > 1)
        }
        targets <- targets[nonsingular]
        target_types <- target_types[nonsingular]
        if (length(targets) == 0) {
            stop("No non-unique targets")
        }

        subsample_df <- inner_join(subsample_df,
            signature_df %>% dplyr::select(sample_id, sig) %>% rename(count = sig),
            by = c("sample_id")
        ) %>% dplyr::select(count, targets, covariates)
        # subsample_df <- sample_df %>% filter(group==grp) %>% dplyr::select(sig, targets, covariates) %>% rename(count=sig)

        results_df <- model_fit(data.frame(subsample_df), covariates, gsub("-", ".", targets), target_types, resampling_method)
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
write_delim(full_results_df, paste0(output_dir, "glm_results.csv"), delim = ",")
