library(tidyverse, MASS)
#, stats)

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
#' data.df <- tibble(intercept=rep(c(1), 10), x1=c(1:10), count=c(21:30))
#' covariates <- c('intercept','x1')
#' results <- negBinGLM(data.df, covariates)
negBinGLM <- function(data.df, covariates, start=NULL, maxit=100, offset=NULL, sample_weights=NULL) {

    # If starting parameters provided
    if (! is.null(start)) {
      start.theta <- start$theta
      start.count <- start$count
    } else {
      start.theta <- 1
      start.count <- NULL
    }

    # If there is an offset column, add to formula
    if (! is.null(offset)) {
      data.df <- data.df %>% mutate(offset=offset)
      formula <- "count ~ 0 + offset(offset) + "
    } else {
      formula <- "count ~ 0 + "
    }

    formula <- as.formula( paste0( formula, paste(covariates$count, collapse=" + ") ) )

    fit_nb <- MASS::glm.nb( formula, data=data.df,
                            start=start.count, init.theta=start.theta,
                            maxit=maxit, weights=sample_weights )

    return(list(betas=fit_nb$coefficients,
                residual=data.df$count-as.numeric(fit_nb$fitted.values),
                covar=vcov(fit_nb),
                theta=fit_nb$theta,
                offset=fit_nb$fitted.values,
                loglike=fit_nb$twologlik/2,
                model='negbin'))

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
#' data.df <- tibble(intercept=rep(c(1), 10), x1=c(1:10), count=c(21:30))
#' covariates <- c('intercept','x1')
#' results <- negBinGLM(data.df, covariates)
genLinMod <- function(data.df, covariates, family='X', start=NULL, maxit=100, offset=NULL, sample_weights=NULL) {

    # For binomial regression with n>1 - need to give a success and fail count column
    if (str_detect(family, regex("[a-z]+[0-9]+"))) {
        glm_family <- str_match(family, regex("([a-z]+)([0-9]+)"))[[2]]
        n <- strtoi(str_match(family, regex("([a-z]+)([0-9]+)"))[[3]])
        data.df$countFail <- n-data.df$count
        formula <- "cbind(count,countFail)"
    } else {
        glm_family <- family
        formula <- "count"
    }

    # If there is an offset column, add to formula
    if (! is.null(offset)) {
      data.df <- data.df %>% mutate(offset=offset)
      formula <- paste(formula, "0 + offset(offset)", sep=" ~ ")
    } else {
      formula <- paste(formula, "0", sep=" ~ ")
    }

    formula <- as.formula( paste(c(formula, covariates$count), collapse=" + ") )

    fit_glm <- stats::glm( formula, data=data.df, family=glm_family,
                            start=start$count,
                            maxit=maxit, weights=sample_weights )

    residual <- data.df$count-as.numeric(fit_glm$fitted.values)
    # Estimate standard error from Gaussian residuals
    # Only meaningful in Gaussian models?
    sigma <- sqrt(sum(residual^2)/(length(residual)-length(covariates$count)))

    return(list(betas=fit_glm$coefficients,
                residual=residual,
                sigma=sigma,
                covar=vcov(fit_glm),
                offset=fit_glm$fitted.values,
                loglike=logLik(fit_glm),
                model=family))

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
#' data.df <- tibble(intercept=rep(c(1), 10), x1=c(1:10), count=c(21:30))
#' covariates <- c('intercept','x1')
#' results <- negBinGLM(data.df, covariates, head(covariates,1))
zeroInflGLM <- function(data.df, covariates, family='X',
                     start=NULL, maxit=100, offset=NULL, sample_weights=NULL) {

    if (length(covariates$zero)==0) {
        if (family=='negbin') {
            return( negBinGLM(data.df, covariates, start=start,
                              maxit=maxit, offset=offset$count, sample_weights=sample_weights) )
        } else {
            return( genLinMod(data.df, covariates, start=start, family=family,
                              maxit=maxit, offset=offset$count, sample_weights=sample_weights) )
        }
    }

    # If there is an offset column, add to formula
    if (! is.null(offset)) {
      data.df <- data.df %>% mutate(offset_count=offset$count, offset_zero=offset$zero)
      formula_count <- "count ~ 0 + offset(offset_count)"
      formula_zero <- "0 + offset(offset_zero)"
    } else {
      formula_count <- "count ~ 0"
      formula_zero <- "0"
    }

    formula <- as.formula( paste0( paste(c(formula_count, covariates$count), collapse=" + "), " | ",
                                   paste(c(formula_zero, covariates$zero), collapse=" + ") ) )

    fit_nb <- pscl::zeroinfl( formula, data=data.df,
                              dist=family, link="logit", start=start,
                              maxit=maxit,  reltol=1e-12, weights=sample_weights)

    return(list(betas=fit_nb$coefficients$count,
                betas_zero=fit_nb$coefficients$zero,
                residual=data.df$count-as.numeric(fit_nb$fitted.values),
                covar=vcov(fit_nb),
                theta=fit_nb$theta,
                offset=fit_nb$fitted.values,
                loglike=fit_nb$loglik,
                model=family))

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
#' data.df <- tibble(intercept=rep(c(1), 10), x1=c(1:10), count=c(21:30))
#' covariates <- c('intercept','x1')
#' results <- glmFit(data.df, covariates, family='poisson')
glmFit <- function(data.df, covariates, family='X', zeroinfl=F,
                    start=NULL, maxit=100, offset=NULL, sample_weights=NULL) {

    if (family=='negbin' & ! zeroinfl) {
        results <- negBinGLM(data.df, covariates, start=start, offset=offset, sample_weights=sample_weights)
    } else if (! zeroinfl ) {
        results <- genLinMod(data.df, covariates, family=family, start=start, offset=offset, sample_weights=sample_weights)
    } else {
        results <- zeroInflGLM(data.df, covariates, family=family, start=start, offset=offset, sample_weights=sample_weights)
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
#' X <- tibble(intercept=rep(c(1), 10), x1=c(1:10))
#' y=c(21:30)
#' results <- glmFit(X, y)
runGlmFit <- function(X, y, covariates=NULL, family="X", zeroinfl=F, start=NULL,
                      maxit=100, offset=NULL, sample_weights=NULL) {

    if (is.null(coariates)) {
        covariates <- list(count=colnames(X), zero=colnames(X))
    }
    data.df <- X %>% mutate(count=y)

    results <- glmFit(data.df, covariates, family=family, zeroinfl=zeroinfl,
                      start=start, maxit=maxit, offset=offset, sample_weights=sample_weights)

    return(results)

}

targetFit <- function(X, covariates, target,
                      family=NULL, zeroinfl=F, maxit=100) {

    # Target may be split into multiple components
    fittarget <- list()
    # Split family into consituent parts if they exist
    families = str_split(family, ",")[[1]]
    if ( length(families)==1 ) {
        # Only one target family
        fittarget <- glmFit(X %>% rename(count=target), covariates,
                            family=family, zeroinfl=F, maxit=100)
    } else {
        # Multiple target families
        for (i in c(2:length(families))) {
            fittarget[[i-1]] <- glmFit(X %>% rename(count=paste(target, i-2, sep="_")), covariates,
                                                            family=families[i], zeroinfl=F, maxit=100)
            print(c(i, fittarget[[i-1]]$betas))
        }
    }

    return(fittarget)

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
#' X <- tibble(intercept=rep(c(1), 10), x1=c(1:10))
#' y=c(21:30)
#' results <- resample_glmFit(X, y, 'x1', nulliter=100, family='negbin')
resampleGLM <- function(X, y, target, covariates=NULL, nulliter=100,
                        family='X', target_family='X',
                        zeroinfl=F, maxit=100, use_offsets=F, sample_weights=NULL,
                        notarget_result=NULL, fittarget_result=NULL, alt_result=NULL,
                        resampling='dCRT') {

    # Covariates excluding the target parameter
    if (is.null(covariates) ) {
        covariates <- list( count=c(colnames(X)[colnames(X)!=target]) )
        if (zeroinfl) { covariates$zero <- covariates$count }
    }
    data.df <- X %>% mutate(count=y)

    # Fit resampling model
    if ( (is.null(fittarget_result)) & (! is.null(covariates$target)) ) {
        fittarget_result <- targetFit(X, list(count=covariates$target), target,
                                      family=target_family, zeroinfl=F, maxit=100)
    }
    # print(c("Fit to target: ", fittarget_result$betas))

    # Fit model without target variable
    if (is.null(notarget_result) ) {
        notarget_result <- glmFit(data.df, covariates, family=family, zeroinfl=zeroinfl, maxit=maxit)
    }
    # print(c("No target: ", notarget_result$betas))

    # Get fixed starting parameters from model fit without target
    if (zeroinfl) {
        start = list(count=c(notarget_result$betas, rnorm(1,0,0)), zero=notarget_result$betas_zero)
        offset = list(count=as.numeric( as.vector(notarget_result$betas) %*% t(data.matrix(data.df[covariates$count]))),
                      zero=as.numeric( as.vector(notarget_result$betas_zero) %*% t(data.matrix(data.df[covariates$zero])) ))
        sample_weights <- 1-expit(offset$zero)
    } else {
        start = list(count=c(notarget_result$betas, rnorm(1,0,0)))
        offset = list(count=as.numeric( as.vector(notarget_result$betas) %*% t(data.matrix(data.df[covariates$count]))))#log(notarget_result$offset)
        sample_weights <- rep(1,length(offset$count))
    }
    if (family=='negbin') { start$theta = notarget_result$theta }

    # Fit model with target variable and hot start
    covariates$count <- c(covariates$count, target)
    if ( is.null(alt_result) ) {
        alt_result <- glmFit(data.df, covariates, family=family, zeroinfl=zeroinfl, maxit=maxit, start=start)
    }

    # # Stratified resample where targets not well distributed from model
    strat_weights <- 1

    # print(c("Alt result: ", alt_result$betas))
    z_star <- scoreTestZ(data.df[target], data.df$count, offset$count,
                         notarget_result$theta, sample_weights=sample_weights*strat_weights, family=family)
    # print(table(data.df[target]))
    counts <- data.frame(index=0)

    # Generate resample-with-replacement of target parameter
    target_resampled <- paste0(target, "_resampled")
    covariates$count <- c(covariates$count[covariates$count!=target], target_resampled)
    # Run resample and fit lots of times
    z_null = c()
    # print(head(data.df))
    for (i in c(1:nulliter)) {

        if (resampling=="bootstrap") {
            # Bootstrap
            resample_idx <- sample(c(1:nrow(data.df)), size=nrow(data.df), replace=T)
            # Get null z value from score test method
            z_null <- c(z_null, scoreTestZ(data.df[resample_idx,][target],
                                           data.df[resample_idx,]$count,
                                           offset$count[resample_idx],
                                           notarget_result$theta,
                                           sample_weights=sample_weights[resample_idx], family=family))
            next

        } else if ( resampling=='dCRT' & is.null(covariates$target)) {
            # Permutation
            data.df[target_resampled] <- sample(data.df[,target], size = length(y), replace = FALSE)
        } else if ( resampling=='dCRT' ) {
            # Distilled CRT
            data.df[target_resampled] <- drawTarget(data.df, covariates$target, fittarget_result, family=target_family)
            # print(table(data.df[target_resampled]))
            # counts <- rbind(counts, pivot_wider(data.frame(table(data.df[target_resampled])), names_from="Var1", values_from="Freq"))
            if (target_family!="gaussian") {
                # print(counts)
                # print(head(data.frame(table(data.df[target_resampled]))))
                counts <- merge(counts, data.frame(table(data.df[target_resampled])), all=T, by.x="index", by.y="Var1")
                colnames(counts)[i+1] <- i
                # print(c("Resampled", table(data.df[target_resampled])))
            }
        } else {
            stop(paste0("resampling method, ", resampling,", unknown"))
        }

        # Get null z value from score test method
        z_null <- c(z_null, scoreTestZ(data.df[target_resampled], data.df$count, offset$count,
                                       notarget_result$theta, sample_weights=sample_weights, family=family))
    }

    if ( resampling=='dCRT' & !is.null(covariates$target) & !(target_family=='gaussian') ) {
        counts[is.na(counts)] <- 0
        means <- rowMeans(counts[,2:ncol(counts)])
        stds <- sqrt(rowSums((counts[,2:ncol(counts)]-means)^2)/(ncol(counts)-1))
        print(c("idx", counts$index))
        print(c("mu", round(means,2)))
        print(c("sd", round(stds,2)))
    }

    # Get p-value target parameter
    if (resampling=="bootstrap") {z_null <- z_null-mean(z_null)}
    pv_score <- nullSamplePV(na.omit(z_null), z_star)

    return(append(alt_result,
                  list(z_star=z_star, z_null=z_null, pv_score=pv_score,
                       null_loglike=notarget_result$loglike)))

}


#'
drawTarget <- function(X, covariates, parameter_sets, family="X") {

    if (length(str_split(family, ",")[[1]])==1) {
        sample <- drawFromGLM(X, covariates, parameter_sets, family=family)
    } else {
        sample <- integer(nrow(X))
        for (i in c(2:length(str_split(family, ",")[[1]]))) {
            sample <- cbind(sample, drawFromGLM(X, covariates, parameter_sets[[i-1]], family=parameter_sets[[i-1]]$model))
            # print(c(i, table(sample[,i])))
        }
        if (str_detect(str_split(family, ",")[[1]][1], regex("[a-z]+[0-9]+"))) {
            n <- strtoi(str_match(str_split(family, ",")[[1]][1], regex("[a-z]+([0-9]+)"))[[2]])
            sample <- rowSums(sample)
            sample[sample>n] <- n
        } else {
            print(paste0("No family ", family, " for multiple resampled GLM"))
        }
    }

    return(sample)

}

#' Draw sample from GLM
#'
#' Draw sample from GLM with given parameters
#' @param X=data.frame. Dataframe of covariates for model.
#' @param covariates=vector of str. Column headers of covariates used for model.
#' @param parameters=list(betas=vector, sigma?). Parameters of GLM fit.
#' @param family = str. GLM distribution family used. Currently binomial and gaussian implemented.
#' @return vector of same length as X. Samples from the model.
drawFromGLM <- function(X, covariates, parameters, family="X") {

    K <- as.matrix(X[covariates]) %*% as.vector(parameters$betas)
    if (family=='binomial') {
        mu <- expit(K)
        sample <- rbinom(dim(X)[1], 1, mu)
    } else if (str_detect(family, regex("binomial[0-9]+"))) {
        n <- strtoi(str_match(family, regex("binomial([0-9]+)"))[[2]])
        mu <- expit(K)
        sample <- rbinom(dim(X)[1], n, mu)
    } else if (str_detect(family, regex("quasibinomial[0-9]+"))) {
        n <- strtoi(str_match(family, regex("quasibinomial([0-9]+)"))[[2]])
        beta_draw = q_fit$betas + matrix(rnorm(N*length(parameters$betas)), nrow=N) %*% chol(parameters$cov)
        K <- rowSums(as.matrix(X) * as.matrix(beta_draw))
        mu <- expit(K)
        sample <- rbinom(dim(X)[1], n, mu)
    } else if (family=='poisson') {
        mu <- exp(K)
        sample <- rpois(dim(X)[1], mu)
    } else if (family=='gaussian') {
        sample <- rnorm(dim(X)[1], K, parameters$sigma)
    } else {
        print(paste0("No family ", family, " for GLM draw"))
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

#' Calculate expit (aka sigmoid) of x
#'
#' Calculate expit (aka sigmoid) of x
#' @param x = numeric type. Value or array of values.
#' @return p. expit(x).
#' @examples
#' x=c(-10,-0.1,0.,0.1,10)
#' p <- expit(x)
expit <- function(x) {
    p <- 1/(1+exp(-x))
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
#' x=c(-10,-0.1,0.,0.1,10)
#' p <- expit(x)
scoreTestZ <-  function(X, count, offset, theta=1, sample_weights=1, family="X") {

    if (family=="negbin") {
        mu <- exp(offset)
        dlogL <- sum( sample_weights * X * ( count - mu*(count+theta)/(mu+theta) ) )
        fisherI <- sum( sample_weights * X^2 * mu/(mu+theta) )
    } else if (family=="poisson") {
        mu <- exp(offset)
        dlogL <- sum( sample_weights * X * ( count - mu ) )
        fisherI <- sum( sample_weights * X^2 * mu )
    } else if (str_detect(family, regex("binomial[0-9]*"))) {
        if (str_detect(family, regex("binomial[0-9]+"))) {
            # Binomial with n>1
            n <- strtoi(str_match(family, regex("binomial([0-9]+)"))[[2]])
        } else { n <- 1 }
        p <- expit(offset)
        dlogL <- sum( sample_weights * X * ( count - n*p ) )
        fisherI <- sum( sample_weights * X^2 * p*(1-p) * n )
    } else {
        print(paste0("No family - ", family, " - for scoreTestZ"))
    }

    return(dlogL/sqrt(fisherI))

}


#' Stratified resampling from binomial distribution
#'
#' Perform stratified resampling from binomial distribution
#' @param k = vector_int(N): value of k for all samples
#' @param p = vector_flt(N): binomial model probability for each sample
#' @param n = int: number of binomial trials
#' @param breaks = seq(0,1,0.1): stratification bins used for p
#' @return subsample = vector_int(N): resampled indices
binomialStratifiedResample <- function(k, p, n, breaks=seq(0,1,0.1), output="samples") {

    # Possible values of k
    k_unique <- seq(0,n,1)

    # Get k probabilities
    p_matrix <- t(replicate(n+1,p))
    p_matrix <- (p_matrix^k_unique * (1-p_matrix)^(n-k_unique)) * choose(n,k_unique)

    # Get number of samples with k=k and p in bin
    breaks <- seq(0,1,0.1)
    break_idx <- findInterval(p, breaks)
    # break_idx <- ecdf(unique(breaks))(break_idx+0.1)
    # Get ID of k in bins (equivalent to k+1)
    k_idx <- findInterval(k, k_unique-0.5)
    k_break_matrix <- table(break_idx, k_idx)

    # Get expected number of samples with k=k and p in bin
    expected_matrix <- rowsum(t(p_matrix), break_idx)

    # Get weights for given bin and k
    weights <- expected_matrix[,unique(k_idx)]/k_break_matrix
    weights[which(!is.finite(weights))]=0
    sample_probs <- weights/max(weights)
    # print(c("k_break_matrix", k_break_matrix))
    # print(c("expected_matrix", expected_matrix))
    # print(c("sample_probs", sample_probs))
    # Get unique positions of break_idx
    break_idx <- match(break_idx, sort(unique(break_idx)))
    sample_probs <- sample_probs[cbind(break_idx, k_idx)]

    # Get resamples - with or without replacement?
    # print(sample_probs)
    if (output=="samples") {
        # subsample <- c(1:length(sample_probs))[rbinom(length(sample_probs),1,sample_probs)==1]
        subsample <- sample(c(1:length(sample_probs)), prob=sample_probs, replace = T)
        # print(c("k", table(k)))
        # print(c("new_k", table(k[subsample])))
        return(sort(subsample))
    } else if (output=="weights") {
        weights <- sample_probs/mean(sample_probs)
        # print(weights[1:10])
        return(weights)
    }

}
