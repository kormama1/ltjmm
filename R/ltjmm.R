RowSums <- function(x){
  if(is.null(dim(x))){
    return(x)
  }else{
    return(rowSums(x))
  }
}

rows_dot_product <- function(X, Y) RowSums(X * Y)
dot_product <- function(x, y) x %*% y
normal <- function(mu, sd) stats::rnorm(length(mu), mu, sd)
rep_vector <- function(x, times) rep(x, times)
rep_matrix <- function(x, nrow, ncol) matrix(x, nrow, ncol)
diag_matrix <- function(x) diag(x)
rep_vector <- function(...) rep(...)
multi_normal <- function(mu, Sigma) mvtnorm::rmvnorm(1, mean = mu, sigma = Sigma)
fabs <- abs
diag_pre_multiply <- function(v, m) diag(v) %*% m

#' LTJMM data setup for fitting with rstan
#'
#' @aliases ltjmm
#' @param formula a \code{\link[Formula]{Formula}} with \code{rhs} denoted the outcome in the 
#' stacked dataset, and \code{rhs} with 4 parts, e.g.: \code{Y ~ variable for observation time | fixed
#' effects | subject id | outcome id}
#' @param data data.frame containing the variables in the model.
#' @param subset an optional vector specifying a subset of observations to be used in the fitting process.
#' @param na.action a function which indicates what should happen when the data contain NAs.
#' @seealso \code{\link[rstan]{stan}}
#' 
#' @export
ltjmm <- function(formula, data, subset, na.action){
  mf <- match.call(expand.dots = FALSE)
  m <- match(c('formula', 'data', 'subset', 'na.action'), names(mf), 0)
  mf <- mf[c(1, m)]
  mf <- eval(mf[[3]], parent.frame())
  f <- Formula::Formula(formula)
  
  subject <- unclass(as.factor(stats::model.frame(stats::formula(stats::terms(f, lhs = 0, rhs = 3)), data = mf)[,1]))
  outcome <- unclass(as.factor(stats::model.frame(stats::formula(stats::terms(f, lhs = 0, rhs = 4)), data = mf)[,1]))
  
  outcome_design <- stats::model.matrix(f, data = mf, rhs = 4)
  outcome_name <- attr(stats::terms(f, lhs = 0, rhs = 4), 'term.labels')
  if(length(outcome_name) > 1) stop('Only one outcome variable is supported. Check the 6th part of the formula.')
  colnames(outcome_design) <- gsub(outcome_name, '', colnames(outcome_design))
  mf <- cbind(mf, outcome_design)
  
  data <- list(
    N_obs = length(subject),
    N_sub = length(unique(subject)),
    N_out = length(unique(outcome)),
    y = as.numeric(stats::model.frame(stats::formula(stats::terms(f, lhs = 1, rhs = 0)), data = mf)[,1]),
    subject = subject,
    outcome = outcome
  )
  
  data$obs_time <- stats::model.frame(stats::formula(stats::terms(f, lhs = 0, rhs = 1)), data = mf)[,1]
  
  # get unique observation times and sort them by subject and time
  unq_time <- data.frame(subject = subject, time = data$obs_time)
  unqs <- !duplicated(paste(unq_time$subject, unq_time$time))
  unq_time <- unq_time[unqs,]
  unq_time <- unq_time[with(unq_time, order(subject, time)), ]
  data$unq_id <- unlist(lapply(1:data$N_obs, function(i) which(unq_time$subject == subject[i] & unq_time$time == data$obs_time[i])))
  mf$unq_id <- data$unq_id
  data$N_unq <- nrow(unq_time)
  data$unq_time_subject <- unq_time$subject
  data$unq_time <- unq_time$time
  # the number of unique obs for each subject in turn
  data$sub_obs_size <- as.numeric(table(data$subject))
  
  # covariates for direct effects on outcome
  data$X <- stats::model.matrix(f, data = mf, rhs = 2)
  data$N_X <- dim(data$X)[2] # N latent-time (fixed) effects
  z <- list(formula = formula, data = data)
  class(z) <- "ltjmm"
  z
}

#' Obtain predictions and confidence intervals for LTJMM models
#' 
#' @aliases predict.ltjmm
#' @param object a \code{\link[ltjmm]{ltjmm}} object
#' @param stanfit.object a corresponding \code{\link[rstan]{stanfit}} object
#' @param newdata new data on which to generate predictions
#' @param link link function (defaults to identity)
#' @param chains number of MCMC chains.
#' @param level for credible intervals
#' @param grouping character indicating whether 'population' or 'subject' level precictions
#' are desired.
#' @param ... additional optional arguments.
#' @seealso \code{\link[ltjmm]{ltjmm}}
#' @method predict ltjmm
#' @export 
predict.ltjmm <- function(object, stanfit.object, newdata, link, chains = NULL, level = 0.95,
                          grouping = c('population', 'subject')[1], ...){
  pd <- setup.predict.ltjmm(object, stanfit.object, newdata, link, chains, level)
  X <- pd$data$X
  obs_time <- pd$data$obs_time
  unq_time <- pd$data$unq_time
  unq_id <- pd$data$unq_id
  outcome <- pd$data$outcome
  sub_obs_size <- pd$data$sub_obs_size
  subject <- as.numeric(as.character(pd$data$subject))
  
  if(grouping == 'subject'){
    outcomes <- do.call(cbind, lapply(1:length(pd$samples), function(i){
      pd$link(
        as.vector(RowSums(X * as.matrix(pd$samples[[i]]$beta[outcome,]))) +
          pd$samples[[i]]$gamma[outcome,] * (obs_time + pd$samples[[i]]$delta[subject,]) +
          pd$samples[[i]]$alpha1[cbind(subject, outcome)] * obs_time + pd$samples[[i]]$alpha0[cbind(subject, outcome)]
      )
    }))
  }
  
  if(grouping == 'population'){
    outcomes <- do.call(cbind, lapply(1:length(pd$samples), function(i){
      pd$link(
        as.vector(RowSums(X * as.matrix(pd$samples[[i]]$beta[outcome,]))) +
          pd$samples[[i]]$gamma[outcome,] * obs_time
      )
    }))
  }
  
  outcomes <- as.data.frame(t(apply(outcomes, 1, function(x){
    c(mean = mean(x),
      se_mean = stats::sd(x)/sqrt(length(x)),
      sd = stats::sd(x),
      lower = stats::quantile(x, probs = (1-level)/2)[[1]],
      median = stats::quantile(x, probs = 0.5)[[1]],
      upper = stats::quantile(x, probs = 1-(1-level)/2)[[1]])
  })))
  
  newdata$.unq_id <- unq_id
  newdata_unq <- newdata[order(unq_id), ]
  newdata_unq <- newdata_unq[!duplicated(newdata_unq$.unq_id),]
  newdata <- newdata[,!colnames(newdata) == '.unq_id']
  
  cbind(newdata, outcomes)
}

#' Create design matrices and extract parameter estimates to obtain predictions and confidence intervals for LTJMMs
#' 
#' @param object a \code{\link[ltjmm]{ltjmm}} object
#' @param stanfit.object a corresponding \code{\link[rstan]{stanfit}} object
#' @param newdata new data on which to generate predictions
#' @param link link function (defaults to identity)
#' @param chains number of MCMC chains.
#' @param level for credible intervals
#' @seealso \code{\link[ltjmm]{ltjmm}}
#' @export
setup.predict.ltjmm <- function(object, stanfit.object, newdata, link, chains = NULL, level = 0.95){
  if (missing(newdata) || is.null(newdata)){
    data <- object$data
  }else{
    data <- ltjmm(formula = object$formula, data = newdata)$data
  }
  if (missing(link) || is.null(link)) {
    link <- function(x) x
  }
  
  # assuming all chains have same structure:
  chain <- 1
  all_names <- names(stanfit.object@sim$samples[[chain]])
  thin <- stanfit.object@stan_args[[chain]]$thin
  s0 <- stanfit.object@stan_args[[chain]]$warmup/thin + 1
  s1 <- stanfit.object@stan_args[[chain]]$iter/thin
  
  beta_names <- all_names[grepl('beta', all_names)]
  gamma_names <- all_names[grepl('gamma', all_names)]
  delta_names <- all_names[grepl('delta', all_names)]
  delta_names <- delta_names[-length(delta_names)]
  alpha0_names <- all_names[(grepl('alpha0', all_names)) & !(grepl('sigma_alpha0', all_names))]
  alpha1_names <- all_names[(grepl('alpha1', all_names)) & !(grepl('sigma_alpha1', all_names))]
  
  samples.array <- as.array(stanfit.object)
  if(is.null(chains)){
    chains <- 1:dim(samples.array)[2]
  }
  samples.matrix <- do.call(rbind, lapply(chains, function(chain) samples.array[,chain,]))
  
  samples <- lapply(1:nrow(samples.matrix), function(samp){
    list(
      beta = matrix(samples.matrix[samp, beta_names], nrow = data$N_out, ncol = data$N_X),
      gamma = matrix(samples.matrix[samp, gamma_names], nrow = data$N_out, ncol = 1),
      delta = matrix(samples.matrix[samp, delta_names], nrow = data$N_sub, ncol = 1),
      alpha0 = matrix(samples.matrix[samp, alpha0_names], nrow = data$N_sub, ncol = data$N_out),
      alpha1 = matrix(samples.matrix[samp, alpha1_names], nrow = data$N_sub, ncol = data$N_out)
    )
  })
  return(list(samples = samples, data = data, link = link))
}


#' Simulate data from LTJMM with multivariate normal distribution for random effects
#' 
#' @param object a \code{\link[ltjmm]{ltjmm}} object
#' @param nsim number of response vectors to simulate. Defaults to 1.
#' @param beta fixed effects for covariates
#' @param gamma latent time slope
#' @param sigma_diag (Cholesky factorization diag_matrix(sigma_diag) * Lcorr ) for random intercepts and slopes
#' @param Lcorr (Cholesky factorization diag_matrix(sigma_diag) * Lcorr) for random intercepts and slopes
#' @param sigma_delta standard deviation for latent time
#' @param delta vector of latent times. If NULL, latent times delta are simulated by normal(0, sigma_delta).
#' @param sigma_y standard deviation for residual variance
#' @param seed random seed
#' @param ... additional optional arguments.
#' @seealso \code{\link[ltjmm]{ltjmm}}
#' @importFrom stats simulate model.frame model.matrix terms
#' @method simulate ltjmm
#' @export
simulate.ltjmm <- function(object,
                           nsim = object$data$N_obs,
                           seed = NULL,
                           beta = array(1, c(object$data$N_out, object$data$N_X)),
                           gamma = rep(1, object$data$N_out),
                           sigma_delta = 1,
                           sigma_y = rep(1, object$data$N_out),
                           sigma_diag = diag(rep(1, 2*object$data$N_out-1)),
                           Lcorr = diag(rep(1, 2*object$data$N_out-1)),
                           delta = NULL, ...){
  if(!is.null(seed)) set.seed(seed)
  N_obs <- object$data[['N_obs']]
  N_sub <- object$data[['N_sub']]
  N_out <- object$data[['N_out']]
  y <- object$data[['y']]
  subject <- object$data[['subject']]
  outcome <- object$data[['outcome']]
  obs_time <- object$data[['obs_time']]
  unq_id <- object$data[['unq_id']]
  N_unq <- object$data[['N_unq']]
  unq_time_subject <- object$data[['unq_time_subject']]
  unq_time <- object$data[['unq_time']]
  sub_obs_size <- object$data[['sub_obs_size']]
  X <- object$data[['X']]
  N_X <- object$data[['N_X']]

  # model
  if(is.null(delta)){
    delta <- rep(NA, N_sub)
    for(i in 1:N_sub){
      delta[i] <- normal(0, sigma_delta)
    }
  }
  
  # transformed parameters
  gamma_aug <- rep(NA, N_obs)
  delta_aug <- rep(NA, N_obs)
  beta_aug <- matrix(NA, nrow = N_obs, ncol = N_X)
  alpha0_aug <- matrix(NA, nrow = N_obs, ncol = N_out)
  alpha1_aug <- matrix(NA, nrow = N_obs, ncol = N_out)
  sigma_y_aug <- rep(NA, N_obs)
  mu <- rep(NA, N_obs)
  
  sigma_L <- diag_pre_multiply(sigma_diag, Lcorr)
  z_alpha <- matrix(stats::rnorm(n=(2*N_out-1)*N_sub, mean=0, sd=1), nrow=(2*N_out-1), ncol=N_sub)
  alpha <- t(sigma_L %*% z_alpha)
	alpha0_raw <- alpha[, 1:(N_out-1)]
	alpha1 <- alpha[, N_out:(2*N_out-1)]
  
	# Sum-to-zero constraint for subject-specific random intercepts over outcomes
	alpha0 <- matrix(NA, N_sub, N_out)
	for(n_sub in 1:N_sub){
		for(n_out in 1:(N_out-1)){
			alpha0[n_sub, n_out] = alpha0_raw[n_sub, n_out];
		}
		alpha0[n_sub, N_out] = -sum(alpha0_raw[n_sub,]);
	}
  
  for(n_obs in 1:N_obs){
		beta_aug[n_obs] <- beta[outcome[n_obs]]
		gamma_aug[n_obs] <- gamma[outcome[n_obs]]
		delta_aug[n_obs] <- delta[subject[n_obs]]
		alpha0_aug[n_obs] <- alpha0[subject[n_obs], outcome[n_obs]] 
		alpha1_aug[n_obs] <- alpha1[subject[n_obs], outcome[n_obs]]    
		sigma_y_aug[n_obs] <- sigma_y[outcome[n_obs]]
		mu[n_obs] <- dot_product(X[n_obs], beta_aug[n_obs]) + gamma_aug[n_obs] * (obs_time[n_obs] + delta_aug[n_obs]) + alpha0_aug[n_obs] + alpha1_aug[n_obs] * obs_time[n_obs]
  }
    
  for(n_obs in 1:N_obs){
    y[n_obs] <- normal(mu[n_obs], sigma_y_aug[n_obs])
  }
  
  list(y = y, delta = delta, alpha0 = alpha0, alpha1 = alpha1)
}
