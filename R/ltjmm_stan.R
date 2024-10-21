#' Latent time joint mixed effect model with univariate normal distributions on random effects fit with stan
#'
#' @export
#' @param formula a \code{\link[Formula]{Formula}} with \code{rhs} denoted the outcome in the 
#' stacked dataset, and \code{rhs} with 4 parts, e.g.: \code{Y ~ variable for observation time | fixed
#' effects | subject id | outcome id}
#' @param lt logical, indicating whether or not latent time effect should be included.
#' @param random_effects character specifying distribution for random intercepts and slopes. Option 'univariate' 
#' specifies that random intercepts and slopes for each outcome follow univariate independent normal distributions.
#' Option 'multivariate' specifies that random intercepts and slopes follow a single mutilvariate normal distribution.)
#' @param data data.frame containing the variables in the model.
#' @param subset an optional vector specifying a subset of observations to be used in the fitting process.
#' @param na.action a function which indicates what should happen when the data contain NAs.
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @seealso \code{\link[ltjmm]{ltjmm}}
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'
ltjmm_stan <- function(formula, lt=TRUE, random_effects='univariate', data, subset, na.action, ...){
  mod <- NULL
  if(lt & random_effects == 'univariate') mod <- stanmodels$ltjmm
  if(lt & random_effects == 'multivariate') mod <- stanmodels$ltjmm_mvnorm_ranef
  if(!lt & random_effects == 'univariate') mod <- stanmodels$mm
  if(!lt & random_effects == 'multivariate') mod <- stanmodels$jmm_mvnorm_ranef
  if(is.null(mod)) stop('Invalid specification for lt and/or random_effects.')
  out <- rstan::sampling(mod, 
    data = ltjmm(formula, data, subset, na.action)$data, ...)
  return(out)
}