#' Latent Time Joint Mixed Effect Models ('ltjmm')
#' 
#' @description Diseases that progress over long periods of time are often studied
#'   by observing cohorts at different stages of disease for shorter periods of
#'   time. We apply MCMC sampling to estimate Latent Time Joint Mixed Effect Models
#'   (LTJMM) from short-term observations with unknown relative observation times.
#' 
#' @docType package
#' @name ltjmm-package
#' @useDynLib ltjmm, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @import rstantools
#' @importFrom rstan sampling
#' @author Michael Donohue <mdonohue@usc.edu>
#' 
#' @references 
#' Stan Development Team (2018). RStan: the R interface to Stan. R package version 2.18.1. http://mc-stan.org
#' 
NULL
