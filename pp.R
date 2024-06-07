# Computing an Approximate Predictive Probability (aPP) From Interim Z-Score
# Berry Consultants
# Joe Marion, Liz Lorenzi, Cora Allen-Savietta, Kert Viele, & Scott Berry
# last update: June 2024


#' Create an approximate predictive probability from an interim Z-score
#' @param z: z-score 
#' @param n: information in data used to compute the z-score
#' @param N: information when the data are complete.
#' @param alpha: the level of the final analysis (one-sided)
#' @return predictive probability scalar in [0, 1]
z_to_pp = function(z, n, N, alpha=0.025){
  stopifnot(n>0)
  stopifnot(N>=n)
  stopifnot(alpha>=0 & alpha<=1)
  r = n/N
  pnorm((z - qnorm(1-alpha)*sqrt(r))/sqrt(1-r))
}

#' Create an approximate predictive probability from a p-value
#' @param p: one-sided p-value from a frequentest analysis. 
#' @param n: information used to compute the z-score
#' @param N: information when the data are complete.
#' @param alpha: the level of the final analysis (one-sided). 
#'               for Bayes analysis use alpha=1-thr.
#' @return predictive probability scalar in [0, 1]
p_to_pp = function(p, n, N, alpha=0.025) {
  z = p_to_z(p)
  z_to_pp(z, n, N, alpha)
}

#' Create a predictive probability from a posterior probability
#'  for Bayesian models
#' @param prob: posterior probability
#' @param n: information used to compute the z-score
#' @param N: information when the data are complete.
#' @param thr: the level of the final analysis (one-sided). 
#'               for Bayes analysis use alpha=1-thr.
#' @return predictive probability scalar in [0, 1]
prob_to_pp = function(prob, n, N, thr=0.975){
  p_to_pp(1-prob, n, N, alpha=1-thr)
}

#' Estimates the amount of information in an analysis when the patients are
#' not randomized 1:1. Information is given in total patients randomized 1:1 
#' For example, if you have 90 treated patients and 60 controls the effective n
#' is 144 total patients worth of information randomized 1:1
#' @param n_t: number of patients assigned to treatment
#' @param n_c: number of patients assigned to control
#' @return scalar, the total size of an equivalent 1:1 trial
effective_n = function(n_t, n_c){
  4 * (n_t * n_c) / (n_t + n_c)
}

#' Estimates the information of an analysis by comparing to a reference analysis 
#' with known information.  
#' @param n_known: information in the reference analysis
#' @param var_known: the variance of the parameter in the reference analysis.
#' @param var_unknown: the variance of the same parameter analysis of inference.
#' @return scalar, the estimated information in the new analysis.
relative_n = function(n_known, var_known, var_unknown){
  var_known / var_unknown * n_known
}

#' Creates a predictive probability for a future trial from a z-score
#' @param z: z-score 
#' @param n: information used to compute the z-score
#' @param NN: information for the future trial
#' @param alpha: the level of the final analysis (one-sided)
#' @return predictive probability scalar in [0, 1]
z_to_pf = function(z, n, NN, alpha=0.025){
  stopifnot(n>0)
  stopifnot(alpha>=0 & alpha<=1)
  r = NN /(n+NN)
  pnorm(z * sqrt(r) + qnorm(1-alpha) * sqrt(1-r))
}

#' Creates a predictive probability for a future trial from a p-value
#' @param p: p-value
#' @param n: information used to compute the z-score
#' @param NN: information for the future trial
#' @param alpha: the level of the final analysis (one-sided)
#' @return predictive probability scalar in [0, 1]
p_to_pf = function(p, n, NN, alpha=0.025) {
  z = p_to_z(p)
  z_to_pf(z, n, NN, alpha)
}

#-------------------------------------------------------------------------------
# Helper Functions
#-------------------------------------------------------------------------------

#' Inverses of the approximation
pp_to_z = function(pp, n, N, alpha=0.025){
  r = n/N
  qnorm(pp)*sqrt(1-r) + qnorm(1-alpha)*sqrt(r)
}


pp_to_p = function(pp, n, N, alpha=0.025){
  z = pp_to_z(pp, n, N, alpha)
  z_to_p(z)
}

pp_to_prob = function(prob, n, N, thr=0.975){
  1 - pp_to_p(prob, n, N, alpha=1-thr)
}

#' Convert p-value or posterior probability to z-score
p_to_z = function(p) qnorm(1-p)

#' Convert z-score to p-value
z_to_p = function(z) pnorm(-z)


z_to_prob = function(z) 1 - z_to_p(z)


