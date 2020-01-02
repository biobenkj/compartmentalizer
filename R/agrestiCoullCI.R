#' Normal alpha/2 quantile
#' 
#' @param q   the quantile at which to extract Z
#' 
#' @return    Z
#'
.z <- function(q) qnorm(1 - ((1-q)/2))

#' n_tilde in AC
#'
#' @param n1  number of successes/ones 
#' @param n0  number of failures/zeroes
#' @param q   quantile for eventual CI (e.g. 0.95 for a 95 percent binomial CI)
#'
#' @return    the effective sample size for smoothed CIs
#'
.n_approx <- function(n1, n0, q) (n0 + n1) + (.z(q)^2)

#' p_tilde in AC 
#'
#' @param n1  number of successes/ones 
#' @param n0  number of failures/zeroes
#' @param q   quantile for eventual CI (e.g. 0.95 for a 95 percent binomial CI)
#'
#' @return    the approximate success probability for a smoothed CIs
#' 
.p_approx <- function(n1, n0, q) (1/.n_approx(n1, n0, q)) * (n1 + ((.z(q)^2)/2))

#' Agresti-Coull confidence interval for a binomial proportion
#' 
#' @param n1  number of successes/ones 
#' @param n0  number of failures/zeroes
#' @param q   quantile for eventual CI (e.g. 0.95 for a 95 percent binomial CI)
#'
#' @return    the approximate (q x 100) percent confidence interval for (p|n1,n0,q)
#' @export
#' 
#' @examples 
#' binom.ci <- agrestiCoullCI(10, 3, 0.95)
#' 
agrestiCoullCI <- function(n1, n0, q) { 
  p_apx <- .p_approx(n1, n0, q)
  n_apx <- .n_approx(n1, n0, q)
  width <- .z(q) * sqrt( (p_apx/n_apx) * (1 - p_apx) )
  res <- c(low=pmax(0, (p_apx - width)),est=p_apx,high=pmin((p_apx + width), 1))
  return(res) 
}
