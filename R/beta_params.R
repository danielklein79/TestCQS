#' Calculate Beta Distribution Parameters for i = 1..(k-1)
#'
#' @description
#' Calculates the parameters `a` and `b` of the beta distribution
#' for the i-th set of eigenvalues (excluding the last one) according to Theorem 3.1.
#'
#' @param i Integer. Index of the eigenvalue set (1 to k-1).
#' @param j Integer. Index of the eigenvalue within the i-th set.
#' @param qq Integer vector. Multiplicities of the eigenvalues (length k).
#' @param qt Numeric vector. (k − 1)-dimensional vector of the counts of the remaining eigenvalues,
##     where the i-th component is the sum of the numbers of eigenvalues from (i + 1) to k.
#' @param n Integer. Sample size.
#'
#' @return Numeric vector of length 2: parameters `a` and `b` for the beta distribution.
#' @keywords internal
#' @noRd
par_ab <- function(i, j, qq, qt, n) {
  #qt_sum <- qt[i]
  qt_sum <- if (length(qt) >= i) qt[i] else 0
  q <- qq[i]
  a <- (n - qt_sum - j) / 2
  b <- ((q + 2) * (j - 1)) / (2 * q) + qt_sum / 2
  return(c(a, b))
}

#' Calculate Beta Distribution Parameters for i = k
#'
#' @description
#' Calculates the parameters `a` and `b` of the beta distribution
#' for the last eigenvalue set (i = k) according to Theorem 3.1.
#'
#' @param j Integer. Index of the eigenvalue within the k-th set.
#' @param qq Integer vector. Multiplicities of the eigenvalues (length k).
#' @param k Integer. Number of distinct eigenvalues.
#' @param n Integer. Sample size.
#'
#' @return Numeric vector of length 2: parameters `a` and `b` for the beta distribution.
#' @keywords internal
#' @noRd
par_ab_k <- function(j, qq, k, n) {
  q <- qq[k]
  a <- (n - 1 - j) / 2
  b <- j * (q + 2) / (2 * q)
  return(c(a, b))
}

#' Characteristic Function of the Exact LRT Distribution
#'
#' @description
#' Computes the characteristic function of the exact distribution
#' of the likelihood ratio test (LRT) statistic as described in Theorem 3.1.
#'
#' @param t Numeric. Argument of the characteristic function.
#' @param qq Integer vector. Multiplicities of the eigenvalues (length k).
#' @param qt Numeric vector. (k − 1)-dimensional vector of the counts of the remaining eigenvalues,
##     where the i-th component is the sum of the numbers of eigenvalues from (i + 1) to k.
#' @param n Integer. Sample size.
#'
#' @return Numeric. Value of the characteristic function at `t`.
#' @keywords internal
#' @noRd
cf_cf <- function(t, qq, qt, n) {
  k <- length(qq)
  cf <- 1

  # Cases for multiple eigenvalues
  if (k > 1) {
  for (i in 1:(k-1)) {
    for (j in 1:qq[i]) {
      pars <- par_ab(i, j, qq, qt, n)
      cf <- cf * CharFunToolR::cf_LogRV_Beta(n * t, pars[1], pars[2], -1)
    }
  }

  # Last component
  if (qq[k] > 1) {
    for (j in 1:(qq[k] - 1)) {
      pars <- par_ab_k(j, qq, k, n)
      cf <- cf * CharFunToolR::cf_LogRV_Beta(n * t, pars[1], pars[2], -1)
    }
  }
  } else if (k == 1){
    # Cases for single eigenvalue, for example structure S
    if (qq[k] > 1) {
      for (j in 1:(qq[k] - 1)) {
        pars <- par_ab_k(j, qq, k, n)
        cf <- cf * CharFunToolR::cf_LogRV_Beta(n * t, pars[1], pars[2], -1)
      }
    }
  }
  return(cf)
}
