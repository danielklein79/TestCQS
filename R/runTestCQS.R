#' @title Testing covariance matrix belonging to commutative quadratic subspace (CQS)
#'
#' @description
#' This function tests the commutative quadratic subspace (CQS) structure under
#' a multivariate or, more generally, elliptically contoured models. It computes the likelihood ratio test (LRT)
#' and Rao score test (RST) statistics, along with their p-values based on
#' chi-squared approximation, exact distribution for LRT, and empirical
#' distribution for RST via simulation. Additionally, the function provides the creation
#' of graphs of the exact distribution of the LRT statistic, the empirical RST and
#' the asymptotic chi-square distribution.
#'
#' @param X Numeric matrix of observations (rows = samples, columns = variables).
#' @param p Integer. Interclass dimension (used only for doubly multivariate structures).
#' @param m Integer. Intraclass dimension (used only for doubly multivariate structures).
#' @param qq Vector of positive integers representing multiplicities of eigenvalues.
#'          If we consider one of the structures S, D, CS, CT, BS_D, BS_CS, BS_CT, BD_S, BD_CS,
#'          BD_CT, BCS_S, BCS_D, BCS_CS, BCS_CT, BCT_S, BCT_D, BCT_CS, BCT_CT, UBS_CS, DCS
#'          the argument should be left empty.
#' @param H Either a string naming a predefined structure (S, D, CS, CT, BS_D, BS_CS, BS_CT,
#'          BD_S, BD_CS, BD_CT, BCS_S, BCS_D, BCS_CS, BCS_CT, BCT_S, BCT_D, BCT_CS, BCT_CT,UBS_CS,DCS)
#'          or a numeric matrix representing the base matrix. If a string is given,
#'          the corresponding `H` matrix and `qq` vector are automatically generated.
#' @param trials Integer. Number of Monte Carlo simulations used to estimate
#'          the empirical p-value of the Rao score test. If \code{trials = 0}
#'          (default), no simulations are performed and the empirical p-value
#'          is returned as \code{NA}. Larger values improve the Monte Carlo
#'          accuracy but may substantially increase computation time.
#' @param show_progress Logical; whether to display a progress indicator.
#'   Defaults to `interactive()`.
#' @param plotsPDF_CDF Logical. If TRUE, plots the PDF/CDF comparison between chi-squared
#'                     and exact LRT distribution (default FALSE).
#' @param set_seed Integer. Seed for random number generation to ensure reproducibility
#'                 (default 123443).
#'
#' @return A list of class \code{TestCQS} with components:
#' \describe{
#'   \item{lr}{Likelihood ratio test statistic.}
#'   \item{rs}{Rao score test statistic.}
#'   \item{pval_lrt}{Chi-squared p-value for LRT.}
#'   \item{pval_rst}{Chi-squared p-value for RST.}
#'   \item{p_lr}{Exact p-value for LRT based on characteristic function inversion.}
#'   \item{p_rs}{Empirical Monte Carlo p-value for the Rao score test.
#'               Returned only when \code{trials > 0}; otherwise \code{NA}.}
#' }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Validates input arguments.
#'   \item Determines the base matrix \code{H} and multiplicities \code{qq}
#'         based on the specified structure.
#'   \item Computes the maximum likelihood estimate of the covariance
#'         matrix under null (H0) and alternative (H1) hypotheses.
#'   \item Calculates LRT and RST statistics and their approximate
#'         chi-squared p-values.
#'   \item Computes exact p-value for LRT using characteristic function inversion.
#'   \item Optionally, estimates the empirical p-value of the Rao score
#'         statistic by Monte Carlo simulation when \code{trials > 0}.
#'         Since this step repeatedly generates multivariate normal
#'         samples and recomputes the test statistic, it may be computationally
#'         intensive for large values of \code{trials}, sample sizes, or dimensions.
#'   \item Optionally, plots the chi-squared and exact PDF of the LRT statistic.
#' }
#'
#' @examples
#' X <- matrix(rnorm(100*5), nrow=100, ncol=5)
#' result <- runTestCQS(X, H="S", trials=100)
#' print(result)
#'
#' @export
runTestCQS <- function(X, p=0, m=0, qq=NULL, H, trials=0, show_progress = interactive(), plotsPDF_CDF=FALSE, set_seed=123443){
  call = match.call()
  package_version <- utils::packageVersion("TestCQS")
  #print("Testing the commutative quadratic subspace structure under a multivariate model")
  #VALIDATION
  ##########
  ##########
  if (!is.numeric(p) || length(p) != 1 || p < 0 || p %% 1 != 0)
    stop("Argument 'p' must be a non-negative integer.")

  if (!is.matrix(X))
    stop("Argument 'X' must be a matrix.")

  if (!is.numeric(X))
    stop("Matrix X must contain numeric values.")

  if (!is.numeric(trials) || length(trials) != 1 || trials < 0 || trials %% 1 != 0)
    stop("Argument 'trials' must be a non-negative integer.")

  if (!is.logical(plotsPDF_CDF) || length(plotsPDF_CDF) != 1)
    stop("Argument 'plotsPDF_CDF' must be TRUE or FALSE.")

  if (!is.numeric(set_seed) || length(set_seed) != 1 || set_seed %% 1 != 0)
    stop("Argument 'set_seed' must be an integer.")
  p = as.numeric(p)
  set.seed(set_seed)
  q = ncol(X)

  # Check H argument and store the user-specified structure for printing
  if (is.character(H)) {

    if (length(H) != 1L || is.na(H)) {
      stop("'H' must be a single character string or a square orthogonal matrix.")
    }

    if (H %in% ordinary_structures) {

      structure_name <- H
      func_name <- paste0("basis_", H)

      fun <- do.call(get(func_name), list(q = q))
      H <- fun$H
      qq <- fun$qq

    } else if (H %in% doubly_structures) {

      if (!is.numeric(m) || length(m) != 1 || m < 0 || m %% 1 != 0)
        stop("Argument 'm' must be a non-negative integer.")
      m = as.numeric(m)

      structure_name <- H
      func_name <- paste0("basis_", H)

      fun <- do.call(get(func_name), list(p = p, m = m))
      H <- fun$H
      qq <- fun$qq

    } else if (H %in% UBS_structure) {

      if (!is.numeric(m) || length(m) <2 || any(m) < 0 || any(m) %% 1 != 0 || anyNA(m))
        stop("Argument 'm' must be a vector of non-negative integers of length at least 2.")
      m = as.numeric(m)

      structure_name <- H
      func_name <- paste0("basis_", H)

      fun <- do.call(get(func_name), list(p = p, m = m))
      H <- fun$H
      qq <- fun$qq

    } else {

      stop(
        "Unknown structure '", H, "'. ",
        "Choose one of the supported structure names or provide ",
        "a square orthogonal matrix."
      )
    }

  } else if (is.matrix(H)) {

    structure_name <- "User-defined"

    if (!is.numeric(H)) {
      stop("When 'H' is a matrix, it must be numeric.")
    }

    if (nrow(H) != ncol(H)) {
      stop("When 'H' is a matrix, it must be square.")
    }

    if (anyNA(H)) {
      stop("When 'H' is a matrix, it cannot contain missing values.")
    }

    # Check orthogonality
    if (!isTRUE(all.equal(
      crossprod(H),
      diag(ncol(H)),
      tolerance = 1e-8
    ))) {
      stop("When 'H' is a matrix, it must be orthogonal.")
    }

    # For a user-supplied H, qq must be supplied
    if (is.null(qq)) {
      stop("'qq' must be supplied when 'H' is a user-defined matrix.")
    }

  } else {

    stop(
      "'H' must be either a supported structure name ",
      "or a square orthogonal numeric matrix."
    )
  }

  # Check of qq vector
  if (!is.numeric(qq) ||
      length(qq) == 0L ||
      anyNA(qq) ||
      any(qq <= 0) ||
      any(qq %% 1 != 0) ||
      sum(qq) != ncol(X)) {
    stop("'qq' must be a vector of positive integers, which sums to model dimension.")
  }

  qq <- as.integer(qq)
  k = length(qq)
  q = sum(qq)
  n = nrow(X)

  mu = rep(0,q)
  S = t(X)%*%Q1(n)%*%X
  a = 1
  b = qq[1]

  # Calculating the MLE of Sigma under H0
  Sig0MLE = matrix(0, nrow = q, ncol = q)

  for (j in 1:k) {
    V=H[,a:b]%*%t(H[,a:b])
    lambda_j = tr(V%*%S)/(n*qq[j])
    Sig0MLE = Sig0MLE + lambda_j*V
    if (j<k){
      a = sum(qq[1:j])+1
      b = sum(qq[1:(j+1)])
    }
  }

  Sig1MLE = S/n

  # Calculating the test statistics
  lr = -2*log((det(Sig1MLE)/det(Sig0MLE))^(n/2))
  rs = n/2*tr((Sig1MLE%*%solve(Sig0MLE)-diag(q))%*%(Sig1MLE%*%solve(Sig0MLE)-diag(q)))

  #Calculating the asymptotical p-values
  df = q*(q+1)/2 - k
  pval_rst = 1 - stats::pchisq(rs,df)
  pval_lrt = 1 - stats::pchisq(lr,df)

  if (length(qq) > 1) {
    qt <- rev(cumsum(rev(qq)))[-1]
  } else {
    qt <- numeric(0)
  }

  # Calculating the exact p-value for LRT
  x_upper <- max(stats::qchisq(0.9999, df), lr * 1.5, 50)
  points  <- seq(0, x_upper, length.out = 2000)
  prob <- seq(0.0005, 0.9995, by = 0.0005)

  options <- list(
    isPlot = plotsPDF_CDF,
    xMin   = 0,
    xMax   = x_upper
  )

  attempt <- 1
  repeat {
    result_try <- try(
      CharFunToolR::cf2DistGP(
        cf = function(t) cf_cf(t, qq, qt, n),
        x  = points,
        prob = prob,
        options = options
      ),
      silent = TRUE
    )

    if (!inherits(result_try, "try-error")) {
      result <- result_try
      break
    } else {
      attempt <- attempt + 1
      if (attempt > 6) {
        warning("cf2DistGP failed after multiple attempts; skipping exact LRT distribution.")
        result <- NA
        break
      }
      options$xMax <- options$xMax * 2
      points <- seq(0, options$xMax, length.out = min(12000, length(points) * 2))
    }
  }

  p_lr <- NA_real_
  if (!is.na(result)[1]) {
    if (lr > max(result$x)) {
      warning(sprintf("LR=%.6f exceeds computed support max=%.6f; using edge CDF (conservative).", lr, max(result$x)))
    }
    p_lr <- 1 - stats::approx(x = result$x, y = result$cdf,
                     xout = min(lr, max(result$x)),
                     rule = 2)$y
  }

  # Calculating the empirical p-value for RST
  p_rs <- elapsed.time <- NA_real_
  if (trials>0){
    message(
      sprintf("Estimating empirical RST p-value using %d Monte Carlo simulations...", trials))
    if (show_progress) {
      pb <- utils::txtProgressBar(min = 0, max = trials, style = 3)
    }
    A = matrix(stats::rnorm(q^2), q, q)
    svd_decomp = svd(A)
    H = svd_decomp$u  #Orthonormal columns

    #Creating a matrix Lambda (eigenvalues)
    eigen_Val = 4*stats::runif(k) + 1
    diagonal = rep(eigen_Val, times=qq)
    Lam = diag(diagonal)

    Sig0 = H%*%Lam%*%t(H)
    Sig0

    RST = c()
    start.time <- Sys.time()
    for (i in 1:trials) {
      X = MASS::mvrnorm(n,mu,Sig0)
      S = t(X)%*%Q1(n)%*%X
      a = 1
      b = qq[1]
      Sig0MLE = matrix(0, nrow = q, ncol = q)

      for (j in 1:k) {
        V = H[,a:b]%*%t(H[,a:b])
        lambda_j = tr(V%*%S)/(n*qq[j])
        Sig0MLE = Sig0MLE+lambda_j*V
        if (j<k){
          a = sum(qq[1:j])+1
          b = sum(qq[1:(j+1)])
        }
      }

      Sig1MLE = S/n

      RST[i] = n/2*tr((Sig1MLE%*%solve(Sig0MLE)-diag(q))%*%(Sig1MLE%*%solve(Sig0MLE)-diag(q)))

      if (show_progress) utils::setTxtProgressBar(pb, i)
    }

    if (show_progress) close(pb)

    elapsed.time <- as.numeric(difftime(Sys.time(), start.time, units = "secs"))

    p_rs = length(RST[RST>rs])/trials
  }

  #Creating the plots
  if (plotsPDF_CDF == TRUE){
    xseq = seq(0, max(points), length.out = 1000)

    dens_chisq = stats::dchisq(xseq, df)
    dens_exact = result$pdf
    ymax = max(c(dens_chisq, dens_exact))

    graphics::plot(xseq, dens_chisq, type = "l", lwd = 2, col = "blue",
         ylim = c(0, 1.1*ymax),
         main = paste0("Chi-squared vs exact distribution of LRT"),
         xlab = "", ylab = "")

    graphics::lines(points, dens_exact, col = "red", lwd = 2)
    graphics::grid(col = "gray", lty = "dotted")
    graphics::legend("topright", legend = c("Chi-squared", "Exact LRT"),
           col = c("blue", "red"), lwd = 2)
  }

  #Creating output
  res <- list(
    # Call and package version
    call = call,
    package_version = package_version,

    # Test statistics
    lr = lr,
    rs = rs,

    # p-values
    pval_lrt = pval_lrt,
    pval_rst = pval_rst,
    p_lr = p_lr,
    p_rs = p_rs,

    # Model information
    n = n,
    q = q,
    k = k,
    df = df,
    qq = qq,

    # Additional information
    trials = trials,
    structure = structure_name,
    elapsed_time = elapsed.time
  )

  class(res) = "TestCQS"
  return(res)
}
