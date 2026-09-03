#------------------------------------------------------------
# Print method
#------------------------------------------------------------

#' @title Print a TestCQS object
#'
#' @description
#' Displays the main results of the likelihood ratio and Rao score tests.
#'
#' @param x Object of class \code{"TestCQS"}.
#' @param ... Further arguments (currently ignored).
#'
#' @export
print.TestCQS <- function(x, ...) {

  cat("Call:\n")
  cat(deparse(x$call), sep = "\n")
  cat("\n")

  cat("\n")
  cat("TestCQS: Testing covariance matrix structure\n")
  cat("============================================\n\n")

  cat(sprintf("Null covariance structure: %s\n\n", paste0(structure_labels[x$structure]," (",x$structure,")")))

  cat("Likelihood ratio test\n")
  cat("---------------------\n")
  cat(sprintf("%-25s %12.4f\n", "Statistic", x$lr))
  cat(sprintf("%-25s %12.4s\n", "Asymptotic p-value", pretty_pvalue(x$pval_lrt)))
  cat(sprintf("%-25s %12.4s\n\n", "Exact p-value", pretty_pvalue(x$p_lr)))

  cat("Rao score test\n")
  cat("--------------\n")
  cat(sprintf("%-25s %12.4f\n", "Statistic", x$rs))
  cat(sprintf("%-25s %12.4s\n", "Asymptotic p-value", pretty_pvalue(x$pval_rst)))

  if (is.na(x$p_rs)) {
    cat(sprintf("%-25s %12s\n", "Empirical p-value", "Not computed"))
    cat("\n")
    cat("Note:\n")
    cat("  Empirical RST p-value was not computed.\n")
    cat("  Set 'trials > 0' to estimate it by Monte Carlo simulation.\n")
  } else {
    cat(sprintf("%-25s %12.4s\n", "Empirical p-value", pretty_pvalue(x$p_rs)))
  }

  cat("\n Use summary() for additional information.\n\n")

  invisible(x)
}


#------------------------------------------------------------
# Summary method
#------------------------------------------------------------

#' @title Summary of a TestCQS object
#'
#' @description
#' Produces a summary of the results from \code{runTestCQS()}.
#'
#' @param object Object of class \code{"TestCQS"}.
#' @param ... Further arguments (currently ignored).
#'
#' @return
#' An object of class \code{"summary.TestCQS"}.
#'
#' @export
summary.TestCQS <- function(object, ...) {

  cat("Call:\n")
  print(object$call)
  cat("\n")

  out <- list(

    structure = object$structure,

    n = object$n,
    q = object$q,
    k = object$k,
    qq = object$qq,
    df = object$df,

    lr = object$lr,
    rs = object$rs,

    pval_lrt = object$pval_lrt,
    pval_rst = object$pval_rst,

    p_lr = object$p_lr,
    p_rs = object$p_rs,

    trials = object$trials,
    elapsed_time = object$elapsed_time

  )

  class(out) <- "summary.TestCQS"

  out
}


#------------------------------------------------------------
# Print summary method
#------------------------------------------------------------

#' @title Print summary of a TestCQS object
#'
#' @param x Object of class \code{"summary.TestCQS"}.
#' @param ... Further arguments (currently ignored).
#'
#' @export
print.summary.TestCQS <- function(x, ...) {

  cat("\n")
  cat("Summary of TestCQS\n")
  cat("==================\n\n")

  cat(sprintf("Null covariance structure: %s\n\n", paste0(structure_labels[x$structure]," (",x$structure,")")))

  cat("Sample information\n")
  cat("------------------\n")
  cat(sprintf("%-30s %8d\n", "Sample size", x$n))
  cat(sprintf("%-30s %8d\n", "Matrix dimension", x$q))
  cat(sprintf("%-30s %8d\n", "Distinct eigenvalues", x$k))
  cat(sprintf("%-30s %8d\n\n", "Degrees of freedom", x$df))

  cat("Likelihood ratio test\n")
  cat("---------------------\n")
  cat(sprintf("%-30s %12.4f\n", "Statistic", x$lr))
  cat(sprintf("%-30s %12.4s\n", "Asymptotic p-value", pretty_pvalue(x$pval_lrt)))
  cat(sprintf("%-30s %12.4s\n\n", "Exact p-value", pretty_pvalue(x$p_lr)))

  cat("Rao score test\n")
  cat("--------------\n")
  cat(sprintf("%-30s %12.4f\n", "Statistic", x$rs))
  cat(sprintf("%-30s %12.4s\n", "Asymptotic p-value", pretty_pvalue(x$pval_rst)))

  if (is.na(x$p_rs)) {
    cat(sprintf("%-30s %12s\n\n",
                "Empirical p-value",
                "Not computed"))
  } else {
    cat(sprintf("%-30s %12.4s\n",
                "Empirical p-value",
                pretty_pvalue(x$p_rs)))
    cat(sprintf("%-30s %12d\n",
                "Monte Carlo trials",
                x$trials))
    cat(sprintf("%-30s %12.2f sec\n\n",
                "Computation time",
                x$elapsed_time))
  }

  cat("Eigenvalue multiplicities\n")
  cat("-------------------------\n")
  cat(paste(x$qq, collapse = " "), "\n\n")

  invisible(x)
}
