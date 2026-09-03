# utils.R
# -----------------------------------------------------------------------------
# Helper functions for the TestCQS package (internal use only)
# -----------------------------------------------------------------------------

#' Check dimension input for a dimension of a matrix
#'
#' Validates that the input dimension \code{q} is a positive integer greater than 1.
#' If the input is invalid, the function raises an informative error.
#'
#' @param q A positive integer representing the dimension.
#' @return Invisibly returns \code{q} as an integer if valid; otherwise, throws an error.
#' @keywords internal
#' @noRd
check <- function(q) {
  if (!is.numeric(q) || length(q) != 1L || q <= 1L || q != as.integer(q)) {
    stop("`q` must be an integer dimension greater than 1.", call. = FALSE)
  }
  invisible(as.integer(q))
}

#' Check dimensions input p and m for doubly multivariate model
#'
#' Validates that interclass dimension \code{p} and intraclass dimension \code{m}
#' are positive integers. The combination \code{p = 1, m = 1} is not allowed.
#'
#' @param p Integer. Interclass dimension.
#' @param m Integer. Intraclass dimension.
#' @return Invisibly returns \code{list(p, m)} as integers if valid; otherwise, throws an error.
#' @keywords internal
#' @noRd
check_pm <- function(p, m) {
  # Check p
  if (!is.numeric(p) || length(p) != 1L || p < 1L || p != as.integer(p)) {
    stop("`p` must be a positive integer.", call. = FALSE)
  }

  # Check m
  if (!is.numeric(m) || length(m) != 1L || m < 1L || m != as.integer(m)) {
    stop("`m` must be a positive integer.", call. = FALSE)
  }

  # Forbid combination p = 1, m = 1
  if (p == 1L && m == 1L) {
    stop("The combination p = 1 and m = 1 is not allowed.", call. = FALSE)
  }

  # Return as integers
  invisible(list(p = as.integer(p), m = as.integer(m)))
}


#' Warn about constructing a large dense matrix
#'
#' Issues a warning when a dense matrix of dimension \code{q x q} is about
#' to be constructed. This is intended to alert users that the resulting
#' allocation may require a substantial amount of memory.
#'
#' @param q A positive integer representing the matrix dimension.
#'
#' @return Invisibly returns \code{NULL}. This function is called for its
#'   side effect of issuing a warning when appropriate.
#'
#' @keywords internal
#' @noRd
.warn_large_dense <- function(q) {
  if (q > 1e6L) {
    warning(
      sprintf(
        "Creating a %d x %d matrix. This may require substantial memory.",
        q, q
      ),
      call. = FALSE
    )
  }

  invisible(NULL)
}



#' Orthogonal projection matrix onto the orthogonal complement of the
#' column space of an n-dimensional vector of ones
#'
#' @param n A positive integer. Dimension of the identity matrix.
#' @return An n x n orthogonal projection matrix.
#' @keywords internal
#' @noRd
Q1 <- function(n){
  check(n)

  diag(n) - matrix(1/n, n, n)
}

#' Trace operator
#'
#' Computes the trace of a square matrix (sum of diagonal elements)
#'
#' @param A A square numeric matrix.
#' @return The trace of a matrix A.
#' @keywords internal
#' @noRd
tr <- function(A) {
  # Check is inout is a square matrix
  if (!is.matrix(A) || nrow(A) != ncol(A)) {
    stop("`A` must be a square matrix.", call. = FALSE)
  }

  sum(diag(A))
}

#' Eigenvectors of circular Toeplitz (as well as compound symmetry) structure
#'
#' Constructs the orthogonal matrix of eigenvectors for a circular Toeplitz
#' (CT) or compound symmetry (CS) covariance structure.
#'
#' @param q A positive integer. Dimension of the matrix.
#' @param CTorder logical; if TRUE, order eigenvectors by eigenvalue multiplicity:
#' (normalized constant vector first (multiplicity 1), followed by pairs of vectors (multiplicity 2).
#' In case of even dimension finished by final single vector (multiplicity 1).
#' If FALSE, eigenvectors are ordered according to the positions of equal elements
#' in the circular Toeplitz structure.
#' @return An n x n matrix of orthonormal eigenvectors.
#' @keywords internal
#' @noRd
CTEigvec <- function(q, CTorder = FALSE) {
  q<-check(q)

  s <- seq.int(0L, q - 1L)

  M <- outer(s, s, `*`) * (2 * pi / q)
  M <- (cos(M) + sin(M)) / sqrt(q)

  if (CTorder) {
    if (q %% 2L == 1L) {  # odd dimension
      k <- (q - 1L) %/% 2L
      if(k > 0L) {
        idx <- c(
          1L,
          as.vector(rbind(
            seq.int(2L, k + 1L),
            seq.int(q, q - k + 1L)
          ))
        )
      } else {
        idx <- c(1L)
      }
    } else {              # even dimension
      k <- q %/% 2L - 1L
      if (k > 0L) {
        idx <- c(
          1L,
          as.vector(rbind(
            seq.int(2L, k + 1L),
            seq.int(q, q - k + 1L)
          )),
          q %/% 2L + 1L
        )
      } else {
        idx <- c(1L, q %/% 2L + 1L)
      }
    }
    M[, idx, drop = FALSE]
  } else {
    M
  }
}


#' @title Format p-values for printing
#'
#' @description
#' Formats p-values for display in printed output. Values smaller than a
#' specified threshold are reported as being less than that threshold,
#' while missing values are returned as \code{"NA"}.
#'
#' @param p A numeric p-value.
#' @param digits Integer specifying the number of decimal places to display
#'   for p-values greater than or equal to \code{eps}. Default is \code{4}.
#' @param eps Numeric threshold below which p-values are displayed as
#'   \code{"<eps"}. Default is \code{1e-4}.
#'
#' @return
#' A character string representing the formatted p-value.
#'
#' @keywords internal
#' @noRd
pretty_pvalue <- function(p, digits = 4, eps = 1e-4) {

  if (is.na(p))
    return("Not computed")

  if (p < eps)
    return(paste0("<", format(eps, scientific = FALSE)))

  formatC(p, format = "f", digits = digits)
}

