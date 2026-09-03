# --------------------------------
# ORDINARY MULTIVARIATE STRUCTURES
# --------------------------------

#' @title Basis Matrix for a sphericity structure
#'
#' @description
#' Generates the basis matrix of eigenvectors and a vector of
#' eigenvalue multiplicities for an ordinary multivariate model with a sphericity
#' covariance structure.
#'
#' @param q Integer, dimension of the matrix (q > 1).
#' @param sparse Logical; if TRUE, returns a sparse identity matrix.
#'
#' @return A list with two components:
#'   \item{H}{Basis matrix containing eigenvectors.}
#'   \item{qq}{Integer vector of eigenvalue multiplicities. Each \code{qq[i]}
#'             specifies how many columns of \code{H} belong to the i-th
#'             distinct eigenvalue.}
#'
#' @details
#' The output defines the spectral decomposition of covariance matrix \eqn{\Sigma}
#' in the form \eqn{\Sigma = H \Lambda H^T} where \eqn{\Lambda} is a diagonal matrix of
#' eigenvalues, with eigenvalue \eqn{\lambda_i} having multiplicity \eqn{qq_i}.
#'
#' @export
basis_S <- function(q, sparse = FALSE) {
  q <- check(q)

  H <- if (sparse) Matrix::Diagonal(q) else diag(q)

  list(
    H  = H,
    qq = q
  )
}

#' @title Basis Matrix for a diagonal structure
#'
#' @description
#' Generates the basis matrix of eigenvectors and a vector of
#' eigenvalue multiplicities for an ordinary multivariate model with a diagonal
#' covariance structure.
#'
#' @param q Integer, dimension of the matrix (q > 1).
#' @param sparse Logical; if TRUE, returns a sparse identity matrix.
#'
#' @return A list with two components:
#'   \item{H}{Basis matrix containing eigenvectors.}
#'   \item{qq}{Integer vector of eigenvalue multiplicities. Each \code{qq[i]}
#'             specifies how many columns of \code{H} belong to the i-th
#'             distinct eigenvalue.}
#'
#' @details
#' The output defines the spectral decomposition of covariance matrix \eqn{\Sigma}
#' in the form \eqn{\Sigma = H \Lambda H^T} where \eqn{\Lambda} is a diagonal matrix of
#' eigenvalues, with eigenvalue \eqn{\lambda_i} having multiplicity \eqn{qq_i}.
#'
#' @export
basis_D <- function(q, sparse = FALSE) {
  q <- check(q)

  H <- if (sparse) Matrix::Diagonal(q) else diag(q)

  list(
    H  = H,
    qq = rep.int(1L, q)
  )
}

#' @title Basis Matrix for a compound symmetry structure
#'
#' @description
#' Generates the basis matrix of eigenvectors and a vector of
#' eigenvalue multiplicities for an ordinary multivariate model with a compound symmetry
#' covariance structure.
#'
#' @param q Integer, dimension of the matrix (q > 1).
#'
#' @return A list with two components:
#'   \item{H}{Basis matrix containing eigenvectors.}
#'   \item{qq}{Integer vector of eigenvalue multiplicities. Each \code{qq[i]}
#'             specifies how many columns of \code{H} belong to the i-th
#'             distinct eigenvalue.}
#'
#' @details
#' The output defines the spectral decomposition of covariance matrix \eqn{\Sigma}
#' in the form \eqn{\Sigma = H \Lambda H^T} where \eqn{\Lambda} is a diagonal matrix of
#' eigenvalues, with eigenvalue \eqn{\lambda_i} having multiplicity \eqn{qq_i}.
#'
#' @export
basis_CS <- function(q) {
  q <- check(q)

  .warn_large_dense(q)

  H <- CTEigvec(q)

  list(
    H  = H,
    qq = c(1L, q - 1L)
  )
}

#' @title Basis Matrix for a circular Toeplitz structure
#'
#' @description
#' Generates the basis matrix of eigenvectors and a vector of
#' eigenvalue multiplicities for an ordinary multivariate model with a circular Toeplitz
#' covariance structure.
#'
#' @param q Integer, dimension of the matrix (q > 1).
#'
#' @return A list with two components:
#'   \item{H}{Basis matrix containing eigenvectors.}
#'   \item{qq}{Integer vector of eigenvalue multiplicities. Each \code{qq[i]}
#'             specifies how many columns of \code{H} belong to the i-th
#'             distinct eigenvalue.}
#'
#' @details
#' The output defines the spectral decomposition of covariance matrix \eqn{\Sigma}
#' in the form \eqn{\Sigma = H \Lambda H^T} where \eqn{\Lambda} is a diagonal matrix of
#' eigenvalues, with eigenvalue \eqn{\lambda_i} having multiplicity \eqn{qq_i}.
#'
#' @export
basis_CT <- function(q) {
  q <- check(q)

  .warn_large_dense(q)

  H <- CTEigvec(q, CTorder = TRUE)

  if (q %% 2L == 1L) {
    k <- (q - 1L) %/% 2L
    qq  <- c(1L, rep.int(2L, k))
  } else {
    k <- q %/% 2L - 1L
    qq  <- c(1L, rep.int(2L, k), 1L)
  }

  list(
    H  = H,
    qq = qq
  )
}

# -------------------------------
# DOUBLY MULTIVARIATE STRUCTURES
# -------------------------------

#' @title Basis Matrix for BS_D structure
#'
#' @description
#' Generates the basis matrix of eigenvectors and a vector of
#' eigenvalue multiplicities for a doubly multivariate model with
#' block-sphericity structure, where each block has a diagonal structure.
#'
#' @param p Integer, between-class dimension (must be >= 1).
#' @param m Integer, within-class dimension (must be >= 1).
#' @param sparse Logical; if TRUE, returns a sparse identity matrix.
#'
#' @return A list with two components:
#'   \item{H}{Basis matrix containing eigenvectors.}
#'   \item{qq}{Integer vector of eigenvalue multiplicities. Each \code{qq[i]}
#'             specifies how many columns of \code{H} belong to the i-th
#'             distinct eigenvalue.}
#'
#' @details
#' The output defines the spectral decomposition of covariance matrix \eqn{\Sigma}
#' in the form \eqn{\Sigma = H \Lambda H^T} where \eqn{\Lambda} is a diagonal matrix of
#' eigenvalues, with eigenvalue \eqn{\lambda_i} having multiplicity \eqn{qq_i}.
#'
#' @export
basis_BS_D <- function(p, m, sparse = FALSE) {
  params <- check_pm(p, m)
  p <- params$p
  m <- params$m

  idx <- m * seq.int(0L, p - 1L) + rep(seq_len(m), each = p)

  H <- if (sparse) Matrix::Diagonal(p * m) else diag(p * m)

  list(
    H = H[, idx, drop = FALSE],
    qq = rep.int(p, m)
  )
}

#' @title Basis Matrix for BS_CS structure
#'
#' @description
#' Generates the basis matrix of eigenvectors and a vector of
#' eigenvalue multiplicities for a doubly multivariate model with
#' block-sphericity structure, where each block has a compound symmetry structure.
#'
#' @param p Integer, between-class dimension (must be >= 1).
#' @param m Integer, within-class dimension (must be >= 1).
#' @param sparse Logical; if TRUE, returns a sparse identity matrix.
#'
#' @return A list with two components:
#'   \item{H}{Basis matrix containing eigenvectors.}
#'   \item{qq}{Integer vector of eigenvalue multiplicities. Each \code{qq[i]}
#'             specifies how many columns of \code{H} belong to the i-th
#'             distinct eigenvalue.}
#'
#' @details
#' The output defines the spectral decomposition of covariance matrix \eqn{\Sigma}
#' in the form \eqn{\Sigma = H \Lambda H^T} where \eqn{\Lambda} is a diagonal matrix of
#' eigenvalues, with eigenvalue \eqn{\lambda_i} having multiplicity \eqn{qq_i}.
#'
#' @export
basis_BS_CS <- function(p, m, sparse = FALSE) {
  params <- check_pm(p, m)
  p <- params$p
  m <- params$m

  L <- if (sparse) Matrix::Diagonal(p) else diag(p)
  M <- CTEigvec(m, CTorder = TRUE)

  H <- kronecker(L, M)

  firsts <- (0L:(p-1L))*m + 1L
  idx <- c(firsts, setdiff(1L:(p*m), firsts))

  list(
    H = H[, idx, drop = FALSE],
    qq = c(p, p*(m - 1L))
  )
}

#' @title Basis Matrix for BS_CT structure
#'
#' @description
#' Generates the basis matrix of eigenvectors and a vector of
#' eigenvalue multiplicities for a doubly multivariate model with
#' block-sphericity structure, where each block has a circular Toeplitz structure.
#'
#' @param p Integer, between-class dimension (must be >= 1).
#' @param m Integer, within-class dimension (must be >= 1).
#' @param sparse Logical; if TRUE, returns a sparse identity matrix.
#'
#' @return A list with two components:
#'   \item{H}{Basis matrix containing eigenvectors.}
#'   \item{qq}{Integer vector of eigenvalue multiplicities. Each \code{qq[i]}
#'             specifies how many columns of \code{H} belong to the i-th
#'             distinct eigenvalue.}
#'
#' @details
#' The output defines the spectral decomposition of covariance matrix \eqn{\Sigma}
#' in the form \eqn{\Sigma = H \Lambda H^T} where \eqn{\Lambda} is a diagonal matrix of
#' eigenvalues, with eigenvalue \eqn{\lambda_i} having multiplicity \eqn{qq_i}.
#'
#' @export
basis_BS_CT <- function(p, m, sparse = FALSE) {
  params <- check_pm(p, m)
  p <- params$p
  m <- params$m

  offsets <- if (p >1L) {
    m * seq.int(0L, p - 1L)
  } else {
    integer(0)
  }

  pair_starts <- if (m > 2L) {
    seq.int(2L, m - 1L, by = 2L)
  } else {
    integer(0)
  }

  idx <- 1L + offsets
  np <- length(pair_starts)

  if (np > 0) {
    rep_offsets <- rep(offsets, times = np)

    B1 <-  rep(pair_starts, each = p) + rep_offsets
    B2 <-  rep(pair_starts + 1L, each = p) + rep_offsets

    rest <- as.vector(rbind(B1, B2))
  } else {
    rest <- integer(0)
  }

  if ((m %% 2L) == 1L) {
    idxp <- c(idx, rest)
    half_m <- (m - 1L) %/% 2L
    qq  <- c(p, rep(2L * p, half_m))
  } else {
    idxp <- c(idx, rest, idx + m - 1L)
    half_m_minus_1 <- (m - 2L) %/% 2L
    qq  <- c(p, rep(2L*p, half_m_minus_1), p)
  }

  L <- if (sparse) Matrix::Diagonal(p) else diag(p)
  M <- CTEigvec(m, CTorder = TRUE)

  H <- kronecker(L, M)

  list(
    H = H[, idxp, drop=FALSE],
    qq = qq
  )
}

#' @title Basis Matrix for BD_S structure
#'
#' @description
#' Generates the basis matrix of eigenvectors and a vector of
#' eigenvalue multiplicities for a doubly multivariate model with
#' block-diagonal structure, where each block has a sphericity structure.
#'
#' @param p Integer, between-class dimension (must be >= 1).
#' @param m Integer, within-class dimension (must be >= 1).
#' @param sparse Logical; if TRUE, returns a sparse identity matrix.
#'
#' @return A list with two components:
#'   \item{H}{Basis matrix containing eigenvectors.}
#'   \item{qq}{Integer vector of eigenvalue multiplicities. Each \code{qq[i]}
#'             specifies how many columns of \code{H} belong to the i-th
#'             distinct eigenvalue.}
#'
#' @details
#' The output defines the spectral decomposition of covariance matrix \eqn{\Sigma}
#' in the form \eqn{\Sigma = H \Lambda H^T} where \eqn{\Lambda} is a diagonal matrix of
#' eigenvalues, with eigenvalue \eqn{\lambda_i} having multiplicity \eqn{qq_i}.
#'
#' @export
basis_BD_S <- function(p, m, sparse = FALSE) {
  params <- check_pm(p, m)
  p <- params$p
  m <- params$m

  H <- if (sparse) Matrix::Diagonal(p * m) else diag(p * m)

  list(
    H = H,
    qq = rep.int(m, p)
  )
}

#' @title Basis Matrix for BD_CS structure
#'
#' @description
#' Generates the basis matrix of eigenvectors and a vector of
#' eigenvalue multiplicities for a doubly multivariate model with
#' block-diagonal structure, where each block has a compound symmetry structure.
#'
#' @param p Integer, between-class dimension (must be >= 1).
#' @param m Integer, within-class dimension (must be >= 1).
#' @param sparse Logical; if TRUE, returns a sparse identity matrix.
#'
#' @return A list with two components:
#'   \item{H}{Basis matrix containing eigenvectors.}
#'   \item{qq}{Integer vector of eigenvalue multiplicities. Each \code{qq[i]}
#'             specifies how many columns of \code{H} belong to the i-th
#'             distinct eigenvalue.}
#'
#' @details
#' The output defines the spectral decomposition of covariance matrix \eqn{\Sigma}
#' in the form \eqn{\Sigma = H \Lambda H^T} where \eqn{\Lambda} is a diagonal matrix of
#' eigenvalues, with eigenvalue \eqn{\lambda_i} having multiplicity \eqn{qq_i}.
#'
#' @export
basis_BD_CS <- function(p, m, sparse = FALSE) {
  params <- check_pm(p, m)
  p <- params$p
  m <- params$m


  L <- if (sparse) Matrix::Diagonal(p) else diag(p)
  M <- CTEigvec(m)

  H <- kronecker(L, M)

  list(
    H  = H,
    qq = rep.int(c(1L, m - 1L), p)
  )
}

#' @title Basis Matrix for BD_CT structure
#'
#' @description
#' Generates the basis matrix of eigenvectors and a vector of
#' eigenvalue multiplicities for a doubly multivariate model with
#' block-diagonal structure, where each block has a circular Toeplitz structure.
#'
#' @param p Integer, between-class dimension (must be >= 1).
#' @param m Integer, within-class dimension (must be >= 1).
#' @param sparse Logical; if TRUE, returns a sparse identity matrix.
#'
#' @return A list with two components:
#'   \item{H}{Basis matrix containing eigenvectors.}
#'   \item{qq}{Integer vector of eigenvalue multiplicities. Each \code{qq[i]}
#'             specifies how many columns of \code{H} belong to the i-th
#'             distinct eigenvalue.}
#'
#' @details
#' The output defines the spectral decomposition of covariance matrix \eqn{\Sigma}
#' in the form \eqn{\Sigma = H \Lambda H^T} where \eqn{\Lambda} is a diagonal matrix of
#' eigenvalues, with eigenvalue \eqn{\lambda_i} having multiplicity \eqn{qq_i}.
#'
#' @export
basis_BD_CT <- function(p, m, sparse=FALSE) {
  params <- check_pm(p, m)
  p <- params$p
  m <- params$m

  if ((m %% 2L) == 1L) {
    qq <- rep(c(1L, rep(2L, (m - 1L) %/% 2L)), p)
  } else {
    qq <- rep(c(1L, rep(2L ,(m - 2L) %/% 2L), 1L), p)
  }

  L <- if (sparse) Matrix::Diagonal(p) else diag(p)
  M <- CTEigvec(m, CTorder = TRUE)

  H <- kronecker(L, M)

  list(
    H = H,
    qq = qq
  )
}

#' @title Basis Matrix for BCS_S structure
#'
#' @description
#' Generates the basis matrix of eigenvectors and a vector of
#' eigenvalue multiplicities for a doubly multivariate model with
#' block compound symmetry structure, where each block has a sphericity structure.
#'
#' @param p Integer, between-class dimension (must be >= 1).
#' @param m Integer, within-class dimension (must be >= 1).
#' @param sparse Logical; if TRUE, returns a sparse identity matrix.
#'
#' @return A list with two components:
#'   \item{H}{Basis matrix containing eigenvectors.}
#'   \item{qq}{Integer vector of eigenvalue multiplicities. Each \code{qq[i]}
#'             specifies how many columns of \code{H} belong to the i-th
#'             distinct eigenvalue.}
#'
#' @details
#' The output defines the spectral decomposition of covariance matrix \eqn{\Sigma}
#' in the form \eqn{\Sigma = H \Lambda H^T} where \eqn{\Lambda} is a diagonal matrix of
#' eigenvalues, with eigenvalue \eqn{\lambda_i} having multiplicity \eqn{qq_i}.
#'
#' @export
basis_BCS_S <- function(p, m, sparse=FALSE) {
  params <- check_pm(p, m)
  p <- params$p
  m <- params$m

  L <- CTEigvec(p)
  M <- if (sparse) Matrix::Diagonal(m) else diag(m)

  H <- kronecker(L, M)

  list(
    H = H,
    qq = c(m, (p - 1L) * m)
  )
}

#' @title Basis Matrix for BCS_D structure
#'
#' @description
#' Generates the basis matrix of eigenvectors and a vector of
#' eigenvalue multiplicities for a doubly multivariate model with
#' block compound symmetry structure, where each block has a diagonal structure.
#'
#' @param p Integer, between-class dimension (must be >= 1).
#' @param m Integer, within-class dimension (must be >= 1).
#' @param sparse Logical; if TRUE, returns a sparse identity matrix.
#'
#' @return A list with two components:
#'   \item{H}{Basis matrix containing eigenvectors.}
#'   \item{qq}{Integer vector of eigenvalue multiplicities. Each \code{qq[i]}
#'             specifies how many columns of \code{H} belong to the i-th
#'             distinct eigenvalue.}
#'
#' @details
#' The output defines the spectral decomposition of covariance matrix \eqn{\Sigma}
#' in the form \eqn{\Sigma = H \Lambda H^T} where \eqn{\Lambda} is a diagonal matrix of
#' eigenvalues, with eigenvalue \eqn{\lambda_i} having multiplicity \eqn{qq_i}.
#'
#' @export
basis_BCS_D <- function(p, m, sparse=FALSE) {
  params <- check_pm(p, m)
  p <- params$p
  m <- params$m

  L <- CTEigvec(p)
  M <- if (sparse) Matrix::Diagonal(m) else diag(m)

  idx <- c(seq_len(m), rep(seq_len(m), each = p - 1L) + m * rep(seq_len(p - 1L) - 1L, times = m) + m)

  H <- kronecker(L, M)

  list(
    H = H[, idx, drop = FALSE],
    qq = c(rep(1L, m), rep(p - 1L, m))
  )
}

#' @title Basis Matrix for BCS_CS structure
#'
#' @description
#' Generates the basis matrix of eigenvectors and a vector of
#' eigenvalue multiplicities for a doubly multivariate model with
#' block compound symmetry structure, where each block has a compound symmetry structure.
#'
#' @param p Integer, between-class dimension (must be >= 1).
#' @param m Integer, within-class dimension (must be >= 1).
#'
#' @return A list with two components:
#'   \item{H}{Basis matrix containing eigenvectors.}
#'   \item{qq}{Integer vector of eigenvalue multiplicities. Each \code{qq[i]}
#'             specifies how many columns of \code{H} belong to the i-th
#'             distinct eigenvalue.}
#'
#' @details
#' The output defines the spectral decomposition of covariance matrix \eqn{\Sigma}
#' in the form \eqn{\Sigma = H \Lambda H^T} where \eqn{\Lambda} is a diagonal matrix of
#' eigenvalues, with eigenvalue \eqn{\lambda_i} having multiplicity \eqn{qq_i}.
#'
#' @export
basis_BCS_CS <- function(p, m) {
  params <- check_pm(p, m)
  p <- params$p
  m <- params$m

  if (as.numeric(p) * as.numeric(m) > 1e6) {  # 1 million elements
    warning(sprintf(
      "Creating a %d x %d matrix. This may require substantial memory.",
      p * m, p * m
    ), call. = FALSE)
  }

  L <- CTEigvec(p)
  M <- CTEigvec(m)

  firsts <- seq_len(p - 2L) * m + 1L
  idx <- c(seq_len(m),
           c(firsts, setdiff(seq_len((p - 1L) * m), firsts)) + m)

  H <- kronecker(L, M)

  list(
    H = H[, idx, drop = FALSE],
    qq = c(1L, m - 1L, p - 1L, (p - 1L) * (m - 1L))
  )
}

#' @title Basis Matrix for BCS_CT structure
#'
#' @description
#' Generates the basis matrix of eigenvectors and a vector of
#' eigenvalue multiplicities for a doubly multivariate model with
#' block compound symmetry structure, where each block has a circular Toeplitz structure.
#'
#' @param p Integer, between-class dimension (must be >= 1).
#' @param m Integer, within-class dimension (must be >= 1).
#'
#' @return A list with two components:
#'   \item{H}{Basis matrix containing eigenvectors.}
#'   \item{qq}{Integer vector of eigenvalue multiplicities. Each \code{qq[i]}
#'             specifies how many columns of \code{H} belong to the i-th
#'             distinct eigenvalue.}
#'
#' @details
#' The output defines the spectral decomposition of covariance matrix \eqn{\Sigma}
#' in the form \eqn{\Sigma = H \Lambda H^T} where \eqn{\Lambda} is a diagonal matrix of
#' eigenvalues, with eigenvalue \eqn{\lambda_i} having multiplicity \eqn{qq_i}.
#'
#' @export
basis_BCS_CT <- function(p, m) {
  params <- check_pm(p, m)
  p <- params$p
  m <- params$m

  if (as.numeric(p) * as.numeric(m) > 1e6) {  # 1 million elements
    warning(sprintf(
      "Creating a %d x %d matrix. This may require substantial memory.",
      p * m, p * m
    ), call. = FALSE)
  }

  L <- CTEigvec(p)
  M <- CTEigvec(m, CTorder = TRUE)
  H <- kronecker(L, M)

  p_minus_1 <- p - 1L

  offsets <- if (p >2L) {
    m * seq.int(0L, p - 2L)
  } else {
    integer(0)
  }
  pair_starts <- if (m > 2L) {
    seq.int(2L, m - 1L, by = 2L)
  } else {
    integer(0)
  }
  idx <- 1L + offsets
  np <- length(pair_starts)

  if (np > 0) {
    rep_offsets <- rep(offsets, times = np)

    B1 <- rep(pair_starts, each = p_minus_1) + rep_offsets
    B2 <- rep(pair_starts + 1L, each = p_minus_1) + rep_offsets

    rest <- as.vector(rbind(B1, B2))
  } else {
    rest <- integer(0)
  }

  if ((m %% 2L) == 1L) {
    idxp <- c(idx, rest)
    half_m <- (m - 1L) %/% 2L
    qq <- c(1L, rep(2L, half_m), p_minus_1, rep(2L * p_minus_1, half_m))
  } else {
    idxp <- c(idx, rest, idx + m - 1L)
    half_m_minus_1 <- (m - 2L) %/% 2L
    qq <- c(1L, rep(2L, half_m_minus_1), 1L, p_minus_1, rep(2L * p_minus_1, half_m_minus_1), p_minus_1)
  }

  idx_final <- c(seq_len(m), idxp + m)

  list(
    H = H[, idx_final, drop = FALSE],
    qq = qq
  )
}

#' @title Basis Matrix for BCT_S structure
#'
#' @description
#' Generates the basis matrix of eigenvectors and a vector of
#' eigenvalue multiplicities for a doubly multivariate model with
#' block circular Toeplitz structure, where each block has a shericity structure.
#'
#' @param p Integer, between-class dimension (must be >= 1).
#' @param m Integer, within-class dimension (must be >= 1).
#' @param sparse Logical; if TRUE, returns a sparse identity matrix.
#'
#' @return A list with two components:
#'   \item{H}{Basis matrix containing eigenvectors.}
#'   \item{qq}{Integer vector of eigenvalue multiplicities. Each \code{qq[i]}
#'             specifies how many columns of \code{H} belong to the i-th
#'             distinct eigenvalue.}
#'
#' @details
#' The output defines the spectral decomposition of covariance matrix \eqn{\Sigma}
#' in the form \eqn{\Sigma = H \Lambda H^T} where \eqn{\Lambda} is a diagonal matrix of
#' eigenvalues, with eigenvalue \eqn{\lambda_i} having multiplicity \eqn{qq_i}.
#'
#' @export
basis_BCT_S <- function(p, m, sparse=FALSE) {
  params <- check_pm(p, m)
  p <- params$p
  m <- params$m

  L <- CTEigvec(p, CTorder = TRUE)
  M <- if (sparse) Matrix::Diagonal(m) else diag(m)

  H <- kronecker(L, M)

  if ((p %% 2L) == 1L) {
    qq <- c(m, rep.int(2L * m, (p - 1L) %/% 2L))
  } else {
    qq <- c(m, rep.int(2L * m, (p - 2L) %/% 2L), m)
  }

  list(
    H = H,
    qq = qq
  )
}

#' @title Basis Matrix for BCT_D structure
#'
#' @description
#' Generates the basis matrix of eigenvectors and a vector of
#' eigenvalue multiplicities for a doubly multivariate model with
#' block circular Toeplitz structure, where each block has a diagonal structure.
#'
#' @param p Integer, between-class dimension (must be >= 1).
#' @param m Integer, within-class dimension (must be >= 1).
#' @param sparse Logical; if TRUE, returns a sparse identity matrix.
#'
#' @return A list with two components:
#'   \item{H}{Basis matrix containing eigenvectors.}
#'   \item{qq}{Integer vector of eigenvalue multiplicities. Each \code{qq[i]}
#'             specifies how many columns of \code{H} belong to the i-th
#'             distinct eigenvalue.}
#'
#' @details
#' The output defines the spectral decomposition of covariance matrix \eqn{\Sigma}
#' in the form \eqn{\Sigma = H \Lambda H^T} where \eqn{\Lambda} is a diagonal matrix of
#' eigenvalues, with eigenvalue \eqn{\lambda_i} having multiplicity \eqn{qq_i}.
#'
#' @export
basis_BCT_D <- function(p, m, sparse=FALSE) {
  params <- check_pm(p, m)
  p <- params$p
  m <- params$m

  L <- CTEigvec(p, CTorder = TRUE)
  M <- if (sparse) Matrix::Diagonal(m) else diag(m)

  H <- kronecker(L, M)

  # indices within each 2m-sized eigenblock
  offsets <- rep(seq_len(m), each = 2L) + c(0L, m)

  if (p %% 2L == 1L) {
    n_blocks <- (p - 1L) %/% 2L
    if (n_blocks > 0L) {
      block_rep <- (seq_len(n_blocks) - 1L) * (2L * m)
      rest <- c(
        matrix(offsets, nrow = 2L * m, ncol = n_blocks) +
           rep.int(block_rep, 2L * m)
      )
    } else {
      rest <- integer(0)
    }
    qq <- c(rep.int(1L, m), rep.int(2L, (p - 1L) %/% 2L *m))
  } else {
    n_blocks <- (p - 2L) %/% 2L
    if (n_blocks > 0L) {
      block_rep <- (seq_len(n_blocks) - 1L) * (2L * m)
      rest <- c(
        matrix(offsets, nrow = 2L * m, ncol = n_blocks) +
          rep.int(block_rep, 2L * m)
      )
    } else {
      rest <- integer(0)
    }
    rest <- c(rest, seq_len(m) + (p - 1L) * m - m)

    qq <- c(rep(1L, m), rep(2L, (p - 2L) %/% 2L *m), rep(1L, m))
  }

  idx <- c(seq_len(m), rest + m)

  list(
    H = H[,idx, drop = FALSE],
    qq = qq
  )
}

#' @title Basis Matrix for BCT_CS structure
#'
#' @description
#' Generates the basis matrix of eigenvectors and a vector of
#' eigenvalue multiplicities for a doubly multivariate model with
#' block circular Toeplitz structure, where each block has a compound symmetry structure.
#'
#' @param p Integer, between-class dimension (must be >= 1).
#' @param m Integer, within-class dimension (must be >= 1).
#'
#' @return A list with two components:
#'   \item{H}{Basis matrix containing eigenvectors.}
#'   \item{qq}{Integer vector of eigenvalue multiplicities. Each \code{qq[i]}
#'             specifies how many columns of \code{H} belong to the i-th
#'             distinct eigenvalue.}
#'
#' @details
#' The output defines the spectral decomposition of covariance matrix \eqn{\Sigma}
#' in the form \eqn{\Sigma = H \Lambda H^T} where \eqn{\Lambda} is a diagonal matrix of
#' eigenvalues, with eigenvalue \eqn{\lambda_i} having multiplicity \eqn{qq_i}.
#'
#' @export
basis_BCT_CS <- function(p, m) {
  params <- check_pm(p, m)
  p <- params$p
  m <- params$m

  L <- CTEigvec(p, CTorder = TRUE)
  M <- CTEigvec(m)

  H <- kronecker(L, M)

  k <- (p + 1L) %% 2L + 1L
  firsts <- seq_len(p - k) * m + 1L - m
  idx <- c(seq_len(m),
           c(firsts, setdiff(seq_len((p - 1L) * m), firsts)) + m
           )

  if (p %% 2L == 1L) {
    qq=c(1L, m - 1L, rep(2L, (p - 1L) %/% 2L), rep(2L * (m - 1L), (p - 1L) %/% 2))
  }
  if (p %% 2L == 0L) {
    qq=c(1L, m - 1L, rep(2L, (p - 2L) %/% 2L), rep(2L * (m - 1L), (p - 2L) %/% 2L), 1L, m - 1L)
  }

  list(
    H = H[,idx, drop = FALSE],
    qq = qq
  )
}

#' @title Basis Matrix for BCT_CT structure
#'
#' @description
#' Generates the basis matrix of eigenvectors and a vector of
#' eigenvalue multiplicities for a doubly multivariate model with
#' block circular Toeplitz structure, where each block has a circular Toeplitz structure.
#'
#' @param p Integer, between-class dimension (must be >= 1).
#' @param m Integer, within-class dimension (must be >= 1).
#'
#' @return A list with two components:
#'   \item{H}{Basis matrix containing eigenvectors.}
#'   \item{qq}{Integer vector of eigenvalue multiplicities. Each \code{qq[i]}
#'             specifies how many columns of \code{H} belong to the i-th
#'             distinct eigenvalue.}
#'
#' @details
#' The output defines the spectral decomposition of covariance matrix \eqn{\Sigma}
#' in the form \eqn{\Sigma = H \Lambda H^T} where \eqn{\Lambda} is a diagonal matrix of
#' eigenvalues, with eigenvalue \eqn{\lambda_i} having multiplicity \eqn{qq_i}.
#'
#' @export
basis_BCT_CT <- function(p, m) {
  params <- check_pm(p, m)
  p <- params$p
  m <- params$m

  L <- CTEigvec(p, CTorder=TRUE)
  M <- CTEigvec(m, CTorder=TRUE)

  H <- kronecker(L, M)

  r <- 2L - (p %% 2L)

  if (m %% 2L == 1L){
    t <- seq_len((m - 1L) %/% 2L)
    block_mat <- matrix(c(
      2L * t,
      2L * t + 1L,
      m + 2L * t,
      m + 2L * t + 1L
    ), nrow = 4, byrow = TRUE)

    B <- c(1L, m + 1L, as.vector(block_mat))

    k_max <- (p - 3L) %/% 2L
    if (k_max < 0L) k_max <- 0L

    offsets <- (2L * m) * (0L:k_max)

    idx <- c(seq_len(m), as.vector(B + rep(offsets, each = length(B))) + m)

    qq=c(1L, rep(2L, (m - 1L) %/% 2L),
             rep(2L, (p - r) %/% 2L),
             rep(4L, (m - 1L) %/% 2L * (p - r) %/% 2L))
  }

  if (m %% 2L == 0L) {
    t <- seq(2L, m - 1L, by = 2L)

    if (length(t) > 0) {
      block_mat <- matrix(c(
        t,
        t + 1L,
        m + t,
        m + t + 1L
      ), nrow = 4, byrow = TRUE)

      mid_part <- as.vector(block_mat)
    } else {
      mid_part <- integer(0)
    }

    if (m %% 2L == 1L) {
      B <- c(1L, m + 1L, mid_part)
    } else {
      B <- c(1L, m + 1L, mid_part, m, 2L * m)
    }

    k_max <- (p - 3L) %/% 2L
    if (k_max < 0) k_max <- 0L
    offsets <- (2L * m) * (0L:k_max)

    idx <- c(seq_len(m), as.vector(B + rep(offsets, each = length(B))) + m)

    qq=c(1L, rep(2L, (m - 2L) %/% 2L), 1L,
             rep(2L,(p - r) %/% 2L),
             rep(4L, (m - 2L) %/% 2L * (p - r) %/% 2L),
             rep(2L, (p - r) %/% 2))
  }

  if (p %% 2L == 0L) {
    idx <- c(idx, seq_len(m) + (p - 1L) * m)

    if (m %% 2L == 0L) {
      qq<-c(qq, 1L, rep(2L, (m - 2L) %/% 2L), 1L)
    } else {
      qq<-c(qq, 1L, rep(2L, (m - 1L) %/% 2L))
    }
  }

  list(
    H = H[,idx, drop = FALSE],
    qq = qq
  )
}

#' @title Basis Matrix for UBS_CS structure
#'
#' @description
#' Generates the basis matrix of eigenvectors and a vector of
#' eigenvalue multiplicities for a doubly multivariate model with
#' uniform block structure, where the block-eigenmatrix has a compound symmetry structure.
#'
#' @param p Integer, between-class dimension (must be >= 1).
#' @param m_vec Integer vector of within-class dimensions, must be of length \code{p}.
#'
#' @return A list with two components:
#'   \item{H}{Basis matrix containing eigenvectors.}
#'   \item{qq}{Integer vector of eigenvalue multiplicities. Each \code{qq[i]}
#'             specifies how many columns of \code{H} belong to the i-th
#'             distinct eigenvalue.}
#'
#' @details
#' The output defines the spectral decomposition of covariance matrix \eqn{\Sigma}
#' in the form \eqn{\Sigma = H \Lambda H^T} where \eqn{\Lambda} is a diagonal matrix of
#' eigenvalues, with eigenvalue \eqn{\lambda_i} having multiplicity \eqn{qq_i}.
#'
#' @export
basis_UBS_CS <- function(p, m_vec) {
  # Check p
  if (!is.numeric(p) || length(p) != 1L || p < 1L || p != as.integer(p)) {
    stop("`p` must be a positive integer.", call. = FALSE)
  }
  p <- as.integer(p)
  # Check m_vec
  if (length(m_vec) != p) stop("length(m_vec) must equal p", call. = FALSE)
  if (any(m_vec < 1L)) stop("all m_vec must be >= 1", call. = FALSE)
    m_vec <- as.integer(m_vec)

  # Dimensions
  m <- sum(m_vec)
  Q_cols <- sum(pmax(m_vec - 1L, 0L))  # Handle m_i = 1 case
  total_cols <- Q_cols + p

  # Preallocate
  H <- matrix(0, m, total_cols)

  # Get L matrix
  L <- CTEigvec(p)

  # Row indices for blocks
  row_starts <- c(1L, cumsum(m_vec[-p]) + 1L)

  # --- Fill Q blocks ---
  Q_current_col <- 1L
  Q_block_cols <- m_vec - 1L
  d_col_starts <- c(Q_current_col, cumsum(Q_block_cols[-p]) + Q_current_col)

  for (i in seq_len(p)) {
    mi <- m_vec[i]
    if (mi > 1L) {
      rows <- row_starts[i]:(row_starts[i] + mi - 1L)
      cols <- d_col_starts[i]:(d_col_starts[i] + mi - 2L)
      H[rows, cols] <- CTEigvec(mi)[, -1L, drop = FALSE]
    }
  }

  # Prepare Pm scaling factors
  pm_scales <- 1 / sqrt(m_vec)

  # Single column: Pm %*% L[, 1]
  pm_col <- Q_cols + 1L
  for (i in seq_len(p)) {
    rows <- row_starts[i]:(row_starts[i] + m_vec[i] - 1L)
    H[rows, pm_col] <- L[i, 1L] * pm_scales[i]
  }

  # Multiple columns: Pm %*% L[, 2:p]
  if (p > 1L) {
    for (j in 2L:p) {
      pm_col <- pm_col + 1L
      for (i in seq_len(p)) {
        rows <- row_starts[i]:(row_starts[i] + m_vec[i] - 1L)
        H[rows, pm_col] <- L[i, j] * pm_scales[i]
      }
    }
  }

  list(
    H = H,
    qq = c(m_vec-1L, 1L, p - 1L)
  )
}

#' @title Basis Matrix for DCS structure
#'
#' @description
#' Generates the basis matrix of eigenvectors and a vector of
#' eigenvalue multiplicities for a doubly multivariate model with
#' double complete symmetry structure.
#'
#' @param p Integer, between-class dimension (must be >= 1).
#' @param m Integer, within-class dimension (must be >= 1).
#'
#' @return A list with two components:
#'   \item{H}{Basis matrix containing eigenvectors.}
#'   \item{qq}{Integer vector of eigenvalue multiplicities. Each \code{qq[i]}
#'             specifies how many columns of \code{H} belong to the i-th
#'             distinct eigenvalue.}
#'
#' @details
#' The output defines the spectral decomposition of covariance matrix \eqn{\Sigma}
#' in the form \eqn{\Sigma = H \Lambda H^T} where \eqn{\Lambda} is a diagonal matrix of
#' eigenvalues, with eigenvalue \eqn{\lambda_i} having multiplicity \eqn{qq_i}.
#'
#' @export
basis_DCS <- function(p, m) {
  params <- check_pm(p, m)
  p <- params$p
  m <- params$m

  L <- CTEigvec(p)
  M <- CTEigvec(m)

  H <- cbind(kronecker(L[,1],M[,1]),kronecker(diag(p),M[,2:m]),kronecker(L[,2:p],M[,1]))

  list(
    H = H,
    qq = c(1L, p*(m - 1L), p - 1L)
  )
}

