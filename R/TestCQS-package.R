#' @title
#' Tools for testing commutative quadratic subspace covariance structures
#'
#' @description
#' Provides tools for testing whether a scale/covariance matrix belongs
#' to a specified commutative quadratic subspace (CQS) structure under
#' multivariate and, more generally, elliptically contoured models.
#'
#' The package provides likelihood ratio and Rao score tests through
#' [runTestCQS()] and functions for constructing orthogonal
#' basis matrices for selected CQS covariance structures.
#'
#' The functions for generating the basis matrices `H` consisted of eigenvectors and
#' vector of multiplicities `qq` of corresponding eigenvalues for selected
#' scale/covariance structures used in multivariate analysis.
#'
#' All basis_*() functions return a list with two components:
#' \describe{
#'   \item{H}{Basis matrix containing eigenvectors.}
#'   \item{qq}{Integer vector of eigenvalue multiplicities. Each `qq[i]`
#'             specifies how many columns of `H` belong to the i-th
#'             distinct eigenvalue.}
#' }
#'
#' For a scale/covariance matrix \eqn{\Sigma}, we have
#' \eqn{\Sigma = H \Lambda H^T} where \eqn{\Lambda} is a diagonal matrix of
#' eigenvalues, with eigenvalue \eqn{\lambda_i} having multiplicity \eqn{qq_i}.
#'
#' @details
#' The package supports the following covariance structures:
#'
#' **Ordinary multivariate model structures:**
#' \itemize{
#'   \item Sphericity (S)
#'   \item Diagonality (D)
#'   \item Compound symmetry (CS)
#'   \item Circular Toeplitz (CT)
#' }
#'
#' **Doubly multivariate model structures:**
#' \itemize{
#'   \item Block sphericity (BS)
#'   \item Block diagonality (BD)
#'   \item Block compound symmetry (BCS)
#'   \item Block circular Toeplitz (BCT)
#'   \item Their combinations
#' }
#'
#' @name TestCQS
"_PACKAGE"

