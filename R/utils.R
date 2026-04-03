#' Construct Companion Matrix for AR(p) Process
#'
#' Builds the companion matrix representation for an AR(p) process, which is
#' used in the state-space formulation of SV(p) models.
#'
#' @param phi Numeric vector of AR coefficients (length p).
#' @param p Integer. Order of the AR process.
#' @param q Integer. Dimension of each block (typically 1 for univariate).
#'
#' @return A \code{(q*p) x (q*p)} companion matrix.
#'
#' @keywords internal
companionMat <- function(phi, p, q) {
  F_tmp <- phi
  if (p > 1) {
    diagmat <- diag(q * (p - 1))
    diagzero <- matrix(0, q * (p - 1), q)
    Mn <- cbind(diagmat, diagzero)
    compMat <- rbind(F_tmp, Mn)
  } else {
    compMat <- matrix(F_tmp, nrow = q, ncol = q)
  }
  return(compMat)
}
