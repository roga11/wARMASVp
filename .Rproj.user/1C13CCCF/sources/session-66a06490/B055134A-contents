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
#' @examples
#' companionMat(c(0.7, 0.2), p = 2, q = 1)
#'
#' @export
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

#' Create Lagged Matrix from a Time Series
#'
#' @param y Numeric vector. The time series.
#' @param lags Integer. Number of lags to include.
#' @param balance Logical. If \code{TRUE}, remove rows with NAs.
#'
#' @return A matrix with columns \code{y(t), y(t-1), ..., y(t-lags)}.
#'
#' @keywords internal
laggedMat <- function(y, lags, balance = TRUE) {
  Tsize <- length(y)
  lagged_matrix <- matrix(NA, nrow = Tsize, ncol = lags + 1)
  for (i in 0:lags) {
    lagged_matrix[(i + 1):Tsize, i + 1] <- y[1:(Tsize - i)]
  }
  colnames(lagged_matrix) <- c("y.lag(0)", paste0("y.lag(", 1:lags, ")"))
  if (balance) {
    lagged_matrix <- lagged_matrix[stats::complete.cases(lagged_matrix), ]
    rownames(lagged_matrix) <- NULL
  }
  return(lagged_matrix)
}
