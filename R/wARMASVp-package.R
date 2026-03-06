#' @title wARMASVp: Winsorized ARMA Estimation for Higher-Order Stochastic
#'   Volatility Models
#'
#' @description Estimation, simulation, hypothesis testing, and forecasting for
#'   univariate higher-order stochastic volatility SV(p) models. Supports
#'   Gaussian, Student-t, and GED innovations with optional leverage effects.
#'
#' @details The main user-facing functions are:
#' \itemize{
#'   \item \code{\link{svp}} -- Estimate SV(p) model with optional leverage
#'   \item \code{\link{sim_svp}} -- Simulate SV(p) processes
#'   \item \code{\link{lmc_lev}}, \code{\link{mmc_lev}} -- Test for leverage
#'   \item \code{\link{kalman_filter}}, \code{\link{forecast_svp}} --
#'     Filter and forecast volatility
#' }
#'
#' @references
#' Ahsan, N. and Dufour, J.-M. (2021). Simple estimators and inference for
#' higher-order stochastic volatility models. \emph{Journal of Econometrics},
#' 224(1), 181-197.
#'
#' Ahsan, N., Dufour, J.-M., and Rodriguez Rondon, G. (2025). Estimation and
#' inference for higher-order stochastic volatility models with leverage.
#' \emph{Journal of Time Series Analysis}.
#'
#' Ahsan, N., Dufour, J.-M., and Rodriguez Rondon, G. (2025). Estimation and
#' inference for stochastic volatility models with heavy-tailed distributions.
#'
#' @useDynLib wARMASVp, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats qnorm quantile var rnorm uniroot
#' @importFrom gsignal poly
#' @keywords internal
"_PACKAGE"
