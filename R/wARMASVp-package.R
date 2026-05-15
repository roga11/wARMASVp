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
#'   \item \code{\link{svpSE}} -- Simulation-based standard errors
#'   \item \code{\link{sim_svp}} -- Simulate SV(p) processes
#'   \item \code{\link{filter_svp}} -- Kalman/mixture/particle filtering
#'   \item \code{\link{forecast_svp}} -- Multi-step volatility forecasts
#'   \item \code{\link{lmc_ar}}, \code{\link{mmc_ar}} -- AR order tests
#'   \item \code{\link{lmc_lev}}, \code{\link{mmc_lev}} -- Leverage tests
#'   \item \code{\link{lmc_t}}, \code{\link{mmc_t}} -- Student-t tail tests
#'   \item \code{\link{lmc_ged}}, \code{\link{mmc_ged}} -- GED tail tests
#' }
#'
#' @references
#' Ahsan, M. N. and Dufour, J.-M. (2021). Simple estimators and inference for
#' higher-order stochastic volatility models. \emph{Journal of Econometrics},
#' 224(1), 181-197. \doi{10.1016/j.jeconom.2021.03.008}
#'
#' Ahsan, M. N., Dufour, J.-M., and Rodriguez-Rondon, G. (2025). Estimation and
#' inference for higher-order stochastic volatility models with leverage.
#' \emph{Journal of Time Series Analysis}, 46(6), 1064-1084.
#' \doi{10.1111/jtsa.12851}
#'
#' Ahsan, M. N., Dufour, J.-M., and Rodriguez-Rondon, G. (2026). Estimation and
#' inference for stochastic volatility models with heavy-tailed distributions.
#' Bank of Canada Staff Working Paper 2026-8. \doi{10.34989/swp-2026-8}
#'
#' @useDynLib wARMASVp, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom gsignal poly
#' @importFrom stats qnorm quantile var rnorm uniroot dnorm pgamma qgamma
#'   rchisq rf rgamma rt runif setNames
#' @keywords internal
"_PACKAGE"
