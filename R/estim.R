#' Estimate a Stochastic Volatility Model
#'
#' Master estimation function for SV(p) models using the Winsorized ARMA-SV
#' (W-ARMA-SV) method. Supports Gaussian, Student-t, and GED error
#' distributions, with optional leverage effects.
#'
#' @details
#' The model is:
#' \deqn{y_t = \sigma_y \exp(w_t / 2) z_t}
#' \deqn{w_t = \phi_1 w_{t-1} + \cdots + \phi_p w_{t-p} + \sigma_v v_t}
#'
#' where \eqn{z_t} follows a distribution specified by \code{errorType}
#' (Gaussian, Student-t, or GED), and \eqn{v_t} is i.i.d. standard normal.
#' When \code{leverage = TRUE}, the correlation between \eqn{z_t} and
#' \eqn{v_t} is estimated as \eqn{\rho}.
#'
#' For Student-t errors with leverage, the correction factor
#' \eqn{C_t(\nu)} from the scale-mixture representation is applied.
#' For GED errors with leverage, the exact implicit equation is solved
#' via 1D root-finding with Gauss-Hermite quadrature.
#'
#' @param y Numeric vector. Observed returns (e.g., de-meaned log returns).
#' @param p Integer. Order of the volatility process. Default is 1.
#' @param J Integer. Winsorizing parameter controlling the number of
#'   autocovariance blocks used. Default is 10.
#' @param leverage Logical. If \code{TRUE}, estimate leverage parameter
#'   \eqn{\rho}. Default is \code{FALSE}.
#' @param errorType Character. Error distribution: \code{"Gaussian"} (default),
#'   \code{"Student-t"}, or \code{"GED"}.
#' @param rho_type Character. Correlation type for leverage estimation. One of
#'   \code{"pearson"} (default), \code{"kendall"}, or \code{"both"}.
#' @param del Numeric. Small constant for log transformation:
#'   \eqn{\log(y_t^2 + \delta)}. Default is \code{1e-10}.
#' @param trunc_lev Logical. If \code{TRUE}, truncate the estimated leverage
#'   parameter to \code{[-0.999, 0.999]}. Default is \code{TRUE}. Setting to
#'   \code{FALSE} can reduce bias in some cases but may yield estimates
#'   outside the parameter space.
#' @param wDecay Logical. Use linearly decaying weights in the WLS estimation.
#'   Default is \code{FALSE}.
#' @param logNu Logical. Solve for \eqn{\nu} in log-space for numerical
#'   stability (Student-t only). Default is \code{FALSE}.
#' @param sigvMethod Character. Method for estimating \eqn{\sigma_v}. One of:
#'   \code{"factored"} (default) — factored-variance estimator (recommended;
#'   dominates RMSE in most settings, see ADRR 2025);
#'   \code{"direct"} — direct variance decomposition;
#'   \code{"hybrid"} — AD2021 closed-form for \eqn{p = 1}, falls back to
#'   \code{"direct"} for \eqn{p \ge 2} (Student-t and GED only).
#' @param winsorize_eps Integer. Number of extreme autocovariance lags to
#'   winsorize (\code{0} = none). Used in Student-t and GED
#'   \eqn{\sigma_\varepsilon^2} estimation. Default \code{0}.
#'
#' @return Depending on \code{errorType}:
#' \itemize{
#'   \item \code{"Gaussian"}: An object of class \code{"svp"} (see below).
#'   \item \code{"Student-t"}: An object of class \code{"svp_t"}.
#'   \item \code{"GED"}: An object of class \code{"svp_ged"}.
#' }
#'
#' The \code{"svp"} class contains:
#' \describe{
#'   \item{mu}{Mean of \eqn{\log(y_t^2 + \delta)}.}
#'   \item{phi}{Numeric vector of AR coefficients.}
#'   \item{sigv}{Standard deviation of volatility innovations.}
#'   \item{sigy}{Unconditional standard deviation.}
#'   \item{rho}{Leverage parameter (if estimated, otherwise \code{NA}).}
#'   \item{y}{The original data.}
#'   \item{p, J}{Model order and winsorizing parameter.}
#'   \item{errorType}{The error distribution used.}
#'   \item{call}{The matched call.}
#' }
#'
#' @examples
#' \donttest{
#' # Gaussian SV(1) without leverage (default)
#' y <- sim_svp(1000, phi = 0.95, sigy = 1, sigv = 0.2)$y
#' fit <- svp(y)
#' summary(fit)
#'
#' # With leverage
#' y2 <- sim_svp(1000, phi = 0.95, sigy = 1, sigv = 0.2, leverage = TRUE, rho = -0.3)$y
#' fit2 <- svp(y2, leverage = TRUE)
#' coef(fit2)
#'
#' # Student-t errors
#' y3 <- sim_svp(1000, phi = 0.9, sigy = 1, sigv = 0.2, errorType = "Student-t", nu = 5)$y
#' fit3 <- svp(y3, errorType = "Student-t")
#' summary(fit3)
#' }
#'
#' @seealso \code{\link{svpSE}} for standard errors.
#'
#' @references
#' Ahsan, N. and Dufour, J.-M. (2021). Simple estimators and inference for
#' higher-order stochastic volatility models. \emph{Journal of Econometrics},
#' 224(1), 181-197.
#'
#' Ahsan, N., Dufour, J.-M., and Rodriguez Rondon, G. (2025). Estimation and
#' inference for stochastic volatility models with heavy-tailed distributions.
#'
#' @export
svp <- function(y, p = 1, J = 10, leverage = FALSE, errorType = "Gaussian",
                rho_type = "pearson", del = 1e-10, trunc_lev = TRUE,
                wDecay = FALSE, logNu = FALSE,
                sigvMethod = "factored", winsorize_eps = 0) {
  cl <- match.call()
  y <- as.numeric(y)
  if (length(y) < 1L) stop("'y' must be a non-empty numeric vector.")
  if (any(!is.finite(y))) stop("'y' must not contain NA, NaN, or Inf values.")
  p <- as.integer(p)
  if (p < 1L) stop("'p' must be a positive integer (>= 1).")
  J <- as.integer(J)
  if (J < 1L) stop("'J' must be a positive integer (>= 1).")
  min_length <- 2L * p + J
  if (length(y) < min_length)
    stop("'y' must have length >= 2*p + J = ", min_length,
         " (got ", length(y), ").")
  if (!is.numeric(del) || del <= 0) stop("'del' must be a positive number.")
  errorType <- match.arg(errorType, c("Gaussian", "Student-t", "GED"))
  sigvMethod <- match.arg(sigvMethod, c("hybrid", "direct", "factored"))
  # Dispatch: estimate base model (without leverage)
  if (errorType == "Gaussian") {
    out <- .svp_gaussian(y, p, J, leverage, rho_type, del, trunc_lev, wDecay,
                         sigvMethod)
  } else if (errorType == "Student-t") {
    out <- .svp_t(y, p, J, del, wDecay, logNu, sigvMethod, winsorize_eps)
    if (leverage) {
      out <- .add_leverage(out, y, as.integer(p), rho_type, del, trunc_lev,
                           wDecay, "Student-t")
    }
  } else if (errorType == "GED") {
    out <- .svp_ged(y, p, J, del, wDecay, sigvMethod, winsorize_eps)
    if (leverage) {
      out <- .add_leverage(out, y, as.integer(p), rho_type, del, trunc_lev,
                           wDecay, "GED")
    }
  }
  out$errorType <- errorType
  out$call <- cl
  return(out)
}


#' Simulation-Based Standard Errors for SV(p) Models
#'
#' Computes standard errors and confidence intervals for estimated parameters
#' by simulating from the fitted model and re-estimating. Supports all model
#' types returned by \code{\link{svp}}: Gaussian (with or without leverage),
#' Student-t, and GED.
#'
#' @param object A fitted model object from \code{\link{svp}}. Can be of class
#'   \code{"svp"}, \code{"svp_t"}, or \code{"svp_ged"}.
#' @param n_sim Integer. Number of Monte Carlo replications. Default 199.
#' @param alpha Numeric. Significance level for confidence intervals. Default 0.05.
#' @param burnin Integer. Burn-in period for simulation. Default 500.
#' @param logNu Logical. Solve for \eqn{\nu} in log-space for numerical
#'   stability (Student-t only). Default is \code{FALSE}.
#'
#' @return A list with:
#' \describe{
#'   \item{CI}{2 x k matrix of confidence intervals (lower, upper).}
#'   \item{SEsim0}{Standard errors relative to true parameter values.}
#'   \item{SEsim}{Standard errors relative to sample mean.}
#'   \item{ISEconservative}{Conservative interval-based standard errors.}
#'   \item{ISEliberal}{Liberal interval-based standard errors.}
#'   \item{thetamat}{Matrix of parameter estimates from simulations.}
#' }
#'
#' @examples
#' \donttest{
#' # Gaussian SV(1)
#' y <- sim_svp(1000, phi = 0.95, sigy = 1, sigv = 0.2)$y
#' fit <- svp(y)
#' se <- svpSE(fit, n_sim = 49)
#' se$CI
#' }
#'
#' @export
svpSE <- function(object, n_sim = 199, alpha = 0.05, burnin = 500,
                  logNu = FALSE) {
  if (!inherits(object, c("svp", "svp_t", "svp_ged"))) {
    stop("object must be of class 'svp', 'svp_t', or 'svp_ged'.")
  }
  if (n_sim < 1) stop("n_sim must be >= 1.")
  if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1)
    stop("'alpha' must be a number in (0, 1).")
  if (!is.numeric(burnin) || burnin < 0)
    stop("'burnin' must be a non-negative integer.")

  if (inherits(object, "svp_t")) {
    .svpSE_t(object, n_sim, alpha, burnin, logNu)
  } else if (inherits(object, "svp_ged")) {
    .svpSE_ged(object, n_sim, alpha, burnin)
  } else {
    .svpSE_gaussian(object, n_sim, alpha, burnin)
  }
}

