#' Simulate from a Stochastic Volatility Model
#'
#' Master simulation function for SV(p) models. Supports Gaussian, Student-t,
#' and GED error distributions, with optional leverage effects. This mirrors
#' the interface of \code{\link{svp}} for estimation.
#'
#' @details
#' The model is:
#' \deqn{y_t = \sigma_y \exp(w_t / 2) z_t}
#' \deqn{w_t = \phi_1 w_{t-1} + \cdots + \phi_p w_{t-p} + \sigma_v v_t}
#'
#' where \eqn{z_t} follows a distribution specified by \code{errorType}
#' (Gaussian, Student-t, or GED), and \eqn{v_t} is i.i.d. standard normal.
#' When \code{leverage = TRUE}, the correlation between \eqn{z_t} and
#' \eqn{v_t} is \eqn{\rho}.
#'
#' For Student-t errors with leverage, the scale-mixture representation
#' \eqn{u_t = z_t \lambda_t^{-1/2}} is used, where leverage operates through
#' the Gaussian component \eqn{z_t}. For GED errors with leverage, a Gaussian
#' copula construction \eqn{u_t = F_{GED}^{-1}(\Phi(z_t))} is used.
#'
#' @param n Integer. Length of the simulated series.
#' @param phi Numeric vector. AR coefficients for log-volatility (length p).
#' @param sigy Numeric. Unconditional standard deviation of returns.
#' @param sigv Numeric. Standard deviation of volatility innovations.
#' @param errorType Character. Error distribution: \code{"Gaussian"} (default),
#'   \code{"Student-t"}, or \code{"GED"}.
#' @param leverage Logical. If \code{TRUE}, simulate with leverage effects
#'   (correlated return and volatility shocks). Default is \code{FALSE}.
#' @param rho Numeric. Leverage parameter (correlation between return and
#'   volatility shocks). Must be in \code{[-1, 1]}. Only used when
#'   \code{leverage = TRUE}. Default is 0.
#' @param nu Numeric. Shape parameter for heavy-tailed distributions.
#'   Degrees of freedom for Student-t (must be > 2) or GED shape (must be > 0).
#'   Required when \code{errorType} is \code{"Student-t"} or \code{"GED"}.
#' @param K Integer. Number of independent series to simulate. Default 1.
#'   Only used for Student-t and GED.
#' @param burnin Integer. Number of initial observations to discard. Default 500.
#'
#' @return Depending on the configuration:
#' \itemize{
#'   \item Gaussian without leverage: a numeric vector of length \code{n}.
#'   \item Gaussian with leverage: a list with components \code{y} (returns),
#'     \code{h} (log-volatility), \code{zeta} (return shocks), \code{veta}
#'     (volatility shocks).
#'   \item Student-t or GED: if \code{K = 1}, a numeric vector of length
#'     \code{n}; otherwise an \code{n x K} matrix.
#' }
#'
#' @examples
#' \donttest{
#' # Gaussian SV(1), no leverage
#' y <- sim_svp(1000, phi = 0.95, sigy = 1, sigv = 0.2)
#' plot(y, type = "l")
#'
#' # Gaussian SV(1) with leverage
#' sim <- sim_svp(1000, phi = 0.95, sigy = 1, sigv = 0.2,
#'                leverage = TRUE, rho = -0.5)
#' plot(sim$y, type = "l")
#'
#' # Student-t SV(1)
#' y_t <- sim_svp(1000, phi = 0.95, sigy = 1, sigv = 0.2,
#'                errorType = "Student-t", nu = 5)
#' plot(y_t, type = "l")
#'
#' # GED SV(1)
#' y_ged <- sim_svp(1000, phi = 0.95, sigy = 1, sigv = 0.2,
#'                  errorType = "GED", nu = 1.5)
#' plot(y_ged, type = "l")
#' }
#'
#' @seealso \code{\link{svp}} for estimation.
#'
#' @export
sim_svp <- function(n, phi, sigy, sigv, errorType = "Gaussian",
                    leverage = FALSE, rho = 0, nu = NULL, K = 1,
                    burnin = 500) {
  phi <- as.numeric(phi)
  p <- length(phi)
  errorType <- match.arg(errorType, c("Gaussian", "Student-t", "GED"))

  if (errorType == "Gaussian") {
    if (leverage) {
      if (abs(rho) > 1) stop("rho must be in [-1, 1].")
      beta <- c(phi, sigy, sigv, rho)
      out <- sim_svp_leverage_norm_cpp(beta, p, n, burnin)
      out$y <- as.numeric(out$y)
      out$h <- as.numeric(out$h)
      out$zeta <- as.numeric(out$zeta)
      out$veta <- as.numeric(out$veta)
      return(out)
    } else {
      beta <- c(phi, sigy, sigv)
      y <- sim_svp_norm_cpp(beta, p, n, burnin)
      return(as.numeric(y))
    }
  } else if (errorType == "Student-t") {
    if (is.null(nu)) stop("nu is required for Student-t errors.")
    if (nu <= 2) stop("nu must be > 2 for Student-t errors.")
    if (leverage) {
      if (abs(rho) > 1) stop("rho must be in [-1, 1].")
      beta <- c(phi, sigy, sigv, nu, rho)
      out <- sim_svp_leverage_t_cpp(beta, p, n, burnin)
      out$y <- as.numeric(out$y)
      out$h <- as.numeric(out$h)
      out$zeta <- as.numeric(out$zeta)
      out$veta <- as.numeric(out$veta)
      return(out)
    } else {
      beta <- c(phi, sigy, sigv, nu)
      y <- sim_sv_t_cpp(beta, p, n, K, burnin)
      if (K == 1) return(as.numeric(y))
      return(y)
    }
  } else if (errorType == "GED") {
    if (is.null(nu)) stop("nu is required for GED errors.")
    if (nu <= 0) stop("nu must be > 0 for GED errors.")
    if (leverage) {
      if (abs(rho) > 1) stop("rho must be in [-1, 1].")
      beta <- c(phi, sigy, sigv, nu, rho)
      out <- sim_svp_leverage_ged_cpp(beta, p, n, burnin)
      out$y <- as.numeric(out$y)
      out$h <- as.numeric(out$h)
      out$zeta <- as.numeric(out$zeta)
      out$veta <- as.numeric(out$veta)
      return(out)
    } else {
      beta <- c(phi, sigy, sigv, nu)
      y <- sim_sv_ged_cpp(beta, p, n, K, burnin)
      if (K == 1) return(as.numeric(y))
      return(y)
    }
  }
}

