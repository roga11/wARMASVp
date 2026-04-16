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
#' \eqn{v_{t+1}} is \eqn{\rho}.
#'
#' For Student-t errors with leverage, the scale-mixture representation
#' \eqn{z_t = \zeta_t \lambda_t^{-1/2}} is used, where leverage operates through
#' the Gaussian component \eqn{\zeta_t}. For GED errors with leverage, a Gaussian
#' copula construction \eqn{z_t = F_{\mathrm{GED}}^{-1}(\Phi(\zeta_t))} is used.
#' In both cases the returned \code{z} is the \emph{effective} return innovation
#' (not the latent \eqn{\zeta_t}), with marginal distribution matching the
#' \code{errorType}.
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
#' @param burnin Integer. Number of initial observations to discard. Default 500.
#'
#' @return A named list of four length-\code{n} numeric vectors:
#' \describe{
#'   \item{\code{y}}{Observed returns \eqn{y_t}.}
#'   \item{\code{h}}{Log-volatility process \eqn{w_t} (equivalently \eqn{h_t}).}
#'   \item{\code{z}}{Return innovation such that
#'     \eqn{y_t = \sigma_y \exp(h_t/2)\, z_t}.
#'     Marginal distribution matches \code{errorType}: N(0,1) for Gaussian,
#'     t(\eqn{\nu}) for Student-t, unit-variance GED(\eqn{\nu}) for GED.}
#'   \item{\code{v}}{Volatility innovation such that
#'     \eqn{h_t - \sum_{j=1}^p \phi_j h_{t-j} = \sigma_v\, v_t}.
#'     Always N(0,1); under leverage, \eqn{v_t = \rho\, \zeta_{t-1} +
#'     \sqrt{1-\rho^2}\, \epsilon_t}.}
#' }
#'
#' @examples
#' \donttest{
#' # Gaussian SV(1), no leverage
#' sim <- sim_svp(1000, phi = 0.95, sigy = 1, sigv = 0.2)
#' plot(sim$y, type = "l")
#'
#' # Gaussian SV(1) with leverage
#' sim_lev <- sim_svp(1000, phi = 0.95, sigy = 1, sigv = 0.2,
#'                    leverage = TRUE, rho = -0.5)
#' plot(sim_lev$y, type = "l")
#'
#' # Student-t SV(1)
#' sim_t <- sim_svp(1000, phi = 0.95, sigy = 1, sigv = 0.2,
#'                  errorType = "Student-t", nu = 5)
#' plot(sim_t$y, type = "l")
#'
#' # GED SV(1)
#' sim_ged <- sim_svp(1000, phi = 0.95, sigy = 1, sigv = 0.2,
#'                    errorType = "GED", nu = 1.5)
#' plot(sim_ged$y, type = "l")
#' }
#'
#' @seealso \code{\link{svp}} for estimation.
#'
#' @export
sim_svp <- function(n, phi, sigy, sigv, errorType = "Gaussian",
                    leverage = FALSE, rho = 0, nu = NULL,
                    burnin = 500) {
  if (!is.numeric(n) || length(n) != 1L || n < 1L || n != as.integer(n))
    stop("'n' must be a positive integer.")
  n <- as.integer(n)
  phi <- as.numeric(phi)
  if (length(phi) < 1L || any(!is.finite(phi)))
    stop("'phi' must be a non-empty numeric vector with finite values.")
  p <- length(phi)
  if (!is.numeric(sigy) || length(sigy) != 1L || sigy <= 0)
    stop("'sigy' must be a positive number.")
  if (!is.numeric(sigv) || length(sigv) != 1L || sigv <= 0)
    stop("'sigv' must be a positive number.")
  if (!is.numeric(burnin) || length(burnin) != 1L || burnin < 0)
    stop("'burnin' must be a non-negative integer.")
  errorType <- match.arg(errorType, c("Gaussian", "Student-t", "GED"))

  if (errorType == "Gaussian") {
    if (leverage) {
      if (abs(rho) > 1) stop("rho must be in [-1, 1].")
      beta <- c(phi, sigy, sigv, rho)
      out <- sim_svp_leverage_norm_cpp(beta, p, n, burnin)
    } else {
      beta <- c(phi, sigy, sigv)
      out <- sim_svp_norm_cpp(beta, p, n, burnin)
    }
  } else if (errorType == "Student-t") {
    if (is.null(nu)) stop("nu is required for Student-t errors.")
    if (nu <= 2) stop("nu must be > 2 for Student-t errors.")
    if (leverage) {
      if (abs(rho) > 1) stop("rho must be in [-1, 1].")
      beta <- c(phi, sigy, sigv, nu, rho)
      out <- sim_svp_leverage_t_cpp(beta, p, n, burnin)
    } else {
      beta <- c(phi, sigy, sigv, nu)
      out <- sim_sv_t_cpp(beta, p, n, burnin)
    }
  } else {  # "GED"
    if (is.null(nu)) stop("nu is required for GED errors.")
    if (nu <= 0) stop("nu must be > 0 for GED errors.")
    if (leverage) {
      if (abs(rho) > 1) stop("rho must be in [-1, 1].")
      beta <- c(phi, sigy, sigv, nu, rho)
      out <- sim_svp_leverage_ged_cpp(beta, p, n, burnin)
    } else {
      beta <- c(phi, sigy, sigv, nu)
      out <- sim_sv_ged_cpp(beta, p, n, burnin)
    }
  }

  # Coerce 1-column matrices from C++ to plain numeric vectors.
  list(
    y = as.numeric(out$y),
    h = as.numeric(out$h),
    z = as.numeric(out$z),
    v = as.numeric(out$v)
  )
}
