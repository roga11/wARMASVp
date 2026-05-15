# =========================================================================== #
# Information criteria for SV(p) AR-order selection
#
# Six criteria computed against an svp() fit:
#   1. AIC_YW       Yule-Walker projection-error AIC on log(y^2+del)
#   2. BIC_YW       Yule-Walker projection-error BIC on log(y^2+del)
#   3. BIC_Whittle  Whittle log-likelihood BIC at SV(p) signal-plus-noise spectrum
#   4. AIC_Kalman   Quasi-likelihood AIC using filter_svp() log-likelihood
#   5. BIC_Kalman   Quasi-likelihood BIC using filter_svp() log-likelihood
#   6. AICc_Kalman  Hurvich-Tsai (1989) small-sample-corrected AIC_Kalman
# =========================================================================== #

# Number of free parameters in an SV(p) model object
.svp_n_params <- function(p, errorType, leverage) {
  p + 2L +
    as.integer(errorType != "Gaussian") +
    as.integer(isTRUE(leverage))
}

# sigma_eps^2 = Var(log z_t^2) by error type. Used in spectral density.
.svp_sigma_eps2 <- function(errorType, nu = NULL) {
  switch(errorType,
    "Gaussian"  = pi^2 / 2,
    "Student-t" = trigamma(0.5) + trigamma(nu / 2),
    "GED"       = (2 / nu)^2 * trigamma(1 / nu),
    stop("Unknown errorType: ", errorType))
}

# Sample autocovariances of y*_t = log(y^2 + del) - mean
# Returns numeric vector of length max_lag + 1: gamma_y*(0..max_lag).
# Divisor 1/T (consistent with FFT periodogram / Whittle convention).
.svp_acov_ystar <- function(y, max_lag, del = 1e-10) {
  ystar <- log(y^2 + del)
  ystar <- ystar - mean(ystar)
  Tn <- length(ystar)
  vapply(0:max_lag,
         function(k) sum(ystar[(1L + k):Tn] * ystar[1:(Tn - k)]) / Tn,
         numeric(1))
}

# Yule-Walker projection-error variance for the AR(p) projection of y*_t
# onto its own past at lags 1..p, evaluated at the W-ARMA-SV phi estimate.
# sigma_pred^2 = gamma_y*(0) - phi' gamma_y*(1:p).
# Diagnostic; consistency in the SV(p) setting is demonstrated by simulation
# (Ahsan, Dufour & Rodriguez-Rondon, 2026) rather than analytically.
.svp_residual_var <- function(fit, y, del = 1e-10) {
  phi <- as.numeric(fit$phi)
  p <- length(phi)
  gv <- .svp_acov_ystar(y, p, del)
  v <- gv[1] - sum(phi * gv[2:(p + 1L)])
  if (!is.finite(v) || v <= 0) NA_real_ else v
}

# Hannan-Rissanen (1982, Biometrika) two-stage ARMA(p,p) residual variance.
# y*_t = log(y^2 + del) - mean has the form h_t + epsilon_t under SV(p), which
# is ARMA(p, p) (Harvey, Ruiz & Shephard 1994).  Stage 1 fits a long AR(L)
# pre-whitening regression; Stage 2 is OLS of y* on its lags 1..p plus the
# stage-1 residuals at lags 1..p.  Residual variance from Stage 2 is the
# Hannan-Rissanen estimator of the ARMA innovation variance.
#
# L_const: scaling for stage-1 AR length; default 1.5 (so L = floor(1.5*T^{1/3})).
# Returns scalar sigma_hr^2, T_eff (rows used in Stage 2), or NA if any guard fails.
.svp_hr_stage1 <- function(ystar, L_const = 1.5) {
  Tn <- length(ystar)
  L  <- max(5L, as.integer(floor(L_const * Tn^(1/3))))
  L  <- min(L, as.integer(floor(Tn / 4L)))
  if (Tn - L < 30L) return(list(eps = NULL, L = L))
  # Build long-AR design at lags 1..L
  X <- matrix(NA_real_, Tn - L, L)
  for (k in 1:L) X[, k] <- ystar[(L - k + 1L):(Tn - k)]
  y <- ystar[(L + 1L):Tn]
  X <- cbind(1, X)
  fit <- tryCatch(qr.coef(qr(X), y), error = function(e) NULL)
  if (is.null(fit) || any(!is.finite(fit))) return(list(eps = NULL, L = L))
  fitted <- as.numeric(X %*% fit)
  eps <- y - fitted
  list(eps = eps, L = L)
}

.svp_hr_stage2 <- function(ystar, eps_hat, L, p) {
  Tn <- length(ystar)
  start <- L + p + 1L  # need ystar lags 1..p AND eps_hat lags 1..p available
  if (start >= Tn) return(NA_real_)
  T_eff <- Tn - start + 1L
  if (T_eff < 2L * p + 5L) return(NA_real_)
  # Build the design: intercept, p AR lags of ystar, p MA lags of eps_hat
  yvec <- ystar[start:Tn]
  Xmat <- matrix(NA_real_, T_eff, 1L + 2L * p)
  Xmat[, 1L] <- 1.0
  for (k in 1:p) Xmat[, 1L + k]      <- ystar[(start - k):(Tn - k)]
  for (k in 1:p) Xmat[, 1L + p + k]  <- eps_hat[(start - L - k):(Tn - L - k)]
  beta <- tryCatch(qr.coef(qr(Xmat), yvec), error = function(e) NULL)
  if (is.null(beta) || any(!is.finite(beta))) return(NA_real_)
  resid <- yvec - as.numeric(Xmat %*% beta)
  ss <- sum(resid * resid)
  df <- T_eff - (2L * p + 1L)
  if (df <= 0L || !is.finite(ss) || ss <= 0) return(NA_real_)
  ss / df
}

# Hannan-Rissanen ARMA(p,p) residual variance — wrapper.
.svp_hr_residual_var <- function(fit, y, del = 1e-10, L_const = 1.5) {
  ystar <- log(y^2 + del); ystar <- ystar - mean(ystar)
  s1 <- .svp_hr_stage1(ystar, L_const)
  if (is.null(s1$eps)) return(list(sigma2 = NA_real_, T_eff = NA_integer_))
  p <- length(fit$phi)
  sigma2 <- .svp_hr_stage2(ystar, s1$eps, s1$L, p)
  if (is.na(sigma2)) return(list(sigma2 = NA_real_, T_eff = NA_integer_))
  T_eff <- length(ystar) - s1$L - p
  list(sigma2 = sigma2, T_eff = as.integer(T_eff))
}

# Whittle log-likelihood of y*_t at the SV(p) signal-plus-noise spectrum.
# Periodogram convention (no 1/(2*pi) factors): for consistency between
# periodogram and spectral density, both omit the 1/(2*pi). This is a
# constant additive shift -m * log(2*pi)/2 from the Brockwell-Davis form
# and does not affect ranking.
#
#   I(omega_j) = |FFT(y*)|^2 / T,   omega_j = 2*pi*j/T
#   f_y*(omega) = sigma_v^2 / |1 - sum phi_j e^{-i j omega}|^2 + sigma_eps^2(nu)
#   ell_W = -0.5 * sum_j [log f(omega_j) + I(omega_j) / f(omega_j)]
.svp_whittle_loglik <- function(fit, y, errorType, del = 1e-10) {
  ystar <- log(y^2 + del)
  ystar <- ystar - mean(ystar)
  Tn <- length(ystar)
  m  <- floor((Tn - 1L) / 2L)
  if (m < 1L) return(NA_real_)
  freqs <- 2 * pi * seq_len(m) / Tn

  Y <- stats::fft(ystar)
  I <- (Mod(Y)^2 / Tn)[2:(m + 1L)]

  phi <- as.numeric(fit$phi)
  sigv2 <- as.numeric(fit$sigv)^2
  nu <- if (errorType == "Gaussian") NA_real_ else as.numeric(fit$v)
  sigeps2 <- .svp_sigma_eps2(errorType, nu)

  # |phi(e^{-iw})|^2 with phi(z) = 1 - sum phi_j z^j
  poly_mag2 <- vapply(freqs, function(om) {
    z <- exp(-1i * om)
    Mod(1 - sum(phi * z^seq_along(phi)))^2
  }, numeric(1))

  f_y <- sigv2 / poly_mag2 + sigeps2

  if (any(!is.finite(f_y)) || any(f_y <= 0)) return(NA_real_)
  -0.5 * sum(log(f_y) + I / f_y)
}


#' Information criteria for SV(p) AR-order selection
#'
#' Computes information criteria for an \code{\link{svp}} fit to support
#' AR-order selection. Eight criteria are computable; \strong{four are
#' returned by default} — \code{BIC_Kalman}, \code{AIC_Kalman},
#' \code{BIC_HR}, \code{AIC_HR}. These span two families (state-space QML
#' and Hannan--Rissanen two-stage ARMA) and two penalty philosophies
#' (Schwarz-consistent BIC / Shibata-efficient AIC), and were selected as
#' the most informative criteria across the simulation grid of the
#' SVHT methodology paper (Ahsan, Dufour and Rodriguez-Rondon 2026).
#' The remaining four are available on request via the \code{criteria}
#' argument.
#'
#' \strong{Default criteria (returned by \code{svp_IC(fit)}):}
#' \itemize{
#'   \item \code{BIC_Kalman}, \code{AIC_Kalman}: \eqn{-2\,\hat\ell_K + k\log T}
#'         and \eqn{-2\,\hat\ell_K + 2k} where \eqn{\hat\ell_K} is the
#'         (quasi-)log-likelihood from \code{\link{filter_svp}}; default
#'         \code{filter_method = "mixture"} uses the Gaussian mixture
#'         Kalman filter (Kim, Shephard and Chib 1998). \code{BIC_Kalman}
#'         is the primary recommended criterion: Schwarz-consistent under
#'         the Bayes-optimal leverage proxy (see \code{proxy} argument)
#'         and strong finite-sample performance across the simulation
#'         grid (Ahsan, Dufour and Rodriguez-Rondon 2026).
#'         \code{AIC_Kalman} is Shibata-efficient
#'         and often selects larger \eqn{p} sooner at \eqn{p_{\mathrm{true}}
#'         \ge 2}.
#'   \item \code{BIC_HR}, \code{AIC_HR}: Hannan--Rissanen (1982) two-stage
#'         ARMA(\eqn{p,p}) criteria. Stage 1: long-AR pre-whitening at
#'         order \eqn{L = \lfloor 1.5\, T^{1/3}\rfloor} produces residuals
#'         \eqn{\hat\varepsilon_t}. Stage 2: OLS regression of \eqn{y_t^*}
#'         on AR lags \eqn{1{:}p} of \eqn{y_t^*} and MA lags \eqn{1{:}p}
#'         of \eqn{\hat\varepsilon_t} gives \eqn{\hat\sigma_u^2}. Then
#'         \eqn{T_{\mathrm{eff}} \log \hat\sigma_u^2 + \{2(2p{+}1),
#'         (2p{+}1)\log T_{\mathrm{eff}}\}}. Filter-free anchor, robust to
#'         mis-specification of the GMKF mixture. \code{BIC_HR} is
#'         Schwarz-consistent for ARMA(\eqn{p,p}) (Hannan & Rissanen
#'         1982; Pötscher 1989).
#' }
#'
#' \strong{Opt-in criteria (request via \code{criteria = ...}):}
#' \itemize{
#'   \item \code{AICc_Kalman}: \code{AIC_Kalman} with the Hurvich--Tsai
#'         (1989) small-sample correction \eqn{2k(k+1)/(T-k-1)}.
#'         Numerically equivalent to \code{AIC_Kalman} at \eqn{T \ge 500};
#'         use when \eqn{T < 500}.
#'   \item \code{BIC_Whittle}: \eqn{-2\,\hat\ell_W + k\log T} where
#'         \eqn{\hat\ell_W} is the Whittle log-likelihood evaluated at the
#'         SV(\eqn{p}) signal-plus-noise spectral density
#'         \eqn{f(\omega) = \sigma_v^2 / |1 - \sum_j \phi_j e^{-ij\omega}|^2
#'                        + \sigma_\varepsilon^2(\nu)}.
#'         Schwarz-consistent but collapses to \eqn{\hat p = 1} in 98--100\%
#'         of cells at \eqn{p_{\mathrm{true}} \ge 2} under near-unit-root
#'         persistence (Ahsan, Dufour and Rodriguez-Rondon 2026). Useful
#'         as a conservative
#'         diagnostic: a Whittle selection of \eqn{p > 1} is strong
#'         evidence against \eqn{p = 1}.
#'   \item \code{AIC_YW}, \code{BIC_YW}: \emph{Legacy / not recommended.}
#'         Yule--Walker projection-error criteria on
#'         \eqn{y_t^* = \log(y_t^2 + \delta) - \mu}, computed as
#'         \eqn{T \log \hat\sigma_{\mathrm{pred}}^2 + \{2k, k\log T\}}
#'         with the AR(\eqn{p}) projection-error variance. Under
#'         SV(\eqn{p}), \eqn{y_t^*} is ARMA(\eqn{p,p}) (not AR(\eqn{p})),
#'         so the AR projection error does not saturate at
#'         \eqn{p_{\mathrm{true}}} and the criteria are
#'         \emph{inconsistent}: the AR(\eqn{p}) projection-error variance
#'         keeps decreasing past \eqn{p_{\mathrm{true}}}, producing
#'         non-monotone (sometimes anti-Schwarz) behaviour in \eqn{T}.
#'         Simulation evidence: 0--29\% correct selection at
#'         \eqn{p_{\mathrm{true}} = 2} across all DGP cells and
#'         \eqn{T \le 10{,}000} (Ahsan, Dufour and Rodriguez-Rondon
#'         2026). Retained for paper-reproducibility of the
#'         documented failure-case results; \strong{use \code{BIC_HR} /
#'         \code{AIC_HR} for theoretically consistent AR-order selection.}
#' }
#'
#' Lower is better; \code{argmin} over a grid of candidate \code{p} (see
#' \code{\link{svp_AR_order}}) selects the AR order.
#'
#' @section Leverage invariance of non-Kalman criteria:
#' Leverage does not affect \code{AIC_YW}, \code{BIC_YW}, or
#' \code{BIC_Whittle}: under the W-ARMA-SV parameterization
#' \eqn{\mathrm{Cov}(v_t, \varepsilon_{t-1}) = 0} for all three error
#' distributions (odd-times-even moment symmetry), so the autocovariance
#' structure of \eqn{y_t^*} is invariant to the leverage parameter. The
#' \code{*_HR} and \code{*_Kalman} criteria do incorporate leverage
#' through the estimated \eqn{\delta_p} and the conditional state
#' innovation variance.
#'
#' @param fit Output of \code{\link{svp}}. Must carry the original \code{y}
#'   series (which \code{svp()} stores by default), \code{errorType}, and
#'   \code{leverage} fields.
#' @param criteria Character vector. Subset of
#'   \code{c("BIC_Kalman", "AIC_Kalman", "AICc_Kalman", "BIC_Whittle",
#'   "BIC_HR", "AIC_HR", "BIC_YW", "AIC_YW")}. Default returns the four
#'   recommended criteria: \code{c("BIC_Kalman", "AIC_Kalman", "BIC_HR",
#'   "AIC_HR")}. See the description for the rationale for each
#'   opt-in criterion.
#' @param filter_method Character. Filter method passed to
#'   \code{\link{filter_svp}} for \code{*_Kalman} criteria. One of
#'   \code{"mixture"} (default, recommended), \code{"corrected"}, or
#'   \code{"particle"}.
#' @param proxy Character. Leverage proxy passed to \code{\link{filter_svp}}.
#'   \code{"bayes_optimal"} (default here, unlike \code{filter_svp})
#'   removes the \eqn{O(T)} log-likelihood bias of the û-proxy under
#'   Student-t leverage and restores Schwarz consistency of \code{BIC_Kalman}.
#'   \code{"u"} reproduces the paper-faithful (Remark 3.5) Kalman likelihood;
#'   set this if you need IC values that match the filter's default behavior.
#' @param K Integer. Number of mixture components for
#'   \code{filter_method = "mixture"}. Default 7 (KSC).
#' @param M Integer. Number of particles for
#'   \code{filter_method = "particle"}. Default 1000. Ignored for other
#'   filter methods.
#' @param seed Integer. Random seed for the bootstrap particle filter.
#'   Default 42. Ignored for non-particle filters.
#' @param del Numeric. Small constant added inside \eqn{\log} to avoid
#'   \eqn{\log 0}. Default \code{1e-10}.
#'
#' @return Named numeric vector of the requested criteria. Lower is better.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' y <- sim_svp(2000, phi = 0.95, sigy = 1, sigv = 0.5)$y
#' fit1 <- svp(y, p = 1)
#' fit2 <- svp(y, p = 2)
#' svp_IC(fit1)
#' svp_IC(fit2)
#' }
#'
#' @references
#' Akaike, H. (1974). A new look at the statistical model identification.
#' \emph{IEEE Transactions on Automatic Control} 19(6), 716--723.
#' \doi{10.1109/TAC.1974.1100705}
#'
#' Schwarz, G. (1978). Estimating the dimension of a model.
#' \emph{Annals of Statistics} 6(2), 461--464.
#' \doi{10.1214/aos/1176344136}
#'
#' Shibata, R. (1976). Selection of the order of an autoregressive model
#' by Akaike's information criterion. \emph{Biometrika} 63(1), 117--126.
#' \doi{10.1093/biomet/63.1.117}
#'
#' Hannan, E. J. (1980). The estimation of the order of an ARMA process.
#' \emph{Annals of Statistics} 8(5), 1071--1081.
#' \doi{10.1214/aos/1176345144}
#'
#' Hannan, E. J., and Rissanen, J. (1982). Recursive estimation of mixed
#' autoregressive-moving average order. \emph{Biometrika} 69(1), 81--94.
#' \doi{10.1093/biomet/69.1.81}
#'
#' Pötscher, B. M. (1989). Model selection under nonstationarity:
#' Autoregressive models and stochastic linear regression models.
#' \emph{Annals of Statistics} 17(3), 1257--1274.
#' \doi{10.1214/aos/1176347267}
#'
#' Whittle, P. (1953). Estimation and information in stationary time series.
#' \emph{Arkiv f\"or Matematik} 2, 423--434. \doi{10.1007/BF02590998}
#'
#' Dunsmuir, W. (1979). A central limit theorem for parameter estimation in
#' stationary vector time series and its application to models for a signal
#' observed with noise. \emph{Annals of Statistics} 7(3), 490--506.
#' \doi{10.1214/aos/1176344671}
#'
#' Hurvich, C. M., and Tsai, C.-L. (1989). Regression and time series model
#' selection in small samples. \emph{Biometrika} 76(2), 297--307.
#' \doi{10.1093/biomet/76.2.297}
#'
#' Kim, S., Shephard, N., and Chib, S. (1998). Stochastic volatility:
#' likelihood inference and comparison with ARCH models.
#' \emph{Review of Economic Studies} 65(3), 361--393.
#' \doi{10.1111/1467-937X.00050}
#'
#' White, H. (1982). Maximum likelihood estimation of misspecified models.
#' \emph{Econometrica} 50(1), 1--25. \doi{10.2307/1912526}
#'
#' Ahsan, M. N., Dufour, J.-M., and Rodriguez-Rondon, G. (2026). Estimation and
#' inference for stochastic volatility models with heavy-tailed distributions.
#' Bank of Canada Staff Working Paper 2026-8. \doi{10.34989/swp-2026-8}
#'
#' @seealso \code{\link{svp_AR_order}}, \code{\link{svp}}, \code{\link{filter_svp}}
#' @export
svp_IC <- function(fit,
                   criteria = c("BIC_Kalman", "AIC_Kalman",
                                "BIC_HR", "AIC_HR"),
                   filter_method = c("mixture", "corrected", "particle"),
                   proxy = c("bayes_optimal", "u"),
                   K = 7L, M = 1000L, seed = 42L,
                   del = 1e-10) {
  if (!inherits(fit, c("svp", "svp_t", "svp_ged")))
    stop("'fit' must be an svp() output (class 'svp', 'svp_t' or 'svp_ged').")

  all_crit <- c("BIC_Kalman", "AIC_Kalman", "AICc_Kalman",
                "BIC_Whittle",
                "BIC_HR", "AIC_HR",
                "BIC_YW", "AIC_YW")
  criteria <- match.arg(criteria, choices = all_crit, several.ok = TRUE)
  filter_method <- match.arg(filter_method)
  proxy <- match.arg(proxy)

  errorType <- fit$errorType
  if (is.null(errorType))
    stop("'fit$errorType' is missing; refit with the current package version.")
  leverage <- isTRUE(fit$leverage)

  y <- as.numeric(fit$y)
  if (length(y) < 4L) stop("Series too short for IC computation.")

  Tn <- length(y)
  p  <- length(fit$phi)
  k  <- .svp_n_params(p, errorType, leverage)

  out <- setNames(rep(NA_real_, length(criteria)), criteria)

  # YW projection-error criteria (legacy — see roxygen note)
  if (any(c("AIC_YW", "BIC_YW") %in% criteria)) {
    sig2 <- .svp_residual_var(fit, y, del)
    if (is.finite(sig2)) {
      if ("AIC_YW" %in% criteria) out["AIC_YW"] <- Tn * log(sig2) + 2 * k
      if ("BIC_YW" %in% criteria) out["BIC_YW"] <- Tn * log(sig2) + k * log(Tn)
    }
  }

  # Hannan-Rissanen ARMA(p,p) criteria — recommended over YW for SV(p)
  if (any(c("AIC_HR", "BIC_HR") %in% criteria)) {
    hr <- .svp_hr_residual_var(fit, y, del)
    if (!is.na(hr$sigma2) && hr$sigma2 > 0 && !is.na(hr$T_eff) && hr$T_eff > 0) {
      Te <- hr$T_eff
      kHR <- 2L * length(fit$phi) + 1L  # AR(p) coefs + MA(p) coefs + intercept
      logS <- log(hr$sigma2)
      if ("AIC_HR" %in% criteria) out["AIC_HR"] <- Te * logS + 2 * kHR
      if ("BIC_HR" %in% criteria) out["BIC_HR"] <- Te * logS + kHR * log(Te)
    }
  }

  # Whittle BIC
  if ("BIC_Whittle" %in% criteria) {
    ell_W <- .svp_whittle_loglik(fit, y, errorType, del)
    if (is.finite(ell_W))
      out["BIC_Whittle"] <- -2 * ell_W + k * log(Tn)
  }

  # Kalman QML criteria
  if (any(c("AIC_Kalman", "BIC_Kalman", "AICc_Kalman") %in% criteria)) {
    filt <- tryCatch(
      filter_svp(fit, method = filter_method, proxy = proxy,
                 K = K, M = M, seed = seed, del = del),
      error = function(e) NULL
    )
    if (!is.null(filt) && is.finite(filt$loglik)) {
      m2L <- -2 * filt$loglik
      if ("AIC_Kalman"  %in% criteria) out["AIC_Kalman"]  <- m2L + 2 * k
      if ("BIC_Kalman"  %in% criteria) out["BIC_Kalman"]  <- m2L + k * log(Tn)
      if ("AICc_Kalman" %in% criteria) {
        if (Tn - k - 1L > 0L)
          out["AICc_Kalman"] <- m2L + 2 * k + 2 * k * (k + 1L) / (Tn - k - 1L)
      }
    }
  }

  out
}


#' AR-order selection sweep for SV(p)
#'
#' Convenience wrapper around \code{\link{svp_IC}}: fits \code{\link{svp}}
#' at each \code{p = 1, ..., pmax} and returns a matrix of information
#' criteria along with the argmin per criterion.
#'
#' @param y Numeric vector. Observed returns.
#' @param pmax Integer. Maximum AR order to consider. Default 6.
#' @param J Integer. Winsorizing parameter passed to \code{\link{svp}}.
#'   Default 10.
#' @param leverage Logical. Whether to estimate leverage. Default
#'   \code{FALSE}.
#' @param errorType Character. \code{"Gaussian"}, \code{"Student-t"}, or
#'   \code{"GED"}. Default \code{"Gaussian"}.
#' @param rho_type,del,trunc_lev,wDecay,logNu,sigvMethod,winsorize_eps Other
#'   arguments passed to \code{\link{svp}}.
#' @param filter_method Character. Filter method for \code{*_Kalman}
#'   criteria. Default \code{"mixture"}.
#' @param proxy Character. Leverage proxy. Default \code{"bayes_optimal"}
#'   for IC consistency under Student-t leverage. See \code{\link{svp_IC}}.
#' @param K,M,seed Filter arguments passed to \code{\link{filter_svp}}.
#' @param criteria Character vector of criteria to compute. Default
#'   returns the four recommended criteria: \code{c("BIC_Kalman",
#'   "AIC_Kalman", "BIC_HR", "AIC_HR")}. See \code{\link{svp_IC}} for the
#'   full set of eight valid names and the rationale for each opt-in
#'   criterion.
#'
#' @return A list with components:
#' \describe{
#'   \item{IC}{Numeric matrix, one row per criterion, one column per
#'     candidate \code{p} in \code{1:pmax}.}
#'   \item{argmin}{Named integer vector, one entry per criterion, giving the
#'     selected \code{p}. \code{NA_integer_} if all entries for that
#'     criterion are \code{NA}.}
#'   \item{fits}{List of length \code{pmax} containing the fitted
#'     \code{svp()} objects (or \code{NULL} if a fit failed).}
#' }
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' y <- sim_svp(2000, phi = 0.95, sigy = 1, sigv = 0.5)$y
#' res <- svp_AR_order(y, pmax = 4)
#' res$IC
#' res$argmin
#' }
#'
#' @seealso \code{\link{svp_IC}}, \code{\link{svp}}, \code{\link{filter_svp}}
#' @export
svp_AR_order <- function(y, pmax = 6L,
                         J = 10L, leverage = FALSE, errorType = "Gaussian",
                         rho_type = "pearson", del = 1e-10, trunc_lev = TRUE,
                         wDecay = FALSE, logNu = FALSE,
                         sigvMethod = "factored", winsorize_eps = 0L,
                         filter_method = "mixture",
                         proxy = c("bayes_optimal", "u"),
                         K = 7L, M = 1000L, seed = 42L,
                         criteria = c("BIC_Kalman", "AIC_Kalman",
                                      "BIC_HR", "AIC_HR")) {
  pmax <- as.integer(pmax)
  if (pmax < 1L) stop("'pmax' must be a positive integer.")

  all_crit <- c("BIC_Kalman", "AIC_Kalman", "AICc_Kalman",
                "BIC_Whittle",
                "BIC_HR", "AIC_HR",
                "BIC_YW", "AIC_YW")
  criteria <- match.arg(criteria, choices = all_crit, several.ok = TRUE)
  filter_method <- match.arg(filter_method, c("mixture", "corrected", "particle"))
  proxy <- match.arg(proxy)

  # Pre-allocate IC_mat in the order requested by the user (criteria arg).
  # Note: row names follow the user-specified ordering.
  IC_mat <- matrix(NA_real_, nrow = length(criteria), ncol = pmax,
                   dimnames = list(criteria, paste0("p=", seq_len(pmax))))
  fits <- vector("list", pmax)

  for (p in seq_len(pmax)) {
    fit <- tryCatch(
      svp(y, p = p, J = J, leverage = leverage, errorType = errorType,
          rho_type = rho_type, del = del, trunc_lev = trunc_lev,
          wDecay = wDecay, logNu = logNu,
          sigvMethod = sigvMethod, winsorize_eps = winsorize_eps),
      error = function(e) NULL
    )
    fits[[p]] <- fit
    if (!is.null(fit)) {
      IC_mat[, p] <- svp_IC(fit, criteria = criteria,
                             filter_method = filter_method,
                             proxy = proxy,
                             K = K, M = M, seed = seed, del = del)
    }
  }

  argmin <- vapply(seq_len(nrow(IC_mat)), function(i) {
    row <- IC_mat[i, ]
    if (all(is.na(row))) NA_integer_ else as.integer(which.min(row))
  }, integer(1))
  names(argmin) <- rownames(IC_mat)

  list(IC = IC_mat, argmin = argmin, fits = fits)
}
