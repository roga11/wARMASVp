# =========================================================================== #
# Hypothesis testing functions for SV(p) models
# =========================================================================== #


# =========================================================================== #
# AR Order Tests
# =========================================================================== #

#' LMC Test for AR Order in SV(p) Models
#'
#' Performs a Local Monte Carlo (LMC) test of the null hypothesis
#' \eqn{H_0: \phi_{p_0+1} = \cdots = \phi_p = 0} (i.e., that an SV(\eqn{p_0})
#' model is sufficient against an SV(\eqn{p}) alternative).
#'
#' @param y Numeric vector. Observed returns.
#' @param p_null Integer. Order under the null hypothesis.
#' @param p_alt Integer. Order under the alternative (\code{p_alt > p_null}).
#' @param J Integer. Winsorizing parameter. Default 10.
#' @param N Integer. Number of Monte Carlo replications. Default 99.
#' @param burnin Integer. Burn-in for simulation. Default 500.
#' @param del Numeric. Small constant for log transformation. Default \code{1e-10}.
#' @param wDecay Logical. Use decaying weights. Default \code{FALSE}.
#' @param Bartlett Logical. If \code{TRUE}, use Bartlett kernel HAC weighting
#'   matrix for a GMM-LRT-type test statistic. If \code{FALSE} (default), use
#'   the sum of squared extra AR coefficients.
#'
#' @return An object of class \code{"svp_test"}, a list containing:
#' \describe{
#'   \item{s0}{Test statistic from observed data.}
#'   \item{sN}{Simulated null distribution (vector of length N).}
#'   \item{pval}{Monte Carlo p-value.}
#'   \item{test_type}{Character string identifying the test.}
#'   \item{null_param}{Name of the parameter(s) tested.}
#'   \item{null_value}{Value(s) under the null hypothesis.}
#'   \item{call}{The matched call.}
#' }
#'
#' @details
#' When \code{Bartlett = FALSE} (default), the test statistic is
#' \eqn{T \sum_{j=p_0+1}^{p} \hat\phi_j^2}, i.e., the sum of squared extra
#' AR coefficients scaled by sample size.
#'
#' When \code{Bartlett = TRUE}, the test statistic is based on the GMM-LRT
#' approach with a Bartlett kernel HAC weighting matrix:
#' \eqn{S = T \times (M_{H_0} - M_{H_1})}, where \eqn{M} denotes the
#' GMM criterion evaluated at the null and alternative estimates. Simulations
#' yielding negative test statistics are discarded and re-drawn.
#'
#' @examples
#' \donttest{
#' y <- sim_svp(1000, phi = 0.95, sigy = 1, sigv = 0.2)
#' test <- lmc_ar(y, p_null = 1, p_alt = 2, J = 10, N = 49)
#' print(test)
#' }
#'
#' @export
lmc_ar <- function(y, p_null, p_alt, J = 10, N = 99, burnin = 500,
                   del = 1e-10, wDecay = FALSE, Bartlett = FALSE) {
  cl <- match.call()
  y_vec <- as.numeric(y)
  Tsize <- length(y_vec)
  if (p_null >= p_alt) stop("p_alt must be greater than p_null.")
  if (p_null < 1) stop("p_null must be >= 1.")
  # Estimate models
  mdl_alt <- svp(y_vec, p = p_alt, J = J, leverage = FALSE, del = del, wDecay = wDecay)
  mdl_null_est <- svp(y_vec, p = p_null, J = J, leverage = FALSE, del = del, wDecay = wDecay)
  # Compute test statistic
  if (isTRUE(Bartlett)) {
    M_null <- LRT_moment_ar_Amat(y_vec, mdl_null_est, del = del, Bartlett = TRUE)
    M_alt  <- LRT_moment_ar_Amat(y_vec, mdl_alt, del = del, Bartlett = TRUE)
    s0 <- Tsize * (M_null - M_alt)
  } else {
    phi_extra <- mdl_alt$phi[(p_null + 1):p_alt]
    s0 <- Tsize * sum(phi_extra^2)
  }
  # Simulate null distribution
  betasim_null <- c(mdl_null_est$phi, mdl_null_est$sigy, mdl_null_est$sigv)
  sN <- .simnull_ar(betasim_null, p_null, p_alt, J, Tsize, N, burnin,
                    del, wDecay, Bartlett)
  pval <- (N + 1 - sum(s0 >= sN)) / (N + 1)
  test_label <- if (isTRUE(Bartlett)) {
    sprintf("LMC AR Order Bartlett (p0=%d vs p=%d)", p_null, p_alt)
  } else {
    sprintf("LMC AR Order (p0=%d vs p=%d)", p_null, p_alt)
  }
  out <- list(s0 = s0, sN = as.numeric(sN), pval = pval,
              test_type = test_label,
              null_param = paste0("phi_", (p_null + 1):p_alt),
              null_value = rep(0, p_alt - p_null),
              call = cl)
  class(out) <- "svp_test"
  return(out)
}


#' MMC Test for AR Order in SV(p) Models
#'
#' Performs a Maximized Monte Carlo (MMC) test of
#' \eqn{H_0: \phi_{p_0+1} = \cdots = \phi_p = 0}
#' by maximizing the MC p-value over nuisance parameters
#' (\eqn{\phi_1, \ldots, \phi_{p_0}, \sigma_y, \sigma_v}).
#'
#' @inheritParams lmc_ar
#' @param eps Numeric vector. Half-width of search region around estimated
#'   nuisance parameters. Default \code{rep(0.3, p_null+2)}.
#' @param threshold Numeric. Target p-value. Default 1.
#' @param method Character. Optimization method: \code{"pso"} or \code{"GenSA"}.
#'   Default \code{"pso"}.
#' @param maxit Integer. Maximum iterations/evaluations. Default depends on method.
#'
#' @return A list with optimization output including \code{value} (maximized p-value)
#'   and \code{par} (nuisance parameters at the maximum).
#'
#' @examples
#' \donttest{
#' y <- sim_svp(1000, phi = 0.95, sigy = 1, sigv = 0.2)
#' mmc <- mmc_ar(y, p_null = 1, p_alt = 2, J = 10, N = 19,
#'               method = "pso", maxit = 10)
#' mmc$value
#' }
#'
#' @export
mmc_ar <- function(y, p_null, p_alt, J = 10, N = 99, burnin = 500,
                   eps = NULL, threshold = 1, method = "pso", maxit = NULL,
                   del = 1e-10, wDecay = FALSE, Bartlett = FALSE) {
  cl <- match.call()
  y_vec <- as.numeric(y)
  Tsize <- length(y_vec)
  if (p_null >= p_alt) stop("p_alt must be greater than p_null.")
  if (p_null < 1) stop("p_null must be >= 1.")
  # Estimate under alternative for test statistic
  mdl_alt <- svp(y_vec, p = p_alt, J = J, leverage = FALSE, del = del, wDecay = wDecay)
  mdl_null_est <- svp(y_vec, p = p_null, J = J, leverage = FALSE, del = del, wDecay = wDecay)
  theta_0 <- c(mdl_null_est$phi, mdl_null_est$sigy, mdl_null_est$sigv)
  n_nuisance <- p_null + 2
  if (is.null(eps)) eps <- rep(0.3, n_nuisance)
  # Bounds for nuisance parameters
  lower <- c(pmax(theta_0[1:p_null] - eps[1:p_null], rep(-0.999, p_null)),
             max(theta_0[p_null + 1] - eps[p_null + 1], 0.01),
             max(theta_0[p_null + 2] - eps[p_null + 2], 0.01))
  upper <- c(pmin(theta_0[1:p_null] + eps[1:p_null], rep(0.999, p_null)),
             theta_0[p_null + 1] + eps[p_null + 1],
             theta_0[p_null + 2] + eps[p_null + 2])
  # Test statistic
  if (isTRUE(Bartlett)) {
    M_null <- LRT_moment_ar_Amat(y_vec, mdl_null_est, del = del, Bartlett = TRUE)
    M_alt  <- LRT_moment_ar_Amat(y_vec, mdl_alt, del = del, Bartlett = TRUE)
    s0 <- Tsize * (M_null - M_alt)
  } else {
    phi_extra <- mdl_alt$phi[(p_null + 1):p_alt]
    s0 <- Tsize * sum(phi_extra^2)
  }
  out <- .run_mmc_optimizer(method, theta_0, .mmc_pval_ar, lower, upper,
                            threshold, maxit,
                            y = y_vec, p_null = p_null, p_alt = p_alt,
                            j = J, N = N, s0 = s0, ini = burnin,
                            del = del, wDecay = wDecay, Bartlett = Bartlett)
  out$value <- -out$value
  out$s0 <- s0
  out$call <- cl
  return(out)
}


# =========================================================================== #
# Leverage Tests
# =========================================================================== #

#' LMC Test for Leverage in SV(p) Models
#'
#' Performs a Local Monte Carlo (LMC) test of the null hypothesis
#' \eqn{H_0: \rho = \rho_0} (typically \eqn{\rho_0 = 0}, i.e., no leverage)
#' using a GMM likelihood-ratio type statistic.
#'
#' @param y Numeric vector. Observed returns.
#' @param p Integer. Order of the volatility process. Default 1.
#' @param J Integer. Winsorizing parameter. Default 10.
#' @param N Integer. Number of Monte Carlo replications. Default 99.
#' @param rho_null Numeric. Value of \eqn{\rho} under the null. Default 0.
#' @param burnin Integer. Burn-in for simulation. Default 500.
#' @param rho_type Character. Correlation type. Default \code{"pearson"}.
#' @param del Numeric. Small constant for log transformation. Default \code{1e-10}.
#' @param trunc_lev Logical. Truncate leverage correlation estimate to
#'   \code{[-1,1]}. Default \code{TRUE}.
#' @param wDecay Logical. Use decaying weights. Default \code{FALSE}.
#' @param Bartlett Logical. If \code{TRUE}, use Bartlett kernel HAC weighting
#'   matrix. If \code{FALSE}, use identity matrix. Default \code{FALSE}.
#'
#' @return An object of class \code{"svp_test"}, a list containing:
#' \describe{
#'   \item{s0}{Test statistic from observed data.}
#'   \item{sN}{Simulated null distribution (vector of length N).}
#'   \item{pval}{Monte Carlo p-value.}
#'   \item{test_type}{Character string identifying the test.}
#'   \item{null_param}{Name of the parameter tested.}
#'   \item{null_value}{Value under the null hypothesis.}
#'   \item{call}{The matched call.}
#' }
#'
#' @examples
#' \donttest{
#' y <- sim_svp(1000, phi = 0.95, sigy = 1, sigv = 0.2, leverage = TRUE, rho = -0.3)$y
#' fit <- svp(y, p = 1, leverage = TRUE)
#' test <- lmc_lev(y, p = 1, J = 10, N = 99)
#' print(test)
#' }
#'
#' @export
lmc_lev <- function(y, p = 1, J = 10, N = 99, rho_null = 0,
                    burnin = 500, rho_type = "pearson", del = 1e-10,
                    trunc_lev = TRUE, wDecay = FALSE,
                    Bartlett = FALSE) {
  cl <- match.call()
  y_vec <- as.numeric(y)
  y_mat <- as.matrix(y_vec)
  Tsize <- length(y_vec)
  # Estimate model under alternative
  mdl_alt <- svp(y_vec, p, J, leverage = TRUE, rho_type = rho_type, del = del,
                 trunc_lev = trunc_lev, wDecay = wDecay)
  if (!is.null(mdl_alt$rho_type) && !is.na(mdl_alt$rho_type) &&
      mdl_alt$rho_type != rho_type) {
    warning("rho_type in estimated model differs from specified rho_type.")
  }
  mdl_null <- mdl_alt
  mdl_null$rho <- rho_null
  # Compute test statistic from observed data
  if (isTRUE(Bartlett)) {
    s0_tmp <- Tsize * (LRT_moment_lev_svp_Amat(y_mat, mdl_null, rho_type, del, TRUE) -
                         LRT_moment_lev_svp_Amat(y_mat, mdl_alt, rho_type, del, TRUE))
  } else {
    Amat <- diag(p + 3)
    s0_tmp <- Tsize * (LRT_moment_lev_svp(y_mat, mdl_null, Amat, rho_type, del) -
                         LRT_moment_lev_svp(y_mat, mdl_alt, Amat, rho_type, del))
  }
  if (is.na(s0_tmp) || s0_tmp < 0) {
    stop("Test statistic is not valid (NA or negative).")
  }
  s0 <- s0_tmp
  # Simulate null distribution
  betasim_null <- c(mdl_alt$phi, mdl_alt$sigy, mdl_alt$sigv, rho_null)
  if (isTRUE(Bartlett)) {
    sN <- .simnull_Amat(betasim_null, rho_null, p, J, Tsize, N, burnin,
                        rho_type, del, TRUE, wDecay = wDecay,
                        trunc_lev = trunc_lev)
  } else {
    sN <- .simnull(betasim_null, rho_null, p, J, Tsize, N, burnin,
                   Amat, rho_type, del, wDecay = wDecay,
                   trunc_lev = trunc_lev)
  }
  # Compute p-value
  pval <- (N + 1 - sum(s0 >= sN)) / (N + 1)
  out <- list(s0 = s0, sN = as.numeric(sN), pval = pval,
              test_type = "LMC Leverage",
              null_param = "rho", null_value = rho_null,
              call = cl)
  class(out) <- "svp_test"
  return(out)
}


#' MMC Test for Leverage in SV(p) Models
#'
#' Performs a Maximized Monte Carlo (MMC) test of the null hypothesis
#' \eqn{H_0: \rho = \rho_0} by maximizing the MC p-value over nuisance
#' parameters (phi, sigma_y, sigma_v).
#'
#' @inheritParams lmc_lev
#' @param eps Numeric vector. Half-width of the search region around the
#'   estimated nuisance parameters. Default \code{rep(0.3, p+2)}.
#' @param threshold Numeric. Target p-value (optimization stops if reached).
#'   Default 1.
#' @param method Character. Optimization method: \code{"pso"} (particle swarm),
#'   \code{"GenSA"} (generalized simulated annealing), or \code{"gridSearch"}.
#'   Default \code{"pso"}.
#' @param maxit Integer or list. Maximum iterations/evaluations for the
#'   optimizer. Default depends on method.
#'
#' @return A list with the optimization output including:
#' \describe{
#'   \item{value}{Maximized p-value.}
#'   \item{par}{Nuisance parameter values at the maximum.}
#' }
#' Additional fields depend on the optimization method used.
#'
#' @examples
#' \donttest{
#' y <- sim_svp(1000, phi = 0.95, sigy = 1, sigv = 0.2, leverage = TRUE, rho = -0.3)$y
#' mmc <- mmc_lev(y, p = 1, J = 10, N = 19, method = "pso", maxit = 10)
#' mmc$value
#' }
#'
#' @export
mmc_lev <- function(y, p = 1, J = 10, N = 99, rho_null = 0,
                    burnin = 500, eps = NULL, threshold = 1,
                    method = "pso", maxit = NULL,
                    rho_type = "pearson", del = 1e-10,
                    trunc_lev = TRUE, wDecay = FALSE,
                    Bartlett = FALSE) {
  cl <- match.call()
  y_vec <- as.numeric(y)
  y_mat <- as.matrix(y_vec)
  Tsize <- length(y_vec)
  # Estimate model under alternative
  mdl_alt <- svp(y_vec, p, J, leverage = TRUE, rho_type = rho_type, del = del,
                 trunc_lev = trunc_lev, wDecay = wDecay)
  theta_0 <- c(mdl_alt$phi, mdl_alt$sigy, mdl_alt$sigv)
  if (is.null(eps)) {
    eps <- rep(0.3, length(theta_0))
  }
  # Define search bounds
  lower <- c(pmax(theta_0[1:p] - eps[1:p], rep(-0.999, p)),
             max(theta_0[p + 1] - eps[p + 1], 0.01),
             max(theta_0[p + 2] - eps[p + 2], 0.01))
  upper <- c(pmin(theta_0[1:p] + eps[1:p], rep(0.999, p)),
             theta_0[p + 1] + eps[p + 1],
             theta_0[p + 2] + eps[p + 2])
  # Compute test statistic from observed data
  mdl_null <- mdl_alt
  mdl_null$rho <- rho_null
  if (isTRUE(Bartlett)) {
    s0_tmp <- Tsize * (LRT_moment_lev_svp_Amat(y_mat, mdl_null, rho_type, del, TRUE) -
                         LRT_moment_lev_svp_Amat(y_mat, mdl_alt, rho_type, del, TRUE))
  } else {
    Amat <- diag(p + 3)
    s0_tmp <- Tsize * (LRT_moment_lev_svp(y_mat, mdl_null, Amat, rho_type, del) -
                         LRT_moment_lev_svp(y_mat, mdl_alt, Amat, rho_type, del))
  }
  if (is.na(s0_tmp) || s0_tmp < 0) {
    stop("Test statistic is not valid (NA or negative).")
  }
  # Choose p-value function based on Bartlett
  if (isTRUE(Bartlett)) {
    pval_fn <- .mmc_pval_lev_Amat
  } else {
    pval_fn <- .mmc_pval_lev
  }
  # Optimize
  if (method == "gridSearch") {
    if (!requireNamespace("NMOF", quietly = TRUE)) {
      stop("Package 'NMOF' required for gridSearch method. Install it or use method='pso'.")
    }
    if (is.null(maxit)) {
      maxit <- as.list(rep(5, p + 2))
    }
    extra_args <- list(y = y_mat, j = J, N = N, mdl_alt = mdl_alt,
                       rho_null = rho_null, ini = burnin,
                       rho_type = rho_type, del = del,
                       wDecay = wDecay, trunc_lev = trunc_lev)
    if (isTRUE(Bartlett)) extra_args$Bartlett <- TRUE
    if (!isTRUE(Bartlett)) extra_args$Amat <- diag(p + 3)
    out <- do.call(NMOF::gridSearch,
                   c(list(fun = pval_fn, lower = lower, upper = upper,
                          npar = length(theta_0), n = maxit,
                          printDetail = TRUE, asList = FALSE),
                     extra_args))
  } else if (method == "GenSA") {
    if (!requireNamespace("GenSA", quietly = TRUE)) {
      stop("Package 'GenSA' required. Install it or use method='pso'.")
    }
    if (is.null(maxit)) maxit <- 100
    extra_args <- list(y = y_mat, j = J, N = N, mdl_alt = mdl_alt,
                       rho_null = rho_null, ini = burnin,
                       rho_type = rho_type, del = del,
                       wDecay = wDecay, trunc_lev = trunc_lev)
    if (isTRUE(Bartlett)) extra_args$Bartlett <- TRUE
    if (!isTRUE(Bartlett)) extra_args$Amat <- diag(p + 3)
    out <- do.call(GenSA::GenSA,
                   c(list(par = theta_0, fn = pval_fn,
                          lower = lower, upper = upper,
                          control = list(threshold.stop = -threshold,
                                         max.call = maxit, verbose = TRUE)),
                     extra_args))
  } else if (method == "pso") {
    if (!requireNamespace("pso", quietly = TRUE)) {
      stop("Package 'pso' required. Install it or use method='GenSA'.")
    }
    if (is.null(maxit)) maxit <- 100
    extra_args <- list(y = y_mat, j = J, N = N, mdl_alt = mdl_alt,
                       rho_null = rho_null, ini = burnin,
                       rho_type = rho_type, del = del,
                       wDecay = wDecay, trunc_lev = trunc_lev)
    if (isTRUE(Bartlett)) extra_args$Bartlett <- TRUE
    if (!isTRUE(Bartlett)) extra_args$Amat <- diag(p + 3)
    out <- do.call(pso::psoptim,
                   c(list(par = theta_0, fn = pval_fn,
                          lower = lower, upper = upper,
                          control = list(abstol = -threshold,
                                         maxf = maxit, trace = 1,
                                         REPORT = 1, trace.stats = TRUE)),
                     extra_args))
  } else {
    stop("method must be 'pso', 'GenSA', or 'gridSearch'.")
  }
  out$value <- -out$value
  out$call <- cl
  return(out)
}


# =========================================================================== #
# Heavy-Tail Tests (Student-t and GED)
# =========================================================================== #

#' LMC Test for Student-t Tail Parameter
#'
#' Performs a Local Monte Carlo (LMC) test of the null hypothesis
#' \eqn{H_0: \nu = \nu_0} for the degrees of freedom parameter in an
#' SV(1) model with Student-t errors. Testing \eqn{\nu_0 = \infty}
#' (or a large value) corresponds to testing for normality.
#'
#' @param y Numeric vector. Observed returns.
#' @param J Integer. Winsorizing parameter. Default 10.
#' @param N Integer. Number of Monte Carlo replications. Default 99.
#' @param nu_null Numeric. Value of \eqn{\nu} under the null hypothesis.
#' @param burnin Integer. Burn-in for simulation. Default 500.
#' @param del Numeric. Small constant for log transformation. Default \code{1e-10}.
#' @param wDecay Logical. Use decaying weights. Default \code{FALSE}.
#' @param Bartlett Logical. Use Bartlett kernel HAC for weighting matrix.
#'   Default \code{TRUE}.
#' @param Amat Weighting matrix specification. \code{NULL} for identity,
#'   \code{"Weighted"} for data-driven HAC, or a 4x4 matrix. Default \code{NULL}.
#' @param logNu Logical. Use log-space for nu estimation. Default \code{TRUE}.
#'
#' @return An object of class \code{"svp_test"}.
#'
#' @examples
#' \donttest{
#' y <- sim_svp(1000, phi = 0.95, sigy = 1, sigv = 0.2, errorType = "Student-t", nu = 5)
#' test <- lmc_t(y, J = 10, N = 49, nu_null = 5)
#' print(test)
#' }
#'
#' @export
lmc_t <- function(y, J = 10, N = 99, nu_null, burnin = 500,
                  del = 1e-10, wDecay = FALSE, Bartlett = TRUE,
                  Amat = NULL, logNu = TRUE) {
  cl <- match.call()
  y_vec <- as.numeric(y)
  y_mat <- as.matrix(y_vec)
  Tsize <- length(y_vec)
  p <- 1L
  # Estimate model under alternative
  mdl_alt <- svp(y_vec, p = 1, J = J, errorType = "Student-t", del = del,
                 logNu = logNu, wDecay = wDecay)
  mdl_null <- mdl_alt
  mdl_null$v <- nu_null
  # Handle Amat specification
  wa <- .parse_Amat(Amat, p)
  Amat_mat <- wa$Amat
  WAmat <- wa$WAmat
  # Compute test statistic
  s0_tmp <- Tsize * (LRT_moment_t(y_mat, mdl_null, Amat_mat, WAmat, del, Bartlett) -
                       LRT_moment_t(y_mat, mdl_alt, Amat_mat, WAmat, del, Bartlett))
  if (is.na(s0_tmp) || s0_tmp < 0) {
    stop("Test statistic is not valid (NA or negative).")
  }
  s0 <- s0_tmp
  # Simulate null distribution
  betasim_null <- c(mdl_alt$phi, mdl_alt$sigy, mdl_alt$sigv, nu_null)
  sN <- .simnull_t(betasim_null, nu_null, J, Tsize, N, burnin,
                   Amat_mat, del, WAmat, Bartlett, logNu,
                   wDecay = wDecay)
  pval <- (N + 1 - sum(s0 >= sN)) / (N + 1)
  out <- list(s0 = s0, sN = as.numeric(sN), pval = pval,
              test_type = "LMC Student-t",
              null_param = "nu", null_value = nu_null,
              call = cl)
  class(out) <- "svp_test"
  return(out)
}


#' LMC Test for GED Shape Parameter
#'
#' Performs a Local Monte Carlo (LMC) test of the null hypothesis
#' \eqn{H_0: \nu = \nu_0} for the shape parameter in an SV(1) model
#' with GED errors. Testing \eqn{\nu_0 = 2} corresponds to testing normality.
#'
#' @inheritParams lmc_t
#'
#' @return An object of class \code{"svp_test"}.
#'
#' @examples
#' \donttest{
#' y <- sim_svp(1000, phi = 0.95, sigy = 1, sigv = 0.2, errorType = "GED", nu = 1.5)
#' test <- lmc_ged(y, J = 10, N = 49, nu_null = 2)
#' print(test)
#' }
#'
#' @export
lmc_ged <- function(y, J = 10, N = 99, nu_null, burnin = 500,
                    del = 1e-10, wDecay = FALSE, Bartlett = TRUE,
                    Amat = NULL) {
  cl <- match.call()
  y_vec <- as.numeric(y)
  y_mat <- as.matrix(y_vec)
  Tsize <- length(y_vec)
  p <- 1L
  mdl_alt <- svp(y_vec, p = 1, J = J, errorType = "GED", del = del,
                 wDecay = wDecay)
  mdl_null <- mdl_alt
  mdl_null$v <- nu_null
  wa <- .parse_Amat(Amat, p)
  Amat_mat <- wa$Amat
  WAmat <- wa$WAmat
  s0_tmp <- Tsize * (LRT_moment_ged(y_mat, mdl_null, Amat_mat, WAmat, del, Bartlett) -
                       LRT_moment_ged(y_mat, mdl_alt, Amat_mat, WAmat, del, Bartlett))
  if (is.na(s0_tmp) || s0_tmp < 0) {
    stop("Test statistic is not valid (NA or negative).")
  }
  s0 <- s0_tmp
  betasim_null <- c(mdl_alt$phi, mdl_alt$sigy, mdl_alt$sigv, nu_null)
  sN <- .simnull_ged(betasim_null, nu_null, J, Tsize, N, burnin,
                     Amat_mat, del, WAmat, Bartlett, wDecay = wDecay)
  pval <- (N + 1 - sum(s0 >= sN)) / (N + 1)
  out <- list(s0 = s0, sN = as.numeric(sN), pval = pval,
              test_type = "LMC GED",
              null_param = "nu", null_value = nu_null,
              call = cl)
  class(out) <- "svp_test"
  return(out)
}


#' MMC Test for Student-t Tail Parameter
#'
#' Performs a Maximized Monte Carlo (MMC) test of \eqn{H_0: \nu = \nu_0}
#' by maximizing the MC p-value over nuisance parameters (phi, sigma_y, sigma_v).
#'
#' @inheritParams lmc_t
#' @param eps Numeric vector. Half-width of search region. Default \code{rep(0.3, 3)}.
#' @param threshold Numeric. Target p-value. Default 1.
#' @param method Character. Optimization method: \code{"pso"} or \code{"GenSA"}.
#'   Default \code{"pso"}.
#' @param maxit Integer. Maximum iterations/evaluations. Default depends on method.
#'
#' @return A list with optimization output including \code{value} (maximized p-value)
#'   and \code{par} (nuisance parameters at the maximum).
#'
#' @examples
#' \donttest{
#' y <- sim_svp(1000, phi = 0.95, sigy = 1, sigv = 0.2, errorType = "Student-t", nu = 5)
#' mmc <- mmc_t(y, J = 10, N = 19, nu_null = 5, method = "pso", maxit = 10)
#' mmc$value
#' }
#'
#' @export
mmc_t <- function(y, J = 10, N = 99, nu_null, burnin = 500,
                  eps = NULL, threshold = 1, method = "pso", maxit = NULL,
                  del = 1e-10, wDecay = FALSE, Bartlett = TRUE,
                  Amat = NULL, logNu = TRUE) {
  cl <- match.call()
  y_vec <- as.numeric(y)
  y_mat <- as.matrix(y_vec)
  Tsize <- length(y_vec)
  mdl_alt <- svp(y_vec, p = 1, J = J, errorType = "Student-t", del = del,
                 logNu = logNu, wDecay = wDecay)
  p <- length(mdl_alt$phi)
  theta_0 <- c(mdl_alt$phi, mdl_alt$sigy, mdl_alt$sigv)
  if (is.null(eps)) eps <- rep(0.3, length(theta_0))
  wa <- .parse_Amat(Amat, p)
  Amat_mat <- wa$Amat
  WAmat <- wa$WAmat
  lower <- c(pmax(theta_0[1:p] - eps[1:p], rep(-0.999, p)),
             max(theta_0[p + 1] - eps[p + 1], 0.01),
             max(theta_0[p + 2] - eps[p + 2], 0.01))
  upper <- c(pmin(theta_0[1:p] + eps[1:p], rep(0.999, p)),
             theta_0[p + 1] + eps[p + 1],
             theta_0[p + 2] + eps[p + 2])
  mdl_null <- mdl_alt
  mdl_null$v <- nu_null
  s0_tmp <- Tsize * (LRT_moment_t(y_mat, mdl_null, Amat_mat, WAmat, del, Bartlett) -
                       LRT_moment_t(y_mat, mdl_alt, Amat_mat, WAmat, del, Bartlett))
  if (!is.finite(s0_tmp) || s0_tmp < 0) {
    stop("Test statistic is not valid (NA or negative).")
  }
  out <- .run_mmc_optimizer(method, theta_0, .mmc_pval_t, lower, upper,
                            threshold, maxit,
                            y = y_mat, j = J, N = N, mdl_alt = mdl_alt,
                            nu_null = nu_null, ini = burnin, Amat = Amat_mat,
                            del = del, Bartlett = Bartlett, logNu = logNu,
                            wDecay = wDecay)
  out$value <- -out$value
  out$s0 <- s0_tmp
  out$call <- cl
  return(out)
}


#' MMC Test for GED Shape Parameter
#'
#' Performs a Maximized Monte Carlo (MMC) test of \eqn{H_0: \nu = \nu_0}
#' for the GED shape parameter.
#'
#' @inheritParams lmc_ged
#' @param eps Numeric vector. Half-width of search region. Default \code{rep(0.3, 3)}.
#' @param threshold Numeric. Target p-value. Default 1.
#' @param method Character. Optimization method: \code{"pso"} or \code{"GenSA"}.
#'   Default \code{"pso"}.
#' @param maxit Integer. Maximum iterations/evaluations. Default depends on method.
#'
#' @return A list with optimization output including \code{value} (maximized p-value)
#'   and \code{par} (nuisance parameters at the maximum).
#'
#' @examples
#' \donttest{
#' y <- sim_svp(1000, phi = 0.95, sigy = 1, sigv = 0.2, errorType = "GED", nu = 1.5)
#' mmc <- mmc_ged(y, J = 10, N = 19, nu_null = 2, method = "pso", maxit = 10)
#' mmc$value
#' }
#'
#' @export
mmc_ged <- function(y, J = 10, N = 99, nu_null, burnin = 500,
                    eps = NULL, threshold = 1, method = "pso", maxit = NULL,
                    del = 1e-10, wDecay = FALSE, Bartlett = TRUE,
                    Amat = NULL) {
  cl <- match.call()
  y_vec <- as.numeric(y)
  y_mat <- as.matrix(y_vec)
  Tsize <- length(y_vec)
  mdl_alt <- svp(y_vec, p = 1, J = J, errorType = "GED", del = del,
                 wDecay = wDecay)
  p <- length(mdl_alt$phi)
  theta_0 <- c(mdl_alt$phi, mdl_alt$sigy, mdl_alt$sigv)
  if (is.null(eps)) eps <- rep(0.3, length(theta_0))
  wa <- .parse_Amat(Amat, p)
  Amat_mat <- wa$Amat
  WAmat <- wa$WAmat
  lower <- c(pmax(theta_0[1:p] - eps[1:p], rep(-0.999, p)),
             max(theta_0[p + 1] - eps[p + 1], 0.01),
             max(theta_0[p + 2] - eps[p + 2], 0.01))
  upper <- c(pmin(theta_0[1:p] + eps[1:p], rep(0.999, p)),
             theta_0[p + 1] + eps[p + 1],
             theta_0[p + 2] + eps[p + 2])
  mdl_null <- mdl_alt
  mdl_null$v <- nu_null
  s0_tmp <- Tsize * (LRT_moment_ged(y_mat, mdl_null, Amat_mat, WAmat, del, Bartlett) -
                       LRT_moment_ged(y_mat, mdl_alt, Amat_mat, WAmat, del, Bartlett))
  if (!is.finite(s0_tmp) || s0_tmp < 0) {
    stop("Test statistic is not valid (NA or negative).")
  }
  out <- .run_mmc_optimizer(method, theta_0, .mmc_pval_ged, lower, upper,
                            threshold, maxit,
                            y = y_mat, j = J, N = N, mdl_alt = mdl_alt,
                            nu_null = nu_null, ini = burnin, Amat = Amat_mat,
                            del = del, Bartlett = Bartlett, wDecay = wDecay)
  out$value <- -out$value
  out$s0 <- s0_tmp
  out$call <- cl
  return(out)
}


