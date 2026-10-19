# RV-like data: log y^2 = c + w_t + eps_t with small Gaussian eps (no log-chi2)
.sim_rv <- function(n, phi = 0.95, sigv = 0.2, s_eps = 0.45, c = -1, seed = 1) {
  set.seed(seed)
  w <- as.numeric(arima.sim(list(ar = phi), n = n, sd = sigv))
  sqrt(exp(c + w + rnorm(n, 0, s_eps)))
}

test_that("estimated mode: fields, identities, and sigma_v from the estimate", {
  y <- .sim_rv(3000)
  fit <- svp(y, p = 1, sigeps = "estimated")
  expect_equal(fit$sigeps_method, "estimated")
  expect_equal(fit$sigeps_fallback, "none")
  expect_equal(fit$sigeps_nu, fit$sigeps)
  expect_equal(fit$mubar, -fit$sigeps^2 / 2)
  expect_equal(fit$mu, log(fit$sigy^2) + fit$mubar, tolerance = 1e-12)
  ly2 <- log(y^2 + fit$del); rho1 <- as.numeric(fit$phi)
  expect_equal(fit$sigv, sqrt((var(ly2) - fit$sigeps^2) * (1 - rho1 * rho1)), tolerance = 1e-12)
  expect_equal(as.numeric(fit$theta), c(as.numeric(fit$phi), fit$sigy, fit$sigv))
})

test_that("estimated mode recovers an RV-type DGP where the theoretical mode cannot", {
  y <- .sim_rv(5000, sigv = 0.2, s_eps = 0.45)
  fe <- svp(y, p = 1, sigeps = "estimated")
  ft <- svp(y, p = 1)
  expect_lt(abs(fe$sigeps - 0.45), 0.08)
  expect_lt(abs(fe$sigv - 0.2), 0.06)
  expect_gt(abs(ft$sigv - 0.2), 0.2)      # abs() artifact: var(log y^2) << pi^2/2
  expect_identical(ft$sigeps_nu, pi / sqrt(2))
})

test_that("theoretical default is unchanged by the new argument", {
  set.seed(2); y <- sim_svp(1500, phi = 0.95, sigy = 1, sigv = 0.5)$y
  f0 <- svp(y); f1 <- svp(y, sigeps = "theoretical")
  for (f in c("mu", "phi", "sigy", "sigv", "theta", "sigeps", "sigeps_nu", "mubar"))
    expect_identical(f0[[f]], f1[[f]], info = f)
  expect_equal(f0$sigeps_method, "theoretical")
  expect_true(is.na(f0$sigeps_fallback))
})

test_that("estimated mode is refused where it is not supported", {
  y <- .sim_rv(1500)
  expect_error(svp(y, errorType = "Student-t", sigeps = "estimated"), "Gaussian")
  expect_error(svp(y, errorType = "GED", sigeps = "estimated"), "Gaussian")
  expect_error(svp(y, leverage = TRUE, sigeps = "estimated"), "leverage")
  expect_error(svp(y, sigvMethod = "hybrid", sigeps = "estimated"), "hybrid")
  fit <- svp(y, sigeps = "estimated")
  expect_error(filter_svp(fit, method = "mixture"), "corrected")
  expect_error(filter_svp(fit, method = "particle"), "corrected")
  expect_error(svpSE(fit, n_sim = 3), "estimated")
  expect_error(svp_IC(fit, filter_method = "mixture"), "corrected")
  expect_true(is.finite(svp_IC(fit, criteria = "BIC_Kalman", filter_method = "corrected")))
})

test_that("the corrected filter uses the estimated sigma_eps^2 and mubar", {
  y <- .sim_rv(2000)
  fit <- svp(y, p = 1, sigeps = "estimated")
  prm <- wARMASVp:::.get_filter_params(fit)
  expect_equal(prm$sig_eps2, fit$sigeps^2)
  expect_equal(prm$mu_intercept, fit$mu, tolerance = 1e-12)
  ft <- svp(y, p = 1)
  expect_false(isTRUE(all.equal(filter_svp(fit)$loglik, filter_svp(ft)$loglik)))
  fc <- forecast_svp(fit, H = 3)
  expect_equal(fc$var_forecast, fit$sigy^2 * exp(fc$log_var_forecast + fc$P_forecast / 2))
})

test_that("fallback: pooled estimate is used when the single-lag one fails, else error", {
  y <- .sim_rv(1500)
  calls <- 0L
  testthat::local_mocked_bindings(
    .sigeps2_moment = function(ys, var_ly2, phi, N, J, winsorize_eps = 0) {
      calls <<- calls + 1L
      if (isTRUE(winsorize_eps)) 0.2 else -1
    }, .package = "wARMASVp")
  fit <- suppressWarnings(svp(y, sigeps = "estimated"))
  expect_equal(fit$sigeps_fallback, "pooled")
  expect_equal(fit$sigeps_nu, sqrt(0.2))
  testthat::local_mocked_bindings(
    .sigeps2_moment = function(...) -1, .package = "wARMASVp")
  expect_error(svp(y, sigeps = "estimated"), "does not support the estimated mode")
})
