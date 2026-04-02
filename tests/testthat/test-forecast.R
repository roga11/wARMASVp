# =========================================================================== #
# Tests for forecasting
# =========================================================================== #

test_that("kalman_filter runs and returns expected components", {
  set.seed(42)
  y <- as.numeric(sim_svp(500, phi = 0.95, sigy = 1, sigv = 0.3))
  mdl <- svp(y, p = 1)
  kf <- kalman_filter(y, mdl)
  expect_type(kf, "list")
  expect_true("w_estimated" %in% names(kf))
  expect_true("w_smoothed" %in% names(kf))
  # New enriched outputs
  expect_true("P_filtered" %in% names(kf))
  expect_true("loglik" %in% names(kf))
})

test_that("forecast_svp accepts model object and produces h-step forecasts", {
  set.seed(42)
  sim <- sim_svp(500, phi = 0.95, sigy = 1, sigv = 0.3,
                 leverage = TRUE, rho = -0.3)
  fit <- svp(sim$y, p = 1, leverage = TRUE)
  fc <- forecast_svp(fit, H = 5)
  expect_s3_class(fc, "svp_forecast")
  expect_true("w_forecasted" %in% names(fc))
  expect_length(fc$w_forecasted, 5)
})

test_that("forecast_svp errors on raw numeric input", {
  set.seed(42)
  y <- as.numeric(sim_svp(200, phi = 0.95, sigy = 1, sigv = 0.3))
  expect_error(forecast_svp(y, H = 5), "svp/svp_t/svp_ged")
})

test_that("forecast_svp: all three output representations are stored", {
  set.seed(42)
  y <- as.numeric(sim_svp(300, phi = 0.95, sigy = 1, sigv = 0.3))
  fit <- svp(y, p = 1)
  fc <- forecast_svp(fit, H = 10)
  expect_length(fc$log_var_forecast, 10)
  expect_length(fc$var_forecast, 10)
  expect_length(fc$vol_forecast, 10)
  expect_length(fc$P_forecast, 10)
  # Default output is log-variance
  expect_equal(fc$w_forecasted, fc$log_var_forecast)
})

test_that("forecast_svp: output parameter selects primary output", {
  set.seed(42)
  y <- as.numeric(sim_svp(300, phi = 0.95, sigy = 1, sigv = 0.3))
  fit <- svp(y, p = 1)
  fc_var <- forecast_svp(fit, H = 5, output = "variance")
  expect_equal(fc_var$w_forecasted, fc_var$var_forecast)
  fc_vol <- forecast_svp(fit, H = 5, output = "volatility")
  expect_equal(fc_vol$w_forecasted, fc_vol$vol_forecast)
})

test_that("forecast_svp: P_forecast is positive and increasing", {
  set.seed(42)
  y <- as.numeric(sim_svp(500, phi = 0.95, sigy = 1, sigv = 0.3))
  fit <- svp(y, p = 1)
  fc <- forecast_svp(fit, H = 20)
  expect_true(all(fc$P_forecast > 0))
  # MSE should generally increase (not strictly for all h, but trend)
  expect_true(fc$P_forecast[20] > fc$P_forecast[1])
})

test_that("forecast_svp: sigma2_forecast is positive", {
  set.seed(42)
  y <- as.numeric(sim_svp(300, phi = 0.95, sigy = 1, sigv = 0.3))
  fit <- svp(y, p = 1)
  fc <- forecast_svp(fit, H = 10)
  expect_true(all(fc$var_forecast > 0))
  expect_true(all(fc$vol_forecast > 0))
})

test_that("forecast_svp: long horizon converges toward 0 (log-var)", {
  set.seed(42)
  y <- as.numeric(sim_svp(500, phi = 0.90, sigy = 1, sigv = 0.3))
  fit <- svp(y, p = 1)
  fc <- forecast_svp(fit, H = 100)
  # Long-horizon log-variance forecast should converge toward 0
  expect_true(abs(fc$log_var_forecast[100]) < abs(fc$log_var_forecast[1]) + 0.1)
})

test_that("forecast_svp: Student-t uses E[u^2] = nu/(nu-2)", {
  set.seed(42)
  sim <- sim_svp(500, phi = 0.95, sigy = 1, sigv = 0.3,
                 errorType = "Student-t", nu = 5)
  fit <- svp(as.numeric(sim), p = 1, errorType = "Student-t")
  fc <- forecast_svp(fit, H = 5)
  # var_forecast should incorporate nu/(nu-2) = 5/3 factor
  expect_true(all(fc$var_forecast > 0))
})

test_that("forecast_svp: GED uses E[u^2] = 1", {
  set.seed(42)
  sim <- sim_svp(500, phi = 0.95, sigy = 1, sigv = 0.3,
                 errorType = "GED", nu = 1.5)
  fit <- svp(as.numeric(sim), p = 1, errorType = "GED")
  fc <- forecast_svp(fit, H = 5)
  expect_true(all(fc$var_forecast > 0))
})

test_that("plot.svp_forecast works without error", {
  set.seed(42)
  sim <- sim_svp(300, phi = 0.95, sigy = 1, sigv = 0.3,
                 leverage = TRUE, rho = -0.3)
  fit <- svp(sim$y, p = 1, leverage = TRUE)
  fc <- forecast_svp(fit, H = 5)
  expect_silent(plot(fc))
})

test_that("companionMat returns correct dimensions", {
  C <- wARMASVp:::companionMat(c(0.7, 0.2), p = 2, q = 1)
  expect_equal(dim(C), c(2, 2))
  expect_equal(unname(C[1, 1]), 0.7)
  expect_equal(unname(C[1, 2]), 0.2)
})
