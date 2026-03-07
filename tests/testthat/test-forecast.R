test_that("kalman_filter runs and returns expected components", {
  set.seed(42)
  y <- sim_svp(500, phi = 0.95, sigy = 1, sigv = 0.3)
  mdl <- svp(y, p = 1)
  kf <- kalman_filter(y, mdl)
  expect_type(kf, "list")
  expect_true("w_estimated" %in% names(kf))
  expect_true("w_smoothed" %in% names(kf))
})

test_that("forecast_svp produces h-step forecasts", {
  set.seed(42)
  y <- sim_svp(500, phi = 0.95, sigy = 1, sigv = 0.3,
               leverage = TRUE, rho = -0.3)$y
  fc <- forecast_svp(y, p = 1, H = 5)
  expect_s3_class(fc, "svp_forecast")
  expect_true("w_forecasted" %in% names(fc))
  expect_length(fc$w_forecasted, 5)
})

test_that("plot.svp_forecast works without error", {
  set.seed(42)
  y <- sim_svp(500, phi = 0.95, sigy = 1, sigv = 0.3,
               leverage = TRUE, rho = -0.3)$y
  fc <- forecast_svp(y, p = 1, H = 5)
  expect_silent(plot(fc))
})

test_that("companionMat returns correct dimensions", {
  C <- companionMat(c(0.7, 0.2), p = 2, q = 1)
  expect_equal(dim(C), c(2, 2))
  expect_equal(unname(C[1, 1]), 0.7)
  expect_equal(unname(C[1, 2]), 0.2)
})
