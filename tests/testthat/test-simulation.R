test_that("sim_svp Gaussian SV(1) returns correct length", {
  set.seed(1)
  y <- sim_svp(500, phi = 0.95, sigy = 1, sigv = 0.3)
  expect_length(y, 500)
  expect_type(y, "double")
})

test_that("sim_svp Gaussian SV(2) returns correct length", {
  set.seed(1)
  y <- sim_svp(500, phi = c(0.20, 0.63), sigy = 1, sigv = 0.5)
  expect_length(y, 500)
})

test_that("sim_svp Student-t returns correct length", {
  set.seed(1)
  y <- sim_svp(500, phi = 0.90, sigy = 1, sigv = 0.3,
               errorType = "Student-t", nu = 5)
  expect_length(y, 500)
})

test_that("sim_svp GED returns correct length", {
  set.seed(1)
  y <- sim_svp(500, phi = 0.90, sigy = 1, sigv = 0.3,
               errorType = "GED", nu = 1.5)
  expect_length(y, 500)
})

test_that("sim_svp leverage returns list with expected components", {
  set.seed(1)
  sim <- sim_svp(500, phi = 0.95, sigy = 1, sigv = 0.3,
                 leverage = TRUE, rho = -0.5)
  expect_type(sim, "list")
  expect_named(sim, c("y", "h", "zeta", "veta"))
  expect_length(sim$y, 500)
})

test_that("sim_svp K > 1 returns matrix", {
  set.seed(1)
  y <- sim_svp(500, phi = 0.90, sigy = 1, sigv = 0.3,
               errorType = "Student-t", nu = 5, K = 3)
  expect_true(is.matrix(y))
  expect_equal(nrow(y), 500)
  expect_equal(ncol(y), 3)
})

test_that("sim_svp rejects invalid inputs", {
  expect_error(sim_svp(500, phi = 0.9, sigy = 1, sigv = 0.3,
                       errorType = "Student-t", nu = 1.5),
               "nu must be > 2")
  expect_error(sim_svp(500, phi = 0.9, sigy = 1, sigv = 0.3,
                       errorType = "GED", nu = -1),
               "nu must be > 0")
  expect_error(sim_svp(500, phi = 0.9, sigy = 1, sigv = 0.3,
                       errorType = "Student-t"),
               "nu is required")
})
