test_that("sigeps fields exist and are consistent for all three families", {
  set.seed(1)
  y <- sim_svp(1500, phi = c(0.5, 0.4), sigy = 1, sigv = 0.6)$y
  for (et in c("Gaussian", "Student-t", "GED")) {
    fit <- suppressWarnings(svp(y, p = 2, errorType = et))
    expect_true(all(c("sigeps", "sigeps_nu", "mubar", "nu_boundary") %in% names(fit)))
    expect_true(is.finite(fit$sigeps) && fit$sigeps > 0)
    expect_true(is.finite(fit$sigeps_nu) && fit$sigeps_nu > 0)
    # mu = log(sigy^2) + mubar, exact for every family and order
    expect_equal(fit$mu, log(fit$sigy^2) + fit$mubar, tolerance = 1e-12)
  }
})

test_that("sigeps_nu and mubar equal the closed forms at nu-hat", {
  set.seed(2)
  y <- sim_svp(1500, phi = 0.9, sigy = 1, sigv = 0.5, errorType = "Student-t", nu = 5)$y
  ft <- svp(y, errorType = "Student-t")
  expect_equal(ft$sigeps_nu^2, psigamma(0.5, 1) + psigamma(ft$v / 2, 1))
  expect_equal(ft$mubar, digamma(0.5) - digamma(ft$v / 2) + log(ft$v))
  fg <- suppressWarnings(svp(y, errorType = "GED"))
  expect_equal(fg$sigeps_nu^2, (2 / fg$v)^2 * psigamma(1 / fg$v, 1))
  expect_equal(fg$mubar, (2 / fg$v) * digamma(1 / fg$v) + lgamma(1 / fg$v) - lgamma(3 / fg$v))
  fgau <- svp(y)
  expect_equal(fgau$sigeps_nu, pi / sqrt(2))
  expect_equal(fgau$mubar, digamma(0.5) + log(2))
  expect_true(is.na(fgau$nu_boundary))
})

test_that("Gaussian sigeps equals the single-lag moment formula", {
  set.seed(3)
  y <- sim_svp(1200, phi = 0.95, sigy = 1, sigv = 0.3)$y
  fit <- svp(y, p = 1, J = 10, del = 1e-10)
  ly2 <- log(y^2 + 1e-10); ys <- ly2 - mean(ly2); N <- length(y)
  gam1 <- sum(ys[2:N] * ys[1:(N - 1)]) / (N - 1)
  se2 <- var(ly2) - gam1 / as.numeric(fit$phi)
  expect_equal(fit$sigeps, sqrt(se2), tolerance = 1e-12)
})

test_that("interior nu-hat gives sigeps == sigeps_nu; the fields survive leverage", {
  set.seed(4)
  y <- sim_svp(2000, phi = 0.9, sigy = 1, sigv = 0.5, errorType = "Student-t", nu = 5,
               leverage = TRUE, rho = -0.5)$y
  fit <- svp(y, errorType = "Student-t", leverage = TRUE)
  if (!fit$nu_boundary) expect_equal(fit$sigeps, fit$sigeps_nu, tolerance = 1e-8)
  expect_true(all(c("sigeps", "sigeps_nu", "mubar") %in% names(fit)))
})

test_that("non-positive sigma_eps^2 estimate gives NA sigeps; warning only on request", {
  expect_silent(v <- wARMASVp:::.sigeps_raw(-0.3))
  expect_true(is.na(v))
  expect_warning(wARMASVp:::.sigeps_raw(-0.3, warn = TRUE), "not positive")
  expect_equal(wARMASVp:::.sigeps_raw(0.25), 0.5)
})

test_that("sigy uses log(y^2 + del) for every family at p = 1 as well", {
  set.seed(5)
  y <- sim_svp(1200, phi = 0.9, sigy = 1, sigv = 0.5)$y
  y[100] <- 0   # an exact zero return must not give sigy = 0 or -Inf
  for (et in c("Gaussian", "Student-t", "GED")) {
    fit <- suppressWarnings(svp(y, p = 1, errorType = et, del = 1e-6))
    expect_true(is.finite(fit$sigy) && fit$sigy > 0)
    expect_equal(fit$mu, log(fit$sigy^2) + fit$mubar, tolerance = 1e-12)
  }
})
