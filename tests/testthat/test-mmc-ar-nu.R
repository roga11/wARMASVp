# mmc_ar with Student-t / GED errors searches over the shape parameter nu.

sim_ht <- function(et, nu, n = 400, seed = 42) {
  set.seed(seed)
  sim_svp(n, phi = 0.95, sigy = 1, sigv = 0.5, errorType = et, nu = nu)$y
}
run_ar <- function(y, et, ..., seed = 1) {
  set.seed(seed)
  suppressWarnings(mmc_ar(y, p_null = 1, p_alt = 2, N = 9, errorType = et,
                          maxit = 4, silence = TRUE, ...))
}

test_that("mmc_ar carries nu as a nuisance parameter for Student-t and GED", {
  skip_on_cran()
  for (et in c("Student-t", "GED")) {
    y <- sim_ht(et, if (et == "Student-t") 5 else 1.5)
    nu_hat <- svp(y, p = 1, errorType = et, logNu = TRUE)$v
    rng <- if (et == "Student-t") c(2.01, 500) else c(0.1, 20)

    # length p_null+3: user-set half-width for nu
    r <- run_ar(y, et, eps = c(0.1, 0, 0.1, 0.05))
    expect_length(r$par, 4)
    expect_true(r$par[4] >= max(rng[1], nu_hat - 0.05) - 1e-12)
    expect_true(r$par[4] <= min(rng[2], nu_hat + 0.05) + 1e-12)
    expect_true(r$pval >= 0 && r$pval <= 1)

    # length p_null+2: nu gets the fixed default half-width 0.3
    r2 <- run_ar(y, et, eps = c(0.1, 0, 0.1))
    expect_length(r2$par, 4)
    expect_true(abs(r2$par[4] - nu_hat) <= 0.3 + 1e-12)

    # default eps
    r3 <- run_ar(y, et)
    expect_length(r3$par, 4)
    expect_true(abs(r3$par[4] - nu_hat) <= 0.3 + 1e-12)

    # zero half-width holds nu at its null estimate
    r0 <- run_ar(y, et, eps = c(0.1, 0, 0.1, 0))
    expect_identical(r0$par[4], nu_hat)

    # both optimizers
    rg <- run_ar(y, et, eps = c(0.1, 0, 0.1, 0.05), method = "GenSA")
    expect_length(rg$par, 4)
    expect_true(rg$pval >= 0 && rg$pval <= 1)

    # wrong lengths
    expect_error(run_ar(y, et, eps = c(0.1, 0)), "p_null\\+2 or p_null\\+3")
    expect_error(run_ar(y, et, eps = rep(0.1, 5)), "p_null\\+2 or p_null\\+3")
  }
})

test_that("mmc_ar with Gaussian errors has no nu coordinate", {
  skip_on_cran()
  set.seed(42)
  y <- sim_svp(400, phi = 0.95, sigy = 1, sigv = 0.5)$y
  r <- run_ar(y, "Gaussian")
  expect_length(r$par, 3)
  expect_error(run_ar(y, "Gaussian", eps = c(0.1, 0, 0.1, 0.1)), "length p_null\\+2 \\(")
})

test_that(".simnull_ar reuses nu-free draws at every nu", {
  skip_on_cran()
  n <- 300 + 100; nd <- 30
  for (et in c("Student-t", "GED")) {
    set.seed(7)
    inn <- list(zeta = matrix(rnorm(n * nd), n), aux = matrix(rnorm(n * nd), n))
    if (et == "Student-t") inn$U_chi2 <- matrix(runif(n * nd), n)
    nus <- if (et == "Student-t") c(4, 8) else c(1.2, 1.8)
    f <- function(nu) suppressWarnings(as.numeric(wARMASVp:::.simnull_ar(
      c(0.9, 1, 0.5), 1, 2, 10, 300, 9, 100, errorType = et, nu_sim = nu,
      innovations = inn)))
    seed_before <- .Random.seed
    a <- f(nus[1])
    expect_identical(.Random.seed, seed_before)   # no fresh draws consumed
    expect_identical(f(nus[1]), a)                # deterministic given the draws
    expect_false(isTRUE(all.equal(f(nus[2]), a))) # and it does depend on nu
  }
})
