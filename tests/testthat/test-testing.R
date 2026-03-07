test_that("lmc_ar runs and returns svp_test object", {
  set.seed(42)
  y <- sim_svp(1000, phi = 0.95, sigy = 1, sigv = 0.3)
  test <- lmc_ar(y, p_null = 1, p_alt = 2, N = 19)
  expect_s3_class(test, "svp_test")
  expect_true("pval" %in% names(test))
  expect_true(test$pval >= 0 && test$pval <= 1)
})

test_that("lmc_lev runs and returns svp_test object", {
  set.seed(42)
  y <- sim_svp(1000, phi = 0.95, sigy = 1, sigv = 0.3)
  test <- lmc_lev(y, p = 1, N = 19)
  expect_s3_class(test, "svp_test")
  expect_true(test$pval >= 0 && test$pval <= 1)
})

test_that("lmc_t runs and returns svp_test object", {
  set.seed(42)
  y <- sim_svp(1000, phi = 0.90, sigy = 1, sigv = 0.3,
               errorType = "Student-t", nu = 5)
  test <- suppressWarnings(lmc_t(y, nu_null = 5, N = 19))
  expect_s3_class(test, "svp_test")
  expect_true(test$pval >= 0 && test$pval <= 1)
})

test_that("lmc_ged runs and returns svp_test object", {
  set.seed(42)
  y <- sim_svp(1000, phi = 0.90, sigy = 1, sigv = 0.3,
               errorType = "GED", nu = 1.5)
  test <- suppressWarnings(lmc_ged(y, nu_null = 1.5, N = 19))
  expect_s3_class(test, "svp_test")
  expect_true(test$pval >= 0 && test$pval <= 1)
})

test_that("print.svp_test works", {
  set.seed(42)
  y <- sim_svp(1000, phi = 0.95, sigy = 1, sigv = 0.3)
  test <- lmc_ar(y, p_null = 1, p_alt = 2, N = 19)
  expect_output(print(test))
})
