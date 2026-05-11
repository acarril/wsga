data(rddsga_synth)

test_that("wsga runs without error (no IPW, no bootstrap)", {
  fit <- wsga(y ~ 1 | sgroup, data = rddsga_synth,
                running = ~ x, bwidth = 0.5,
                noipsw = TRUE, bootstrap = FALSE)
  expect_s3_class(fit, "wsga")
  expect_length(fit$coefficients, 3)
})

test_that("per-subgroup estimates are close to true values (2 and 4)", {
  fit <- wsga(y ~ 1 | sgroup, data = rddsga_synth,
                running = ~ x, bwidth = 0.5,
                noipsw = TRUE, bootstrap = FALSE)
  expect_equal(fit$coefficients$g0, 2, tolerance = 0.2)
  expect_equal(fit$coefficients$g1, 4, tolerance = 0.2)
  expect_equal(fit$coefficients$diff, 2, tolerance = 0.3)
})

test_that("wsga with IPW runs and balance tables are returned", {
  fit <- wsga(y ~ m | sgroup, data = rddsga_synth,
                running = ~ x, bwidth = 0.5,
                noipsw = FALSE, bootstrap = FALSE)
  expect_false(is.null(fit$balance))
  expect_false(is.null(fit$balance$weighted$aggregate))
  expect_named(fit$balance$unweighted$aggregate$table,
               c("mean_G0", "mean_G1", "std_diff", "p_value"))
})

test_that("weighted balance reduces std_diff relative to unweighted", {
  fit <- wsga(y ~ m | sgroup, data = rddsga_synth,
                running = ~ x, bwidth = 0.5,
                noipsw = FALSE, bootstrap = FALSE)
  expect_lt(fit$balance$weighted$aggregate$avgdiff,
            fit$balance$unweighted$aggregate$avgdiff)
})

test_that("coef(), vcov(), confint(), nobs() S3 methods work", {
  fit <- wsga(y ~ 1 | sgroup, data = rddsga_synth,
                running = ~ x, bwidth = 0.5,
                noipsw = TRUE, bootstrap = FALSE)
  expect_length(coef(fit), 3)
  expect_equal(dim(vcov(fit)), c(2, 2))
  ci <- confint(fit)
  expect_equal(dim(ci), c(3, 2))
  expect_type(nobs(fit), "integer")
})

test_that("polynomial order p=2 runs", {
  fit <- wsga(y ~ 1 | sgroup, data = rddsga_synth,
                running = ~ x, bwidth = 0.5, p = 2,
                noipsw = TRUE, bootstrap = FALSE)
  expect_s3_class(fit, "wsga")
})

test_that("triangular kernel runs", {
  fit <- wsga(y ~ 1 | sgroup, data = rddsga_synth,
                running = ~ x, bwidth = 0.5, kernel = "triangular",
                noipsw = TRUE, bootstrap = FALSE)
  expect_s3_class(fit, "wsga")
})

test_that("fuzzy IV model runs", {
  fit <- wsga(y ~ 1 | sgroup, data = rddsga_synth,
                running = ~ x, bwidth = 0.5,
                model = "iv", fuzzy = ~ D,
                noipsw = TRUE, bootstrap = FALSE)
  expect_s3_class(fit, "wsga")
  # IV estimates should be larger in magnitude than RF (compliance < 1)
  rf <- wsga(y ~ 1 | sgroup, data = rddsga_synth,
             running = ~ x, bwidth = 0.5, noipsw = TRUE, bootstrap = FALSE)
  expect_gt(abs(fit$coefficients$g0), abs(rf$coefficients$g0) * 0.9)
})
