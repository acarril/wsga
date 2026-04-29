data(rddsga_synth)

test_that("rddsga runs without error (no IPW, no bootstrap)", {
  fit <- rddsga(y ~ 1 | sgroup, data = rddsga_synth,
                running = ~ x, bwidth = 0.5,
                noipsw = TRUE, bootstrap = FALSE)
  expect_s3_class(fit, "rddsga")
  expect_length(fit$coefficients, 3)
})

test_that("per-subgroup estimates are close to true values (2 and 4)", {
  fit <- rddsga(y ~ 1 | sgroup, data = rddsga_synth,
                running = ~ x, bwidth = 0.5,
                noipsw = TRUE, bootstrap = FALSE)
  expect_equal(fit$coefficients$g0, 2, tolerance = 0.2)
  expect_equal(fit$coefficients$g1, 4, tolerance = 0.2)
  expect_equal(fit$coefficients$diff, 2, tolerance = 0.3)
})

test_that("rddsga with IPW runs and balance tables are returned", {
  fit <- rddsga(y ~ m | sgroup, data = rddsga_synth,
                running = ~ x, bwidth = 0.5,
                noipsw = FALSE, bootstrap = FALSE)
  expect_false(is.null(fit$balance))
  expect_false(is.null(fit$balance$weighted))
  expect_named(fit$balance$unweighted$table,
               c("mean_G0", "mean_G1", "std_diff", "p_value"))
})

test_that("weighted balance reduces std_diff relative to unweighted", {
  fit <- rddsga(y ~ m | sgroup, data = rddsga_synth,
                running = ~ x, bwidth = 0.5,
                noipsw = FALSE, bootstrap = FALSE)
  expect_lt(fit$balance$weighted$avgdiff,
            fit$balance$unweighted$avgdiff)
})

test_that("coef(), vcov(), confint(), nobs() S3 methods work", {
  fit <- rddsga(y ~ 1 | sgroup, data = rddsga_synth,
                running = ~ x, bwidth = 0.5,
                noipsw = TRUE, bootstrap = FALSE)
  expect_length(coef(fit), 3)
  expect_equal(dim(vcov(fit)), c(2, 2))
  ci <- confint(fit)
  expect_equal(dim(ci), c(3, 2))
  expect_type(nobs(fit), "integer")
})

test_that("polynomial order p=2 runs", {
  fit <- rddsga(y ~ 1 | sgroup, data = rddsga_synth,
                running = ~ x, bwidth = 0.5, p = 2,
                noipsw = TRUE, bootstrap = FALSE)
  expect_s3_class(fit, "rddsga")
})

test_that("triangular kernel runs", {
  fit <- rddsga(y ~ 1 | sgroup, data = rddsga_synth,
                running = ~ x, bwidth = 0.5, kernel = "triangular",
                noipsw = TRUE, bootstrap = FALSE)
  expect_s3_class(fit, "rddsga")
})

test_that("fuzzy IV model runs", {
  fit <- rddsga(y ~ 1 | sgroup, data = rddsga_synth,
                running = ~ x, bwidth = 0.5,
                model = "iv", fuzzy = ~ D,
                noipsw = TRUE, bootstrap = FALSE)
  expect_s3_class(fit, "rddsga")
  # IV estimates should be larger in magnitude than RF (compliance < 1)
  rf <- rddsga(y ~ 1 | sgroup, data = rddsga_synth,
               running = ~ x, bwidth = 0.5, noipsw = TRUE, bootstrap = FALSE)
  expect_gt(abs(fit$coefficients$g0), abs(rf$coefficients$g0) * 0.9)
})
