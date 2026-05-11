data(rddsga_synth)

test_that("default inference is 'empirical' when bootstrap = TRUE", {
  fit <- wsga_rdd(y ~ 1 | sgroup, data = rddsga_synth,
              running = ~ x, bwidth = 0.5,
              noipsw = TRUE, bootstrap = TRUE, bsreps = 20, seed = 1)
  expect_equal(fit$inference, "empirical")
  # CIs should match the bootstrap percentile CIs exactly
  expect_equal(fit$ci$g0[["lb"]],
               unname(quantile(fit$bootstrap$draws[, "g0"], 0.025)))
  expect_equal(fit$ci$g0[["ub"]],
               unname(quantile(fit$bootstrap$draws[, "g0"], 0.975)))
})

test_that("default inference is 'analytical' when bootstrap = FALSE", {
  fit <- wsga_rdd(y ~ 1 | sgroup, data = rddsga_synth,
              running = ~ x, bwidth = 0.5,
              noipsw = TRUE, bootstrap = FALSE)
  expect_equal(fit$inference, "analytical")
  expect_null(fit$bootstrap)
})

test_that("inference = 'normal' yields symmetric CIs from the bootstrap SE", {
  fit <- wsga_rdd(y ~ 1 | sgroup, data = rddsga_synth,
              running = ~ x, bwidth = 0.5,
              noipsw = TRUE, bootstrap = TRUE, bsreps = 50, seed = 7,
              inference = "normal")
  expect_equal(fit$inference, "normal")
  # CI should be exactly b ± qnorm(0.975) * SE
  z <- qnorm(0.975)
  expect_equal(fit$ci$g0[["lb"]], fit$coefficients$g0 - z * fit$se$g0)
  expect_equal(fit$ci$g0[["ub"]], fit$coefficients$g0 + z * fit$se$g0)
  expect_equal(fit$ci$diff[["lb"]], fit$coefficients$diff - z * fit$se$diff)
})

test_that("inference = 'normal' and 'empirical' differ on the same bootstrap", {
  fit_emp  <- wsga_rdd(y ~ 1 | sgroup, data = rddsga_synth,
                   running = ~ x, bwidth = 0.5,
                   noipsw = TRUE, bootstrap = TRUE, bsreps = 100, seed = 42,
                   inference = "empirical")
  fit_norm <- wsga_rdd(y ~ 1 | sgroup, data = rddsga_synth,
                   running = ~ x, bwidth = 0.5,
                   noipsw = TRUE, bootstrap = TRUE, bsreps = 100, seed = 42,
                   inference = "normal")
  # Same point estimates and SEs; different CIs
  expect_equal(fit_emp$coefficients, fit_norm$coefficients)
  expect_equal(fit_emp$se, fit_norm$se)
  expect_false(isTRUE(all.equal(fit_emp$ci$g0, fit_norm$ci$g0)))
})

test_that("invalid combinations error", {
  expect_error(
    wsga_rdd(y ~ 1 | sgroup, data = rddsga_synth,
             running = ~ x, bwidth = 0.5,
             bootstrap = TRUE, inference = "analytical"),
    "requires `bootstrap = FALSE`"
  )
  expect_error(
    wsga_rdd(y ~ 1 | sgroup, data = rddsga_synth,
             running = ~ x, bwidth = 0.5,
             bootstrap = FALSE, inference = "empirical"),
    "requires `bootstrap = TRUE`"
  )
  expect_error(
    wsga_rdd(y ~ 1 | sgroup, data = rddsga_synth,
             running = ~ x, bwidth = 0.5,
             bootstrap = FALSE, inference = "normal"),
    "requires `bootstrap = TRUE`"
  )
})

test_that("bootstrap result exposes failed count", {
  fit <- wsga_rdd(y ~ 1 | sgroup, data = rddsga_synth,
              running = ~ x, bwidth = 0.5,
              noipsw = TRUE, bootstrap = TRUE, bsreps = 30, seed = 3)
  expect_true("failed" %in% names(fit$bootstrap))
  expect_equal(fit$bootstrap$failed + fit$bootstrap$B_ok, 30)
})
