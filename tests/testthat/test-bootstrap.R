data(rddsga_synth)

test_that("bootstrap runs and returns B draws", {
  fit <- wsga(y ~ 1 | sgroup, data = rddsga_synth,
                running = ~ x, bwidth = 0.5,
                noipsw = TRUE, bootstrap = TRUE, bsreps = 20, seed = 1)
  expect_false(is.null(fit$bootstrap))
  expect_equal(nrow(fit$bootstrap$draws), 20)
  expect_equal(ncol(fit$bootstrap$draws), 3)
  expect_equal(colnames(fit$bootstrap$draws), c("g0", "g1", "diff"))
})

test_that("bootstrap SEs are positive and plausible", {
  fit <- wsga(y ~ 1 | sgroup, data = rddsga_synth,
                running = ~ x, bwidth = 0.5,
                noipsw = TRUE, bootstrap = TRUE, bsreps = 50, seed = 2)
  expect_gt(fit$se$g0, 0)
  expect_gt(fit$se$g1, 0)
  expect_gt(fit$se$diff, 0)
})

test_that("empirical p-values are in (0, 1]", {
  fit <- wsga(y ~ 1 | sgroup, data = rddsga_synth,
                running = ~ x, bwidth = 0.5,
                noipsw = TRUE, bootstrap = TRUE, bsreps = 30, seed = 3)
  expect_true(fit$pval$g0   > 0 && fit$pval$g0   <= 1)
  expect_true(fit$pval$g1   > 0 && fit$pval$g1   <= 1)
  expect_true(fit$pval$diff > 0 && fit$pval$diff  <= 1)
})

test_that("diff CI uses empirical percentiles, not normal approximation", {
  fit <- wsga(y ~ 1 | sgroup, data = rddsga_synth,
                running = ~ x, bwidth = 0.5,
                noipsw = TRUE, bootstrap = TRUE, bsreps = 50, seed = 4)
  diff_draws <- fit$bootstrap$draws[, "diff"]
  expect_equal(fit$ci$diff[["lb"]], unname(quantile(diff_draws, 0.025)))
  expect_equal(fit$ci$diff[["ub"]], unname(quantile(diff_draws, 0.975)))
})

test_that("seed makes bootstrap reproducible", {
  fit_a <- wsga(y ~ 1 | sgroup, data = rddsga_synth,
                running = ~ x, bwidth = 0.5,
                noipsw = TRUE, bootstrap = TRUE, bsreps = 20, seed = 99)
  fit_b <- wsga(y ~ 1 | sgroup, data = rddsga_synth,
                running = ~ x, bwidth = 0.5,
                noipsw = TRUE, bootstrap = TRUE, bsreps = 20, seed = 99)
  expect_equal(fit_a$bootstrap$draws, fit_b$bootstrap$draws)
})
