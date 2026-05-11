data(rddsga_synth)

test_that("bootstrap runs and returns B draws", {
  fit <- wsga_rdd(y ~ 1 | sgroup, data = rddsga_synth,
                running = ~ x, bwidth = 0.5,
                noipsw = TRUE, bootstrap = TRUE, bsreps = 20, seed = 1)
  expect_false(is.null(fit$bootstrap))
  expect_equal(nrow(fit$bootstrap$draws), 20)
  expect_equal(ncol(fit$bootstrap$draws), 3)
  expect_equal(colnames(fit$bootstrap$draws), c("g0", "g1", "diff"))
})

test_that("bootstrap SEs are positive and plausible", {
  fit <- wsga_rdd(y ~ 1 | sgroup, data = rddsga_synth,
                running = ~ x, bwidth = 0.5,
                noipsw = TRUE, bootstrap = TRUE, bsreps = 50, seed = 2)
  expect_gt(fit$se$g0, 0)
  expect_gt(fit$se$g1, 0)
  expect_gt(fit$se$diff, 0)
})

test_that("empirical p-values are in (0, 1]", {
  fit <- wsga_rdd(y ~ 1 | sgroup, data = rddsga_synth,
                running = ~ x, bwidth = 0.5,
                noipsw = TRUE, bootstrap = TRUE, bsreps = 30, seed = 3)
  expect_true(fit$pval$g0   > 0 && fit$pval$g0   <= 1)
  expect_true(fit$pval$g1   > 0 && fit$pval$g1   <= 1)
  expect_true(fit$pval$diff > 0 && fit$pval$diff  <= 1)
})

test_that("diff CI uses empirical percentiles, not normal approximation", {
  fit <- wsga_rdd(y ~ 1 | sgroup, data = rddsga_synth,
                running = ~ x, bwidth = 0.5,
                noipsw = TRUE, bootstrap = TRUE, bsreps = 50, seed = 4)
  diff_draws <- fit$bootstrap$draws[, "diff"]
  expect_equal(fit$ci$diff[["lb"]], unname(quantile(diff_draws, 0.025)))
  expect_equal(fit$ci$diff[["ub"]], unname(quantile(diff_draws, 0.975)))
})

test_that("empirical p-values use recentered formula (H0: coef=0)", {
  fit <- wsga_rdd(y ~ 1 | sgroup, data = rddsga_synth,
              running = ~ x, bwidth = 0.5,
              noipsw = TRUE, bootstrap = TRUE, bsreps = 50, seed = 5)
  draws  <- fit$bootstrap$draws
  est    <- c(fit$coefficients$g0, fit$coefficients$g1, fit$coefficients$diff)
  B_ok   <- fit$bootstrap$B_ok
  expected <- sapply(seq_len(3), function(g) {
    cnt <- sum(abs(draws[, g] - est[g]) >= abs(est[g]))
    (1 + cnt) / (B_ok + 1)
  })
  expect_equal(fit$pval$g0,   expected[1])
  expect_equal(fit$pval$g1,   expected[2])
  expect_equal(fit$pval$diff, expected[3])
})

test_that("empirical p-values are small for clearly significant estimates", {
  # g0≈2, g1≈4 with bwidth=5 — large effects, p should be much less than 0.05
  fit <- wsga_rdd(y ~ 1 | sgroup, data = rddsga_synth,
              running = ~ x, bwidth = 5,
              noipsw = TRUE, bootstrap = TRUE, bsreps = 100, seed = 6)
  expect_lt(fit$pval$g0,   0.05)
  expect_lt(fit$pval$g1,   0.05)
  expect_lt(fit$pval$diff, 0.05)
})

test_that("seed makes bootstrap reproducible", {
  fit_a <- wsga_rdd(y ~ 1 | sgroup, data = rddsga_synth,
                running = ~ x, bwidth = 0.5,
                noipsw = TRUE, bootstrap = TRUE, bsreps = 20, seed = 99)
  fit_b <- wsga_rdd(y ~ 1 | sgroup, data = rddsga_synth,
                running = ~ x, bwidth = 0.5,
                noipsw = TRUE, bootstrap = TRUE, bsreps = 20, seed = 99)
  expect_equal(fit_a$bootstrap$draws, fit_b$bootstrap$draws)
})
