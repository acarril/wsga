data(rddsga_synth)

test_that("fixed_ps default is FALSE; output records it", {
  fit_default <- wsga(y ~ m | sgroup, data = rddsga_synth,
                      running = ~ x, bwidth = 0.5,
                      bootstrap = TRUE, bsreps = 20, seed = 1)
  expect_false(fit_default$fixed_ps)
})

test_that("fixed_ps = TRUE flips the recorded flag and runs cleanly", {
  fit_fixed <- wsga(y ~ m | sgroup, data = rddsga_synth,
                    running = ~ x, bwidth = 0.5,
                    bootstrap = TRUE, bsreps = 20, seed = 1,
                    fixed_ps = TRUE)
  expect_true(fit_fixed$fixed_ps)
  expect_equal(nrow(fit_fixed$bootstrap$draws), 20L)
})

test_that("fixed_ps changes the bootstrap distribution vs. refit-per-rep", {
  fit_refit <- wsga(y ~ m | sgroup, data = rddsga_synth,
                    running = ~ x, bwidth = 0.5,
                    bootstrap = TRUE, bsreps = 50, seed = 7,
                    fixed_ps = FALSE)
  fit_fixed <- wsga(y ~ m | sgroup, data = rddsga_synth,
                    running = ~ x, bwidth = 0.5,
                    bootstrap = TRUE, bsreps = 50, seed = 7,
                    fixed_ps = TRUE)
  # Point estimates are identical (PS is computed once on the original sample
  # in both cases, and the bootstrap inherits that fit when fixed_ps = FALSE
  # only inside the loop)
  expect_equal(fit_refit$coefficients, fit_fixed$coefficients)
  # Draws differ — fixed_ps doesn't refit the logit on each replicate
  expect_false(isTRUE(all.equal(fit_refit$bootstrap$draws,
                                fit_fixed$bootstrap$draws)))
})

test_that("fixed_ps is reproducible under seed", {
  fit_a <- wsga(y ~ m | sgroup, data = rddsga_synth,
                running = ~ x, bwidth = 0.5,
                bootstrap = TRUE, bsreps = 20, seed = 42,
                fixed_ps = TRUE)
  fit_b <- wsga(y ~ m | sgroup, data = rddsga_synth,
                running = ~ x, bwidth = 0.5,
                bootstrap = TRUE, bsreps = 20, seed = 42,
                fixed_ps = TRUE)
  expect_equal(fit_a$bootstrap$draws, fit_b$bootstrap$draws)
})

test_that("fixed_ps works alongside cluster bootstrap", {
  set.seed(1)
  d <- rddsga_synth
  d$school <- sample.int(50, nrow(d), replace = TRUE)
  fit <- wsga(y ~ m | sgroup, data = d,
              running = ~ x, bwidth = 0.5,
              bootstrap = TRUE, bsreps = 20, seed = 5,
              cluster_var = "school", fixed_ps = TRUE)
  expect_true(fit$fixed_ps)
  expect_equal(fit$bootstrap$N_clusters, 50L)
})

test_that("print method tags [PS fixed] when fixed_ps = TRUE", {
  fit <- wsga(y ~ m | sgroup, data = rddsga_synth,
              running = ~ x, bwidth = 0.5,
              bootstrap = TRUE, bsreps = 10, seed = 1,
              fixed_ps = TRUE)
  out <- capture.output(print(fit))
  expect_true(any(grepl("\\[PS fixed\\]", out)))
})
