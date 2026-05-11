data(rddsga_synth)

test_that("rddsga() is a deprecated alias for wsga_rdd()", {
  expect_warning(
    fit_old <- rddsga(y ~ 1 | sgroup, data = rddsga_synth,
                      running = ~ x, bwidth = 0.5,
                      noipsw = TRUE, bootstrap = FALSE),
    class = "deprecatedWarning"
  )
  fit_new <- wsga_rdd(y ~ 1 | sgroup, data = rddsga_synth,
                      running = ~ x, bwidth = 0.5,
                      noipsw = TRUE, bootstrap = FALSE)
  # The alias must produce identical estimates and class
  expect_s3_class(fit_old, "wsga")
  expect_equal(fit_old$coefficients, fit_new$coefficients)
})
