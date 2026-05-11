# Synthetic DiD panel: 2 periods, n_units balanced; sgroup × treat × M
# constant within unit.  True per-subgroup effects: tau_0 = 1, tau_1 = 3.
make_did_panel <- function(n_units = 200L, seed = 1L) {
  set.seed(seed)
  unit   <- rep(seq_len(n_units), each = 2L)
  time   <- rep(c(0L, 1L), times = n_units)
  sgroup <- rep(rbinom(n_units, 1L, 0.5), each = 2L)
  M      <- rep(rnorm(n_units, mean = 0.4 * sgroup), each = 2L)
  D      <- rep(rbinom(n_units, 1L, 0.5), each = 2L)
  alpha  <- rep(rnorm(n_units), each = 2L)
  post   <- as.integer(time == 1L)
  tau    <- ifelse(sgroup == 1L, 3, 1)
  y      <- alpha + 0.5 * post + tau * D * post + 0.3 * M + rnorm(length(unit), sd = 0.5)
  data.frame(unit = unit, time = time, sgroup = sgroup,
             m = M, D = D, y = y)
}

test_that("DiD design errors when required args are missing", {
  d <- make_did_panel()
  expect_error(
    wsga_did(y ~ 1 | sgroup, data = d),
    "requires `unit`, `time`, and `treat`"
  )
})

test_that("DiD requires exactly 2 unique time values", {
  d <- make_did_panel()
  d$time[1] <- 99L  # introduces a third time value
  expect_error(
    wsga_did(y ~ 1 | sgroup, data = d,
             unit = "unit", time = "time", treat = "D", bootstrap = FALSE),
    "exactly 2 unique"
  )
})

test_that("DiD errors when `treat` varies within unit", {
  d <- make_did_panel()
  d$D[d$unit == 1L][1] <- 1L - d$D[d$unit == 1L][1]  # flip one observation
  expect_error(
    wsga_did(y ~ 1 | sgroup, data = d,
             unit = "unit", time = "time", treat = "D", bootstrap = FALSE),
    "varies within unit"
  )
})

test_that("DiD errors when `sgroup` varies within unit", {
  d <- make_did_panel()
  d$sgroup[d$unit == 1L][1] <- 1L - d$sgroup[d$unit == 1L][1]
  expect_error(
    wsga_did(y ~ 1 | sgroup, data = d,
             unit = "unit", time = "time", treat = "D", bootstrap = FALSE),
    "varies within unit"
  )
})

# wsga_did() does not expose a `model` argument: fuzzy DiD is impossible
# by construction (the function signature excludes it).

test_that("DiD recovers per-subgroup treatment effects on the synthetic panel", {
  d <- make_did_panel(n_units = 600L, seed = 2L)
  fit <- wsga_did(y ~ 1 | sgroup, data = d,
                  unit = "unit", time = "time", treat = "D",
                  noipsw = TRUE, bootstrap = FALSE)
  # True effects: tau_0 = 1, tau_1 = 3
  expect_equal(fit$coefficients$g0, 1, tolerance = 0.25)
  expect_equal(fit$coefficients$g1, 3, tolerance = 0.25)
  expect_equal(fit$coefficients$diff, 2, tolerance = 0.3)
})

test_that("DiD output records design + unit/time/treat names + post_value", {
  d <- make_did_panel()
  fit <- wsga_did(y ~ 1 | sgroup, data = d,
                  unit = "unit", time = "time", treat = "D",
                  noipsw = TRUE, bootstrap = FALSE)
  expect_equal(fit$design, "did")
  expect_equal(fit$unit_name, "unit")
  expect_equal(fit$time_name, "time")
  expect_equal(fit$treat_name, "D")
  expect_equal(fit$post_value, 1L)
})

test_that("DiD cluster_var defaults to unit when bootstrap is on", {
  d <- make_did_panel(n_units = 60L, seed = 3L)
  fit <- wsga_did(y ~ 1 | sgroup, data = d,
                  unit = "unit", time = "time", treat = "D",
                  noipsw = TRUE, bootstrap = TRUE, bsreps = 20, seed = 1)
  expect_equal(fit$cluster_var_name, "unit")
  expect_equal(fit$bootstrap$N_clusters, 60L)
})

test_that("DiD bootstrap with IPW runs and produces positive SEs", {
  d <- make_did_panel(n_units = 100L, seed = 4L)
  fit <- wsga_did(y ~ m | sgroup, data = d,
                  unit = "unit", time = "time", treat = "D",
                  bootstrap = TRUE, bsreps = 20, seed = 5)
  expect_gt(fit$se$g0, 0)
  expect_gt(fit$se$g1, 0)
  expect_gt(fit$se$diff, 0)
})

test_that("DiD produces aggregate AND treated-only balance tables; RDD only aggregate", {
  # DiD path
  d <- make_did_panel(n_units = 200L, seed = 6L)
  fit_did <- wsga_did(y ~ m | sgroup, data = d,
                      unit = "unit", time = "time", treat = "D",
                      bootstrap = FALSE)
  expect_false(is.null(fit_did$balance$unweighted$aggregate))
  expect_false(is.null(fit_did$balance$unweighted$treated))
  expect_false(is.null(fit_did$balance$weighted$aggregate))
  expect_false(is.null(fit_did$balance$weighted$treated))

  # RDD path (using existing example data)
  data(rddsga_synth)
  fit_rdd <- wsga_rdd(y ~ m | sgroup, data = rddsga_synth,
                      running = ~ x, bwidth = 0.5, bootstrap = FALSE)
  expect_false(is.null(fit_rdd$balance$unweighted$aggregate))
  expect_null(fit_rdd$balance$unweighted$treated)
  expect_false(is.null(fit_rdd$balance$weighted$aggregate))
  expect_null(fit_rdd$balance$weighted$treated)
})

test_that("DiD with post_value override picks the right period", {
  d <- make_did_panel()
  fit <- wsga_did(y ~ 1 | sgroup, data = d,
                  unit = "unit", time = "time", treat = "D",
                  post_value = 1L,
                  noipsw = TRUE, bootstrap = FALSE)
  expect_equal(fit$post_value, 1L)

  expect_error(
    wsga_did(y ~ 1 | sgroup, data = d,
             unit = "unit", time = "time", treat = "D",
             post_value = 99L,
             noipsw = TRUE, bootstrap = FALSE),
    "not in the unique values of `time`"
  )
})
