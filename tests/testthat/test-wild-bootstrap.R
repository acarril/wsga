# Wild cluster bootstrap (boot_type = "wild")
# CGM 2008 motivates this for small-G; tests cover the happy path + the
# validation errors that distinguish WCB from pairs.

data(rddsga_synth)
data(wsga_did_synth)

make_clustered_rdd <- function(n_clusters = 20L, seed = 1L) {
  set.seed(seed)
  d <- rddsga_synth
  d$school <- sample.int(n_clusters, nrow(d), replace = TRUE)
  d
}

# Small DiD panel: 20 units, 2 periods.  Small enough that WCB matters.
make_small_did <- function(n_units = 20L, seed = 1L) {
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

# ── RDD: happy path ──────────────────────────────────────────────────────────

test_that("WCB runs on sharp RDD with explicit cluster_var", {
  d <- make_clustered_rdd(n_clusters = 20L)
  fit <- wsga_rdd(y ~ 1 | sgroup, data = d, running = ~ x, bwidth = 0.5,
                  noipsw = TRUE, bootstrap = TRUE, bsreps = 50, seed = 1,
                  cluster_var = "school", boot_type = "wild")
  expect_equal(fit$boot_type, "wild")
  expect_equal(fit$bootstrap$N_clusters, 20L)
  expect_equal(nrow(fit$bootstrap$draws), 50L)
  expect_true(all(c("g0", "g1", "diff") %in% colnames(fit$bootstrap$draws)))
  expect_gt(fit$se$g0, 0)
  expect_gt(fit$se$g1, 0)
  expect_true(fit$pval$g0 >= 0 && fit$pval$g0 <= 1)
})

test_that("WCB seed makes results reproducible (sharp RDD)", {
  d <- make_clustered_rdd(n_clusters = 20L)
  a <- wsga_rdd(y ~ 1 | sgroup, data = d, running = ~ x, bwidth = 0.5,
                noipsw = TRUE, bootstrap = TRUE, bsreps = 30, seed = 42,
                cluster_var = "school", boot_type = "wild")
  b <- wsga_rdd(y ~ 1 | sgroup, data = d, running = ~ x, bwidth = 0.5,
                noipsw = TRUE, bootstrap = TRUE, bsreps = 30, seed = 42,
                cluster_var = "school", boot_type = "wild")
  expect_equal(a$bootstrap$draws, b$bootstrap$draws)
})

# ── DiD: happy path (cluster_var defaults to unit) ───────────────────────────

test_that("WCB runs on DiD using the default unit clustering", {
  d <- make_small_did(n_units = 20L)
  fit <- wsga_did(y ~ 1 | sgroup, data = d,
                  unit = "unit", time = "time", treat = "D",
                  noipsw = TRUE, bootstrap = TRUE, bsreps = 50, seed = 1,
                  boot_type = "wild")
  expect_equal(fit$boot_type, "wild")
  expect_equal(fit$bootstrap$N_clusters, 20L)
  expect_equal(fit$cluster_var_name, "unit")
  expect_gt(fit$se$diff, 0)
})

# ── Validation errors ────────────────────────────────────────────────────────

test_that("WCB errors when cluster_var is NULL (RDD case)", {
  expect_error(
    wsga_rdd(y ~ 1 | sgroup, data = rddsga_synth, running = ~ x, bwidth = 0.5,
             noipsw = TRUE, bootstrap = TRUE, bsreps = 10,
             boot_type = "wild"),
    "requires a clustering variable"
  )
})

test_that("WCB errors with model = \"iv\" (fuzzy RDD)", {
  d <- make_clustered_rdd(n_clusters = 20L)
  expect_error(
    wsga_rdd(y ~ 1 | sgroup, data = d, running = ~ x, bwidth = 0.5,
             fuzzy = ~ D, model = "iv",
             noipsw = TRUE, bootstrap = TRUE, bsreps = 10,
             cluster_var = "school", boot_type = "wild"),
    "not supported with `model = \"iv\"`"
  )
})

test_that("WCB errors when bootstrap = FALSE", {
  d <- make_clustered_rdd(n_clusters = 20L)
  expect_error(
    wsga_rdd(y ~ 1 | sgroup, data = d, running = ~ x, bwidth = 0.5,
             noipsw = TRUE, bootstrap = FALSE,
             cluster_var = "school", boot_type = "wild"),
    "requires `bootstrap = TRUE`"
  )
})

# ── Advisory behavior ────────────────────────────────────────────────────────

test_that("pairs at G<30 warns and recommends WCB", {
  d <- make_clustered_rdd(n_clusters = 12L, seed = 2L)
  expect_warning(
    wsga_rdd(y ~ 1 | sgroup, data = d, running = ~ x, bwidth = 0.5,
             noipsw = TRUE, bootstrap = TRUE, bsreps = 5, seed = 1,
             cluster_var = "school", boot_type = "pairs"),
    "consider `boot_type = \"wild\"`"
  )
})

test_that("wild at G<30 does NOT emit the pairs advisory", {
  d <- make_clustered_rdd(n_clusters = 12L, seed = 2L)
  # No warning expected — WCB is the recommended fix
  expect_warning(
    wsga_rdd(y ~ 1 | sgroup, data = d, running = ~ x, bwidth = 0.5,
             noipsw = TRUE, bootstrap = TRUE, bsreps = 5, seed = 1,
             cluster_var = "school", boot_type = "wild"),
    NA   # i.e., no warning matching this call
  )
})

# ── Sanity vs analytical cluster-robust SE ───────────────────────────────────
# WCB SE should be in the same ballpark as the analytical cluster-robust SE
# on a non-pathological DGP.  Allows a factor of 3 in either direction — this
# is a sanity floor, not a tight comparison.
test_that("WCB SE is in the same order of magnitude as analytical cluster SE", {
  d <- make_clustered_rdd(n_clusters = 40L)
  fa <- wsga_rdd(y ~ 1 | sgroup, data = d, running = ~ x, bwidth = 0.5,
                 noipsw = TRUE, bootstrap = FALSE,
                 cluster_var = "school", vce = "cluster")
  fw <- wsga_rdd(y ~ 1 | sgroup, data = d, running = ~ x, bwidth = 0.5,
                 noipsw = TRUE, bootstrap = TRUE, bsreps = 200, seed = 1,
                 cluster_var = "school", boot_type = "wild")
  ratio_g0   <- fw$se$g0   / fa$se$g0
  ratio_diff <- fw$se$diff / fa$se$diff
  expect_gt(ratio_g0,   1/3)
  expect_lt(ratio_g0,   3)
  expect_gt(ratio_diff, 1/3)
  expect_lt(ratio_diff, 3)
})

# ── Defensive errors ─────────────────────────────────────────────────────────

test_that("WCB errors when cluster_var has NAs in the active sample", {
  d <- make_clustered_rdd(n_clusters = 20L)
  d$school[1:50] <- NA   # NA-ize some cluster values within bandwidth
  expect_error(
    wsga_rdd(y ~ 1 | sgroup, data = d, running = ~ x, bwidth = 0.5,
             noipsw = TRUE, bootstrap = TRUE, bsreps = 10, seed = 1,
             cluster_var = "school", boot_type = "wild"),
    "cluster variable has .* NA"
  )
})
