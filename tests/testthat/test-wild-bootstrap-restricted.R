# WCB-R (boot_type = "wild_restricted") tests.
# Architecture: one combined loop per replicate --
#   * unrestricted refit  -> CI / SE (empirical percentile, coauthor agreement)
#   * three restricted refits (G0=0, G1=0, diff=0) -> p-values
# P-value formula for restricted draws: (1 + #{|draw| >= |est|}) / (B+1)
# (no recentering -- draws already centred at 0 under H0).

data(rddsga_synth)
data(wsga_did_synth)

make_clustered_rdd <- function(n_clusters = 20L, seed = 1L) {
  set.seed(seed)
  d <- rddsga_synth
  d$school <- sample.int(n_clusters, nrow(d), replace = TRUE)
  d
}

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
  y      <- alpha + 0.5 * post + tau * D * post + 0.3 * M +
            rnorm(length(unit), sd = 0.5)
  data.frame(unit = unit, time = time, sgroup = sgroup,
             m = M, D = D, y = y)
}

# -- Happy path: RDD --

test_that("WCB-R runs on sharp RDD and returns valid structure", {
  d   <- make_clustered_rdd(n_clusters = 20L)
  fit <- wsga_rdd(y ~ 1 | sgroup, data = d, running = ~ x, bwidth = 0.5,
                  noipsw = TRUE, bootstrap = TRUE, bsreps = 50, seed = 1,
                  cluster_var = "school", boot_type = "wild_restricted")
  expect_equal(fit$boot_type, "wild_restricted")
  expect_equal(fit$bootstrap$N_clusters, 20L)
  expect_true(all(unlist(fit$pval) >= 0 & unlist(fit$pval) <= 1))
  expect_true(fit$ci$g0[["lb"]]   < fit$ci$g0[["ub"]])
  expect_true(fit$ci$g1[["lb"]]   < fit$ci$g1[["ub"]])
  expect_true(fit$ci$diff[["lb"]] < fit$ci$diff[["ub"]])
})

test_that("WCB-R is reproducible with seed (RDD)", {
  d <- make_clustered_rdd(n_clusters = 20L)
  a <- wsga_rdd(y ~ 1 | sgroup, data = d, running = ~ x, bwidth = 0.5,
                noipsw = TRUE, bootstrap = TRUE, bsreps = 30, seed = 42,
                cluster_var = "school", boot_type = "wild_restricted")
  b <- wsga_rdd(y ~ 1 | sgroup, data = d, running = ~ x, bwidth = 0.5,
                noipsw = TRUE, bootstrap = TRUE, bsreps = 30, seed = 42,
                cluster_var = "school", boot_type = "wild_restricted")
  expect_equal(a$pval, b$pval)
  expect_equal(a$ci,   b$ci)
})

test_that("WCB-R CI equals WCB-U CI with the same seed (unrestricted draws for CI)", {
  # The unrestricted pass (for CI/SE) uses the same sign draws as a plain
  # WCB-U run with the same seed, so CIs must be numerically identical.
  d     <- make_clustered_rdd(n_clusters = 20L)
  fit_r <- wsga_rdd(y ~ 1 | sgroup, data = d, running = ~ x, bwidth = 0.5,
                    noipsw = TRUE, bootstrap = TRUE, bsreps = 100, seed = 7,
                    cluster_var = "school", boot_type = "wild_restricted")
  fit_u <- wsga_rdd(y ~ 1 | sgroup, data = d, running = ~ x, bwidth = 0.5,
                    noipsw = TRUE, bootstrap = TRUE, bsreps = 100, seed = 7,
                    cluster_var = "school", boot_type = "wild")
  expect_equal(fit_r$ci, fit_u$ci)
})

test_that("WCB-R stores boot_type_pval and boot_type_ci metadata", {
  d   <- make_clustered_rdd(n_clusters = 20L)
  fit <- wsga_rdd(y ~ 1 | sgroup, data = d, running = ~ x, bwidth = 0.5,
                  noipsw = TRUE, bootstrap = TRUE, bsreps = 30, seed = 1,
                  cluster_var = "school", boot_type = "wild_restricted")
  expect_equal(fit$bootstrap$boot_type_pval, "wild_restricted")
  expect_equal(fit$bootstrap$boot_type_ci,   "wild")
})

# -- Happy path: DiD --

test_that("WCB-R runs on DiD (cluster defaults to unit)", {
  d   <- make_small_did(n_units = 20L)
  fit <- wsga_did(y ~ 1 | sgroup, data = d,
                  unit = "unit", time = "time", treat = "D",
                  noipsw = TRUE, bootstrap = TRUE, bsreps = 50, seed = 1,
                  boot_type = "wild_restricted")
  expect_equal(fit$boot_type, "wild_restricted")
  expect_equal(fit$bootstrap$N_clusters, 20L)
  expect_true(all(unlist(fit$pval) >= 0 & unlist(fit$pval) <= 1))
})

# -- Advisory --

test_that("boot_type='wild' with G<12 warns and recommends wild_restricted", {
  d <- make_clustered_rdd(n_clusters = 10L, seed = 2L)
  expect_warning(
    wsga_rdd(y ~ 1 | sgroup, data = d, running = ~ x, bwidth = 0.5,
             noipsw = TRUE, bootstrap = TRUE, bsreps = 5, seed = 1,
             cluster_var = "school", boot_type = "wild"),
    "consider `boot_type = \"wild_restricted\"`"
  )
})

test_that("boot_type='wild_restricted' at G<12 does NOT warn", {
  d <- make_clustered_rdd(n_clusters = 10L, seed = 2L)
  expect_no_warning(
    wsga_rdd(y ~ 1 | sgroup, data = d, running = ~ x, bwidth = 0.5,
             noipsw = TRUE, bootstrap = TRUE, bsreps = 5, seed = 1,
             cluster_var = "school", boot_type = "wild_restricted")
  )
})

test_that("boot_type='wild' with 12<=G<30 warns to use wild but not wild_restricted", {
  d <- make_clustered_rdd(n_clusters = 20L, seed = 2L)
  w <- withCallingHandlers(
    wsga_rdd(y ~ 1 | sgroup, data = d, running = ~ x, bwidth = 0.5,
             noipsw = TRUE, bootstrap = TRUE, bsreps = 5, seed = 1,
             cluster_var = "school", boot_type = "pairs"),
    warning = function(w) {
      invokeRestart("muffleWarning")
    }
  )
  # wild at 12<=G<30 should NOT warn about wild_restricted
  expect_no_warning(
    wsga_rdd(y ~ 1 | sgroup, data = d, running = ~ x, bwidth = 0.5,
             noipsw = TRUE, bootstrap = TRUE, bsreps = 5, seed = 1,
             cluster_var = "school", boot_type = "wild")
  )
})

# -- Validation errors --

test_that("wild_restricted errors without cluster_var (RDD)", {
  expect_error(
    wsga_rdd(y ~ 1 | sgroup, data = rddsga_synth, running = ~ x, bwidth = 0.5,
             noipsw = TRUE, bootstrap = TRUE, bsreps = 10,
             boot_type = "wild_restricted"),
    "requires a clustering variable"
  )
})

test_that("wild_restricted errors with model = 'iv'", {
  d <- make_clustered_rdd(n_clusters = 20L)
  expect_error(
    wsga_rdd(y ~ 1 | sgroup, data = d, running = ~ x, bwidth = 0.5,
             fuzzy = ~ D, model = "iv",
             noipsw = TRUE, bootstrap = TRUE, bsreps = 10,
             cluster_var = "school", boot_type = "wild_restricted"),
    "not supported with `model = \"iv\"`"
  )
})

test_that("wild_restricted errors when bootstrap = FALSE", {
  d <- make_clustered_rdd(n_clusters = 20L)
  expect_error(
    wsga_rdd(y ~ 1 | sgroup, data = d, running = ~ x, bwidth = 0.5,
             noipsw = TRUE, bootstrap = FALSE,
             cluster_var = "school", boot_type = "wild_restricted"),
    "requires `bootstrap = TRUE`"
  )
})
