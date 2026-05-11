data(rddsga_synth)

# Add a synthetic cluster column to the example data.  Clusters are larger
# than 30 by default; specific tests narrow this when needed.
make_clustered <- function(n_clusters = 50L, seed = 1L) {
  set.seed(seed)
  d <- rddsga_synth
  d$school <- sample.int(n_clusters, nrow(d), replace = TRUE)
  d
}

test_that("cluster bootstrap exposes N_clusters and uses pairs-cluster path", {
  d <- make_clustered(n_clusters = 50L)
  fit <- wsga_rdd(y ~ 1 | sgroup, data = d, running = ~ x, bwidth = 0.5,
              noipsw = TRUE, bootstrap = TRUE, bsreps = 20, seed = 1,
              cluster_var = "school")
  expect_true("N_clusters" %in% names(fit$bootstrap))
  expect_equal(fit$bootstrap$N_clusters, 50L)
  expect_equal(fit$cluster_var_name, "school")
})

test_that("row-level bootstrap (cluster_var = NULL) does not expose N_clusters", {
  fit <- wsga_rdd(y ~ 1 | sgroup, data = rddsga_synth,
                 running = ~ x, bwidth = 0.5,
                 noipsw = TRUE, bootstrap = TRUE, bsreps = 20, seed = 1)
  expect_false("N_clusters" %in% names(fit$bootstrap))
  expect_null(fit$cluster_var_name)
})

test_that("build_cluster_rep_data assigns fresh IDs per draw and returns orig_idx", {
  d <- data.frame(school = c("A", "A", "B", "B", "C"),
                  y      = 1:5)
  out <- build_cluster_rep_data(d, "school", d$school,
                                drawn = c("A", "A", "B"))
  # Returns a list with data + orig_idx
  expect_named(out, c("data", "orig_idx"))
  # A drawn twice → two distinct fresh IDs; B drawn once
  expect_equal(length(unique(out$data$school)), 3L)
  expect_true(all(grepl("__d[0-9]+$", out$data$school)))
  # Original outcomes preserved (each A row appears twice in the resample)
  expect_equal(sum(out$data$y %in% c(1, 2)), 4L)  # 2 A rows × 2 draws
  expect_equal(sum(out$data$y %in% c(3, 4)), 2L)  # 2 B rows × 1 draw
  # orig_idx maps resampled rows back to original row indices
  expect_equal(length(out$orig_idx), nrow(out$data))
  expect_equal(out$orig_idx, c(1L, 2L, 1L, 2L, 3L, 4L))
})

test_that("cluster_resample with strata respects strata sizes", {
  set.seed(99)
  uc <- letters[1:10]
  strata <- rep(c("x", "y"), each = 5)
  drawn <- cluster_resample(uc, strata)
  expect_equal(length(drawn), 10L)
  # Each stratum contributes its size in the draw, but order is x-then-y
  expect_true(all(drawn[1:5] %in% uc[1:5]))
  expect_true(all(drawn[6:10] %in% uc[6:10]))
})

test_that("block_var must be constant within cluster_var", {
  d <- make_clustered(n_clusters = 50L)
  # Make block_var vary within school -> should error.
  d$bad <- sample(c("east", "west"), nrow(d), replace = TRUE)
  expect_error(
    wsga_rdd(y ~ 1 | sgroup, data = d, running = ~ x, bwidth = 0.5,
             noipsw = TRUE, bootstrap = TRUE, bsreps = 5, seed = 1,
             cluster_var = "school", block_var = "bad"),
    "must be constant within `cluster_var`"
  )
})

test_that("block_var constant within cluster works", {
  d <- make_clustered(n_clusters = 50L)
  d$region <- ifelse(d$school <= 25, "east", "west")
  fit <- wsga_rdd(y ~ 1 | sgroup, data = d, running = ~ x, bwidth = 0.5,
              noipsw = TRUE, bootstrap = TRUE, bsreps = 20, seed = 1,
              cluster_var = "school", block_var = "region")
  expect_equal(fit$bootstrap$N_clusters, 50L)
  expect_gt(fit$se$g0, 0)
})

test_that("few-clusters (<30) emits the SE-understatement warning", {
  d <- make_clustered(n_clusters = 12L, seed = 2L)
  expect_warning(
    wsga_rdd(y ~ 1 | sgroup, data = d, running = ~ x, bwidth = 0.5,
             noipsw = TRUE, bootstrap = TRUE, bsreps = 5, seed = 1,
             cluster_var = "school"),
    "12 clusters"
  )
})

test_that("seed makes cluster bootstrap reproducible", {
  d <- make_clustered(n_clusters = 50L)
  fit_a <- wsga_rdd(y ~ 1 | sgroup, data = d, running = ~ x, bwidth = 0.5,
                   noipsw = TRUE, bootstrap = TRUE, bsreps = 20, seed = 42,
                   cluster_var = "school")
  fit_b <- wsga_rdd(y ~ 1 | sgroup, data = d, running = ~ x, bwidth = 0.5,
                   noipsw = TRUE, bootstrap = TRUE, bsreps = 20, seed = 42,
                   cluster_var = "school")
  expect_equal(fit_a$bootstrap$draws, fit_b$bootstrap$draws)
})
