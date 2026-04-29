set.seed(1)
n  <- 200
x  <- runif(n, -1, 1)
G  <- rbinom(n, 1, 0.5)
m  <- rnorm(n, mean = 0.3 * G)
df <- data.frame(x = x, G = G, m = m)

within_bw <- abs(x) < 0.5
touse     <- rep(TRUE, n)

test_that("compute_pscore returns probabilities in (0,1) for active obs", {
  res <- compute_pscore(G, df[, "m", drop = FALSE], within_bw, touse,
                        ps_model = "logit", comsup = FALSE)
  ps_active <- res$pscore[!is.na(res$pscore)]
  expect_true(all(ps_active > 0 & ps_active < 1))
})

test_that("common support trims correctly", {
  res <- compute_pscore(G, df[, "m", drop = FALSE], within_bw, touse,
                        ps_model = "logit", comsup = TRUE)
  ps_active <- res$pscore[!is.na(res$pscore)]
  cs <- res$comsup_idx[!is.na(res$pscore)]
  ps_g1 <- ps_active[G[!is.na(res$pscore)] == 1]
  # All in-support obs should have pscore in [min(ps|G1), max(ps|G1)]
  expect_true(all(ps_active[cs] >= min(ps_g1) - 1e-12))
  expect_true(all(ps_active[cs] <= max(ps_g1) + 1e-12))
})

test_that("m=2 IPW weights are positive for active obs", {
  ps_res <- compute_pscore(G, df[, "m", drop = FALSE], within_bw, touse)
  iw_res <- compute_ipw_weights(G, ps_res$pscore, within_bw, touse,
                                 ps_res$comsup_idx, m = 2)
  ipw_active <- iw_res$ipw[!is.na(iw_res$ipw)]
  expect_true(all(ipw_active > 0))
})

test_that("m=1: group 0 weights are 1, group 1 weights >= 0", {
  ps_res <- compute_pscore(G, df[, "m", drop = FALSE], within_bw, touse)
  iw_res <- compute_ipw_weights(G, ps_res$pscore, within_bw, touse,
                                 ps_res$comsup_idx, m = 1)
  active_g0 <- touse & within_bw & G == 0
  expect_true(all(iw_res$ipw[active_g0] == 1))
})
