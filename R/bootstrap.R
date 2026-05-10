#' Stratified bootstrap resample indices
#'
#' @param n Total number of observations.
#' @param strata Optional factor/integer vector for stratification.
#' @return Integer vector of row indices (length = n, with replacement within
#'   each stratum).
#' @noRd
stratified_resample <- function(n, strata = NULL) {
  if (is.null(strata)) return(sample.int(n, replace = TRUE))
  idx <- seq_len(n)
  out <- integer(n)
  for (s in unique(strata)) {
    in_s <- which(strata == s)
    out[in_s] <- sample(in_s, replace = TRUE)
  }
  out
}


#' Bootstrap variance-covariance matrix and inference
#'
#' Mirrors the `myboo` subroutine in `wsga.ado`. Each bootstrap replication
#' resamples the full dataset (optionally stratified), refits the IPW
#' propensity score and main regression, and records the two per-subgroup
#' coefficients. The returned variance matrix is used for SE calculation and
#' empirical CIs/p-values are derived from the bootstrap distribution.
#'
#' @param run_one_rep Function `(data_b) -> c(b_g0, b_g1)` that fits one
#'   bootstrap replicate and returns the two coefficient estimates.
#' @param data Data frame to resample.
#' @param B Integer: number of bootstrap replications.
#' @param est Numeric length-2: the point estimates (b_g0, b_g1) from the
#'   original data.
#' @param block_var Optional character: column name to use as bootstrap strata.
#' @param seed Optional integer RNG seed.
#' @param fixed_fs Logical: fixed first-stage bootstrap for IV (handled by
#'   `run_one_rep` closure).
#' @return A list:
#'   - `draws` (B×2 matrix of bootstrap coefficient draws).
#'   - `vcov` (2×2 variance matrix).
#'   - `pval` (length-2 empirical p-values, formula `(1 + count)/(B + 1)`).
#'   - `ci` (2×2 matrix: rows = g0/g1, cols = lb/ub, empirical 2.5/97.5).
#' @noRd
run_bootstrap <- function(run_one_rep, data, B, est,
                          block_var = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  n      <- nrow(data)
  strata <- if (!is.null(block_var)) data[[block_var]] else NULL
  draws  <- matrix(NA_real_, nrow = B, ncol = 2,
                   dimnames = list(NULL, c("g0", "g1")))

  message(sprintf("Bootstrap replications (%d):", B))
  for (i in seq_len(B)) {
    idx <- stratified_resample(n, strata)
    tryCatch({
      draws[i, ] <- run_one_rep(data[idx, ])
    }, error = function(e) {
      # leave row as NA; warn once at end
    })
    if (i %% 10 == 0 || i == B) cat(sprintf("\r  %d/%d", i, B))
  }
  cat("\n")

  # Drop failed replicates
  ok <- complete.cases(draws)
  if (any(!ok)) warning(sprintf("%d bootstrap replicates failed and were dropped.",
                                sum(!ok)))
  draws_ok <- draws[ok, , drop = FALSE]
  B_ok     <- nrow(draws_ok)

  # Add diff = g1 - g0 as third column
  draws_ok <- cbind(draws_ok, diff = draws_ok[, "g1"] - draws_ok[, "g0"])

  # Variance-covariance from bootstrap distribution (g0, g1 only)
  V_boot <- var(draws_ok[, c("g0", "g1")])

  # Point estimates for all three parameters (needed for empirical p-values)
  est_all <- c(est, diff = est[2] - est[1])

  # Empirical p-values: (1 + #{|draw| >= |est|}) / (B_ok + 1)
  pval <- sapply(seq_len(3), function(g) {
    cnt <- sum(abs(draws_ok[, g]) >= abs(est_all[g]))
    (1 + cnt) / (B_ok + 1)
  })
  names(pval) <- c("g0", "g1", "diff")

  # Empirical CIs: 2.5 / 97.5 percentiles for all three parameters
  ci <- apply(draws_ok, 2, quantile, probs = c(0.025, 0.975))
  rownames(ci) <- c("lb", "ub")

  list(draws = draws_ok, vcov = V_boot, pval = pval, ci = ci,
       B_ok = B_ok)
}
