#' Stratified row-level bootstrap resample indices
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


#' Pairs cluster bootstrap: draw clusters with replacement.
#'
#' @param unique_clusters Vector of unique cluster identifiers (any type).
#' @param cluster_strata Optional vector aligned with `unique_clusters`: each
#'   cluster's stratum. If supplied, clusters are sampled with replacement
#'   *within* each stratum.
#' @return A vector of drawn cluster IDs, same length and type as
#'   `unique_clusters`. Duplicates are expected — that's the whole point.
#' @noRd
cluster_resample <- function(unique_clusters, cluster_strata = NULL) {
  if (is.null(cluster_strata)) return(sample(unique_clusters, replace = TRUE))
  out <- unique_clusters[integer(0)]   # zero-length, preserves type
  for (s in unique(cluster_strata)) {
    mask <- cluster_strata == s
    out <- c(out, sample(unique_clusters[mask], sum(mask), replace = TRUE))
  }
  out
}


#' Build a resampled dataset by replicating each drawn cluster's rows.
#'
#' Fresh-IDs recipe (Cameron–Gelbach–Miller): the cluster identifier column
#' is rewritten as `<orig_id>__d<j>` for the j-th draw, so a cluster drawn N
#' times yields N distinct copies. Required for any downstream specification
#' that uses cluster-level fixed effects (DiD); a no-op when there are no FE.
#'
#' @param data Original data frame.
#' @param cluster_var Character: name of the cluster column.
#' @param cluster_ids The cluster column values, length nrow(data).
#' @param drawn Drawn cluster IDs (output of `cluster_resample`).
#' @param unit_var Optional character: name of a "nested" identifier column
#'   (e.g. the DiD `unit` column when `cluster_var` is something coarser
#'   like "state"). If supplied and different from `cluster_var`, fresh IDs
#'   are also propagated to that column so unit FE remain identified when
#'   two draws of the same cluster bring in overlapping units.
#' @return A list with two elements:
#'   - `data`: data frame with `length(drawn) × rows-per-cluster` rows.
#'   - `orig_idx`: integer vector of the same length, where
#'     `orig_idx[k]` is the row in the original `data` that the k-th
#'     row of the resampled data came from. Used by the `fixed_ps`
#'     code path to look up the original-sample propensity score.
#' @noRd
build_cluster_rep_data <- function(data, cluster_var, cluster_ids, drawn,
                                   unit_var = NULL) {
  parts          <- vector("list", length(drawn))
  orig_idx_parts <- vector("list", length(drawn))
  for (j in seq_along(drawn)) {
    rows <- which(cluster_ids == drawn[j])
    sub <- data[rows, , drop = FALSE]
    sub[[cluster_var]] <- paste0(as.character(drawn[j]), "__d", j)
    if (!is.null(unit_var) && unit_var != cluster_var) {
      sub[[unit_var]] <- paste0(as.character(sub[[unit_var]]), "__d", j)
    }
    parts[[j]]          <- sub
    orig_idx_parts[[j]] <- rows
  }
  out <- do.call(rbind, parts)
  rownames(out) <- NULL
  list(data = out, orig_idx = unlist(orig_idx_parts))
}


#' Bootstrap variance-covariance matrix and inference
#'
#' Mirrors the `myboo` subroutine in `wsga.ado`. Each replicate resamples the
#' data (rows or whole clusters), refits the IPW propensity score and main
#' regression, and records the two per-subgroup coefficients.
#'
#' Two regimes:
#' - `cluster_var = NULL`: row-level bootstrap (optionally stratified by
#'   `block_var`).
#' - `cluster_var` set: pairs cluster bootstrap. Whole clusters are drawn
#'   with replacement; each drawn cluster contributes all its rows; cluster
#'   IDs are made unique per draw (Cameron–Gelbach–Miller). `block_var`, if
#'   supplied, must be constant within `cluster_var` and is used to stratify
#'   at the cluster level.
#'
#' @param run_one_rep Function `(data_b) -> c(b_g0, b_g1)`.
#' @param data Data frame to resample.
#' @param B Integer: number of bootstrap replications.
#' @param est Numeric length-2: point estimates (b_g0, b_g1).
#' @param cluster_var Optional character: column name whose values define
#'   clusters. When set, runs a pairs cluster bootstrap.
#' @param block_var Optional character: stratification variable. Semantics
#'   adapt: strata of rows when `cluster_var` is NULL, strata of clusters
#'   otherwise.
#' @param seed Optional integer RNG seed.
#' @return A list with `draws`, `vcov`, `pval`, `ci`, `B_ok`, `failed`, and —
#'   when `cluster_var` is set — `N_clusters`.
#' @noRd
run_bootstrap <- function(run_one_rep, data, B, est,
                          cluster_var = NULL,
                          block_var = NULL,
                          unit_var = NULL,
                          seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  n <- nrow(data)

  if (is.null(cluster_var)) {
    strata <- if (!is.null(block_var)) data[[block_var]] else NULL
    resample <- function() {
      idx <- stratified_resample(n, strata)
      list(data = data[idx, , drop = FALSE], orig_idx = idx)
    }
    N_clusters <- NULL
  } else {
    cluster_ids <- data[[cluster_var]]
    unique_clusters <- unique(cluster_ids[!is.na(cluster_ids)])
    N_clusters <- length(unique_clusters)

    cluster_strata <- NULL
    if (!is.null(block_var)) {
      strata_vec <- data[[block_var]]
      # Each cluster must map to exactly one stratum value
      strata_per_cluster <- tapply(strata_vec, cluster_ids,
                                   function(s) length(unique(s)))
      bad <- names(strata_per_cluster)[strata_per_cluster > 1]
      if (length(bad) > 0) {
        stop(sprintf(
          "`block_var` must be constant within `cluster_var`. Violating clusters: %s%s",
          paste(head(bad, 3), collapse = ", "),
          if (length(bad) > 3) sprintf(" (and %d more)", length(bad) - 3) else ""
        ))
      }
      cluster_strata <- vapply(
        unique_clusters,
        function(cc) strata_vec[which(cluster_ids == cc)[1]],
        FUN.VALUE = strata_vec[1]
      )
    }

    resample <- function() {
      drawn <- cluster_resample(unique_clusters, cluster_strata)
      build_cluster_rep_data(data, cluster_var, cluster_ids, drawn,
                             unit_var = unit_var)
    }
    # build_cluster_rep_data already returns list(data = ..., orig_idx = ...)
  }

  draws <- matrix(NA_real_, nrow = B, ncol = 2,
                  dimnames = list(NULL, c("g0", "g1")))

  message(sprintf("Bootstrap replications (%d):", B))
  for (i in seq_len(B)) {
    r <- resample()
    tryCatch({
      draws[i, ] <- run_one_rep(r$data, r$orig_idx)
    }, error = function(e) {
      # leave row as NA; warn once at end
    })
    if (i %% 10 == 0 || i == B) cat(sprintf("\r  %d/%d", i, B))
  }
  cat("\n")

  ok <- complete.cases(draws)
  if (any(!ok)) {
    msg <- sprintf("%d of %d bootstrap replicates failed and were dropped",
                   sum(!ok), B)
    if (!is.null(cluster_var)) {
      msg <- paste0(msg, " (often caused by drawn samples lacking variation in `sgroup` or a moderator; consider increasing `bsreps`).")
    } else {
      msg <- paste0(msg, ".")
    }
    warning(msg)
  }
  draws_ok <- draws[ok, , drop = FALSE]
  B_ok     <- nrow(draws_ok)

  draws_ok <- cbind(draws_ok, diff = draws_ok[, "g1"] - draws_ok[, "g0"])
  V_boot   <- var(draws_ok[, c("g0", "g1")])
  est_all  <- c(est, diff = est[2] - est[1])

  pval <- sapply(seq_len(3), function(g) {
    cnt <- sum(abs(draws_ok[, g] - est_all[g]) >= abs(est_all[g]))
    (1 + cnt) / (B_ok + 1)
  })
  names(pval) <- c("g0", "g1", "diff")

  ci <- apply(draws_ok, 2, quantile, probs = c(0.025, 0.975))
  rownames(ci) <- c("lb", "ub")

  result <- list(draws = draws_ok, vcov = V_boot, pval = pval, ci = ci,
                 B_ok = B_ok, failed = sum(!ok))
  if (!is.null(cluster_var)) result$N_clusters <- N_clusters
  result
}


#' Wild cluster bootstrap (unrestricted, Rademacher signs)
#'
#' WCB-U for weighted least squares. The propensity score and IPW weights are
#' fixed at their original-sample values -- WCB conditions on the data and only
#' sign-flips residuals at the cluster level, so it does **not** propagate
#' uncertainty from the propensity-score estimation step. Use pairs bootstrap
#' if propagating IPW uncertainty matters more than small-G size control.
#'
#' Per replicate `b`:
#'   1. Draw cluster signs `v_g in {-1, +1}` (Rademacher), one per cluster.
#'   2. Form `y*_i = yhat_i + v_{g(i)} * e_i`, where `yhat` and `e` are the
#'      original-sample fitted values and raw residuals.
#'   3. Refit the weighted OLS on `(X, y*)` with the original weights.
#'   4. Record the two subgroup coefficients.
#'
#' The fit uses a cached QR decomposition (X is fixed across replicates), so
#' each replicate is O(np), not O(np^2).
#'
#' @param fit The original `lm` fit returned by `run_model()`.
#' @param X Full design matrix (n x p) from `build_*_design_matrix()`.
#' @param y Outcome vector (length n).
#' @param w Weights vector (length n; kernel * IPW).
#' @param active Logical mask: rows with `final_wt > 0` used in estimation.
#' @param cluster_vec Cluster identifier (length n; any type).
#' @param coef_g0_name,coef_g1_name Column names for the two coefficients of
#'   interest (always `"G0_Z"` / `"G1_Z"` in this package).
#' @param B Number of bootstrap replications.
#' @param seed Optional RNG seed.
#' @return Same shape as `run_bootstrap()`: list with `draws`, `vcov`, `pval`,
#'   `ci`, `B_ok`, `failed`, `N_clusters`.
#' @noRd
run_wild_bootstrap <- function(fit, X, y, w, active, cluster_vec,
                               coef_g0_name, coef_g1_name,
                               B, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # Restrict to active rows (the only rows that influenced the fit)
  X_a <- X[active, , drop = FALSE]
  y_a <- y[active]
  w_a <- w[active]
  c_a <- cluster_vec[active]

  # Drop aliased columns (NA coefficients) so the system is full rank
  beta_hat <- coef(fit)
  non_aliased <- names(beta_hat)[!is.na(beta_hat)]
  beta_clean <- beta_hat[non_aliased]
  X_clean    <- X_a[, non_aliased, drop = FALSE]

  # Guard: a coefficient of interest aliased away (e.g. perfect Z-by-G collinearity
  # on the active sample) would silently NA every draw.  Error early instead.
  missing_coefs <- setdiff(c(coef_g0_name, coef_g1_name), non_aliased)
  if (length(missing_coefs) > 0) {
    stop(sprintf(
      "Wild cluster bootstrap: coefficient(s) %s aliased away by `lm()` -- likely perfect collinearity between the subgroup and treatment indicators on the active sample.",
      paste(shQuote(missing_coefs), collapse = ", ")))
  }

  # Fitted values and raw residuals (in the un-transformed scale)
  y_hat <- as.numeric(X_clean %*% beta_clean)
  e     <- y_a - y_hat

  # Cluster bookkeeping
  if (any(is.na(c_a))) {
    stop(sprintf(
      "Wild cluster bootstrap: cluster variable has %d NA value(s) in the active sample.",
      sum(is.na(c_a))))
  }
  unique_clusters <- unique(c_a)
  N_clusters      <- length(unique_clusters)
  if (N_clusters < 2L) {
    stop(sprintf(
      "Wild cluster bootstrap requires at least 2 clusters; found %d in the active sample.",
      N_clusters))
  }
  cluster_idx <- match(as.character(c_a), as.character(unique_clusters))

  # Cache QR of the weighted design -- refits are then a single qr.coef() call.
  sqrt_w <- sqrt(w_a)
  X_w    <- sqrt_w * X_clean
  qr_X   <- qr(X_w)
  if (qr_X$rank < ncol(X_clean)) {
    stop("Wild cluster bootstrap: weighted design matrix is rank-deficient.")
  }

  draws <- matrix(NA_real_, nrow = B, ncol = 2,
                  dimnames = list(NULL, c("g0", "g1")))

  message(sprintf("Wild cluster bootstrap (%d):", B))
  for (i in seq_len(B)) {
    signs   <- sample(c(-1, 1), N_clusters, replace = TRUE)
    v       <- signs[cluster_idx]
    y_star  <- y_hat + v * e
    y_w     <- sqrt_w * y_star
    beta_b  <- tryCatch(qr.coef(qr_X, y_w), error = function(e) NULL)
    if (!is.null(beta_b)) {
      names(beta_b) <- non_aliased
      draws[i, "g0"] <- beta_b[[coef_g0_name]]
      draws[i, "g1"] <- beta_b[[coef_g1_name]]
    }
    if (i %% 10 == 0 || i == B) cat(sprintf("\r  %d/%d", i, B))
  }
  cat("\n")

  ok       <- complete.cases(draws)
  if (any(!ok)) {
    warning(sprintf("%d of %d wild bootstrap replicates failed and were dropped.",
                    sum(!ok), B))
  }
  draws_ok <- draws[ok, , drop = FALSE]
  B_ok     <- nrow(draws_ok)

  draws_ok <- cbind(draws_ok, diff = draws_ok[, "g1"] - draws_ok[, "g0"])
  V_boot   <- var(draws_ok[, c("g0", "g1")])
  est_all  <- c(g0 = beta_clean[[coef_g0_name]],
                g1 = beta_clean[[coef_g1_name]])
  est_all  <- c(est_all, diff = est_all[["g1"]] - est_all[["g0"]])

  pval <- sapply(seq_len(3), function(g) {
    cnt <- sum(abs(draws_ok[, g] - est_all[g]) >= abs(est_all[g]))
    (1 + cnt) / (B_ok + 1)
  })
  names(pval) <- c("g0", "g1", "diff")

  ci <- apply(draws_ok, 2, quantile, probs = c(0.025, 0.975))
  rownames(ci) <- c("lb", "ub")

  list(draws = draws_ok, vcov = V_boot, pval = pval, ci = ci,
       B_ok = B_ok, failed = sum(!ok), N_clusters = N_clusters)
}
