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
#' @return A data frame with N_drawn × (rows-per-cluster) rows.
#' @noRd
build_cluster_rep_data <- function(data, cluster_var, cluster_ids, drawn) {
  parts <- vector("list", length(drawn))
  for (j in seq_along(drawn)) {
    rows <- which(cluster_ids == drawn[j])
    sub <- data[rows, , drop = FALSE]
    sub[[cluster_var]] <- paste0(as.character(drawn[j]), "__d", j)
    parts[[j]] <- sub
  }
  out <- do.call(rbind, parts)
  rownames(out) <- NULL
  out
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
                          seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  n <- nrow(data)

  if (is.null(cluster_var)) {
    strata <- if (!is.null(block_var)) data[[block_var]] else NULL
    resample <- function() {
      idx <- stratified_resample(n, strata)
      data[idx, , drop = FALSE]
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
      build_cluster_rep_data(data, cluster_var, cluster_ids, drawn)
    }
  }

  draws <- matrix(NA_real_, nrow = B, ncol = 2,
                  dimnames = list(NULL, c("g0", "g1")))

  message(sprintf("Bootstrap replications (%d):", B))
  for (i in seq_len(B)) {
    tryCatch({
      draws[i, ] <- run_one_rep(resample())
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
    cnt <- sum(abs(draws_ok[, g]) >= abs(est_all[g]))
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
