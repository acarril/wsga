#' Compute balance table and global balance statistics
#'
#' Mirrors the `balancematrix` subroutine in `rddsga.ado`. Can be called for
#' both the unweighted and IPW-weighted cases.
#'
#' @param balance_vars Character vector of balance variable names.
#' @param data Data frame containing all variables.
#' @param running Character: name of the running variable (centered at cutoff).
#' @param G Integer 0/1 vector: subgroup indicator.
#' @param sgroup0 Integer 0/1 vector: complement of G.
#' @param within_bw Logical: within bandwidth.
#' @param touse Logical: non-missing observations.
#' @param weights Numeric vector of observation weights (kernel * IPW or
#'   kernel-only). Used as `iweights` in balance regressions.
#' @param comsup_idx Logical: common support indicator (all TRUE for unweighted).
#' @param rbalance Integer: 0 = mean at cutoff, 1 = mean in sample.
#' @param N_G0 Pre-computed group 0 count (passed in to avoid recomputing).
#' @param N_G1 Pre-computed group 1 count.
#' @return A list with elements:
#'   - `table`: data frame with columns mean_G0, mean_G1, std_diff, p_value.
#'   - `N_G0`, `N_G1`: integer counts.
#'   - `avgdiff`: mean absolute standardized difference.
#'   - `Fstat`: global F-statistic.
#'   - `pval_global`: global p-value.
#' @noRd
compute_balance <- function(balance_vars, data, running, G, sgroup0,
                            within_bw, touse, weights, comsup_idx,
                            rbalance = 0, N_G0, N_G1) {
  active   <- touse & within_bw & comsup_idx
  x        <- data[[running]]
  w_active <- weights[active]

  n_bal <- length(balance_vars)
  mean_G0   <- numeric(n_bal)
  mean_G1   <- numeric(n_bal)
  std_diff  <- numeric(n_bal)
  pval      <- numeric(n_bal)

  for (j in seq_along(balance_vars)) {
    var <- data[[balance_vars[j]]][active]
    g   <- G[active]
    g0  <- sgroup0[active]
    xa  <- x[active]
    wa  <- w_active

    if (rbalance == 1) {
      # Mean in sample: regress var ~ sgroup0 + sgroup (no intercept) with iw
      df_nc <- data.frame(var = var, g0 = g0, g1 = g)
      fit_nc <- lm(var ~ 0 + g0 + g1, data = df_nc, weights = wa)
      mean_G0[j] <- coef(fit_nc)[["g0"]]
      mean_G1[j] <- coef(fit_nc)[["g1"]]

      # Difference and p-value: regress var ~ sgroup0 (with intercept)
      fit_diff <- lm(var ~ g0, weights = wa)
      ct <- summary(fit_diff)$coefficients
      diff_j <- ct["g0", "Estimate"]
      pval[j] <- ct["g0", "Pr(>|t|)"]
    } else {
      # Mean at cutoff: regress var ~ running_var within each group separately
      idx0 <- g == 0
      idx1 <- g == 1
      fit0 <- lm(var[idx0] ~ xa[idx0], weights = wa[idx0])
      fit1 <- lm(var[idx1] ~ xa[idx1], weights = wa[idx1])
      mean_G0[j] <- coef(fit0)[["(Intercept)"]]
      mean_G1[j] <- coef(fit1)[["(Intercept)"]]

      # Diff and p-value: fully interacted with sgroup0
      # var ~ running_var * sgroup0  → coef on sgroup0 is the intercept difference
      df_int <- data.frame(var = var, xa = xa, g0 = g0)
      fit_int <- lm(var ~ xa * g0, data = df_int, weights = wa)
      ct <- summary(fit_int)$coefficients
      diff_j <- ct["g0", "Estimate"]
      pval[j] <- ct["g0", "Pr(>|t|)"]
    }

    # Standardized mean difference
    sd_pooled <- sd(var)
    std_diff[j] <- if (sd_pooled > 0) diff_j / sd_pooled else 0
  }

  avgdiff <- mean(abs(std_diff))

  # Global F-statistic and p-value
  g_active  <- G[active]
  xa_active <- x[active]

  if (rbalance == 1) {
    bal_df <- as.data.frame(lapply(balance_vars, \(v) data[[v]][active]))
    names(bal_df) <- balance_vars
    bal_df$.G <- g_active
    fit_f <- lm(.G ~ ., data = bal_df, weights = w_active)
    sm    <- summary(fit_f)
    Fstat <- sm$fstatistic[["value"]]
    pval_global <- 1 - pf(Fstat, sm$fstatistic[["numdf"]], sm$fstatistic[["dendf"]])
  } else {
    # Frisch-Waugh style (Stata lines 982-1002):
    # 1. reg G ~ aux_j* + running (no constant) → residuals
    # 2. reg residuals ~ balance_vars → F-stat
    bal_mat <- as.data.frame(lapply(balance_vars, \(v) data[[v]][active]))
    names(bal_mat) <- balance_vars
    # build aux_j = var_j * running
    aux_mat <- as.data.frame(lapply(balance_vars, \(v) data[[v]][active] * xa_active))
    names(aux_mat) <- paste0("aux_", balance_vars)
    aux_mat$.x   <- xa_active
    aux_mat$.G   <- g_active
    fit_res <- lm(.G ~ 0 + ., data = aux_mat, weights = w_active)
    resids  <- residuals(fit_res)
    df_res  <- bal_mat
    df_res$.resid <- resids
    fit_f   <- lm(.resid ~ ., data = df_res, weights = w_active)
    sm      <- summary(fit_f)
    Fstat   <- sm$fstatistic[["value"]]
    pval_global <- 1 - pf(Fstat, sm$fstatistic[["numdf"]], sm$fstatistic[["dendf"]])
  }

  table_out <- data.frame(
    mean_G0  = mean_G0,
    mean_G1  = mean_G1,
    std_diff = std_diff,
    p_value  = pval,
    row.names = balance_vars
  )

  list(
    table       = table_out,
    N_G0        = N_G0,
    N_G1        = N_G1,
    avgdiff     = avgdiff,
    Fstat       = Fstat,
    pval_global = pval_global
  )
}
