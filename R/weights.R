#' Compute propensity scores for subgroup membership
#'
#' Fits a binary response model (logit or probit) for the subgroup indicator
#' within the bandwidth, optionally trimmed to common support.
#'
#' @param G Integer 0/1 vector: subgroup indicator.
#' @param balance_mat Numeric matrix/data frame of balance covariates.
#' @param within_bw Logical vector: observations within bandwidth.
#' @param touse Logical vector: non-missing observations to include.
#' @param obs_weights Optional numeric vector of observation weights.
#' @param ps_model Character: `"logit"` (default) or `"probit"`.
#' @param comsup Logical: restrict to common support? Default `FALSE`.
#' @return A list with elements:
#'   - `pscore`: numeric vector (length = `length(G)`), NA outside `within_bw & touse`.
#'   - `comsup_idx`: logical vector marking observations in common support.
#' @noRd
compute_pscore <- function(G, balance_mat, within_bw, touse,
                           obs_weights = NULL, ps_model = c("logit", "probit"),
                           comsup = FALSE) {
  ps_model <- match.arg(ps_model)
  link <- if (ps_model == "logit") "logit" else "probit"

  active <- which(touse & within_bw & !is.na(G))
  df_fit <- as.data.frame(balance_mat[active, , drop = FALSE])
  df_fit$.G <- G[active]

  w_fit <- if (!is.null(obs_weights)) obs_weights[active] else NULL

  fit <- glm(.G ~ ., data = df_fit, family = binomial(link = link),
             weights = w_fit)

  pscore <- rep(NA_real_, length(G))
  pscore[active] <- fitted(fit)

  # Common support: [min(ps | G==1), max(ps | G==1)]
  ps_active <- pscore[active]
  g_active  <- G[active]
  if (comsup) {
    ps_g1    <- ps_active[g_active == 1]
    cs_min   <- min(ps_g1)
    cs_max   <- max(ps_g1)
    comsup_idx <- rep(FALSE, length(G))
    comsup_idx[active] <- ps_active >= cs_min & ps_active <= cs_max
  } else {
    comsup_idx <- rep(FALSE, length(G))
    comsup_idx[active] <- TRUE
  }

  list(pscore = pscore, comsup_idx = comsup_idx, fit = fit)
}


#' Compute IPW weights for subgroup balance
#'
#' Three modes following the Stata implementation:
#' - `m = 2`: balance both groups toward the pooled distribution.
#' - `m = 1`: balance group 1 toward group 0 (ATT for G=0).
#' - `m = 0`: balance group 0 toward group 1 (ATT for G=1).
#'
#' @param G Integer 0/1 vector: subgroup indicator.
#' @param pscore Numeric vector of propensity scores (NA outside active set).
#' @param within_bw Logical vector: within bandwidth.
#' @param touse Logical vector: non-missing observations.
#' @param comsup_idx Logical vector: in common support.
#' @param m Integer in `{0, 1, 2}`: weighting mode.
#' @param obs_weights Optional pre-existing observation weights.
#' @return A list:
#'   - `ipw`: numeric vector of IPW weights (NA where not applicable).
#'   - `N_G0`, `N_G1`: observation counts per group (given active + common support).
#' @noRd
compute_ipw_weights <- function(G, pscore, within_bw, touse, comsup_idx, m = 2,
                                obs_weights = NULL) {
  # Counts conditional on m (mirrors Stata lines 828-845)
  if (m == 2) {
    N_G0 <- sum(touse & within_bw & comsup_idx & G == 0 & !is.na(pscore))
    N_G1 <- sum(touse & within_bw & comsup_idx & G == 1 & !is.na(pscore))
  } else if (m == 1) {
    N_G0 <- sum(touse & within_bw & comsup_idx & G == 0)
    N_G1 <- sum(touse & within_bw & comsup_idx & G == 1 & !is.na(pscore))
  } else { # m == 0
    N_G0 <- sum(touse & within_bw & comsup_idx & G == 0 & !is.na(pscore))
    N_G1 <- sum(touse & within_bw & comsup_idx & G == 1)
  }

  p <- pscore
  ipw <- rep(NA_real_, length(G))

  if (m == 2) {
    active0 <- touse & within_bw & comsup_idx & G == 0 & !is.na(p)
    active1 <- touse & within_bw & comsup_idx & G == 1 & !is.na(p)
    ipw[active1] <- (N_G1 / (N_G1 + N_G0)) / p[active1]
    ipw[active0] <- (N_G0 / (N_G1 + N_G0)) / (1 - p[active0])
  } else if (m == 1) {
    active0 <- touse & within_bw & comsup_idx & G == 0
    active1 <- touse & within_bw & comsup_idx & G == 1 & !is.na(p)
    ipw[active1] <- (1 - p[active1]) / p[active1]
    ipw[active0] <- 1
  } else { # m == 0
    active0 <- touse & within_bw & comsup_idx & G == 0 & !is.na(p)
    active1 <- touse & within_bw & comsup_idx & G == 1
    ipw[active0] <- p[active0] / (1 - p[active0])
    ipw[active1] <- 1
  }

  # nweights = ipw * obs_weights (for balance regressions)
  nweights <- ipw
  if (!is.null(obs_weights)) nweights <- nweights * obs_weights

  list(ipw = ipw, nweights = nweights, N_G0 = N_G0, N_G1 = N_G1)
}
