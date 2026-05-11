#' @importFrom stats lm glm binomial coef residuals fitted sd var pf pnorm pt
#'   qnorm qt quantile as.formula complete.cases formula nobs
#' @importFrom Formula Formula
NULL

#' Weighted Subgroup Analysis in Regression Discontinuity Designs
#'
#' Estimates per-subgroup RD treatment effects and their difference using
#' inverse probability weighting (IPW) to balance observed moderators across
#' subgroups at the cutoff. Supports sharp and fuzzy RD designs, multiple
#' kernel types, polynomial order selection, three IPW weighting modes, common
#' support trimming, and bootstrap inference.
#'
#' @param formula A two-part `Formula` of the form
#'   `outcome ~ covariates | subgroup`. The left of `|` contains covariates
#'   included in the outcome regression; the right of `|` is the binary (0/1)
#'   subgroup indicator. If there are no covariates, write `outcome ~ 1 | subgroup`.
#' @param data A data frame.
#' @param running A one-sided formula (e.g. `~ score`) or a character string
#'   naming the running variable column in `data`.
#' @param cutoff Numeric cutoff value. The running variable is internally
#'   centered: `x_c = running - cutoff`. Default `0`.
#' @param bwidth Positive numeric half-bandwidth. Required.
#' @param model Character: estimation model.
#'   - `"rf"` (default): reduced form (sharp RD or reduced form of fuzzy RD).
#'   - `"fs"`: first stage of fuzzy RD (outcome is the treatment indicator).
#'   - `"iv"`: fuzzy RD via 2SLS.
#' @param fuzzy A one-sided formula (`~ treat`) or character string naming the
#'   fuzzy treatment indicator. Required when `model` is `"fs"` or `"iv"`.
#' @param p Integer polynomial order for the local polynomial regression.
#'   Default `1` (local linear).
#' @param kernel Character kernel type: `"uniform"` (default), `"triangular"`,
#'   or `"epanechnikov"`.
#' @param balance A one-sided formula (e.g. `~ m1 + m2`) specifying balance
#'   variables for the propensity score model. Defaults to the covariates in
#'   `formula` if not supplied.
#' @param ps_model Propensity score model: `"logit"` (default) or `"probit"`.
#' @param noipsw Logical. If `TRUE`, skips IPW entirely (subgroup analysis
#'   without reweighting). Default `FALSE`.
#' @param m Integer weighting mode when `noipsw = FALSE`:
#'   - `2` (default): balance both groups toward the pooled distribution.
#'   - `1`: balance group 1 toward group 0 distribution.
#'   - `0`: balance group 0 toward group 1 distribution.
#' @param comsup Logical. Restrict to common propensity score support?
#'   Default `FALSE`.
#' @param show_balance Logical. Print balance tables to the console?
#'   Default `FALSE`.
#' @param rbalance Integer. How to compute conditional means for the balance
#'   table: `0` = mean at the cutoff (local linear extrapolation, default),
#'   `1` = mean in the full bandwidth sample.
#' @param bootstrap Logical. Run the bootstrap loop (refits the full pipeline
#'   per replicate)? Default `TRUE`.
#' @param bsreps Positive integer: number of bootstrap replications. Default `200`.
#' @param inference Character. How standard errors, CIs, and p-values are
#'   computed and reported. One of:
#'   - `"empirical"`: bootstrap SE; empirical (percentile) CIs at 2.5/97.5;
#'     `(1 + count)/(B + 1)` p-values. Requires `bootstrap = TRUE`. **Default
#'     when `bootstrap = TRUE`.**
#'   - `"normal"`: bootstrap SE; normal-approximation CIs `b ± z·SE` and
#'     p-values `2·(1 − Φ(|t|))`. Requires `bootstrap = TRUE`.
#'   - `"analytical"`: sandwich (or cluster-robust) SE; CIs and p-values from
#'     the t distribution with the regression residual df. Requires
#'     `bootstrap = FALSE`. **Default when `bootstrap = FALSE`.**
#' @param block_var Character or `NULL`. Column name for stratified bootstrap
#'   resampling. Default `NULL` (unstratified). Semantics adapt to
#'   `cluster_var`: when `cluster_var` is `NULL`, `block_var` defines strata
#'   of rows; when `cluster_var` is set, it defines strata of clusters and
#'   must be constant within each cluster (validated at runtime).
#' @param fixed_fs Logical. Fixed first-stage bootstrap for IV: bootstrap the
#'   reduced form and divide by the point-estimate first-stage coefficients.
#'   Default `FALSE`.
#' @param fixed_ps Logical. Diagnostic flag. When `TRUE`, the propensity score
#'   model is **not** refit inside the bootstrap loop; instead, each
#'   resampled row inherits the propensity score (and common-support
#'   membership) from its original-sample row. This isolates the
#'   variance contribution of the outcome model from that of the
#'   weight-estimation step. **It understates the SE relative to the
#'   paper's intended interpretation**, which is why the default is
#'   `FALSE`. Useful for sanity-checking and for fast iteration on the
#'   outcome specification when the PS fit is expensive. No effect when
#'   `noipsw = TRUE` or when there are no balance variables.
#' @param seed Integer or `NULL`. RNG seed for reproducibility. Default `NULL`.
#' @param vce Character: heteroskedasticity-consistent SE type passed to
#'   `sandwich::vcovHC`. Common choices: `"HC1"` (default), `"HC0"`, `"HC2"`,
#'   `"HC3"`. Use `"cluster"` together with `cluster_var` for cluster-robust SEs.
#' @param cluster_var Character or `NULL`. Column name defining clusters.
#'   Drives two things: when `bootstrap = TRUE`, switches the bootstrap from
#'   row-level to pairs cluster (whole clusters are resampled with
#'   replacement; cluster IDs are made unique per draw so unit fixed effects
#'   remain identified — the Cameron–Gelbach–Miller recipe). When
#'   `inference = "analytical"`, also requires `vce = "cluster"` to deliver
#'   the cluster-robust sandwich SE. With fewer than ~30 unique clusters, a
#'   one-time warning is emitted about likely SE understatement; see
#'   `fwildclusterboot::boottest` for a wild-cluster alternative.
#' @param weights Character or `NULL`. Column name of pre-existing observation
#'   weights. These are multiplied into the kernel weights and (if `!noipsw`)
#'   into the IPW weights.
#'
#' @return An S3 object of class `"wsga"` with the following elements:
#' \describe{
#'   \item{`coefficients`}{Named list `(g0, g1, diff)` of point estimates.}
#'   \item{`se`}{Named list of standard errors. Bootstrap-based for
#'     `inference` `"empirical"` or `"normal"`; sandwich/cluster-robust for
#'     `"analytical"`.}
#'   \item{`pval`}{Named list of two-sided p-values, computed according to
#'     `inference`.}
#'   \item{`ci`}{Named list of 95% confidence intervals, computed according to
#'     `inference`.}
#'   \item{`vcov`}{2×2 variance-covariance matrix for (g0, g1).}
#'   \item{`bootstrap`}{List with `draws` (B×3 matrix), `vcov`, `pval`, `ci`,
#'     `B_ok` (surviving reps), `failed` (dropped reps). `NULL` if
#'     `bootstrap = FALSE`.}
#'   \item{`inference`}{Character: the inference mode in effect.}
#'   \item{`balance`}{List with `unweighted` and `weighted` balance table data
#'     frames (each with global stats). `NULL` if no balance variables.}
#'   \item{`pscore`}{Numeric vector of propensity scores (NA outside active
#'     set). `NULL` if `noipsw = TRUE`.}
#'   \item{`ipw_weights`}{Numeric vector of IPW weights. `NULL` if
#'     `noipsw = TRUE`.}
#'   \item{`kernel_weights`}{Numeric vector of kernel weights.}
#'   \item{`model_fit`}{The underlying `lm` or `ivreg` object.}
#'   \item{`call`}{The matched call.}
#'   \item{`nobs`}{Named integer: `total`, `in_bw`, `G0`, `G1`.}
#'   \item{`bwidth`, `cutoff`, `kernel`, `p`, `m`, `model`}{Configuration.}
#' }
#'
#' @examples
#' data(rddsga_synth)
#' # Sharp RD, no IPW
#' fit <- wsga(y ~ 1 | sgroup, data = rddsga_synth,
#'             running = ~ x, bwidth = 0.5, noipsw = TRUE,
#'             bootstrap = FALSE)
#' print(fit)
#'
#' # Sharp RD with IPW balancing moderator m
#' fit2 <- wsga(y ~ m | sgroup, data = rddsga_synth,
#'              running = ~ x, bwidth = 0.5,
#'              bootstrap = FALSE)
#' summary(fit2)
#'
#' @export
wsga <- function(formula,
                   data,
                   running,
                   cutoff      = 0,
                   bwidth,
                   model       = c("rf", "fs", "iv"),
                   fuzzy       = NULL,
                   p           = 1L,
                   kernel      = c("uniform", "triangular", "epanechnikov"),
                   balance     = NULL,
                   ps_model    = c("logit", "probit"),
                   noipsw      = FALSE,
                   m           = 2L,
                   comsup      = FALSE,
                   show_balance = FALSE,
                   rbalance    = 0L,
                   bootstrap   = TRUE,
                   bsreps      = 200L,
                   inference   = NULL,
                   block_var   = NULL,
                   fixed_fs    = FALSE,
                   fixed_ps    = FALSE,
                   seed        = NULL,
                   vce         = "HC1",
                   cluster_var = NULL,
                   weights     = NULL) {

  cl <- match.call()

  # ── 1. Argument validation ──────────────────────────────────────────────────
  model    <- match.arg(model)
  kernel   <- match.arg(kernel)
  ps_model <- match.arg(ps_model)
  if (missing(bwidth) || !is.numeric(bwidth) || bwidth <= 0)
    stop("`bwidth` must be a positive number.")
  if (!m %in% c(0L, 1L, 2L))
    stop("`m` must be 0, 1, or 2.")
  if (model %in% c("fs", "iv") && is.null(fuzzy))
    stop("`fuzzy` must be specified when model = '", model, "'.")
  if (!is.data.frame(data))
    stop("`data` must be a data frame.")

  # Resolve `inference`. Three modes; default depends on `bootstrap`.
  inference <- if (is.null(inference)) {
    if (bootstrap) "empirical" else "analytical"
  } else {
    match.arg(inference, c("empirical", "normal", "analytical"))
  }
  if (bootstrap && inference == "analytical")
    stop("`inference = 'analytical'` requires `bootstrap = FALSE` (analytical SEs do not use the bootstrap loop).")
  if (!bootstrap && inference %in% c("empirical", "normal"))
    stop(sprintf("`inference = '%s'` requires `bootstrap = TRUE`.", inference))

  # ── 2. Parse formula ────────────────────────────────────────────────────────
  fml <- Formula::Formula(formula)
  if (length(fml)[2] != 2)
    stop("`formula` must be two-part: outcome ~ covariates | subgroup")

  outcome_name  <- all.vars(formula(fml, lhs = 1, rhs = 0))
  covariate_names <- all.vars(formula(fml, lhs = 0, rhs = 1))
  if (identical(covariate_names, "1")) covariate_names <- character(0)
  sgroup_name   <- all.vars(formula(fml, lhs = 0, rhs = 2))
  if (length(sgroup_name) != 1)
    stop("The right-hand side of `|` must name exactly one subgroup variable.")

  # Running variable
  running_name <- if (inherits(running, "formula")) {
    all.vars(running)[[1]]
  } else {
    as.character(running)
  }

  # Fuzzy treatment
  fuzzy_name <- if (!is.null(fuzzy)) {
    if (inherits(fuzzy, "formula")) all.vars(fuzzy)[[1]] else as.character(fuzzy)
  } else NULL

  # Balance variables
  if (!is.null(balance)) {
    balance_names <- all.vars(balance)
  } else {
    balance_names <- covariate_names
  }

  # User-supplied weights column
  obs_weights <- if (!is.null(weights)) {
    if (is.character(weights)) data[[weights]] else weights
  } else NULL

  # Cluster variable
  clust <- if (!is.null(cluster_var)) data[[cluster_var]] else NULL

  # ── 3. Extract vectors ───────────────────────────────────────────────────────
  y  <- data[[outcome_name]]
  x  <- data[[running_name]] - cutoff   # center at cutoff
  G  <- data[[sgroup_name]]
  D  <- if (!is.null(fuzzy_name)) data[[fuzzy_name]] else NULL

  if (!all(G %in% c(0, 1, NA)))
    stop("Subgroup variable must be binary (0/1).")

  sgroup0 <- as.integer(G == 0)
  Z       <- as.integer(x > 0)           # treatment assignment

  # ── 4. Sample masks ─────────────────────────────────────────────────────────
  touse    <- !is.na(y) & !is.na(x) & !is.na(G)
  within_bw <- abs(x) < bwidth

  # ── 5. Kernel weights ───────────────────────────────────────────────────────
  kwt <- kernel_weights(x, bwidth, kernel, obs_weights)

  # ── 6. Propensity score + IPW ────────────────────────────────────────────────
  if (!noipsw && length(balance_names) > 0) {
    bal_mat <- data[, balance_names, drop = FALSE]

    ps_res <- compute_pscore(G, bal_mat, within_bw, touse,
                             obs_weights, ps_model, comsup)
    pscore     <- ps_res$pscore
    comsup_idx <- ps_res$comsup_idx

    ipw_res  <- compute_ipw_weights(G, pscore, within_bw, touse,
                                    comsup_idx, m, obs_weights)
    ipw      <- ipw_res$ipw
    nweights <- ipw_res$nweights
    N_G0_w   <- ipw_res$N_G0
    N_G1_w   <- ipw_res$N_G1

    # Final regression weight: kernel * IPW
    final_wt <- ipw * kwt
    final_wt[is.na(final_wt)] <- 0
  } else {
    pscore     <- NULL
    comsup_idx <- touse & within_bw   # trivially all active
    ipw        <- NULL
    nweights   <- obs_weights
    final_wt   <- kwt
    N_G0_w <- sum(touse & within_bw & G == 0)
    N_G1_w <- sum(touse & within_bw & G == 1)
  }

  # Active observations for estimation (non-zero weight and within bw)
  active <- touse & within_bw & final_wt > 0

  # Unweighted counts (for balance table)
  N_G0_unw <- sum(touse & within_bw & G == 0)
  N_G1_unw <- sum(touse & within_bw & G == 1)

  # ── 7. Balance tables ────────────────────────────────────────────────────────
  balance_result <- NULL
  if (length(balance_names) > 0) {
    # Unweighted balance (kernel weights only)
    bal_unw <- compute_balance(
      balance_names, data, running_name, G, sgroup0,
      within_bw, touse,
      weights    = replace(kwt * coerce_obs_wt(obs_weights, nrow(data)), !touse | !within_bw, 0),
      comsup_idx = (touse & within_bw),
      rbalance   = rbalance,
      N_G0       = N_G0_unw,
      N_G1       = N_G1_unw
    )

    if (show_balance) {
      cat("\nUnweighted balance:\n")
      print(bal_unw$table)
      cat(sprintf("N (G=0): %d   N (G=1): %d   Mean |std diff|: %.4f   F: %.3f   p: %.4f\n\n",
                  bal_unw$N_G0, bal_unw$N_G1, bal_unw$avgdiff,
                  bal_unw$Fstat, bal_unw$pval_global))
    }

    bal_ipsw <- NULL
    if (!noipsw) {
      bal_ipsw <- compute_balance(
        balance_names, data, running_name, G, sgroup0,
        within_bw, touse,
        weights    = replace(nweights * kwt / kwt,  # nweights alone for iw
                             is.na(nweights) | !touse | !within_bw | !comsup_idx, 0),
        comsup_idx = comsup_idx,
        rbalance   = rbalance,
        N_G0       = N_G0_w,
        N_G1       = N_G1_w
      )

      if (show_balance) {
        cat("IPW-weighted balance:\n")
        print(bal_ipsw$table)
        cat(sprintf("N (G=0): %d   N (G=1): %d   Mean |std diff|: %.4f   F: %.3f   p: %.4f\n\n",
                    bal_ipsw$N_G0, bal_ipsw$N_G1, bal_ipsw$avgdiff,
                    bal_ipsw$Fstat, bal_ipsw$pval_global))
      }
    }

    balance_result <- list(unweighted = bal_unw, weighted = bal_ipsw)
  }

  # ── 8. Design matrix + estimation ────────────────────────────────────────────
  # For model = "fs", outcome is the fuzzy treatment variable
  outcome_for_model <- if (model == "fs") fuzzy_name else outcome_name

  # Temporarily replace the outcome column if needed
  data$.wsga_x <- x   # centered running variable

  dm <- build_design_matrix(
    data       = data,
    outcome    = outcome_for_model,
    running    = ".wsga_x",
    Z          = Z,
    G          = G,
    covariates = covariate_names,
    p          = p
  )

  fuzzy_data <- if (model == "iv") {
    list(D = D, G = G, Z = Z)
  } else NULL

  # For fixed_fs bootstrap, save first-stage coefficients now
  fs_coefs <- NULL
  if (model == "iv" && fixed_fs) {
    dm_fs <- build_design_matrix(
      data = data, outcome = fuzzy_name, running = ".wsga_x",
      Z = Z, G = G, covariates = covariate_names, p = p
    )
    fs_fit <- run_model(dm_fs, final_wt, active, model = "rf",
                        vce = vce, cluster_var = clust)
    fs_coefs <- c(
      coef(fs_fit$fit)[fs_fit$coef_g0_name],
      coef(fs_fit$fit)[fs_fit$coef_g1_name]
    )
  }

  model_result <- run_model(dm, final_wt, active, model,
                            fuzzy_data = fuzzy_data,
                            vce = vce, cluster_var = clust)

  est <- extract_estimates(model_result,
                           df_resid   = if (model != "iv") model_result$fit$df.residual else NULL,
                           use_normal = FALSE)

  # ── 9. Bootstrap ─────────────────────────────────────────────────────────────
  boot_result <- NULL
  if (bootstrap) {
    point_est <- c(est$b_g0, est$b_g1)

    # Closure capturing all configuration for one bootstrap replicate
    run_one_rep <- make_boot_rep(
      outcome_name     = outcome_for_model,
      running_name     = running_name,
      cutoff           = cutoff,
      bwidth           = bwidth,
      G_name           = sgroup_name,
      covariate_names  = covariate_names,
      balance_names    = balance_names,
      obs_wt_name      = if (!is.null(weights) && is.character(weights)) weights else NULL,
      obs_weights_vec  = obs_weights,
      kernel           = kernel,
      noipsw           = noipsw,
      ps_model         = ps_model,
      comsup           = comsup,
      m                = m,
      model            = model,
      fuzzy_name       = fuzzy_name,
      p                = p,
      vce              = vce,
      cluster_var_name = cluster_var,
      fixed_fs         = fixed_fs,
      fs_coefs         = fs_coefs,
      fixed_ps            = fixed_ps,
      original_pscore     = pscore,
      original_comsup_idx = comsup_idx
    )

    # Few-clusters advisory when clustering is active
    if (!is.null(cluster_var)) {
      n_clust <- length(unique(data[[cluster_var]][!is.na(data[[cluster_var]])]))
      if (n_clust < 30L) {
        warning(sprintf(
          "Cluster bootstrap with %d clusters. With fewer than ~30 clusters, pairs-cluster bootstrap can substantially understate standard errors. Consider also reporting analytical cluster-robust SEs (set `bootstrap = FALSE`, `vce = \"cluster\"`) and/or a wild cluster bootstrap-t (not implemented; see fwildclusterboot::boottest in R or boottest in Stata).",
          n_clust
        ))
      }
    }

    boot_result <- run_bootstrap(run_one_rep, data, bsreps, point_est,
                                 cluster_var = cluster_var,
                                 block_var = block_var, seed = seed)

    # SE and t-statistics always come from the bootstrap when it runs.
    V_boot  <- boot_result$vcov
    se_g0_b <- sqrt(V_boot[1, 1])
    se_g1_b <- sqrt(V_boot[2, 2])
    cov_b   <- V_boot[1, 2]
    se_diff_b <- sqrt(se_g0_b^2 + se_g1_b^2 - 2 * cov_b)

    t_g0_b   <- est$b_g0   / se_g0_b
    t_g1_b   <- est$b_g1   / se_g1_b
    t_diff_b <- est$b_diff / se_diff_b

    est$se_g0 <- se_g0_b; est$se_g1 <- se_g1_b; est$se_diff <- se_diff_b
    est$t_g0  <- t_g0_b;  est$t_g1  <- t_g1_b;  est$t_diff  <- t_diff_b
    est$vcov_2x2 <- V_boot

    # CIs and p-values: branch on inference mode.
    if (inference == "empirical") {
      est$p_g0   <- boot_result$pval[["g0"]]
      est$p_g1   <- boot_result$pval[["g1"]]
      est$p_diff <- boot_result$pval[["diff"]]
      est$ci_g0   <- c(lb = boot_result$ci["lb", "g0"],   ub = boot_result$ci["ub", "g0"])
      est$ci_g1   <- c(lb = boot_result$ci["lb", "g1"],   ub = boot_result$ci["ub", "g1"])
      est$ci_diff <- c(lb = boot_result$ci["lb", "diff"], ub = boot_result$ci["ub", "diff"])
    } else {  # inference == "normal"
      z <- qnorm(0.975)
      est$p_g0   <- 2 * (1 - pnorm(abs(t_g0_b)))
      est$p_g1   <- 2 * (1 - pnorm(abs(t_g1_b)))
      est$p_diff <- 2 * (1 - pnorm(abs(t_diff_b)))
      est$ci_g0   <- c(lb = est$b_g0   - z * se_g0_b,   ub = est$b_g0   + z * se_g0_b)
      est$ci_g1   <- c(lb = est$b_g1   - z * se_g1_b,   ub = est$b_g1   + z * se_g1_b)
      est$ci_diff <- c(lb = est$b_diff - z * se_diff_b, ub = est$b_diff + z * se_diff_b)
    }
  }

  # ── 10. Assemble output ──────────────────────────────────────────────────────
  data$.wsga_x <- NULL  # clean up

  nobs <- c(
    total = sum(touse),
    in_bw = sum(touse & within_bw),
    G0    = N_G0_w,
    G1    = N_G1_w
  )

  structure(
    list(
      coefficients  = list(g0 = est$b_g0,   g1 = est$b_g1,   diff = est$b_diff),
      se            = list(g0 = est$se_g0,  g1 = est$se_g1,  diff = est$se_diff),
      tstat         = list(g0 = est$t_g0,   g1 = est$t_g1,   diff = est$t_diff),
      pval          = list(g0 = est$p_g0,   g1 = est$p_g1,   diff = est$p_diff),
      ci            = list(g0 = est$ci_g0,  g1 = est$ci_g1,  diff = est$ci_diff),
      vcov          = est$vcov_2x2,
      bootstrap     = boot_result,
      balance       = balance_result,
      pscore        = pscore,
      ipw_weights   = ipw,
      kernel_weights = kwt,
      model_fit     = model_result$fit,
      call          = cl,
      nobs          = nobs,
      outcome       = outcome_name,
      bwidth        = bwidth,
      cutoff        = cutoff,
      kernel        = kernel,
      p             = p,
      m             = m,
      model         = model,
      noipsw        = noipsw,
      inference     = inference,
      cluster_var_name = if (is.null(cluster_var)) NULL else cluster_var,
      fixed_ps      = fixed_ps
    ),
    class = "wsga"
  )
}


# ── Helper: build bootstrap-replicate closure ──────────────────────────────────
make_boot_rep <- function(outcome_name, running_name, cutoff, bwidth,
                          G_name, covariate_names, balance_names,
                          obs_wt_name, obs_weights_vec,
                          kernel, noipsw, ps_model, comsup, m,
                          model, fuzzy_name, p, vce, cluster_var_name,
                          fixed_fs, fs_coefs,
                          fixed_ps = FALSE,
                          original_pscore = NULL,
                          original_comsup_idx = NULL) {
  function(data_b, orig_idx) {
    x_b  <- data_b[[running_name]] - cutoff
    G_b  <- data_b[[G_name]]
    y_b  <- data_b[[outcome_name]]
    Z_b  <- as.integer(x_b > 0)

    touse_b    <- !is.na(y_b) & !is.na(x_b) & !is.na(G_b)
    within_b   <- abs(x_b) < bwidth
    obs_wt_b   <- if (!is.null(obs_wt_name)) data_b[[obs_wt_name]] else obs_weights_vec

    kwt_b <- kernel_weights(x_b, bwidth, kernel, obs_wt_b)

    if (!noipsw && length(balance_names) > 0) {
      if (fixed_ps) {
        # Look up original-sample propensity score and common-support mask
        # via the orig_idx mapping (each resampled row carries its original
        # row index).
        ps_b <- list(
          pscore     = original_pscore[orig_idx],
          comsup_idx = original_comsup_idx[orig_idx]
        )
      } else {
        bal_b <- data_b[, balance_names, drop = FALSE]
        ps_b  <- compute_pscore(G_b, bal_b, within_b, touse_b,
                                obs_wt_b, ps_model, comsup)
      }
      ipw_b   <- compute_ipw_weights(G_b, ps_b$pscore, within_b, touse_b,
                                     ps_b$comsup_idx, m, obs_wt_b)
      fw_b    <- ipw_b$ipw * kwt_b
      fw_b[is.na(fw_b)] <- 0
    } else {
      fw_b <- kwt_b
    }

    active_b <- touse_b & within_b & fw_b > 0
    data_b$.wsga_x <- x_b

    if (model == "iv" && fixed_fs) {
      # Bootstrap reduced form and divide by saved FS coefficients
      dm_b <- build_design_matrix(data_b, fuzzy_name, ".wsga_x",
                                  Z_b, G_b, covariate_names, p)
      res_b <- run_model(dm_b, fw_b, active_b, "rf", vce = vce)
      bc    <- coef(res_b$fit)
      return(c(
        bc[res_b$coef_g0_name] / fs_coefs[1],
        bc[res_b$coef_g1_name] / fs_coefs[2]
      ))
    }

    outcome_b  <- if (model == "fs") fuzzy_name else outcome_name
    fuzzy_b    <- if (model == "iv") list(D = data_b[[fuzzy_name]], G = G_b, Z = Z_b) else NULL
    dm_b       <- build_design_matrix(data_b, outcome_b, ".wsga_x",
                                      Z_b, G_b, covariate_names, p)
    clust_b    <- if (!is.null(cluster_var_name)) data_b[[cluster_var_name]] else NULL
    res_b      <- run_model(dm_b, fw_b, active_b, model,
                            fuzzy_data = fuzzy_b, vce = vce, cluster_var = clust_b)
    bc         <- coef(res_b$fit)
    c(bc[res_b$coef_g0_name], bc[res_b$coef_g1_name])
  }
}


# ── Helper: coerce obs_weights to a length-n vector ───────────────────────────
coerce_obs_wt <- function(obs_weights, n) {
  if (is.null(obs_weights)) rep(1, n) else obs_weights
}


# ── S3 methods ────────────────────────────────────────────────────────────────

#' @export
print.wsga <- function(x, ...) {
  use_boot <- !is.null(x$bootstrap)
  ci_label <- switch(x$inference,
                     empirical  = "[95% CI (emp.)]",
                     normal     = "[95% CI (norm., boot)]",
                     analytical = "[95% CI (norm.)]")
  z_label  <- if (x$inference == "analytical") "t" else "z"
  p_label  <- if (x$inference == "analytical") "P>|t|" else "P>|z|"

  cat(sprintf("\nRD subgroup analysis  |  model: %s  |  bwidth: %g  |  kernel: %s\n\n",
              x$model, x$bwidth, x$kernel))

  hdr <- sprintf("%-12s | %9s  %9s  %6s  %7s  %s\n",
                 x$outcome, "Coef.", "Std. Err.", z_label, p_label, ci_label)
  sep <- paste0(strrep("-", 13), "+", strrep("-", 64), "\n")

  cat(sep)
  cat(hdr)
  cat(sep)
  cat("Subgroup     |\n")

  fmt_row <- function(label, b, se, t, p, ci) {
    p_str <- if (p < 0.001) "<.001" else sprintf("%.3f", p)
    sprintf("  %-10s | %9.4f  %9.4f  %6.2f  %7s  %9.4f  %9.4f\n",
            label, b, se, t, p_str, ci[["lb"]], ci[["ub"]])
  }

  cat(fmt_row("0", x$coefficients$g0, x$se$g0, x$tstat$g0, x$pval$g0, x$ci$g0))
  cat(fmt_row("1", x$coefficients$g1, x$se$g1, x$tstat$g1, x$pval$g1, x$ci$g1))
  cat(sep)
  cat(fmt_row("Difference", x$coefficients$diff, x$se$diff,
              x$tstat$diff, x$pval$diff, x$ci$diff))
  cat(sep)

  if (use_boot) {
    ps_tag <- if (isTRUE(x$fixed_ps)) " [PS fixed]" else ""
    if (!is.null(x$bootstrap$N_clusters)) {
      total_reps <- x$bootstrap$B_ok + x$bootstrap$failed
      cat(sprintf("Cluster bootstrap: %d/%d reps, clustered on '%s' (%d clusters)%s\n",
                  x$bootstrap$B_ok, total_reps,
                  x$cluster_var_name, x$bootstrap$N_clusters, ps_tag))
    } else {
      cat(sprintf("Bootstrap replications: %d%s\n", x$bootstrap$B_ok, ps_tag))
    }
  }

  cat(sprintf("N (G=0): %d   N (G=1): %d\n", x$nobs[["G0"]], x$nobs[["G1"]]))
  invisible(x)
}


#' @export
summary.wsga <- function(object, ...) {
  print(object, ...)

  if (!is.null(object$balance)) {
    cat("\n--- Unweighted balance ---\n")
    print(round(object$balance$unweighted$table, 4))
    cat(sprintf("F = %.3f  (p = %.4f)  Mean |std diff| = %.4f\n",
                object$balance$unweighted$Fstat,
                object$balance$unweighted$pval_global,
                object$balance$unweighted$avgdiff))

    if (!is.null(object$balance$weighted)) {
      cat("\n--- IPW-weighted balance ---\n")
      print(round(object$balance$weighted$table, 4))
      cat(sprintf("F = %.3f  (p = %.4f)  Mean |std diff| = %.4f\n",
                  object$balance$weighted$Fstat,
                  object$balance$weighted$pval_global,
                  object$balance$weighted$avgdiff))
    }
  }

  invisible(object)
}


#' @export
coef.wsga <- function(object, ...) {
  c(g0 = object$coefficients$g0,
    g1 = object$coefficients$g1,
    diff = object$coefficients$diff)
}


#' @export
vcov.wsga <- function(object, ...) object$vcov


#' @export
confint.wsga <- function(object, parm = NULL, level = 0.95, ...) {
  ci <- rbind(
    g0   = object$ci$g0,
    g1   = object$ci$g1,
    diff = object$ci$diff
  )
  colnames(ci) <- c("lb", "ub")
  if (!is.null(parm)) ci <- ci[parm, , drop = FALSE]
  ci
}


#' @export
nobs.wsga <- function(object, ...) object$nobs[["in_bw"]]
