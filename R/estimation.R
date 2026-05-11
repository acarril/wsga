#' Build design matrix for subgroup RD regression
#'
#' Replicates the Stata specification:
#'   y ~ i.sgroup#1._cutoff + i.sgroup + i.sgroup#(covariates x x*Z [poly])
#' Column names use underscores (e.g. "G0_Z") so they are safe in formulas.
#'
#' @param data Data frame with all variables.
#' @param outcome Character: outcome variable name.
#' @param running Character: running variable name (already centered at cutoff).
#' @param Z Integer 0/1 vector: cutoff indicator.
#' @param G Integer 0/1 vector: subgroup indicator.
#' @param covariates Character vector of covariate names.
#' @param p Integer: polynomial order.
#' @return A list: `X` (design matrix), `y` (outcome), `coef_g0_name`,
#'   `coef_g1_name`.
#' @noRd
build_rdd_design_matrix <- function(data, outcome, running, Z, G, covariates, p) {
  n <- nrow(data)
  x <- data[[running]]
  y <- data[[outcome]]

  G0 <- as.integer(G == 0)
  G1 <- as.integer(G == 1)

  G0Z <- G0 * Z
  G1Z <- G1 * Z

  poly_cols  <- list()
  poly_names <- character(0)
  for (k in seq_len(p)) {
    xk <- x^k
    poly_cols  <- c(poly_cols,
                    list(G0 * xk, G1 * xk, G0 * xk * Z, G1 * xk * Z))
    poly_names <- c(poly_names,
                    sprintf("G0_x%d", k), sprintf("G1_x%d", k),
                    sprintf("G0_x%dZ", k), sprintf("G1_x%dZ", k))
  }

  cov_cols  <- list()
  cov_names <- character(0)
  for (v in covariates) {
    cv <- data[[v]]
    # sanitize covariate name for formula safety
    vn <- make.names(v)
    cov_cols  <- c(cov_cols, list(G0 * cv, G1 * cv))
    cov_names <- c(cov_names, paste0("G0_", vn), paste0("G1_", vn))
  }

  X <- cbind(
    G0_Z = G0Z,
    G1_Z = G1Z,
    G0   = G0,
    G1   = G1,
    do.call(cbind, poly_cols),
    do.call(cbind, if (length(cov_cols) > 0) cov_cols else list(matrix(nrow = n, ncol = 0)))
  )
  colnames(X) <- c("G0_Z", "G1_Z", "G0", "G1", poly_names, cov_names)

  list(X = X, y = y, coef_g0_name = "G0_Z", coef_g1_name = "G1_Z")
}


#' Run the main RD regression (sharp: lm, fuzzy: ivreg)
#'
#' @param dm List from `build_design_matrix()`.
#' @param weights Numeric vector of observation weights.
#' @param active Logical index of rows to include.
#' @param model Character: `"rf"`, `"fs"`, or `"iv"`.
#' @param fuzzy_data Named list: `D`, `G`, `Z`. Only used for `model = "iv"`.
#' @param vce Character: HC type or `"cluster"`.
#' @param cluster_var Optional vector for cluster SEs.
#' @return List: `fit`, `vcov_fit`, `coef_g0_name`, `coef_g1_name`,
#'   `coef_names`.
#' @noRd
run_model <- function(dm, weights, active, model = c("rf", "fs", "iv"),
                      fuzzy_data = NULL, vce = "HC1", cluster_var = NULL) {
  model <- match.arg(model)

  X <- dm$X[active, , drop = FALSE]
  y <- dm$y[active]
  w <- weights[active]

  if (model %in% c("rf", "fs")) {
    fit <- lm(y ~ 0 + X, weights = w)
    # Strip the "X" prefix lm adds to column names
    names(fit$coefficients) <- colnames(X)

    if (vce == "cluster" && !is.null(cluster_var)) {
      vcov_fit <- sandwich::vcovCL(fit, cluster = cluster_var[active])
    } else {
      vcov_fit <- sandwich::vcovHC(fit, type = vce)
    }
    rownames(vcov_fit) <- colnames(vcov_fit) <- colnames(X)

    return(list(fit = fit, vcov_fit = vcov_fit,
                coef_g0_name = dm$coef_g0_name,
                coef_g1_name = dm$coef_g1_name,
                coef_names   = colnames(X)))
  }

  # IV (2SLS) ─────────────────────────────────────────────────────────────────
  D <- fuzzy_data$D
  G <- fuzzy_data$G
  Z_full <- fuzzy_data$Z

  G0 <- as.integer(G == 0)[active]
  G1 <- as.integer(G == 1)[active]
  D_a <- D[active]
  Z_a <- Z_full[active]

  G0D <- G0 * D_a   # endogenous
  G1D <- G1 * D_a   # endogenous
  G0Z <- G0 * Z_a   # instrument
  G1Z <- G1 * Z_a   # instrument

  exo_cols <- setdiff(colnames(X), c("G0_Z", "G1_Z"))
  X_exo    <- X[, exo_cols, drop = FALSE]

  # Build a plain data frame (all columns already have formula-safe names)
  df_iv <- as.data.frame(cbind(
    y = y, G0D = G0D, G1D = G1D, X_exo, G0Z = G0Z, G1Z = G1Z
  ))

  rhs_exo   <- paste(exo_cols, collapse = " + ")
  fmla_str  <- if (nchar(rhs_exo) > 0) {
    sprintf("y ~ 0 + G0D + G1D + %s | 0 + G0Z + G1Z + %s", rhs_exo, rhs_exo)
  } else {
    "y ~ 0 + G0D + G1D | 0 + G0Z + G1Z"
  }
  fmla <- as.formula(fmla_str)

  fit <- ivreg::ivreg(fmla, data = df_iv, weights = w)

  # Rename endogenous coefficients to match G0_Z / G1_Z convention
  all_names <- names(fit$coefficients)
  all_names[all_names == "G0D"] <- "G0_Z"
  all_names[all_names == "G1D"] <- "G1_Z"
  names(fit$coefficients) <- all_names

  if (vce == "cluster" && !is.null(cluster_var)) {
    vcov_fit <- sandwich::vcovCL(fit, cluster = cluster_var[active])
  } else {
    vcov_fit <- sandwich::vcovHC(fit, type = vce)
  }
  dimnames(vcov_fit) <- list(all_names, all_names)

  list(fit = fit, vcov_fit = vcov_fit,
       coef_g0_name = "G0_Z",
       coef_g1_name = "G1_Z",
       coef_names   = all_names)
}


#' Extract per-subgroup and difference estimates
#'
#' @param model_result List from `run_model()`.
#' @param df_resid Residual degrees of freedom (NULL = normal approximation).
#' @param use_normal Logical: force normal distribution?
#' @return Named list of scalars and vectors.
#' @noRd
extract_estimates <- function(model_result, df_resid = NULL, use_normal = FALSE) {
  b <- model_result$fit$coefficients
  names(b) <- model_result$coef_names
  V <- model_result$vcov_fit
  dimnames(V) <- list(model_result$coef_names, model_result$coef_names)

  g0 <- model_result$coef_g0_name
  g1 <- model_result$coef_g1_name

  b0  <- b[[g0]];   b1  <- b[[g1]]
  se0 <- sqrt(V[g0, g0]); se1 <- sqrt(V[g1, g1])
  cov01 <- V[g0, g1]

  b_diff  <- b1 - b0
  se_diff <- sqrt(se0^2 + se1^2 - 2 * cov01)

  t0 <- b0 / se0;  t1 <- b1 / se1;  t_diff <- b_diff / se_diff

  if (use_normal || is.null(df_resid)) {
    p0     <- 2 * (1 - pnorm(abs(t0)))
    p1     <- 2 * (1 - pnorm(abs(t1)))
    p_diff <- 2 * (1 - pnorm(abs(t_diff)))
    qt_fn  <- \(alpha) qnorm(1 - alpha / 2)
  } else {
    p0     <- 2 * pt(-abs(t0),     df = df_resid)
    p1     <- 2 * pt(-abs(t1),     df = df_resid)
    p_diff <- 2 * pt(-abs(t_diff), df = df_resid)
    qt_fn  <- \(alpha) qt(1 - alpha / 2, df = df_resid)
  }

  q <- qt_fn(0.05)

  list(
    b_g0    = b0,  b_g1    = b1,  b_diff    = b_diff,
    se_g0   = se0, se_g1   = se1, se_diff   = se_diff,
    t_g0    = t0,  t_g1    = t1,  t_diff    = t_diff,
    p_g0    = p0,  p_g1    = p1,  p_diff    = p_diff,
    ci_g0   = c(lb = b0 - q * se0,        ub = b0 + q * se0),
    ci_g1   = c(lb = b1 - q * se1,        ub = b1 + q * se1),
    ci_diff = c(lb = b_diff - q * se_diff, ub = b_diff + q * se_diff),
    vcov_2x2 = V[c(g0, g1), c(g0, g1)]
  )
}


#' Build design matrix for subgroup DiD regression (long-form TWFE).
#'
#' Encodes:
#'   Y_it = sum_i alpha_i I(unit=i)
#'        + gamma_0 * post_t * (1 - G_i) + gamma_1 * post_t * G_i
#'        + delta_0 * D_i * post_t * (1 - G_i)
#'        + delta_1 * D_i * post_t * G_i
#'        + sum_k beta_{0k} X_{kt} * (1 - G_i)
#'        + sum_k beta_{1k} X_{kt} * G_i
#'        + eps_it
#'
#' Coefficients of interest are `delta_0` and `delta_1`, exposed under the
#' same `G0_Z` / `G1_Z` names used in the RDD design matrix so that
#' downstream `run_model()` / `extract_estimates()` need no changes.
#'
#' @param data Data frame in long format (one row per unit-time).
#' @param outcome Character: outcome variable name.
#' @param unit Character: unit identifier column. Used to construct unit FE.
#' @param time Character: time variable column.
#' @param post_value Scalar: value of `time` denoting the post period.
#' @param treat Character: treatment-group indicator (constant within unit;
#'   0/1).
#' @param G Integer 0/1 vector: subgroup indicator (constant within unit).
#' @param covariates Character vector of covariate names. Time-invariant
#'   covariates will be perfectly collinear with the unit FE and dropped by
#'   `lm()` (a NA coefficient); harmless for the coefficients of interest.
#' @return A list: `X`, `y`, `coef_g0_name`, `coef_g1_name`.
#' @noRd
build_did_design_matrix <- function(data, outcome, unit, time, post_value,
                                    treat, G, covariates) {
  n <- nrow(data)
  y <- data[[outcome]]
  G0 <- as.integer(G == 0)
  G1 <- as.integer(G == 1)
  D  <- data[[treat]]
  post <- as.integer(data[[time]] == post_value)

  # Treatment-effect coefficients (delta_0 = G0_Z, delta_1 = G1_Z)
  G0_DP <- G0 * D * post
  G1_DP <- G1 * D * post

  # Time-by-subgroup interactions (gamma_0, gamma_1 in 2-period DiD)
  G0_post <- G0 * post
  G1_post <- G1 * post

  # Unit fixed effects: one column per unit.  `as.factor` orders the
  # levels lexicographically; column names are `unit_<level>`.
  unit_factor <- as.factor(data[[unit]])
  unit_fe_mat <- stats::model.matrix(~ unit_factor - 1)
  colnames(unit_fe_mat) <- paste0("unit_", levels(unit_factor))

  # Covariates × G (same convention as RDD)
  cov_cols  <- list()
  cov_names <- character(0)
  for (v in covariates) {
    cv <- data[[v]]
    vn <- make.names(v)
    cov_cols  <- c(cov_cols, list(G0 * cv, G1 * cv))
    cov_names <- c(cov_names, paste0("G0_", vn), paste0("G1_", vn))
  }

  X <- cbind(
    G0_Z    = G0_DP,
    G1_Z    = G1_DP,
    G0_post = G0_post,
    G1_post = G1_post,
    do.call(cbind, if (length(cov_cols) > 0) cov_cols
                   else list(matrix(nrow = n, ncol = 0))),
    unit_fe_mat
  )
  colnames(X) <- c("G0_Z", "G1_Z", "G0_post", "G1_post",
                   cov_names, colnames(unit_fe_mat))

  list(X = X, y = y, coef_g0_name = "G0_Z", coef_g1_name = "G1_Z")
}
