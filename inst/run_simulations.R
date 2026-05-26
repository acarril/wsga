# inst/run_simulations.R
# Monte Carlo simulation study for wsga (RDD path)
#
# Four scenarios validating the core claims in the paper:
#   S1  Bias removal    -- IPW removes bias from imbalanced moderators
#   S2  Size control    -- WSGA maintains nominal 5% under the null; naive over-rejects
#   S3  Power vs K      -- WSGA retains power as K grows; full saturation collapses
#   S4  Unobservables   -- residual bias as the untestable IPW assumption fails
# Plus a bootstrap coverage sub-study nested inside S1.
#
# Usage (from rddsga-repo/):
#   Rscript inst/run_simulations.R           # full run (~10 min on 10 cores)
#   Rscript inst/run_simulations.R --quick   # smoke run (20 reps, small grid)
#
# Output: inst/simulation_results/*.csv
# -----------------------------------------------------------------------------

# ---- 0. Setup ----------------------------------------------------------------

args  <- commandArgs(trailingOnly = TRUE)
QUICK <- "--quick" %in% args

suppressPackageStartupMessages(devtools::load_all(quiet = TRUE))
library(parallel)

N_REPS   <- if (QUICK)  20L  else 1000L   # reps for S1-S4
N_REPS_B <- if (QUICK)  10L  else  500L   # reps for bootstrap coverage
BSREPS   <- if (QUICK)  20L  else  100L   # bootstrap replications per sim rep
NCORES   <- max(1L, detectCores() - 1L)
BW       <- 0.5                            # half-bandwidth (fixed throughout)
TRUE_DIFF <- 2                             # true differential RD effect for S1/S2/S4
S3_DIFF   <- 0.5                           # smaller effect for S3 power curves
OUTDIR   <- "inst/simulation_results"

dir.create(OUTDIR, showWarnings = FALSE)
RNGkind("L'Ecuyer-CMRG")
set.seed(42)
message(sprintf(
  "wsga Monte Carlo: %s | %d reps | %d bootstrap reps | %d cores",
  if (QUICK) "QUICK MODE" else "full run", N_REPS, BSREPS, NCORES))


# ---- 1. DGP ------------------------------------------------------------------
#
# Sharp RD: x ~ U(-1, 1), Z = 1[x >= 0], G ~ Bern(0.5)
# M = delta * G + eps_M,  eps_M ~ N(0, 1)   [moderator, imbalanced by delta]
#
# M enters as a treatment-response amplifier, NOT a level shifter.
# A level shifter (beta_M * M) would be absorbed by the local linear fit
# within each subgroup and would NOT bias the RD jump estimate.
# An amplifier ((b + beta_M*M)*Z*G) does bias the jump: the estimated
# within-group RD for G=1 picks up beta_M * E[M | G=1] = beta_M * delta,
# while the G=0 estimate is correct (E[M | G=0] = 0).
#
# DGP:
#   y = (2 + beta_M*M)*Z*(1-G) + (2 + beta_G + beta_M*M)*Z*G
#       + 0.3*x*(1-G) + 0.5*x*G + eps
#
# True per-group RD effects (at reference M=0):   b_g0 = 2, b_g1 = 2 + beta_G
# Naive estimate (no IPW):   b_g1_hat = 2 + beta_G + beta_M * delta
#                            diff_hat = beta_G + beta_M * delta
#                            bias     = beta_M * delta
# WSGA estimate:             diff_hat ~ beta_G  (bias removed by reweighting)
#
dgp_rdd <- function(n, delta, beta_M, beta_G = TRUE_DIFF) {
  x <- runif(n, -1, 1)
  G <- rbinom(n, 1, 0.5)
  M <- rnorm(n, mean = delta * G, sd = 1)
  Z <- as.integer(x >= 0)
  y <- (2 + beta_M*M)*Z*(1-G) + (2 + beta_G + beta_M*M)*Z*G +
       0.3*x*(1-G) + 0.5*x*G + rnorm(n, sd = 0.5)
  data.frame(y=y, x=x, G=G, M=M)
}

# K-moderator variant for Scenario 3 (power vs K).
#
# Each M_k ~ N(delta_k * G, 1) with delta_k = 0.4 (continuous, moderate
# imbalance).  M_k has NO direct effect on y -- it only predicts G.
# Continuous moderators avoid the extreme propensity scores produced by
# many imbalanced binary moderators (which create large IPW weight variance
# and would artificially suppress WSGA's power regardless of K).
#
# Three comparators:
#   naive      -- y ~ 1 | G, noipsw=TRUE  (no adjustment; unbiased here
#                  since M has no direct y-effect, but loses power with K
#                  because it can't tell wsga_rdd there are no confounders)
#   wsga       -- y ~ 1 | G, balance ~ M1+...+MK  (IPW; lean outcome)
#   full_sat   -- y ~ M1+...+MK | G, noipsw=TRUE  (linear covariates; adds
#                  2K df to the outcome; a proxy for a covariate-adjustment
#                  approach that is NOT the Calonico et al. cell-splitting
#                  which would require 2^K interaction terms and collapses
#                  even faster -- cell-splitting is infeasible to implement
#                  in wsga_rdd's formula interface)
#
# Expected story: wsga power stays close to naive (outcome stays lean);
# full_sat declines slowly from 2K df overhead; both far better than
# the theoretical cell-splitting curve (2^K parameters, not shown).
dgp_rdd_multi <- function(n, K, beta_G = S3_DIFF, delta_m = 0.4) {
  x   <- runif(n, -1, 1)
  G   <- rbinom(n, 1, 0.5)
  Z   <- as.integer(x >= 0)
  mdf <- as.data.frame(
    sapply(seq_len(max(K, 1L)), function(k)
      rnorm(n, mean = delta_m * G, sd = 1)))
  colnames(mdf) <- paste0("M", seq_len(max(K, 1L)))
  y <- 2*Z*(1-G) + (2 + beta_G)*Z*G +
       0.3*x*(1-G) + 0.5*x*G + rnorm(n, sd = 0.5)
  cbind(data.frame(y=y, x=x, G=G), mdf)
}


# ---- 2. Single-rep helpers ---------------------------------------------------

# Suppress wsga's per-replicate progress messages
hush <- function(expr) suppressWarnings(suppressMessages(expr))

# S1 / S2: analytical inference, naive vs WSGA
# Naive: no adjustment for M at all (y ~ 1 | G, noipsw=TRUE).
#   When delta > 0 and beta_M > 0, the differential subgroup estimate is
#   biased by approximately beta_M * delta (M is correlated with G and Y).
# WSGA: M absorbed via IPW propensity score, lean outcome regression.
one_rep <- function(n, delta, beta_M, beta_G = TRUE_DIFF) {
  d <- dgp_rdd(n, delta, beta_M, beta_G)

  f_naive <- tryCatch(hush(wsga_rdd(
    y ~ 1 | G, data = d, running = ~ x, bwidth = BW,
    noipsw = TRUE, bootstrap = FALSE)), error = function(e) NULL)

  f_wsga  <- tryCatch(hush(wsga_rdd(
    y ~ 1 | G, data = d, running = ~ x, bwidth = BW,
    balance = ~ M, noipsw = FALSE, bootstrap = FALSE)), error = function(e) NULL)

  if (is.null(f_naive) || is.null(f_wsga)) return(NULL)

  c(naive_est  = unname(coef(f_naive)["diff"]),
    naive_pval = f_naive$pval$diff,
    wsga_est   = unname(coef(f_wsga)["diff"]),
    wsga_pval  = f_wsga$pval$diff)
}

# S1 coverage: bootstrap inference, WSGA only
one_rep_boot <- function(n, delta, beta_M, beta_G = TRUE_DIFF, bsreps = BSREPS) {
  d <- dgp_rdd(n, delta, beta_M, beta_G)
  f <- tryCatch(hush(wsga_rdd(
    y ~ 1 | G, data = d, running = ~ x, bwidth = BW,
    balance = ~ M, noipsw = FALSE,
    bootstrap = TRUE, bsreps = bsreps, inference = "empirical")),
    error = function(e) NULL)
  if (is.null(f)) return(NULL)
  lb <- f$ci$diff["lb"];  ub <- f$ci$diff["ub"]
  c(est    = unname(coef(f)["diff"]),
    ci_lb  = unname(lb),
    ci_ub  = unname(ub),
    covers = as.integer(lb <= beta_G && ub >= beta_G))
}

# S3: power vs K moderators (uses S3_DIFF = 0.5, a smaller effect so the
# power curves separate -- at TRUE_DIFF = 2 all methods hit 100% power)
one_rep_multi <- function(n, K, beta_G = S3_DIFF) {
  d <- dgp_rdd_multi(n, K, beta_G)

  # Naive: no adjustment at all (y ~ 1 | G, noipsw=TRUE); unbiased here
  # since M has no direct y-effect, but included for baseline comparison
  f_naive <- tryCatch(hush(wsga_rdd(
    y ~ 1 | G, data = d, running = ~ x, bwidth = BW,
    noipsw = TRUE, bootstrap = FALSE)), error = function(e) NULL)

  # WSGA: lean outcome regression (y ~ 1 | G), all K in propensity score
  bal <- if (K > 0)
    as.formula(paste("~", paste0("M", seq_len(K), collapse = "+")))
  else NULL
  f_wsga <- tryCatch(hush(wsga_rdd(
    y ~ 1 | G, data = d, running = ~ x, bwidth = BW,
    balance = bal, noipsw = (K == 0), bootstrap = FALSE)), error = function(e) NULL)

  # Full saturation: all K as linear outcome covariates, no IPW (adds 2K df)
  rhs    <- if (K > 0) paste0("M", seq_len(K), collapse = "+") else "1"
  f_full <- tryCatch(hush(wsga_rdd(
    as.formula(paste("y ~", rhs, "| G")),
    data = d, running = ~ x, bwidth = BW,
    noipsw = TRUE, bootstrap = FALSE)), error = function(e) NULL)

  if (is.null(f_naive) || is.null(f_wsga) || is.null(f_full)) return(NULL)
  c(naive_pval = f_naive$pval$diff,
    wsga_pval  = f_wsga$pval$diff,
    full_pval  = f_full$pval$diff)
}

# S4: sensitivity to unobservable confounding.
# U = gamma * G + eps_U is unobservable (not in M), amplifies the treatment
# response via U*Z.  WSGA controls for M but not U, so residual bias ~ gamma.
# At gamma=0: WSGA is unbiased.  At gamma>0: bias grows with gamma.
one_rep_unobs <- function(n, gamma, delta = 1, beta_M = 0.3, beta_G = TRUE_DIFF) {
  d <- dgp_rdd(n, delta, beta_M, beta_G)
  U <- rnorm(nrow(d), mean = gamma * d$G, sd = 1)   # unobservable, correlated with G
  Z <- as.integer(d$x >= 0)
  d$y <- d$y + U * Z                                 # U amplifies treatment response
  f <- tryCatch(hush(wsga_rdd(
    y ~ 1 | G, data = d, running = ~ x, bwidth = BW,
    balance = ~ M, noipsw = FALSE, bootstrap = FALSE)), error = function(e) NULL)
  if (is.null(f)) return(NULL)
  c(wsga_est = unname(coef(f)["diff"]))
}


# ---- 3. Scenario runner ------------------------------------------------------

run_cells <- function(grid, rep_fn, n_reps, label) {
  message(sprintf("\n[%s] %d cells x %d reps", label, nrow(grid), n_reps))
  per_rep_seeds <- sample.int(1e6L, n_reps)
  do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
    p    <- grid[i, , drop = FALSE]
    args <- as.list(p)
    reps <- mclapply(per_rep_seeds, function(s) {
      set.seed(s)
      do.call(rep_fn, args)
    }, mc.cores = NCORES)
    ok    <- Filter(Negate(is.null), reps)
    draws <- as.data.frame(do.call(rbind, ok))
    message(sprintf("  cell %d/%d %s  n_ok=%d/%d",
                    i, nrow(grid),
                    paste(names(p), unlist(p), sep = "=", collapse = " "),
                    nrow(draws), n_reps))
    cbind(p[rep(1L, nrow(draws)), , drop = FALSE], draws,
          row.names = NULL)
  }))
}


# ---- 4. Scenario 1: Bias removal (analytical) --------------------------------

g1 <- if (QUICK) {
  expand.grid(n = 600L, delta = c(0, 1), beta_M = c(0, 0.3),
              stringsAsFactors = FALSE)
} else {
  expand.grid(n = c(300L, 600L, 1500L),
              delta = c(0, 0.5, 1, 2),
              beta_M = c(0, 0.3, 0.5),
              stringsAsFactors = FALSE)
}

reps1 <- run_cells(g1, one_rep, N_REPS, "S1 bias/RMSE")

s1 <- do.call(rbind, lapply(split(reps1, interaction(reps1$n, reps1$delta, reps1$beta_M)), function(x) {
  data.frame(
    n       = x$n[1], delta = x$delta[1], beta_M = x$beta_M[1],
    n_ok    = nrow(x),
    bias_naive  = mean(x$naive_est) - TRUE_DIFF,
    rmse_naive  = sqrt(mean((x$naive_est - TRUE_DIFF)^2)),
    rej_naive   = mean(x$naive_pval < 0.05, na.rm = TRUE),
    bias_wsga   = mean(x$wsga_est)  - TRUE_DIFF,
    rmse_wsga   = sqrt(mean((x$wsga_est  - TRUE_DIFF)^2)),
    rej_wsga    = mean(x$wsga_pval  < 0.05, na.rm = TRUE)
  )
}))
write.csv(s1, file.path(OUTDIR, "s1_bias_rmse.csv"), row.names = FALSE)
message("  -> saved s1_bias_rmse.csv")


# ---- 5. Scenario 1 coverage: bootstrap WSGA ----------------------------------

g1b <- if (QUICK) {
  expand.grid(n = 600L, delta = c(0, 1), beta_M = 0.3,
              stringsAsFactors = FALSE)
} else {
  expand.grid(n = c(300L, 600L, 1500L),
              delta = c(0, 1),
              beta_M = 0.3,
              stringsAsFactors = FALSE)
}

reps1b <- run_cells(g1b, one_rep_boot, N_REPS_B, "S1 coverage")

s1b <- do.call(rbind, lapply(split(reps1b, interaction(reps1b$n, reps1b$delta, reps1b$beta_M)), function(x) {
  data.frame(
    n = x$n[1], delta = x$delta[1], beta_M = x$beta_M[1],
    n_ok    = nrow(x),
    coverage = mean(x$covers, na.rm = TRUE),
    ci_width = mean(x$ci_ub - x$ci_lb, na.rm = TRUE)
  )
}))
write.csv(s1b, file.path(OUTDIR, "s1_coverage.csv"), row.names = FALSE)
message("  -> saved s1_coverage.csv")


# ---- 6. Scenario 2: Size control (null: beta_G = 0) -------------------------

g2 <- if (QUICK) {
  expand.grid(n = 600L, delta = c(0, 1), beta_M = c(0, 0.3),
              stringsAsFactors = FALSE)
} else {
  expand.grid(n = c(600L, 1500L),
              delta = c(0, 0.5, 1, 2),
              beta_M = c(0, 0.3),
              stringsAsFactors = FALSE)
}

one_rep_null <- function(n, delta, beta_M)
  one_rep(n, delta, beta_M, beta_G = 0)

reps2 <- run_cells(g2, one_rep_null, N_REPS, "S2 size")

s2 <- do.call(rbind, lapply(split(reps2, interaction(reps2$n, reps2$delta, reps2$beta_M)), function(x) {
  data.frame(
    n = x$n[1], delta = x$delta[1], beta_M = x$beta_M[1],
    n_ok       = nrow(x),
    size_naive = mean(x$naive_pval < 0.05, na.rm = TRUE),
    size_wsga  = mean(x$wsga_pval  < 0.05, na.rm = TRUE)
  )
}))
write.csv(s2, file.path(OUTDIR, "s2_size.csv"), row.names = FALSE)
message("  -> saved s2_size.csv")


# ---- 7. Scenario 3: Power vs number of moderators K -------------------------

g3 <- if (QUICK) {
  expand.grid(n = 600L, K = c(0L, 2L, 4L),
              stringsAsFactors = FALSE)
} else {
  expand.grid(n = c(600L, 1500L),
              K = c(0L, 1L, 2L, 4L, 6L, 8L),
              stringsAsFactors = FALSE)
}

reps3 <- run_cells(g3, one_rep_multi, N_REPS, "S3 power vs K")

s3 <- do.call(rbind, lapply(split(reps3, interaction(reps3$n, reps3$K)), function(x) {
  data.frame(
    n = x$n[1], K = x$K[1],
    n_ok        = nrow(x),
    power_naive = mean(x$naive_pval < 0.05, na.rm = TRUE),
    power_wsga  = mean(x$wsga_pval  < 0.05, na.rm = TRUE),
    power_full  = mean(x$full_pval  < 0.05, na.rm = TRUE)
  )
}))
write.csv(s3, file.path(OUTDIR, "s3_power_vs_K.csv"), row.names = FALSE)
message("  -> saved s3_power_vs_K.csv")


# ---- 8. Scenario 4: Unobservable confounding sensitivity --------------------

g4 <- if (QUICK) {
  expand.grid(n = 600L, gamma = c(0, 0.5, 1),
              stringsAsFactors = FALSE)
} else {
  expand.grid(n = c(600L, 1500L),
              gamma = c(0, 0.25, 0.5, 1, 2),
              stringsAsFactors = FALSE)
}

reps4 <- run_cells(g4, one_rep_unobs, N_REPS, "S4 unobservables")

s4 <- do.call(rbind, lapply(split(reps4, interaction(reps4$n, reps4$gamma)), function(x) {
  data.frame(
    n = x$n[1], gamma = x$gamma[1],
    n_ok      = nrow(x),
    bias_wsga = mean(x$wsga_est) - TRUE_DIFF
  )
}))
write.csv(s4, file.path(OUTDIR, "s4_unobservable.csv"), row.names = FALSE)
message("  -> saved s4_unobservable.csv")


# ---- 9. Print summary tables -------------------------------------------------

cat("\n\n==== S1: Bias and RMSE (true diff = 2) ====\n")
print(s1[order(s1$n, s1$delta, s1$beta_M), ], row.names = FALSE, digits = 3)

cat("\n==== S1 coverage: 95% CI coverage (nominal = 0.95) ====\n")
print(s1b[order(s1b$n, s1b$delta), ], row.names = FALSE, digits = 3)

cat("\n==== S2: Size control (true diff = 0; nominal alpha = 0.05) ====\n")
print(s2[order(s2$n, s2$delta, s2$beta_M), ], row.names = FALSE, digits = 3)

cat("\n==== S3: Power vs K moderators ====\n")
print(s3[order(s3$n, s3$K), ], row.names = FALSE, digits = 3)

cat("\n==== S4: Residual bias from unobservable confounding ====\n")
print(s4[order(s4$n, s4$gamma), ], row.names = FALSE, digits = 3)

message("\nDone. Results in ", OUTDIR)
