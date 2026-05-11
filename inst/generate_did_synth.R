# Generates the wsga_did_synth balanced 2-period panel dataset and saves it to
# data/wsga_did_synth.rda and stata/wsga_did_synth.dta.
# Run from the package root: Rscript inst/generate_did_synth.R

set.seed(7)

n_units <- 500L

# Unit-level variables
unit_id <- seq_len(n_units)

# Binary subgroup indicator (roughly balanced), constant within unit
sgroup_u <- rbinom(n_units, 1, 0.5)

# Continuous moderator: mean differs across subgroups (drives the IPW need)
m_u <- rnorm(n_units, mean = 0.5 * sgroup_u, sd = 1)

# Binary treatment indicator, independent of sgroup (~50/50)
D_u <- rbinom(n_units, 1, 0.5)

# Unit fixed effects
alpha_u <- rnorm(n_units, mean = 0, sd = 1)

# True per-subgroup effects: tau_0 = 1, tau_1 = 3 (difference = 2)
tau_u <- ifelse(sgroup_u == 1L, 3, 1)

# Long-form balanced panel: 2 periods per unit
unit <- rep(unit_id, each = 2L)
time <- rep(c(0L, 1L), times = n_units)
post <- as.integer(time == 1L)

sgroup <- rep(sgroup_u, each = 2L)
m      <- rep(m_u,      each = 2L)
D      <- rep(D_u,      each = 2L)
alpha  <- rep(alpha_u,  each = 2L)
tau    <- rep(tau_u,    each = 2L)

eps <- rnorm(length(unit), mean = 0, sd = 0.5)

# DGP: y_it = alpha_i + 0.5 * post_t + tau_{sgroup_i} * D_i * post_t + 0.3 * m_i + eps_it
y <- alpha + 0.5 * post + tau * D * post + 0.3 * m + eps

wsga_did_synth <- data.frame(
  unit   = unit,
  time   = time,
  sgroup = sgroup,
  m      = m,
  D      = D,
  y      = y
)

save(wsga_did_synth, file = "data/wsga_did_synth.rda", compress = "bzip2")

# Stata-format export for cross-language testing
if (requireNamespace("haven", quietly = TRUE)) {
  haven::write_dta(wsga_did_synth, "stata/wsga_did_synth.dta")
} else if (requireNamespace("foreign", quietly = TRUE)) {
  foreign::write.dta(wsga_did_synth, "stata/wsga_did_synth.dta")
} else {
  stop("Neither `haven` nor `foreign` is installed; cannot write .dta file.")
}

message(
  "Saved data/wsga_did_synth.rda and stata/wsga_did_synth.dta  ",
  "(n_rows = ", nrow(wsga_did_synth), ", n_units = ", n_units, "). ",
  "True per-subgroup effects: tau_0 = 1, tau_1 = 3 (diff = 2)."
)
