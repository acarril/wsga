# wsga (R package) NEWS

## wsga 1.0.3 (2026-05-12)

### New features

- **R**: New `boot_type = c("pairs", "wild")` argument on `wsga_rdd()` and
  `wsga_did()`. When `"wild"`, runs an unrestricted wild cluster bootstrap
  (WCB-U) with Rademacher signs at the cluster level. Recommended at small
  G (< ~30 clusters), where pairs over-rejects under H0 (Cameron, Gelbach &
  Miller 2008). Requires `cluster_var` to be set; not supported with
  `model = "iv"` (#30).
- **Stata**: New `wildcluster` option on `wsga did`. Same WCB-U scheme as
  the R implementation: predicts the fitted values (`xbu`) and idiosyncratic
  residuals (`e`) from the main `xtreg, fe` fit, then sign-flips residuals
  at the unit level per replicate and refits the weighted regression. New
  `e(boot_type)` macro records which scheme ran (#30).
- **Both**: the G < 30 advisory now recommends the wild-cluster option
  directly (previously pointed to external tools).

### Notes

- WCB conditions on the data and does not refit the propensity score, so it
  does not propagate IPW estimation uncertainty. This is a deliberate
  tradeoff for better size control at small G; use the pairs bootstrap
  (default) if IPW uncertainty propagation matters more.
- Stata RDD WCB is deferred to a follow-up: `_wsga_rdd_myboo` does not
  currently expose a `cluster()` path, which is a prerequisite refactor.
- Bump also unifies pre-existing version drift in `stata/wsga.pkg`,
  `stata/rddsga.pkg`, and `stata/stata.toc`.

---

## wsga 1.0.2 (2026-05-11)

### Bug fixes

- **Stata**: `wsga did` now correctly wires up `ipsweight()` and `pscore()` as
  named output variables in the dataset (previously accepted but silently
  ignored) (#25).
- **Stata**: `wsga did` now implements `comsup` — units outside the G=1
  propensity score range are excluded from estimation and a `comsup` variable
  is created in the dataset, matching RDD behavior. Common support is
  re-evaluated per bootstrap replicate (#25).
- **Stata**: `wsga did` now implements `blockbootstrap(varname)` as a
  stratified unit-level cluster resample (`bsample, strata(varname)`); the
  variable is validated for unit-constancy (#25).
- **Stata**: align version numbers across all Stata files to unified `1.0.x`
  scheme (`wsga.sthlp`, `wsga_rdd.sthlp`, `rddsga.*` previously lagged).

---

## wsga 1.0.1 (2026-05-11)

### Bug fixes

- **Stata**: `wsga did` now calls `ereturn post` at the end of `_wsga_did`,
  trimming `e(b)` and `e(V)` to the two treatment-effect columns (`G0_Z`,
  `G1_Z`). Previously the full `xtreg, fe` result was left in `e()`, so
  `e(b)` / `e(V)` were unreliable in downstream code. Bootstrap path now also
  computes the covariance of the two coefficient draws so that `e(V)` is fully
  bootstrap-derived when bootstrap is on (#24).
- **Stata**: align `*!` version line with unified `1.0.x` versioning.

---

## wsga 1.0.0 (2026-05-11)

### Breaking changes

- `wsga()` is **removed**. Replace all calls with `wsga_rdd()` (RD designs) or
  `wsga_did()` (DiD designs). Each function accepts only the arguments relevant
  to its design; the `design =` argument no longer exists.
- `rddsga()` (deprecated alias) now forwards to `wsga_rdd()`.

---

## wsga 0.7.0 (2026-05-11)

### New features

- **`design = "did"`**: sharp 2-period DiD-SGA pipeline. New required arguments `unit`, `time`, `treat`; optional `post_value` (defaults to `max(time)`). Runs long-form TWFE with subgroup × post interactions, IPW reweighting, and pairs cluster bootstrap over units (Cameron-Gelbach-Miller) by default.
- **`inference` argument**: three-mode inference — `"empirical"` (default with bootstrap, percentile CIs), `"normal"` (normal-approx from bootstrap SE), `"analytical"` (sandwich/cluster-robust, no bootstrap).
- **`fixed_ps`**: opt-in flag to hold the propensity score fixed at the original-sample fit across bootstrap replicates. Useful for variance decomposition. Default `FALSE`.
- **`cluster_var`**: drives both analytical SEs and the bootstrap. DiD defaults to clustering on `unit`. One-time warning when fewer than 30 unique clusters are present.
- **Bundled dataset**: `wsga_did_synth` — balanced 2-period panel of 500 units (seed = 7). True effects `tau_0 = 1`, `tau_1 = 3`.
- **DiD balance tables**: aggregate and treated-only (D = 1) balance, each unweighted and IPW-weighted.
- **Vignette**: `vignette("wsga-did-intro")` walks through the full DiD workflow.
- **`bsreps` default** bumped 50 → 200.

### Breaking changes

- Balance accessor: `fit$balance$unweighted$table` → `fit$balance$unweighted$aggregate$table` (nested structure to accommodate the DiD treated-only block).
- Default `bsreps` changed from 50 to 200.

### Package rename (from rddsga)

- Package renamed `rddsga` → `wsga`. The old function `rddsga()` was retained as a deprecated alias (now forwards to `wsga_rdd()`). Plan to remove in v2.

---

## rddsga 0.3.0 → wsga 0.6.x (internal milestones)

- Umbrella refactor: single `wsga()` entry point dispatching both RD and DiD paths (replaced by `wsga_rdd()` / `wsga_did()` in 1.0.0).
- S3 methods: `print`, `summary`, `coef`, `vcov`, `confint`, `nobs`.
- `block_var` stratification for the bootstrap.
- `seed` argument for reproducible bootstrap draws.
- `fixed_fs` for fuzzy RD: hold first-stage estimate fixed across bootstrap reps.
- Three kernel types: uniform, triangular, Epanechnikov.

---

## rddsga 0.2.x / 0.3.0 (2018–2023)

Initial CRAN-adjacent R implementation of the RDD-SGA estimator accompanying
the working paper "Weighted Subgroup Analysis in Regression Discontinuity
Designs" (Carril, Cazor, Gerardino, Litschig, Pomeranz).
