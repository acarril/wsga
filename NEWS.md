# wsga (R package) NEWS

## wsga 0.8.0 (2026-05-11)

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

- Umbrella refactor: single `wsga()` entry point dispatching both RD and DiD paths (replaced by `wsga_rdd()` / `wsga_did()` in 0.8.0).
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
