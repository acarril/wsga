# wsga (R package) NEWS

## wsga 1.3.1 (2026-05-26)

### Internal

- Add Monte Carlo simulation study (`inst/run_simulations.R`) validating
  four core estimator properties: bias removal (S1), size control under
  the null (S2), power as a function of K moderators (S3), and sensitivity
  to unobservable confounding (S4). Full results (1000 reps, 10 cores)
  saved to `inst/simulation_results/*.csv`. Closes #29.

---

## wsga 1.3.0 (2026-05-26)

### New features

- **R + Stata**: add restricted wild cluster bootstrap (WCB-R) for better
  size control at very small cluster counts (G < ~12), following MacKinnon
  and Webb (2017, 2018). In R: `boot_type = "wild_restricted"`; in Stata:
  `wcbrestricted` option on both `wsga rdd` and `wsga did`.
- **Architecture**: WCB-R uses a combined loop -- unrestricted draws for
  empirical percentile CIs (per coauthor agreement) and three restricted
  draws (H0: G0=0, H0: G1=0, H0: G0=G1) for p-values. P-value formula
  for restricted draws does not recenter (`|draw| >= |est|`), matching
  MacKinnon-Webb eq. 8.
- **Advisories updated**: G < 30 with `pairs` recommends `boot_type =
  "wild"` (unchanged); G < 12 with `boot_type = "wild"` now recommends
  `boot_type = "wild_restricted"` (#36).

---

## wsga 1.2.2 (2026-05-26)

### Bug fixes

- **Stata** (`wsga rdd`, WCB path): the bootstrap loop no longer mutates
  the user's outcome variable in place. Previously `_wsga_rdd_myboo` did
  `replace <depvar> = <ystar>` inside `preserve/restore`, which was correct
  in the happy path but left the outcome corrupted in memory if a replicate
  errored between `replace` and `restore`. WCB now parses the cmdline once
  via `gettoken` and refits as `<cmd> <ystar_tempvar> <rhs>` per replicate,
  never touching the dataset's copy of the outcome (#37).
- **Stata** (`wsga rdd`): missing cluster IDs in the estimation sample now
  error out (`exit 198`) instead of being silently lumped into one implicit
  group by `by cluster`. Matches R behavior. DiD is unaffected because
  `markout touse unit time treat sgroup` already drops missing-unit rows
  (#37).

### Internal

- New regression tests in `stata/tests/smoke_rdd_wild.do` (2 added, 10 total):
  (a) `Y` is byte-for-byte unchanged after a WCB run, (b) missing cluster
  IDs cause `wsga rdd` to error with `rc != 0`.

---

## wsga 1.2.1 (2026-05-26)

### Documentation

- Update the `wsga-version-bump` skill to include a post-merge Step 8
  for tagging the release commit and publishing it via `gh release create`.
  Codifies the release-per-version-bump policy adopted retroactively when
  v1.0.1 through v1.2.0 were backfilled.

---

## wsga 1.2.0 (2026-05-26)

### New features

- **Stata**: `wsga rdd` gains two new bootstrap options matching `wsga did`:
  - `cluster(varname)` -- pairs-cluster bootstrap. Whole clusters are
    resampled with replacement (`bsample, cluster()`). Composes cleanly
    with `blockbootstrap()` for stratified cluster resampling. No silent
    default: clustering must be specified explicitly.
  - `wildcluster` -- unrestricted wild cluster bootstrap (WCB-U) with
    Rademacher signs at the cluster level. Conditions on the data, sign-flips
    the residuals from the main fit, refits on `y_star = xb + sign * resid`.
    Recommended at small G (<~30 clusters), where pairs over-rejects under
    H0 (Cameron, Gelbach & Miller 2008). Requires `cluster()`. Not supported
    with `ivregress` (fuzzy RD via 2SLS); WCB on a 2SLS fit needs a different
    recipe (#30).
- **Stata**: `wsga rdd` now emits a G<30 advisory when `cluster()` is set
  and there are fewer than 30 unique clusters, recommending `wildcluster`.
  Mirrors the existing DiD advisory.
- **Stata**: `wsga rdd` posts new e-class results: `e(B_ok)` (successful
  bootstrap reps), `e(N_clust)` (cluster count, when `cluster()` set),
  `e(boot_type)` (`pairs` or `wild`), `e(clustvar)` (clustering variable
  name).
- **Stata**: `wsga rdd` display label now distinguishes
  `Bootstrap replications` / `Cluster bootstrap` / `Wild cluster bootstrap`
  to match what was actually run.

### Documentation

- `wsga_rdd.sthlp` gains a "Stored results" section mirroring
  `wsga_did.sthlp`, plus examples for `cluster()` and `wildcluster`.

### Internal

- New `stata/tests/smoke_rdd_wild.do`: 8 checks covering pairs-cluster
  bookkeeping, wildcluster runs and bookkeeping, pairs no-cluster
  unchanged, validation errors (wildcluster+nobootstrap,
  wildcluster-without-cluster, wildcluster+ivregress), seed reproducibility,
  and the G<30 advisory path.

---

## wsga 1.1.0 (2026-05-25)

### Breaking changes

- **Stata**: `wsga rdd` e-class return names standardized to match `wsga did`.
  RDD previously posted `e(lb_g0)`, `e(ub_g0)`, `e(lb_g1)`, `e(ub_g1)`,
  `e(lb_diff)`, `e(ub_diff)`, `e(pval0)`, `e(pval1)`, `e(pval_diff)`. These
  are now `e(ci_lb_g0)`, `e(ci_ub_g0)`, `e(ci_lb_g1)`, `e(ci_ub_g1)`,
  `e(ci_lb_diff)`, `e(ci_ub_diff)`, `e(p_g0)`, `e(p_g1)`, `e(p_diff)`,
  matching what `wsga did` already returned. Downstream code reading these
  values from `wsga rdd` must be updated.

### Bug fixes

- **Stata**: `wsga rdd` now overwrites `e(p_g0)`, `e(p_g1)`, `e(p_diff)`,
  `e(ci_lb_g0)`, `e(ci_ub_g0)`, `e(ci_lb_g1)`, `e(ci_ub_g1)`, `e(ci_lb_diff)`,
  `e(ci_ub_diff)` with the mode-aware locals after the display block.
  Previously `_wsga_rdd_myboo` posted these unconditionally as empirical
  percentiles, so when `normal` was specified the displayed CI columns and
  p-values were correct but the e-class returns silently kept the empirical
  values. Empirical mode is a no-op (#32).
- **Stata**: stripped non-ASCII characters from `wsga.ado`, `wsga_rdd.sthlp`,
  `wsga_did.sthlp`, and `NEWS.md` to satisfy the repository's ASCII-only
  policy. No semantic content changed.

### Internal

- New `stata/tests/smoke_inference_modes.do`: 18 checks (9 RDD + 9 DiD)
  asserting `e()` CI bounds and p-values differ between empirical and normal
  modes (regression guard for #32). Now uses a unified key list across
  both designs, courtesy of the rename above.

---

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
- **Stata**: `wsga did` now implements `comsup` -- units outside the G=1
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

- **`design = "did"`**: sharp 2-period DiD-SGA pipeline. New required arguments `unit`, `time`, `treat`; optional `post_value` (defaults to `max(time)`). Runs long-form TWFE with subgroup x post interactions, IPW reweighting, and pairs cluster bootstrap over units (Cameron-Gelbach-Miller) by default.
- **`inference` argument**: three-mode inference -- `"empirical"` (default with bootstrap, percentile CIs), `"normal"` (normal-approx from bootstrap SE), `"analytical"` (sandwich/cluster-robust, no bootstrap).
- **`fixed_ps`**: opt-in flag to hold the propensity score fixed at the original-sample fit across bootstrap replicates. Useful for variance decomposition. Default `FALSE`.
- **`cluster_var`**: drives both analytical SEs and the bootstrap. DiD defaults to clustering on `unit`. One-time warning when fewer than 30 unique clusters are present.
- **Bundled dataset**: `wsga_did_synth` -- balanced 2-period panel of 500 units (seed = 7). True effects `tau_0 = 1`, `tau_1 = 3`.
- **DiD balance tables**: aggregate and treated-only (D = 1) balance, each unweighted and IPW-weighted.
- **Vignette**: `vignette("wsga-did-intro")` walks through the full DiD workflow.
- **`bsreps` default** bumped 50 -> 200.

### Breaking changes

- Balance accessor: `fit$balance$unweighted$table` -> `fit$balance$unweighted$aggregate$table` (nested structure to accommodate the DiD treated-only block).
- Default `bsreps` changed from 50 to 200.

### Package rename (from rddsga)

- Package renamed `rddsga` -> `wsga`. The old function `rddsga()` was retained as a deprecated alias (now forwards to `wsga_rdd()`). Plan to remove in v2.

---

## rddsga 0.3.0 -> wsga 0.6.x (internal milestones)

- Umbrella refactor: single `wsga()` entry point dispatching both RD and DiD paths (replaced by `wsga_rdd()` / `wsga_did()` in 1.0.0).
- S3 methods: `print`, `summary`, `coef`, `vcov`, `confint`, `nobs`.
- `block_var` stratification for the bootstrap.
- `seed` argument for reproducible bootstrap draws.
- `fixed_fs` for fuzzy RD: hold first-stage estimate fixed across bootstrap reps.
- Three kernel types: uniform, triangular, Epanechnikov.

---

## rddsga 0.2.x / 0.3.0 (2018-2023)

Initial CRAN-adjacent R implementation of the RDD-SGA estimator accompanying
the working paper "Weighted Subgroup Analysis in Regression Discontinuity
Designs" (Carril, Cazor, Gerardino, Litschig, Pomeranz).
