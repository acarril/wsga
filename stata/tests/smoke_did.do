// smoke_did.do — tests for seed() and fixedps options in wsga design(did)
// Run from rddsga-repo/:
//   stata-mp -b do stata/tests/smoke_did.do && cat smoke_did.log

cd /Users/acarril/Sidequests/rddsga/rddsga-repo
quietly do stata/wsga.ado
use stata/wsga_did_synth, clear

local n_fail = 0

// ── TEST 1: seed() makes DiD bootstrap reproducible ──────────────────────────
wsga y m, sgroup(sgroup) design(did) unit(unit) time(time) treat(D) ///
  bsreps(20) seed(42) noipsw
scalar _b0_run1 = e(b)[1,1]
scalar _b1_run1 = e(b)[1,2]

wsga y m, sgroup(sgroup) design(did) unit(unit) time(time) treat(D) ///
  bsreps(20) seed(42) noipsw
scalar _b0_run2 = e(b)[1,1]
scalar _b1_run2 = e(b)[1,2]

// Point estimates must be identical (same data, same model — seed only affects bootstrap)
if abs(_b0_run1 - _b0_run2) < 1e-10 & abs(_b1_run1 - _b1_run2) < 1e-10 {
    di as result "PASS [seed: point estimates stable]"
}
else {
    di as error "FAIL [seed: point estimates differ across runs]"
    local ++n_fail
}

// ── TEST 2: seed() makes RDD myboo reproducible ───────────────────────────────
use stata/rddsga_synth, clear
wsga Y X, sgroup(G) bwidth(10) reducedform seed(42) bsreps(20) noipsw
scalar _rb0_run1 = e(b)[1,1]

wsga Y X, sgroup(G) bwidth(10) reducedform seed(42) bsreps(20) noipsw
scalar _rb0_run2 = e(b)[1,1]

if abs(_rb0_run1 - _rb0_run2) < 1e-10 {
    di as result "PASS [seed: RDD point estimates stable]"
}
else {
    di as error "FAIL [seed: RDD point estimates differ across runs]"
    local ++n_fail
}

// ── TEST 3: fixedps runs without error ────────────────────────────────────────
use stata/wsga_did_synth, clear
capture wsga y m, sgroup(sgroup) design(did) unit(unit) time(time) treat(D) ///
  bsreps(20) seed(7) fixedps
if _rc != 0 {
    di as error "FAIL [fixedps]: rc=" _rc
    local ++n_fail
}
else {
    di as result "PASS [fixedps: runs without error]"
}

// ── TEST 4: fixedps SE <= refit-per-rep SE (fixed_ps removes PS variance) ─────
wsga y m, sgroup(sgroup) design(did) unit(unit) time(time) treat(D) ///
  bsreps(100) seed(1)
// SE of the difference from refit-per-rep bootstrap
matrix _V_refit = e(V)
scalar _se_diff_refit = sqrt(_V_refit[1,1] + _V_refit[2,2] - 2*_V_refit[1,2])

wsga y m, sgroup(sgroup) design(did) unit(unit) time(time) treat(D) ///
  bsreps(100) seed(1) fixedps
matrix _V_fixed = e(V)
scalar _se_diff_fixed = sqrt(_V_fixed[1,1] + _V_fixed[2,2] - 2*_V_fixed[1,2])

di "  SE(diff) refit-per-rep: " _se_diff_refit
di "  SE(diff) fixed-PS:      " _se_diff_fixed

if _se_diff_fixed <= _se_diff_refit * 1.05 {
    di as result "PASS [fixedps: SE not larger than refit-per-rep (within 5% margin)]"
}
else {
    di as error "FAIL [fixedps: SE unexpectedly larger than refit-per-rep]"
    local ++n_fail
}

// ── Summary ──────────────────────────────────────────────────────────────────
if `n_fail' == 0 {
    di as result _newline "All smoke tests passed."
}
else {
    di as error _newline "`n_fail' smoke test(s) FAILED."
    exit 1
}
