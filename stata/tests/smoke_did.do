// smoke_did.do — tests for wsga did
// Run from rddsga-repo/:
//   stata-mp -b do stata/tests/smoke_did.do && cat smoke_did.log

quietly do stata/wsga.ado
use stata/wsga_did_synth, clear

local n_fail = 0

// ── TEST 1: seed() makes DiD bootstrap reproducible ──────────────────────────
wsga did y m, sgroup(sgroup) unit(unit) time(time) treat(D) ///
  bsreps(20) seed(42) noipsw
scalar _b0_run1 = e(b)[1,1]
scalar _b1_run1 = e(b)[1,2]

wsga did y m, sgroup(sgroup) unit(unit) time(time) treat(D) ///
  bsreps(20) seed(42) noipsw
scalar _b0_run2 = e(b)[1,1]
scalar _b1_run2 = e(b)[1,2]

if abs(_b0_run1 - _b0_run2) < 1e-10 & abs(_b1_run1 - _b1_run2) < 1e-10 {
    di as result "PASS [seed: point estimates stable]"
}
else {
    di as error "FAIL [seed: point estimates differ across runs]"
    local ++n_fail
}

// ── TEST 2: fixedps runs without error ────────────────────────────────────────
capture wsga did y m, sgroup(sgroup) unit(unit) time(time) treat(D) ///
  bsreps(20) seed(7) fixedps
if _rc != 0 {
    di as error "FAIL [fixedps]: rc=" _rc
    local ++n_fail
}
else {
    di as result "PASS [fixedps: runs without error]"
}

// ── TEST 3: fixedps SE <= refit-per-rep SE ────────────────────────────────────
wsga did y m, sgroup(sgroup) unit(unit) time(time) treat(D) ///
  bsreps(100) seed(1)
matrix _V_refit = e(V)
scalar _se_diff_refit = sqrt(_V_refit[1,1] + _V_refit[2,2] - 2*_V_refit[1,2])

wsga did y m, sgroup(sgroup) unit(unit) time(time) treat(D) ///
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

// ── TEST 4: nobootstrap returns analytical SEs ────────────────────────────────
capture wsga did y m, sgroup(sgroup) unit(unit) time(time) treat(D) nobootstrap
if _rc != 0 {
    di as error "FAIL [nobootstrap]: rc=" _rc
    local ++n_fail
}
else {
    di as result "PASS [nobootstrap: runs without error]"
}

if `n_fail' == 0 {
    di as result _newline "All smoke tests passed."
}
else {
    di as error _newline "`n_fail' smoke test(s) FAILED."
    exit 1
}
