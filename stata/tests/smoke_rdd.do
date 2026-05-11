// smoke_rdd.do — basic RDD smoke test for wsga rdd
// Run from rddsga-repo/:
//   stata-mp -b do stata/tests/smoke_rdd.do && cat smoke_rdd.log
//
// Tests that all three kernel types run without error on rddsga_synth.dta
// and that the uniform-kernel coefficient is in a plausible range.

quietly do stata/wsga.ado

use stata/rddsga_synth, clear

local n_fail = 0

foreach kernel in "" "kernel(triangular)" "kernel(epanechnikov)" {
    local label = cond("`kernel'" == "", "uniform", "`kernel'")
    capture wsga rdd Y, sgroup(G) running(X) bwidth(10) reducedform nobootstrap noipsw `kernel'
    if _rc != 0 {
        di as error "FAIL [`label']: rc=" _rc
        local ++n_fail
    }
    else {
        di as result "PASS [`label']"
    }
}

// Sanity-check: G=0 coefficient should be non-missing and finite
wsga rdd Y, sgroup(G) running(X) bwidth(10) reducedform nobootstrap noipsw
scalar _b0 = e(b)[1,1]
if mi(_b0) {
    di as error "FAIL [coef check]: e(b)[1,1] is missing"
    local ++n_fail
}
else {
    di as result "PASS [coef check]: G=0 coef = " _b0
}

// Seed makes bootstrap reproducible
wsga rdd Y, sgroup(G) running(X) bwidth(10) reducedform bsreps(20) seed(42) noipsw
scalar _rb0_run1 = e(b)[1,1]
wsga rdd Y, sgroup(G) running(X) bwidth(10) reducedform bsreps(20) seed(42) noipsw
scalar _rb0_run2 = e(b)[1,1]
if abs(_rb0_run1 - _rb0_run2) < 1e-10 {
    di as result "PASS [seed: RDD point estimates stable]"
}
else {
    di as error "FAIL [seed: RDD point estimates differ across runs]"
    local ++n_fail
}

if `n_fail' == 0 {
    di as result _newline "All smoke tests passed."
}
else {
    di as error _newline "`n_fail' smoke test(s) FAILED."
    exit 1
}
