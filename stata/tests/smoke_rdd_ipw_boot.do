// smoke_rdd_ipw_boot.do -- IPW reduced-form bootstrap regression test (issue #43)
//
// Before the fix, the per-replicate propensity-score logit clobbered e(cmdline),
// so the pairs-bootstrap loop re-ran the logit instead of the outcome model and
// aborted with r(111) "[0.G#1._cutoff] not found" on the FIRST IPW replicate.
// Every IPW (m()) RDD bootstrap was affected; noipsw was not. This test runs an
// IPW bootstrap (which used to crash) and a noipsw bootstrap (regression guard).
//
// Run from rddsga-repo/:
//   stata-mp -b do stata/tests/smoke_rdd_ipw_boot.do && cat smoke_rdd_ipw_boot.log

quietly do stata/wsga.ado

local n_fail = 0

// --- Test 1: IPW reduced-form bootstrap must run (was aborting with r(111)) ---
use stata/rddsga_synth, clear
capture wsga rdd Y, sgroup(G) running(X) bwidth(20) reducedform ///
    balance(W1 W2) m(2) bsreps(40) seed(42) rbalance(0)
if _rc != 0 {
    di as error "FAIL [IPW bootstrap]: rc=" _rc " (expected 0)"
    local ++n_fail
}
else if mi(e(b)[1,1]) | mi(e(b)[1,2]) {
    di as error "FAIL [IPW estimates]: G0=" e(b)[1,1] " G1=" e(b)[1,2]
    local ++n_fail
}
else {
    di as result "PASS [IPW bootstrap]: G0=" e(b)[1,1] " G1=" e(b)[1,2]
}

// --- Test 2: noipsw pairs bootstrap still works (regression guard) ---
use stata/rddsga_synth, clear
capture wsga rdd Y, sgroup(G) running(X) bwidth(20) reducedform noipsw ///
    bsreps(40) seed(42) rbalance(0)
if _rc == 0 & !mi(e(b)[1,1]) {
    di as result "PASS [noipsw bootstrap]"
}
else {
    di as error "FAIL [noipsw bootstrap]: rc=" _rc
    local ++n_fail
}

// --- Test 3: IPW fuzzy IV + first-stage bootstrap (same shared loop) ---
foreach spec in "ivregress" "firststage" {
    use stata/rddsga_synth, clear
    capture wsga rdd Y, sgroup(G) running(X) bwidth(20) `spec' fuzzy(D) ///
        balance(W1 W2) m(2) bsreps(30) seed(42) rbalance(0)
    if _rc == 0 & !mi(e(b)[1,1]) {
        di as result "PASS [IPW `spec' bootstrap]"
    }
    else {
        di as error "FAIL [IPW `spec' bootstrap]: rc=" _rc
        local ++n_fail
    }
}

if `n_fail' == 0 {
    di as result _newline "All IPW-bootstrap smoke tests passed."
}
else {
    di as error _newline "`n_fail' IPW-bootstrap smoke test(s) FAILED."
    exit 1
}
