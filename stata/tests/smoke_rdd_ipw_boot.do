// smoke_rdd_ipw_boot.do -- IPW bootstrap + degenerate-replicate handling (issue #43)
//
// Regression test for the bug where the propensity-score logit clobbered
// e(cmdline), so the pairs-bootstrap loop re-ran the logit instead of the
// outcome model and aborted with r(111) "[0.G#1._cutoff] not found" on the
// FIRST IPW replicate. Also covers the graceful handling of genuinely
// degenerate replicates (thin subgroup x cutoff cells) and the fail-fast
// guard for structurally unidentified subgroups.
//
// Run from rddsga-repo/:
//   stata-mp -b do stata/tests/smoke_rdd_ipw_boot.do && cat smoke_rdd_ipw_boot.log

quietly do stata/wsga.ado

local n_fail = 0

// --- Test 1: healthy IPW bootstrap must run (was aborting with r(111)) ---
use stata/rddsga_synth, clear
capture wsga rdd Y, sgroup(G) running(X) bwidth(20) reducedform ///
    balance(W1 W2) m(2) bsreps(40) seed(42) rbalance(0)
if _rc != 0 {
    di as error "FAIL [healthy IPW bootstrap]: rc=" _rc " (expected 0)"
    local ++n_fail
}
else if e(B_ok) != e(N_reps) {
    di as error "FAIL [healthy IPW B_ok]: B_ok=" e(B_ok) " != N_reps=" e(N_reps)
    local ++n_fail
}
else if mi(e(b)[1,1]) | mi(e(b)[1,2]) {
    di as error "FAIL [healthy IPW estimates]: G0=" e(b)[1,1] " G1=" e(b)[1,2]
    local ++n_fail
}
else {
    di as result "PASS [healthy IPW bootstrap]: B_ok=" e(B_ok) "/" e(N_reps)
}

// --- Test 2: noipsw pairs bootstrap still works (regression guard) ---
use stata/rddsga_synth, clear
capture wsga rdd Y, sgroup(G) running(X) bwidth(20) reducedform noipsw ///
    bsreps(40) seed(42) rbalance(0)
if _rc == 0 & e(B_ok) == e(N_reps) {
    di as result "PASS [noipsw bootstrap]: B_ok=" e(B_ok) "/" e(N_reps)
}
else {
    di as error "FAIL [noipsw bootstrap]: rc=" _rc " B_ok=" e(B_ok) " N_reps=" e(N_reps)
    local ++n_fail
}

// --- Test 3: degenerate replicates (thin cell) are dropped, not silently kept ---
// G0 has only 2 below-cutoff obs, so some resamples empty that cell.
use stata/rddsga_synth, clear
set seed 55
keep if abs(X) < 12
gen byte _above = X > 0
gen byte sg = 1
replace sg = 0 if _above==1 & runiform() < 0.5
gen double _r = runiform() if _above==0
sort _r
gen long _rk = _n if _above==0
replace sg = 0 if _above==0 & _rk <= 2
replace sg = 1 if _above==0 & _rk >  2
capture wsga rdd Y, sgroup(sg) running(X) bwidth(12) reducedform ///
    balance(W1 W2) m(2) bsreps(100) seed(42) rbalance(0)
if _rc != 0 {
    di as error "FAIL [thin-cell run]: rc=" _rc " (expected 0)"
    local ++n_fail
}
else if e(B_ok) >= 1 & e(B_ok) < e(N_reps) {
    di as result "PASS [thin-cell partial drop]: B_ok=" e(B_ok) "/" e(N_reps) " (some reps dropped)"
}
else {
    di as error "FAIL [thin-cell partial drop]: B_ok=" e(B_ok) " N_reps=" e(N_reps) " (expected 1 <= B_ok < N_reps)"
    local ++n_fail
}
drop _above _r _rk

// --- Test 4: fail-fast on a structurally unidentified subgroup ---
// One subgroup lives entirely on one side of the cutoff -> abort with r(111).
use stata/rddsga_synth, clear
gen byte sg1 = (X < 0)        // sg1==1 iff below cutoff: no below-cutoff obs for sg1==0
capture wsga rdd Y, sgroup(sg1) running(X) bwidth(20) reducedform ///
    noipsw bsreps(10) seed(1)
if _rc == 111 {
    di as result "PASS [fail-fast]: unidentified subgroup correctly aborted with r(111)"
}
else {
    di as error "FAIL [fail-fast]: expected r(111), got rc=" _rc
    local ++n_fail
}

// --- Test 5: IPW fuzzy ivregress + first-stage bootstrap (same shared loop) ---
foreach spec in "ivregress" "firststage" {
    use stata/rddsga_synth, clear
    capture wsga rdd Y, sgroup(G) running(X) bwidth(20) `spec' fuzzy(D) ///
        balance(W1 W2) m(2) bsreps(30) seed(42) rbalance(0)
    if _rc == 0 & e(B_ok) == e(N_reps) & !mi(e(b)[1,1]) {
        di as result "PASS [IPW `spec' bootstrap]: B_ok=" e(B_ok) "/" e(N_reps)
    }
    else {
        di as error "FAIL [IPW `spec' bootstrap]: rc=" _rc " B_ok=" e(B_ok)
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
