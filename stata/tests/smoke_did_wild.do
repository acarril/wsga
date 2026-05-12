// smoke_did_wild.do — tests for `wsga did, wildcluster`
// Run from rddsga-repo/:
//   stata-mp -b do stata/tests/smoke_did_wild.do && cat smoke_did_wild.log

quietly do stata/wsga.ado
use stata/wsga_did_synth, clear

local n_fail = 0

// ── TEST 1: wildcluster runs and posts non-trivial inference ─────────────────
capture wsga did y m, sgroup(sgroup) unit(unit) time(time) treat(D) ///
  bsreps(50) seed(1) noipsw wildcluster
if _rc != 0 {
    di as error "FAIL [wildcluster: runs]: rc=" _rc
    local ++n_fail
}
else {
    di as result "PASS [wildcluster: runs without error]"
    di "  B_ok = " e(B_ok) "  N_clust = " e(N_clust)
    di "  b_g0 = " e(b_g0) "   se_g0 = " e(se_g0) "   p_g0 = " e(p_g0)
    di "  b_g1 = " e(b_g1) "   se_g1 = " e(se_g1) "   p_g1 = " e(p_g1)
    di "  b_diff = " e(b_diff) "  se_diff = " e(se_diff) "  p_diff = " e(p_diff)
    di "  boot_type = " e(boot_type)
}

// ── TEST 2: boot_type ereturn is "wild" ──────────────────────────────────────
wsga did y m, sgroup(sgroup) unit(unit) time(time) treat(D) ///
  bsreps(20) seed(7) noipsw wildcluster
if "`e(boot_type)'" == "wild" {
    di as result "PASS [wildcluster: e(boot_type) == 'wild']"
}
else {
    di as error "FAIL [wildcluster: e(boot_type) = '`e(boot_type)''; expected 'wild']"
    local ++n_fail
}

// ── TEST 3: pairs path still tags boot_type = "pairs" ────────────────────────
wsga did y m, sgroup(sgroup) unit(unit) time(time) treat(D) ///
  bsreps(20) seed(7) noipsw
if "`e(boot_type)'" == "pairs" {
    di as result "PASS [pairs: e(boot_type) == 'pairs']"
}
else {
    di as error "FAIL [pairs: e(boot_type) = '`e(boot_type)''; expected 'pairs']"
    local ++n_fail
}

// ── TEST 4: seed makes wildcluster reproducible ──────────────────────────────
wsga did y m, sgroup(sgroup) unit(unit) time(time) treat(D) ///
  bsreps(30) seed(42) noipsw wildcluster
scalar _se0_run1   = e(se_g0)
scalar _sediff_run1 = e(se_diff)

wsga did y m, sgroup(sgroup) unit(unit) time(time) treat(D) ///
  bsreps(30) seed(42) noipsw wildcluster
scalar _se0_run2   = e(se_g0)
scalar _sediff_run2 = e(se_diff)

if abs(_se0_run1 - _se0_run2) < 1e-10 & abs(_sediff_run1 - _sediff_run2) < 1e-10 {
    di as result "PASS [seed: wildcluster SEs reproducible across runs]"
}
else {
    di as error "FAIL [seed: wildcluster SEs differ across identical-seed runs]"
    local ++n_fail
}

// ── TEST 5: wildcluster + nobootstrap errors ─────────────────────────────────
capture wsga did y m, sgroup(sgroup) unit(unit) time(time) treat(D) ///
  nobootstrap wildcluster
if _rc != 0 {
    di as result "PASS [wildcluster+nobootstrap rejected, rc=" _rc "]"
}
else {
    di as error "FAIL [wildcluster+nobootstrap should have errored]"
    local ++n_fail
}

// ── TEST 6: WCB SE in same order of magnitude as analytical cluster SE ───────
// On the bundled DiD synth (500 units, plenty of clusters), wild SE and
// analytical cluster-robust SE should be within a factor of 3.
wsga did y m, sgroup(sgroup) unit(unit) time(time) treat(D) noipsw nobootstrap
scalar _se_g0_analytical = e(se_g0)

wsga did y m, sgroup(sgroup) unit(unit) time(time) treat(D) ///
  noipsw bsreps(200) seed(1) wildcluster
scalar _se_g0_wcb = e(se_g0)
scalar _ratio = _se_g0_wcb / _se_g0_analytical

di "  analytical cluster SE(g0): " _se_g0_analytical
di "  WCB SE(g0):                " _se_g0_wcb
di "  ratio:                     " _ratio

if _ratio > 1/3 & _ratio < 3 {
    di as result "PASS [wildcluster SE in same order of magnitude as analytical]"
}
else {
    di as error "FAIL [wildcluster SE ratio = " _ratio "; outside (1/3, 3)]"
    local ++n_fail
}

if `n_fail' == 0 {
    di as result _newline "All wildcluster smoke tests passed."
}
else {
    di as error _newline "`n_fail' wildcluster smoke test(s) FAILED."
    exit 1
}
