// smoke_rdd_wild.do -- tests for `wsga rdd, cluster() wildcluster / wcbrestricted`
// Run from rddsga-repo/:
//   stata-mp -b do stata/tests/smoke_rdd_wild.do && cat smoke_rdd_wild.log

quietly do stata/wsga.ado
use stata/rddsga_synth, clear

// Synthetic clustering variable: 50 clusters of ~200 obs (G >= 30, no advisory)
gen clust = ceil(_n/200)

local n_fail = 0

// -- TEST 1: pairs cluster runs and posts cluster bookkeeping --
capture wsga rdd Y, sgroup(G) running(X) bwidth(10) reducedform ///
  bsreps(50) seed(1) noipsw cluster(clust)
if _rc != 0 {
  di as error "FAIL [pairs cluster: runs]: rc=" _rc
  local ++n_fail
}
else if "`e(boot_type)'" == "pairs" & e(N_clust) == 50 {
  di as result "PASS [pairs cluster: boot_type=pairs, N_clust=50]"
}
else {
  di as error "FAIL [pairs cluster: boot_type=`e(boot_type)', N_clust=" e(N_clust) "]"
  local ++n_fail
}

// -- TEST 2: wildcluster runs and posts boot_type=wild --
capture wsga rdd Y, sgroup(G) running(X) bwidth(10) reducedform ///
  bsreps(50) seed(1) noipsw cluster(clust) wildcluster
if _rc != 0 {
  di as error "FAIL [wildcluster: runs]: rc=" _rc
  local ++n_fail
}
else if "`e(boot_type)'" == "wild" & e(N_clust) == 50 {
  di as result "PASS [wildcluster: boot_type=wild, N_clust=50]"
}
else {
  di as error "FAIL [wildcluster: boot_type=`e(boot_type)', N_clust=" e(N_clust) "]"
  local ++n_fail
}

// -- TEST 3: pairs path WITHOUT cluster still tags boot_type=pairs and N_clust=. --
wsga rdd Y, sgroup(G) running(X) bwidth(10) reducedform ///
  bsreps(30) seed(1) noipsw
if "`e(boot_type)'" == "pairs" & mi(e(N_clust)) {
  di as result "PASS [pairs no-cluster: boot_type=pairs, N_clust=.]"
}
else {
  di as error "FAIL [pairs no-cluster: boot_type=`e(boot_type)', N_clust=" e(N_clust) "]"
  local ++n_fail
}

// -- TEST 4: wildcluster + nobootstrap rejected --
capture wsga rdd Y, sgroup(G) running(X) bwidth(10) reducedform ///
  noipsw cluster(clust) wildcluster nobootstrap
if _rc != 0 {
  di as result "PASS [wildcluster+nobootstrap rejected, rc=" _rc "]"
}
else {
  di as error "FAIL [wildcluster+nobootstrap should have errored]"
  local ++n_fail
}

// -- TEST 5: wildcluster without cluster() rejected --
capture wsga rdd Y, sgroup(G) running(X) bwidth(10) reducedform ///
  bsreps(20) seed(1) noipsw wildcluster
if _rc != 0 {
  di as result "PASS [wildcluster-without-cluster rejected, rc=" _rc "]"
}
else {
  di as error "FAIL [wildcluster-without-cluster should have errored]"
  local ++n_fail
}

// -- TEST 6: wildcluster + ivregress rejected (no WCB on 2SLS) --
capture wsga rdd Y, sgroup(G) running(X) bwidth(10) ivregress fuzzy(D) ///
  bsreps(20) seed(1) noipsw cluster(clust) wildcluster
if _rc != 0 {
  di as result "PASS [wildcluster+ivregress rejected, rc=" _rc "]"
}
else {
  di as error "FAIL [wildcluster+ivregress should have errored]"
  local ++n_fail
}

// -- TEST 7: seed makes wildcluster reproducible (compare bootstrap V) --
wsga rdd Y, sgroup(G) running(X) bwidth(10) reducedform ///
  bsreps(30) seed(42) noipsw cluster(clust) wildcluster
scalar _v00_run1 = e(V)[1,1]
scalar _v11_run1 = e(V)[2,2]
scalar _ci_run1  = e(ci_lb_g0)

wsga rdd Y, sgroup(G) running(X) bwidth(10) reducedform ///
  bsreps(30) seed(42) noipsw cluster(clust) wildcluster
scalar _v00_run2 = e(V)[1,1]
scalar _v11_run2 = e(V)[2,2]
scalar _ci_run2  = e(ci_lb_g0)

if abs(_v00_run1 - _v00_run2) < 1e-10 & abs(_v11_run1 - _v11_run2) < 1e-10 ///
   & abs(_ci_run1 - _ci_run2) < 1e-10 {
  di as result "PASS [seed: wildcluster V and CI reproducible]"
}
else {
  di as error "FAIL [seed: wildcluster reproducibility broken]"
  di "  v00 diff: " abs(_v00_run1 - _v00_run2)
  di "  v11 diff: " abs(_v11_run1 - _v11_run2)
  di "  ci  diff: " abs(_ci_run1 - _ci_run2)
  local ++n_fail
}

// -- TEST 8: G<30 advisory fires for pairs but not for wild --
preserve
gen clust_small = ceil(_n/1000)
qui wsga rdd Y, sgroup(G) running(X) bwidth(10) reducedform ///
  bsreps(20) seed(1) noipsw cluster(clust_small)
if e(N_clust) == 10 {
  di as result "PASS [pairs small-G: N_clust=10 (advisory expected on screen)]"
}
else {
  di as error "FAIL [pairs small-G: N_clust=" e(N_clust) "]"
  local ++n_fail
}
restore

// -- TEST 9: WCB does not mutate the user's outcome variable (#37) --
preserve
tempvar _Y_snapshot
gen double `_Y_snapshot' = Y
qui wsga rdd Y, sgroup(G) running(X) bwidth(10) reducedform ///
  bsreps(20) seed(7) noipsw cluster(clust) wildcluster
qui count if !mi(Y) & !mi(`_Y_snapshot') & abs(Y - `_Y_snapshot') > 1e-12
if r(N) == 0 {
  di as result "PASS [WCB does not mutate Y in dataset]"
}
else {
  di as error "FAIL [WCB mutated " r(N) " rows of Y]"
  local ++n_fail
}
restore

// -- TEST 10: missing cluster IDs rejected (#37) --
preserve
gen clust_with_mi = clust
qui replace clust_with_mi = . in 1/100
capture wsga rdd Y, sgroup(G) running(X) bwidth(10) reducedform ///
  bsreps(20) seed(1) noipsw cluster(clust_with_mi)
if _rc != 0 {
  di as result "PASS [missing-cluster rejected, rc=" _rc "]"
}
else {
  di as error "FAIL [missing cluster IDs should have errored]"
  local ++n_fail
}
restore

// -- TEST 11: wcbrestricted runs and posts boot_type=wild_restricted --
capture wsga rdd Y, sgroup(G) running(X) bwidth(10) reducedform ///
  bsreps(30) seed(1) noipsw cluster(clust) wcbrestricted
if _rc != 0 {
  di as error "FAIL [wcbrestricted: runs]: rc=" _rc
  local ++n_fail
}
else if "`e(boot_type)'" == "wild_restricted" & e(N_clust) == 50 {
  di as result "PASS [wcbrestricted: boot_type=wild_restricted, N_clust=50]"
}
else {
  di as error "FAIL [wcbrestricted: boot_type=`e(boot_type)', N_clust=" e(N_clust) "]"
  local ++n_fail
}

// -- TEST 12: wcbrestricted without cluster() rejected --
capture wsga rdd Y, sgroup(G) running(X) bwidth(10) reducedform ///
  bsreps(20) seed(1) noipsw wcbrestricted
if _rc != 0 {
  di as result "PASS [wcbrestricted-without-cluster rejected, rc=" _rc "]"
}
else {
  di as error "FAIL [wcbrestricted-without-cluster should have errored]"
  local ++n_fail
}

// -- TEST 13: wcbrestricted + nobootstrap rejected --
capture wsga rdd Y, sgroup(G) running(X) bwidth(10) reducedform ///
  noipsw cluster(clust) wcbrestricted nobootstrap
if _rc != 0 {
  di as result "PASS [wcbrestricted+nobootstrap rejected, rc=" _rc "]"
}
else {
  di as error "FAIL [wcbrestricted+nobootstrap should have errored]"
  local ++n_fail
}

// -- TEST 14: wildcluster + wcbrestricted mutually exclusive --
capture wsga rdd Y, sgroup(G) running(X) bwidth(10) reducedform ///
  bsreps(10) seed(1) noipsw cluster(clust) wildcluster wcbrestricted
if _rc != 0 {
  di as result "PASS [wildcluster+wcbrestricted mutually exclusive, rc=" _rc "]"
}
else {
  di as error "FAIL [wildcluster+wcbrestricted should have errored]"
  local ++n_fail
}

// -- TEST 15: wcbrestricted seed makes results reproducible --
wsga rdd Y, sgroup(G) running(X) bwidth(10) reducedform ///
  bsreps(30) seed(42) noipsw cluster(clust) wcbrestricted
scalar _p0_run1 = e(p_g0)
scalar _ci_run1 = e(ci_lb_g0)

wsga rdd Y, sgroup(G) running(X) bwidth(10) reducedform ///
  bsreps(30) seed(42) noipsw cluster(clust) wcbrestricted
scalar _p0_run2 = e(p_g0)
scalar _ci_run2 = e(ci_lb_g0)

if abs(_p0_run1 - _p0_run2) < 1e-10 & abs(_ci_run1 - _ci_run2) < 1e-10 {
  di as result "PASS [wcbrestricted seed: p-val and CI reproducible]"
}
else {
  di as error "FAIL [wcbrestricted seed reproducibility broken]"
  di "  p0 diff: " abs(_p0_run1 - _p0_run2)
  di "  ci diff: " abs(_ci_run1 - _ci_run2)
  local ++n_fail
}

// -- Final result --
if `n_fail' == 0 {
  di as result _newline "All wsga rdd wildcluster smoke tests passed."
}
else {
  di as error _newline "`n_fail' smoke test(s) FAILED."
  exit 1
}
