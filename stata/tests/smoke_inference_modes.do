// smoke_inference_modes.do -- e-class CI bounds and p-values must respect
// the `normal' option, for both RDD and DiD (issue #32).
//
// Run from rddsga-repo/:
//   stata-mp -b do stata/tests/smoke_inference_modes.do && cat smoke_inference_modes.log

quietly do stata/wsga.ado

local n_fail = 0
local tol = 1e-10
local keys ci_lb_g0 ci_ub_g0 ci_lb_g1 ci_ub_g1 ci_lb_diff ci_ub_diff p_g0 p_g1 p_diff

// -- RDD
use stata/rddsga_synth, clear
wsga rdd Y, sgroup(G) running(X) bwidth(10) reducedform bsreps(200) seed(42) noipsw
foreach q of local keys {
  scalar rdd_`q'_emp = e(`q')
}
wsga rdd Y, sgroup(G) running(X) bwidth(10) reducedform bsreps(200) seed(42) noipsw normal
foreach q of local keys {
  scalar d = abs(rdd_`q'_emp - e(`q'))
  if d < `tol' {
    di as error "FAIL [RDD e(`q')]: identical across modes (emp=" rdd_`q'_emp ", nrm=" e(`q') ")"
    local ++n_fail
  }
  else {
    di as result "PASS [RDD e(`q')]: differs by " d
  }
}

// -- DiD (regression guard; was already correct prior to #32)
use stata/wsga_did_synth, clear
wsga did y m, sgroup(sgroup) unit(unit) time(time) treat(D) bsreps(50) seed(42) noipsw
foreach q of local keys {
  scalar did_`q'_emp = e(`q')
}
wsga did y m, sgroup(sgroup) unit(unit) time(time) treat(D) bsreps(50) seed(42) noipsw normal
foreach q of local keys {
  scalar d = abs(did_`q'_emp - e(`q'))
  if d < `tol' {
    di as error "FAIL [DiD e(`q')]: identical across modes (emp=" did_`q'_emp ", nrm=" e(`q') ")"
    local ++n_fail
  }
  else {
    di as result "PASS [DiD e(`q')]: differs by " d
  }
}

if `n_fail' == 0 di as result _newline "All inference-mode e-return checks passed."
else {
  di as error _newline "`n_fail' check(s) FAILED."
  exit 1
}
