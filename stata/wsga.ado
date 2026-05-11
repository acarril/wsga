*! 1.0.2 Alvaro Carril 2026-05-11

// ── Dispatcher ────────────────────────────────────────────────────────────────
program define wsga
  version 11.1
  gettoken sub rest : 0, parse(" ")
  local sub = strlower("`sub'")
  if "`sub'" == "rdd" {
    _wsga_rdd `rest'
  }
  else if "`sub'" == "did" {
    _wsga_did `rest'
  }
  else {
    di as error `"wsga: unknown subcommand "`sub'". Valid subcommands: rdd, did."'
    di as text "Type {help wsga} for details."
    exit 198
  }
end

// ── RDD implementation ────────────────────────────────────────────────────────
program define _wsga_rdd, eclass
version 11.1
syntax varlist(min=1 numeric fv) [if] [in], ///
  SGroup(name) RUnning(varname) ///
  [ BWidth(real -1) Cutoff(real 0) Kernel(string) fuzzy(name) ///
    IPSWeight(name) PSCore(name) comsup ///
    BALance(varlist numeric) DIBALance probit ///
    IVregress REDUCEDform FIRSTstage vce(string) p(int 1) m(int 2) rbalance(int 0) ///
    noBOOTstrap bsreps(real 200) FIXEDbootstrap FIXEDps BLOCKbootstrap(string) NORMal noipsw weights(string) ///
    Seed(string) ]

// ─── Validate required options ────────────────────────────────────────────────
if `bwidth' < 0 {
  di as error "bwidth() is required and must be positive."
  exit 198
}
if ("`ivregress'" != "" | "`firststage'" != "") & "`fuzzy'" == "" {
  di as error "fuzzy() must be specified with ivregress or firststage."
  exit 198
}
// ──────────────────────────────────────────────────────────────────────────────

*-------------------------------------------------------------------------------
* Check inputs
*-------------------------------------------------------------------------------

// Check that depvar is not a factor variable
local fvops = "`s(fvops)'" == "true" | _caller() >= 11
if `fvops' {
  gettoken first rest : varlist
  _fv_check_depvar `first'
}


if ("`weights'"~="") {
tempvar oweights
g `oweights'=`weights'
}
else {
tempvar oweights
g `oweights'=1
}


/*
// Check that fuzzy() is specified if ivregress is specified
if "`ivregress'" != "" & "`fuzzy'" == "" {
  di as error "fuzzy() must be specified with ivregress"
  exit 198
}*/

// ipsweight(): define new propensity score weighting variable or use a tempvar
if "`ipsweight'" != "" confirm new variable `ipsweight'
else tempvar ipsweight


// comsup(): define tempvar for no common support (default option)
if "`comsup'" == "" {
tempvar NOCOMSUP
g `NOCOMSUP'=1
}
// pscore(): define new propensity score variable or use a tempvar
if "`pscore'" != "" confirm new variable `pscore'
else tempvar pscore


// Issue warning if no covariates and no vars in balance when propensity score weighting is used:
if ("`ipsw'"!="noipsw")  {
if `: list sizeof varlist'<=2 & `: list sizeof balance'==0 {
  di as error "either {it:indepvars} or {bf:balance()} must be specified"
  exit 198
}
}

//  When Propensity score weighting is not used, request covariates only if dibalance is specified:
else {
if "`dibalance'" != ""  & `: list sizeof varlist'<=2 & `: list sizeof balance'==0 {
  di as error "either {it:indepvars} or {bf:balance()} must be specified"
  exit 198
}
}



// Issue warning if options bsreps and normal are specified along with nobootstrap
if "`bootstrap'" == "nobootstrap" & (`bsreps' != 200 | "`normal'" != "") {
  di as text "Warning: options " as result "bsreps" as text " and " as result "normal" ///
    as text " are irrelevant if " as result "nobootstrap" as text " is specified"
}




*-------------------------------------------------------------------------------
* Process inputs
*-------------------------------------------------------------------------------

// Mark observations to be used
marksample touse, novarlist

// Extract outcome variable
local depvar : word 1 of `varlist'

// Running variable (from named option)
local assignvar `running'

// Define covariates list
local covariates : list varlist - depvar

// Add c. stub to continuous covariates in balance for factor interactions
foreach var in `balance' {
  capture _fv_check_depvar `var'
  if _rc != 0 local fv_balance `fv_balance' `var'
  else local fv_balance `fv_balance' c.`var'
}


// Add c. stub to continuous covariates for factor interactions
foreach var in `covariates' {
  capture _fv_check_depvar `var'
  if _rc != 0 local fv_covariates `fv_covariates' `var'
  else local fv_covariates `fv_covariates' c.`var'
}

// Create complementary sgroup var
tempvar sgroup0
qui gen `sgroup0' = (`sgroup' == 0) if !mi(`sgroup')

// Extract balance variables
if "`balance'" == "" & "`ipsw'" != "noipsw"  local balance `covariates'
local n_balance `: word count `balance''

if "`balance'" == "" & "`dibalance'" != "" & "`ipsw'" == "noipsw"  local balance `covariates'
local n_balance `: word count `balance''

// Define model to fit (logit is default)
if "`probit'" != "" local binarymodel probit
else local binarymodel logit

// Create bandwidth condition 
local bwidthtab `bwidth'
local bwidth abs(`assignvar') < `bwidth'


// Create condition whether only group 1 or group 2 is weighted or both groups are weighted. 

if `m'==2 /* both groups weighted*/ {
local mc=`m'
}
if `m'==1 /*only group 1*/ {
local mc=`m'
}
if `m'==0 /*only group 0*/ {
local mc=`m'
}

if `rbalance'==1 {
local rbalance=1 // mean in sample
}
if `rbalance'==0 {
local rbalance=0 // mean in cutoff
}



// Create indicator cutoff variable
*tempvar cutoffvar
*gen _cutoff = (`assignvar'>`cutoff')
*lab var _cutoff "fuzzy"
cap drop _cutoff
gen _cutoff = (`assignvar'>`cutoff')

// Compute spline options
forval i=1/`p' {
tempvar assignvar_`i'
qui gen `assignvar_`i''=`assignvar'^`i'
tempvar gassignvar_`i'
qui gen `gassignvar_`i''=`assignvar_`i''*_cutoff
local polynm `polynm' i.`sgroup'#c.`assignvar_`i'' i.`sgroup'#c.`gassignvar_`i''
}


tempvar kwt
// create weights for kernel 

* default 
if "`kernel'"=="" {
local kernel = "uni" 
}

  if ("`kernel'"=="epanechnikov" | "`kernel'"=="epa") {
    local kernel_type = "Epanechnikov"
    qui g double `kwt'=max(0,3/4*(`bwidthtab'^2-abs(`assignvar')^2))*`oweights'
  }
  else if ("`kernel'"=="triangular" | "`kernel'"=="tri") {
      local kernel_type = "Triangular"
      qui g double `kwt'=max(0,`bwidthtab'-abs(`assignvar'))*`oweights'
  }
  else {
    local kernel_type = "Uniform"
    qui g double `kwt'=(-`bwidthtab'<=(`assignvar') & `assignvar'<`bwidthtab')*`oweights'
  }
       

// Create weight local for balance
if "`ipsw'" == "noipsw" local weight = ""
else local weight "[pw=`ipsweight']"


*** Count number of observations when noipsw is specified and there are not variables to balance:

if `: list sizeof balance'==0 {

 qui count if `touse' & `bwidth' & `sgroup'==0
 local N_G0 = `r(N)'
 qui count if `touse' & `bwidth' & `sgroup'==1
 local N_G1 = `r(N)'
 
scalar unw_N_G1 = `N_G1'
scalar unw_N_G0 = `N_G0' 


}


*-------------------------------------------------------------------------------
* Compute balance table matrices
*-------------------------------------------------------------------------------

* Original balance
*-------------------------------------------------------------------------------

if `: list sizeof balance'!=0 {

// Compute balanace matrix 
_wsga_rdd_balancematrix, matname(unw)  ///
 weights(`oweights') assignvar(`assignvar') touse(`touse') bwidth(`bwidth') balance(`balance') m(`m') ///
  sgroup(`sgroup') sgroup0(`sgroup0') n_balance(`n_balance') rbalance(`rbalance')
  

// Store balance matrix and computed balance stats
matrix unw = e(unw)
foreach s in unw_N_G0 unw_N_G1 unw_pvalue unw_Fstat unw_avgdiff {
  scalar `s' = e(`s')
}
// Display balance matrix and balance stats
if "`dibalance'" != "" {
  di _newline as result "Unweighted"
  matlist unw, ///
    border(rows) format(%9.3g) noblank
  di "Obs. in subgroup 0: " unw_N_G0
  di "Obs. in subgroup 1: " unw_N_G1
  di "Mean abs(std_diff): " unw_avgdiff
  di "F-statistic: " unw_Fstat
  di "Global p-value: " unw_pval_global
}


if "`ipsw'" != "noipsw" {

* Propensity Score Weighting balance
*-------------------------------------------------------------------------------
// Compute balanace matrix 

_wsga_rdd_balancematrix, matname(ipsw)  ///
  psw ipsweight(`ipsweight') weights(`oweights') touse(`touse') assignvar(`assignvar') bwidth(`bwidth') balance(`balance') m(`m') ///
  pscore(`pscore') comsup nocomsup(`NOCOMSUP') binarymodel(`binarymodel') ///
	sgroup(`sgroup') sgroup0(`sgroup0') n_balance(`n_balance')  rbalance(`rbalance') 
	
// Store balance matrix and computed balance stats
matrix ipsw = e(ipsw)
foreach s in ipsw_N_G0 ipswN_G1 ipsw_pvalue ipsw_Fstat ipsw_avgdiff {
  scalar `s' = e(`s')
}
// Display balance matrix and balance stats
if "`dibalance'" != "" {
  di _newline as result "Inverse Propensity Score Weighting"
  matlist ipsw, ///
  border(rows) format(%9.3g) noblank
  di "Obs. in subgroup 0: " ipsw_N_G0
  di "Obs. in subgroup 1: " ipsw_N_G1
  di "Mean abs(std_diff): " ipsw_avgdiff
  di "F-statistic: " ipsw_Fstat
  di "Global p-value: " ipsw_pval_global
}
}
}



*-------------------------------------------------------------------------------
* Estimation
*-------------------------------------------------------------------------------


tempvar kernelipsw


// Create weight local for regression
if "`ipsw'" == "noipsw" { 
local weight = "[pw=`kwt']"
}
else { 
qui g double `kernelipsw'=`ipsweight'*`kwt'
local weight "[pw=`kernelipsw']"
}




// Build psw options for myboo bootstrap calls
local psw_boot_opts
if "`ipsw'" != "noipsw" & `: list sizeof balance' > 0 {
  local psw_boot_opts psw ipsweight(`ipsweight') kernelipsw(`kernelipsw') ///
    kwt(`kwt') touse(`touse') bwidth(`bwidth') balance(`balance') ///
    oweights(`oweights') m(`m') binarymodel(`binarymodel') pscore(`pscore')
  if "`comsup'" != "" local psw_boot_opts `psw_boot_opts' comsup
}
if "`seed'" != "" local psw_boot_opts `psw_boot_opts' seed(`seed')

* First stage
*-------------------------------------------------------------------------------
if "`firststage'" != "" {
    mat IndIV=[0]
  // Regression
    qui reg `fuzzy' i.`sgroup'#1._cutoff i.`sgroup' ///
    i.`sgroup'#(`fv_covariates' c.`assignvar' c.`assignvar'#_cutoff)  `polynm' ///
    `weight' if `touse' & `bwidth', vce(`vce')
   
  ** Escalar used as a indicator to compute bootstrap. 
    mat fixed=[0] 

if "`blockbootstrap'"!="" {
ereturn local blockbootstrap `blockbootstrap'
}
    
  // Compute bootstrapped variance-covariance matrix and post results
  if "`bootstrap'" != "nobootstrap" _wsga_rdd_myboo `sgroup' _cutoff `bsreps', `psw_boot_opts'
  // If no bootstrap, trim b and V to show only RD estimates
  else _wsga_rdd_epost `sgroup' _cutoff 
}

* Reduced form
*-------------------------------------------------------------------------------
if "`reducedform'" != "" {
  // Regression
    mat IndIV=[0]
	
	** this is only used to report R2 with covariates in balance	
	qui reg `depvar' i.`sgroup'#1._cutoff i.`sgroup' ///
    i.`sgroup'#(`fv_balance' c.`assignvar' c.`assignvar'#_cutoff) `polynm' ///
    `weight' if `touse' & `bwidth', vce(`vce')

	local r2_covbal= e(r2)
	local r2_acovbal= e(r2_a)

    ** Main Regression
	qui reg `depvar' i.`sgroup'#1._cutoff i.`sgroup' ///
    i.`sgroup'#(`fv_covariates' c.`assignvar' c.`assignvar'#_cutoff) `polynm' ///
    `weight' if `touse' & `bwidth', vce(`vce')
	
	ereturn scalar r2_covbal=`r2_covbal'
	ereturn scalar r2_acovbal=`r2_acovbal'
  
  ** Escalar used as a indicator to compute bootstrap. 
    mat fixed=[0]
   
if "`blockbootstrap'"!="" {
ereturn local blockbootstrap `blockbootstrap'
}

  // Compute bootstrapped variance-covariance matrix and post results
  if "`bootstrap'" != "nobootstrap" _wsga_rdd_myboo `sgroup' _cutoff `bsreps', `psw_boot_opts'
  // If no bootstrap, trim b and V to show only RD estimates
  else _wsga_rdd_epost `sgroup' _cutoff 
}

* Instrumental variables
*-------------------------------------------------------------------------------
if "`ivregress'" != "" {
  // Regression

  mat IndIV=[1]
 
 qui reg `fuzzy' i.`sgroup'#1._cutoff i.`sgroup' ///
   i.`sgroup'#(`fv_covariates' c.`assignvar' c.`assignvar'#_cutoff) `polynm' ///
   `weight' if `touse' & `bwidth', vce(`vce')
   
   
 local coeffFSg0: di _b[0.`sgroup'#1._cutoff]
 local coeffFSg1: di _b[1.`sgroup'#1._cutoff]
 
 ** If RDD is fuzzy and we use a fixed bootstrap, so we save the first stage coefficient 
   mat FS=[`coeffFSg0',`coeffFSg1']

   qui reg `depvar' i.`sgroup'#1._cutoff i.`sgroup' ///
    i.`sgroup'#(`fv_covariates' c.`assignvar' c.`assignvar'#_cutoff) `polynm' ///
    `weight' if `touse' & `bwidth', vce(`vce')
	
local RFline `e(cmdline)'
  
 qui ivregress 2sls `depvar' i.`sgroup' ///
   i.`sgroup'#(`fv_covariates' c.`assignvar' c.`assignvar'#_cutoff) `polynm' ///
    (i.`sgroup'#1.`fuzzy' = i.`sgroup'#1._cutoff) ///
    `weight' if `touse' & `bwidth', vce(`vce')

** Escalar used as a indicator to compute bootstrap. 0 if bootstrap is computed using IV estimation
mat fixed=[0]

if "`blockbootstrap'" != "" {
ereturn local blockbootstrap `blockbootstrap'
}


if "`fixedbootstrap'" != "" {  

** If RDD is fuzzy and we use fixed bootstrap, so we use reduced form estimation to compute bootstrap in the IV
ereturn local cmdline `RFline'  

** Escalar used as a indicator to compute bootstrap. 1 if bootstrap is computed using reduced form estimation
** and keeping first stage fixed. 
 mat fixed=[1]
    }

  // Compute bootstrapped variance-covariance matrix and post results
  if "`bootstrap'" != "nobootstrap" _wsga_rdd_myboo `sgroup' `fuzzy' `bsreps', `psw_boot_opts'

  // If no bootstrap, trim b and V to show only RD estimates
  else _wsga_rdd_epost `sgroup' `fuzzy'
*  mat cumulative = e(cumulative)
*  ereturn matrix cumulative = cumulative

}


* Post balance results
*-------------------------------------------------------------------------------
// Post global balance stats
if `: list sizeof balance'!=0 & "`ipsw'" != "noipsw"   {

foreach w in unw ipsw {
  foreach s in N_G0 N_G1 pvalue Fstat avgdiff {
    ereturn scalar `w'_`s' = `w'_`s'
  }

}
// Post balance matrices
ereturn matrix ipsw ipsw
ereturn matrix unw unw

}
if "`ipsw'" == "noipsw" & `: list sizeof balance'!=0 {

local w unw
  foreach s in N_G0 N_G1 pvalue Fstat avgdiff {
    ereturn scalar `w'_`s' = `w'_`s'
  }

// Post balance matrices
ereturn matrix unw unw
}
if "`ipsw'" == "noipsw" & `: list sizeof balance'==0  {

ereturn scalar unw_N_G1 = `N_G1'
ereturn scalar unw_N_G0 = `N_G0'

}

*-------------------------------------------------------------------------------
* Results
*-------------------------------------------------------------------------------

* Post and display estimation results
*-------------------------------------------------------------------------------
if "`ivregress'" != "" | "`reducedform'" != "" | "`firststage'" != "" {

  // Post abridged b and V matrices
  ereturn repost b=b V=V, resize
  // Display estimates by subgroup
*  di as result "Subgroup estimates"
*  ereturn display
  // Display difference of subgroup estimates 
*  di _newline as result "Difference estimate"
  if "`ivregress'" == "" {
*   di as text "_nl_1 = _b[1.`sgroup'#1._cutoff] - _b[0.`sgroup'#1._cutoff]" _continue
    qui nlcom _b[1.`sgroup'#1._cutoff] - _b[0.`sgroup'#1._cutoff]
  }
  else {
*   di as text "_nl_1 = _b[1.`sgroup'#1.`fuzzy'] - _b[0.`sgroup'#1.`fuzzy']" _continue
    qui nlcom _b[1.`sgroup'#1.`fuzzy'] - _b[0.`sgroup'#1.`fuzzy']
  } 

  * Compute and store subgroup estimates 
  *-------------------------------------------------------------------------------
  if "`ivregress'" == "" scalar df = e(df_r)
  else scalar df = e(df_m)

// Estimates by subgroup
  forvalues g = 0/1 {
    // Coefficient
    matrix e_b = e(b)
    scalar b_g`g' = e_b[1,`=`g'+1']
    // Standard error
    matrix e_V = e(V)
    scalar se_g`g' = sqrt(e_V[`=`g'+1',`=`g'+1'])
    // t-stat
    scalar t_g`g' = b_g`g'/se_g`g'
    // P-value and confidence interval
    if "`bootstrap'" != "nobootstrap" {
      if "`normal'" != "" {
        // Normal approximation using bootstrap SE
        scalar p_g`g'     = 2*(1 - normal(abs(t_g`g')))
        scalar ci_lb_g`g' = b_g`g' - invnormal(0.975)*se_g`g'
        scalar ci_ub_g`g' = b_g`g' + invnormal(0.975)*se_g`g'
      }
      else {
        // Empirical (percentile) bootstrap CIs and (1+count)/(B+1) p-values
        scalar p_g`g'     = e(pval`g')
        scalar ci_lb_g`g' = e(lb_g`g')
        scalar ci_ub_g`g' = e(ub_g`g')
      }
    }
    else {
      scalar p_g`g'     = ttail(df, abs(t_g`g'))*2
      scalar ci_lb_g`g' = b_g`g' + invttail(df, 0.975)*se_g`g'
      scalar ci_ub_g`g' = b_g`g' + invttail(df, 0.025)*se_g`g'
    }
  }

  * Compute and store difference estimates 
  *-------------------------------------------------------------------------------
  // Coefficient
  matrix e_b_diff = r(b)
  scalar b_diff = e_b_diff[1,1]
  // Standard error
  matrix e_V_diff = r(V)
  scalar se_diff = sqrt(e_V_diff[1,1])
  // t-stat
  scalar t_diff = b_diff/se_diff
  // P-value and confidence interval
  if "`bootstrap'" != "nobootstrap" {
    if "`normal'" != "" {
      // Normal approximation using bootstrap SE
      scalar p_diff     = 2*(1 - normal(abs(t_diff)))
      scalar ci_lb_diff = b_diff - invnormal(0.975)*se_diff
      scalar ci_ub_diff = b_diff + invnormal(0.975)*se_diff
    }
    else {
      // Empirical (percentile) bootstrap CIs and (1+count)/(B+1) p-values
      scalar p_diff     = e(pval_diff)
      scalar ci_lb_diff = e(lb_diff)
      scalar ci_ub_diff = e(ub_diff)
    }
  }
  else {
    scalar p_diff     = ttail(df, abs(t_diff))*2
    scalar ci_lb_diff = b_diff + invttail(df, 0.975)*se_diff
    scalar ci_ub_diff = b_diff + invttail(df, 0.025)*se_diff
  }


  * Display estimation results
  *-------------------------------------------------------------------------------
  // Normal based
  if "`normal'" != "" {
    di as text "{hline 13}{c TT}{hline 64}"
    di as text %12s abbrev("`depvar'",12) " {c |}" ///
      _col(15) "{ralign 11:Coef.}" ///
      _col(26) "{ralign 12:Std. Err.}" ///
      _col(38) "{ralign 8:t }" /// notice extra space
      _col(46) "{ralign 8:P>|t|}" ///
      _col(54) "{ralign 25:[95% Conf. Interval]}" 
    di as text "{hline 13}{c +}{hline 64}"
    di as text "Subgroup" _col(14) "{c |}"
    forvalues g = 0/1 {
      display as text %12s abbrev("`g'",12) " {c |}" ///
        as result ///
        "  " %9.0g b_g`g' ///
        "  " %9.0g se_g`g' ///
        "    " %5.2f t_g`g' ///
        "   " %5.3f p_g`g' ///
        "    " %9.0g ci_lb_g`g' ///
        "   " %9.0g ci_ub_g`g'
    }
    di as text "{hline 13}{c +}{hline 64}"
    display as text "Difference   {c |}" ///
      as result ///
      "  " %9.0g b_diff ///
      "  " %9.0g se_diff ///
      "    " %5.2f t_diff ///
      "   " %5.3f p_diff ///
      "    " %9.0g ci_lb_diff ///
      "   " %9.0g ci_ub_diff
    di as text "{hline 13}{c BT}{hline 64}"
  }
  // Empirical 
  else {
    di as text "{hline 13}{c TT}{hline 64}"
    di as text %12s abbrev("`depvar'",12) " {c |}" ///
      _col(15) "{ralign 11:Coef.}" ///
      _col(26) "{ralign 12:Std. Err.}" ///
      _col(38) "{ralign 8:z }" /// notice extra space
      _col(46) "{ralign 8:P>|z|}" ///
      _col(58) "{ralign 25:[95% Conf. Interval] (P)}" 
    di as text "{hline 13}{c +}{hline 64}"
    di as text "Subgroup" _col(14) "{c |}"
    forvalues g = 0/1 {
      display as text %12s abbrev("`g'",12) " {c |}" ///
        as result ///
        "  " %9.0g b_g`g' ///
        "  " %9.0g se_g`g' ///
        "    " %5.2f t_g`g' ///
        "   " %5.3f p_g`g' ///
        "    " %9.0g ci_lb_g`g' ///
        "   " %9.0g ci_ub_g`g'
    }
    di as text "{hline 13}{c +}{hline 64}"
    display as text "Difference   {c |}" ///
      as result ///
      "  " %9.0g b_diff ///
      "  " %9.0g se_diff ///
      "    " %5.2f t_diff ///
      "   " %5.3f p_diff ///
      "    " %9.0g ci_lb_diff ///
      "   " %9.0g ci_ub_diff
    di as text "{hline 13}{c BT}{hline 64}"
  }
}

* End
*-------------------------------------------------------------------------------
cap drop _cutoff
end

*===============================================================================
* Define auxiliary subroutines
*===============================================================================

*-------------------------------------------------------------------------------
* epost: post matrices in e(b) and e(V); leave other ereturn results unchanged
*-------------------------------------------------------------------------------
program _wsga_rdd_epost, eclass
  // Store results: scalars
  local scalars: e(scalars)
  foreach scalar of local scalars {
    local `scalar' = e(`scalar')
  }
  // Store results: macros
  local macros: e(macros)
  foreach macro of local macros {
    local `macro' = e(`macro')
  }
  // Store results: matrices (drop V_modelbased; b and V are computed below)
  local matrices: e(matrices)
  // Store results: functions
  tempvar esample
  gen `esample' = e(sample)
  // b and V matrices
  matrix b = e(b)
  matrix V = e(V)
  matrix b = b[1, "0.`1'#1.`2'".."1.`1'#1.`2'"]
  matrix V = V["0.`1'#1.`2'".."1.`1'#1.`2'", "0.`1'#1.`2'".."1.`1'#1.`2'"]
  ereturn post, esample(`esample')
  // Post results: scalars
  foreach scalar of local scalars {
    ereturn scalar `scalar' = ``scalar''
  }
  // Post results: macros
  foreach macro of local macros {
    ereturn local `macro' ``macro''
  }
end

*-------------------------------------------------------------------------------
* myboo: compute bootstrapped variance-covariance matrix & adjust ereturn results
*-------------------------------------------------------------------------------
program define _wsga_rdd_myboo, eclass
  syntax anything [, PSW IPSWeight(name) KERNELipsw(name) KWT(name) ///
    TOuse(name) BWIDTH(string) BALance(varlist) OWeights(name) ///
    M(integer 2) BINarymodel(string) PSCore(name) COmsup Seed(string)]

  local svar : word 1 of `anything'
  local cvar : word 2 of `anything'
  local B    : word 3 of `anything'

  // Store results: scalars
  local scalars: e(scalars)
  foreach scalar of local scalars {
    local `scalar' = e(`scalar')
  }
  // Store results: macros
  local macros: e(macros)
  foreach macro of local macros {
    local `macro' = e(`macro')
  }

  // Store results: functions
  tempvar esample
  gen `esample' = e(sample)
  // Extract b submatrix with subgroup coefficients
  matrix b = e(b)
  matrix b = b[1, "0.`svar'#1.`cvar'".."1.`svar'#1.`cvar'"]
  matrix colnames b = 0.`svar'#1.`cvar' 1.`svar'#1.`cvar'

  // Start bootstrap
  di ""
  _dots 0, title(Bootstrap replications) reps(`B')
  cap mat drop cumulative
  if "`seed'" != "" set seed `seed'
  tempvar COMSUP_b
  forvalues i=1/`B' {
    preserve
    bsample, strata(`e(blockbootstrap)')

    // Re-estimate propensity score on bootstrap sample
    if "`psw'" != "" {
      // Re-fit binary model and predict new pscore
      qui drop `pscore'
      qui `binarymodel' `svar' `balance' [pw=`oweights'] if `touse' & `bwidth'
      qui predict double `pscore' if `touse' & `bwidth' & !mi(`svar')
      // Re-evaluate common support
      cap drop `COMSUP_b'
      if "`comsup'" != "" {
        qui sum `pscore' if `svar' == 1 & `touse' & `bwidth'
        qui gen byte `COMSUP_b' = (`pscore' >= r(min) & `pscore' <= r(max)) ///
          if `touse' & `bwidth' & !mi(`svar')
      }
      else {
        qui gen byte `COMSUP_b' = 1 if `touse' & `bwidth' & !mi(`svar')
      }
      // Recompute N_G0, N_G1 on bootstrap sample
      if `m' == 2 {
        qui count if `touse' & `bwidth' & `COMSUP_b' & `svar'==0 & !mi(`pscore')
        local N_G0_b = r(N)
        qui count if `touse' & `bwidth' & `COMSUP_b' & `svar'==1 & !mi(`pscore')
        local N_G1_b = r(N)
      }
      else if `m' == 1 {
        qui count if `touse' & `bwidth' & `COMSUP_b' & `svar'==0
        local N_G0_b = r(N)
        qui count if `touse' & `bwidth' & `COMSUP_b' & `svar'==1 & !mi(`pscore')
        local N_G1_b = r(N)
      }
      else {
        qui count if `touse' & `bwidth' & `COMSUP_b' & `svar'==0 & !mi(`pscore')
        local N_G0_b = r(N)
        qui count if `touse' & `bwidth' & `COMSUP_b' & `svar'==1
        local N_G1_b = r(N)
      }
      // Recompute ipsweight: clear first, then set per mode
      qui replace `ipsweight' = . if `touse' & `bwidth'
      if `m' == 2 {
        qui replace `ipsweight' = ///
          (`N_G1_b'/(`N_G1_b'+`N_G0_b')/`pscore'*(`svar'==1) + ///
           `N_G0_b'/(`N_G1_b'+`N_G0_b')/(1-`pscore')*(`svar'==0)) ///
          if `touse' & `bwidth' & `COMSUP_b' & !mi(`svar')
      }
      else if `m' == 1 {
        qui replace `ipsweight' = (1-`pscore')/`pscore' ///
          if `touse' & `bwidth' & `COMSUP_b' & `svar'==1
        qui replace `ipsweight' = 1 if `touse' & `bwidth' & `svar'==0
      }
      else {
        qui replace `ipsweight' = `pscore'/(1-`pscore') ///
          if `touse' & `bwidth' & `COMSUP_b' & `svar'==0
        qui replace `ipsweight' = 1 if `touse' & `bwidth' & `svar'==1
      }
      // Recompute kernelipsw = ipsweight * kwt
      qui replace `kernelipsw' = `ipsweight' * `kwt'
    }

    qui `e(cmdline)'
    tempname this_run
    // Non-IV or bootstrap-both-stages
    if IndIV[1,1]==0 | fixed[1,1]==0 {
      local b_g0_i = _b[0.`svar'#1.`cvar']
      local b_g1_i = _b[1.`svar'#1.`cvar']
      mat `this_run' = (`b_g0_i', `b_g1_i', `b_g1_i' - `b_g0_i')
    }
    // IV with fixed first stage: bootstrap reduced form, divide by saved FS coefs
    if IndIV[1,1]==1 & fixed[1,1]==1 {
      local CoeffIVg0_`i' = _b[0.`svar'#1._cutoff]/FS[1,1]
      local CoeffIVg1_`i' = _b[1.`svar'#1._cutoff]/FS[1,2]
      mat `this_run' = (`CoeffIVg0_`i'', `CoeffIVg1_`i'', ///
        `CoeffIVg1_`i'' - `CoeffIVg0_`i'')
    }
    mat cumulative = nullmat(cumulative) \ `this_run'
    restore
    _dots `i' 0
  }

  di _newline
  // Compute 2x2 VCV from first two columns (g0, g1)
  cap mat drop V
  mata: cumulative = st_matrix("cumulative")
  mata: st_matrix("V", variance(cumulative[., 1..2]))
  // Add names to V
  mat rownames V = 0.`svar'#1.`cvar' 1.`svar'#1.`cvar'
  mat colnames V = 0.`svar'#1.`cvar' 1.`svar'#1.`cvar'
  // Return
  ereturn post, esample(`esample')
  // Post results: scalars
  foreach scalar of local scalars {
    ereturn scalar `scalar' = ``scalar''
  }
  ereturn scalar N_reps = `B'
  ereturn scalar level = 95
  // Empirical p-values for subgroups
  // Recentered: count draws where |draw - est| >= |est|, testing H0: coef=0
  cap scalar drop bscoef
  forvalues g = 0/1 {
    local count = 0
    forvalues i = 1/`B' {
      scalar bscoef = cumulative[`i',`=`g'+1']
      if abs(bscoef - b[1,`=`g'+1']) >= abs(b[1,`=`g'+1']) local count = `count'+1
    }
    scalar pval`g' = (1+`count') / (`B' + 1)
    ereturn scalar pval`g' = pval`g'
  }
  // Empirical p-value for diff (column 3)
  scalar orig_diff = b[1,2] - b[1,1]
  local count_diff = 0
  forvalues i = 1/`B' {
    scalar bscoef = cumulative[`i',3]
    if abs(bscoef - orig_diff) >= abs(orig_diff) local count_diff = `count_diff' + 1
  }
  scalar pval_diff = (1 + `count_diff') / (`B' + 1)
  ereturn scalar pval_diff = pval_diff
  // Empirical confidence intervals
  svmat cumulative, names(_subgroup)
  forvalues g = 0/1 {
    qui centile _subgroup`=`g'+1', centile(2.5 97.5)
    drop _subgroup`=`g'+1'
    scalar lb_g`g' = r(c_1)
    ereturn scalar lb_g`g' = lb_g`g'
    scalar ub_g`g' = r(c_2)
    ereturn scalar ub_g`g' = ub_g`g'
  }
  // Empirical CI for diff (column 3)
  qui centile _subgroup3, centile(2.5 97.5)
  drop _subgroup3
  scalar lb_diff = r(c_1)
  ereturn scalar lb_diff = lb_diff
  scalar ub_diff = r(c_2)
  ereturn scalar ub_diff = ub_diff
  // Post results: macros
  foreach macro of local macros {
    if "`macro'" == "clustvar" continue
    ereturn local `macro' ``macro''
  }
  ereturn local vcetype "Bootstrap"
  ereturn local vce "bootstrap"
  ereturn local prefix "bootstrap"
end


*-------------------------------------------------------------------------------
* _wsga_rdd_balancematrix: compute balance table matrices and other statistics
*-------------------------------------------------------------------------------
program define _wsga_rdd_balancematrix, eclass
syntax, matname(string) /// important inputs, differ by call
  touse(name)  weights(string) bwidth(string) balance(varlist) m(int) /// unchanging inputs
  [psw ipsweight(name)  pscore(name) comsup nocomsup(name) binarymodel(string)] /// only needed for PSW balance
  sgroup(name) sgroup0(name) n_balance(int) rbalance(int) assignvar(name) // todo: eliminate these? can be computed by subroutine at low cost


* Create variables specific to PSW matrix
*-------------------------------------------------------------------------------
if "`psw'" != "" { // if psw
  // Fit binary response model
 qui cap drop comsup
 qui `binarymodel' `sgroup' `balance' [pw=`weights'] if `touse' & `bwidth'


  // Generate pscore variable and clear stored results
  qui predict double `pscore' if `touse' & `bwidth' & !mi(`sgroup')
  ereturn clear

  // No compute common support area as default (create a aux variable)
if "`nocomsup'" != "" {
  tempvar COMSUP
  qui gen `COMSUP' = 1 if `touse' & `bwidth' & !mi(`sgroup')
  qui gen comsup=`COMSUP'
}  
else {
qui sum `pscore' if `sgroup' == 1 /* todo: check why this is like that */
tempvar COMSUP
    qui gen `COMSUP' = ///
      (`pscore' >= r(min) & ///
       `pscore' <= r(max))
     
   label var `COMSUP' "Dummy for obs. in common support"
	
   qui g comsup = `COMSUP'

   label var comsup "Dummy for obs. in common support"
}
  // Count observations in each fuzzy group
  
	if `m'==2 {
  qui count if `touse' & `bwidth' & `COMSUP' & `sgroup'==0 & !mi(`pscore')
  local N_G0 = `r(N)'
  qui count if `touse' & `bwidth' & `COMSUP' & `sgroup'==1 & !mi(`pscore')
  local N_G1 = `r(N)'
	}
	if `m'==1 {
  qui count if `touse' & `bwidth' & `COMSUP' & `sgroup'==0
  local N_G0 = `r(N)'
  qui count if `touse' & `bwidth' & `COMSUP' & `sgroup'==1 & !mi(`pscore')
  local N_G1 = `r(N)'
	}	
	if `m'==0 {
  qui count if `touse' & `bwidth' & `COMSUP' & `sgroup'==0 & !mi(`pscore')
  local N_G0 = `r(N)'
  qui count if `touse' & `bwidth' & `COMSUP' & `sgroup'==1 
  local N_G1 = `r(N)'
	}	
  	
 
  
  // Compute propensity score weighting vector
  cap drop `ipsweight'
	if `m'==2 {
  qui gen `ipsweight' = ///
   ( `N_G1'/(`N_G1'+`N_G0')/`pscore'*(`sgroup'==1) + ///
    `N_G0'/(`N_G1'+`N_G0')/(1-`pscore')*(`sgroup'==0)) ///
    if `touse' & `bwidth' & `COMSUP' & !mi(`sgroup')
	}
	if `m'==1 {
  qui gen `ipsweight' = ///
   (1-`pscore')/`pscore' ///
    if `touse' & `bwidth' & `COMSUP' & `sgroup'==1 	
	qui replace `ipsweight'=1 ///
	if `touse' & `bwidth' & `sgroup'==0 	
	}
if `m'==0 {
  qui gen `ipsweight' = ///
   `pscore'/(1-`pscore') ///
    if `touse' & `bwidth' & `COMSUP' & `sgroup'==0
	qui replace `ipsweight'=1 ///
	if `touse' & `bwidth' & `sgroup'==1 

	}
	

	tempvar nweights
	qui gen `nweights'=`ipsweight'*`weights'	
    
} // end if psw


* Count obs. in each fuzzy group if not PSW matrix
*-------------------------------------------------------------------------------
else { // if nopsw
  qui count if `touse' & `bwidth' & `sgroup'==0
  local N_G0 = `r(N)'
  qui count if `touse' & `bwidth' & `sgroup'==1
  local N_G1 = `r(N)'
} // end if nopsw

* Compute stats specific for each covariate 
*-------------------------------------------------------------------------------
local j = 0
foreach var of varlist `balance' {
  local ++j
  // Compute and store conditional expectations
  if `rbalance'==1 {
  if "`psw'" == ""   qui reg `var' `sgroup0' `sgroup' [iw=`weights'] if `touse' & `bwidth', noconstant /* */
  else  qui reg `var' `sgroup0' `sgroup' [iw=`nweights'] if `touse' & `bwidth' & `COMSUP', noconstant
  local coef`j'_G0 = _b[`sgroup0']
  local coef`j'_G1 = _b[`sgroup']
  
  // Compute and store mean differences and their p-values
  if "`psw'" == "" qui reg `var' `sgroup0' [iw=`weights'] if `touse' & `bwidth'
  else qui reg `var' `sgroup0' [iw=`nweights'] if `touse' & `bwidth' & `COMSUP'
  matrix m = r(table)
  scalar diff`j'=m[1,1] // mean difference
  local pval`j' = m[4,1] // p-value 

  // Standardized mean difference
  if "`psw'" == "" qui summ `var' if `touse' & `bwidth' & !mi(`sgroup')
  else qui summ `var' if `touse' & `bwidth' & `COMSUP' & !mi(`sgroup')
  local stddiff`j' = (diff`j')/r(sd)
  }
  
    if `rbalance'==0 {
 
  
  // Compute and store mean differences and their p-values
 if "`psw'" == ""  { 
    qui reg `var' `assignvar' [iw=`weights'] if `touse' & `bwidth' & `sgroup0'
	qui estimates store model1
    local coef`j'_G0 = _b[_cons]
	local se`j'_G0 =_se[_cons]
	qui reg `var' `assignvar'  [iw=`weights'] if `touse' & `bwidth' &  `sgroup'  
	qui estimates store model2
	  local coef`j'_G1 =  _b[_cons]
	local se`j'_G1 =_se[_cons]
	
	qui reg `var' c.`assignvar'##`sgroup0'   [iw=`weights'] if `touse' & `bwidth'  
matrix m=r(table)

  scalar diff`j'=m[1,3] // mean difference
 
  local pval`j' = m[4,3] // p-value 
  }
  else  { 
    qui reg `var' `assignvar' [iw=`nweights'] if `touse' & `bwidth' & `sgroup0'
	qui estimates store model1
    local coef`j'_G0 = _b[_cons]

	qui reg `var' `assignvar'  [iw=`nweights'] if `touse' & `bwidth' &  `sgroup'  
	qui estimates store model2
	  local coef`j'_G1 =  _b[_cons]
  	qui reg `var' c.`assignvar'##`sgroup0'   [iw=`nweights'] if `touse' & `bwidth'  
matrix m=r(table)

  scalar diff`j'=m[1,3] // mean difference
 
  local pval`j' = m[4,3] // p-value 
  }

  // Standardized mean difference
  if "`psw'" == "" qui summ `var' if `touse' & `bwidth' & !mi(`sgroup')
  else qui summ `var' if `touse' & `bwidth' & `COMSUP' & !mi(`sgroup')
  local stddiff`j' = (diff`j')/r(sd)
  }
}



* Compute global stats
*-------------------------------------------------------------------------------
// Mean of absolute standardized mean differences (ie. stddiff + ... + stddiff`k')
/* todo: this begs to be vectorized */

local avgdiff = 0
forvalues j = 1/`n_balance' {
  local avgdiff = abs(`stddiff`j'') + `avgdiff' // sum over `j' (balance)
}
local avgdiff = `avgdiff'/`n_balance' // compute mean 

// F-statistic and global p-value
if  `rbalance'==1 {
if "`psw'" == "" qui reg `sgroup' `balance' [iw=`weights']  if `touse' & `bwidth'
else qui reg `sgroup' `balance' [iw=`nweights'] if `touse' & `bwidth' & `COMSUP' 
local Fstat = e(F)
local pval_global = 1-F(e(df_m),e(df_r),e(F))
}


if  `rbalance'==0 {
local j = 0 
foreach var of varlist `balance' {
local j=`j'+1
qui gen aux`j'=`var' *`assignvar'
}

if "`psw'" == "" {

  qui reg `sgroup' aux* `assignvar' [iw=`weights']  if `touse' & `bwidth', noconstant
  predict res_aux, res
  qui reg res_aux `balance' [iw=`weights']  if `touse' & `bwidth' 
  cap drop res_aux
 }
else  { 
  qui reg `sgroup' aux* `assignvar'  [iw=`nweights']  if `touse' & `bwidth', noconstant
  predict res_aux, res
   qui reg res_aux `balance' [iw=`nweights']  if `touse' & `bwidth' 
     cap drop res_aux
  }
local Fstat = e(F)
local pval_global = 1-F(e(df_m),e(df_r),e(F))
cap drop aux*
}





* Create balance matrix
*-------------------------------------------------------------------------------
// Matrix parameters
matrix `matname' = J(`n_balance', 4, .)
matrix colnames `matname' = mean_G0 mean_G1 std_diff p-value
matrix rownames `matname' = `balance'

// Add per-covariate values 
forvalues j = 1/`n_balance' {
  matrix `matname'[`j',1] = `coef`j'_G0'
  matrix `matname'[`j',2] = `coef`j'_G1'
  matrix `matname'[`j',3] = `stddiff`j''
  matrix `matname'[`j',4] = `pval`j''
}



// Return matrix and other scalars
scalar `matname'_N_G0 = `N_G0'
scalar `matname'_N_G1 = `N_G1'
scalar `matname'_avgdiff = `avgdiff'
scalar `matname'_Fstat = `Fstat'
scalar `matname'_pval_global = `pval_global'

ereturn matrix `matname' = `matname', copy
ereturn scalar `matname'_avgdiff = `avgdiff'
ereturn scalar `matname'_Fstat = `Fstat'
ereturn scalar `matname'_pvalue = `pval_global'
ereturn scalar `matname'_N_G1 = `N_G1'
ereturn scalar `matname'_N_G0 = `N_G0'




end

********************************************************************************
* _wsga_did: sharp 2-period DiD pipeline
********************************************************************************

program define _wsga_did, eclass
version 11.1
syntax varlist(min=1 numeric fv) [if] [in], ///
  SGroup(name) UNIt(varname) TIMe(varname) TReat(varname) ///
  [ POST_value(string) ///
    BALance(varlist numeric) DIBALance probit ///
    IPSWeight(name) PSCore(name) comsup ///
    vce(string) m(int 2) ///
    noBOOTstrap bsreps(real 200) FIXEDbootstrap FIXEDps BLOCKbootstrap(string) ///
    NORMal noipsw weights(string) Seed(string) ]

  // ── Sample mask
  marksample touse, novarlist
  markout `touse' `unit' `time' `treat' `sgroup'

  // ── Outcome and covariates
  local depvar : word 1 of `varlist'
  local covariates : list varlist - depvar

  // ipsweight()/pscore(): user-named output variables or tempvars (mirrors RDD)
  if "`ipsweight'" != "" confirm new variable `ipsweight'
  else tempvar ipsweight
  if "`pscore'" != "" confirm new variable `pscore'
  else tempvar pscore

  // ── Validate: time has exactly 2 unique non-missing values
  qui levelsof `time' if `touse', local(t_vals)
  local n_tvals : word count `t_vals'
  if `n_tvals' != 2 {
    di as error ///
"design(did) requires `time' to have exactly 2 unique non-missing values; found `n_tvals'. Staggered adoption is not currently supported."
    exit 198
  }

  // ── Resolve post_value
  if "`post_value'" == "" {
    qui summarize `time' if `touse', meanonly
    local post_value = r(max)
  }
  else {
    local pv_ok = 0
    foreach tv of local t_vals {
      if "`tv'" == "`post_value'" local pv_ok = 1
    }
    if `pv_ok' == 0 {
      di as error ///
"post_value(`post_value') is not one of the unique values of `time' (`t_vals')."
      exit 198
    }
  }

  // ── Unit-constancy checks (treat, sgroup, balance moderators, blockbootstrap)
  foreach v in `treat' `sgroup' `balance' `blockbootstrap' {
    tempvar _dist
    qui by `unit' (`v'), sort: gen byte `_dist' = (`v'[1] != `v'[_N])
    qui summarize `_dist' if `touse', meanonly
    if r(max) > 0 {
      di as error ///
"`v' varies within `unit' (the DiD design requires `v' to be unit-level)."
      exit 198
    }
    drop `_dist'
  }

  // ── Build design variables
  tempvar G0 G1 post G0_Z G1_Z G0_post G1_post
  qui gen byte `G0'      = (`sgroup' == 0)
  qui gen byte `G1'      = (`sgroup' == 1)
  qui gen byte `post'    = (`time' == `post_value')
  qui gen byte `G0_Z'    = `G0' * `treat' * `post'    // delta_0
  qui gen byte `G1_Z'    = `G1' * `treat' * `post'    // delta_1
  qui gen byte `G0_post' = `G0' * `post'
  qui gen byte `G1_post' = `G1' * `post'

  // ── Covariate × G interactions
  local g0_covs ""
  local g1_covs ""
  foreach v of local covariates {
    tempvar g0_`v' g1_`v'
    qui gen double `g0_`v'' = `G0' * `v'
    qui gen double `g1_`v'' = `G1' * `v'
    local g0_covs "`g0_covs' `g0_`v''"
    local g1_covs "`g1_covs' `g1_`v''"
  }

  // ── IPW: fit logit/probit on (G, balance) and compute weights.
  // Because G and the moderators are unit-constant (validated above), fitting
  // on the long panel produces unit-constant pscore values — the per-row
  // values are correct as-is and no broadcast is required.
  // When comsup is specified, units outside the G=1 pscore range are excluded.
  qui gen double `pscore'    = .
  qui gen double `ipsweight' = 1

  if "`ipsw'" != "noipsw" & `: list sizeof balance' > 0 {
    if "`probit'" != "" {
      qui probit `sgroup' `balance' if `touse'
    }
    else {
      qui logit `sgroup' `balance' if `touse'
    }
    qui predict double _wsga_ps if `touse', pr
    qui replace `pscore' = _wsga_ps
    qui drop _wsga_ps

    // Unit-level counts for normalization (count one row per unit)
    tempvar unit_first
    qui by `unit' (`touse'), sort: gen byte `unit_first' = (_n == 1) & `touse'

    // Common support restriction: use G=1 pscore range (mirrors RDD convention)
    tempvar _comsup
    if "`comsup'" != "" {
      qui sum `pscore' if `sgroup' == 1 & `unit_first' & `touse', meanonly
      cap drop comsup
      qui gen byte comsup    = (`pscore' >= r(min) & `pscore' <= r(max)) if `touse' & !mi(`sgroup')
      qui gen byte `_comsup' = comsup
    }
    else {
      qui gen byte `_comsup' = 1 if `touse' & !mi(`sgroup')
    }

    qui count if `unit_first' & `sgroup' == 0 & !mi(`pscore') & `_comsup'
    local N_G0 = r(N)
    qui count if `unit_first' & `sgroup' == 1 & !mi(`pscore') & `_comsup'
    local N_G1 = r(N)
    local N_tot = `N_G0' + `N_G1'

    // m=2: balance both groups toward pooled (the default and paper-preferred)
    // m=1: balance group 1 toward group 0  (ATT for G=0)
    // m=0: balance group 0 toward group 1  (ATT for G=1)
    qui replace `ipsweight' = .
    if `m' == 2 {
      qui replace `ipsweight' = (`N_G1' / `N_tot') / `pscore'         if `sgroup' == 1 & !mi(`pscore') & `_comsup'
      qui replace `ipsweight' = (`N_G0' / `N_tot') / (1 - `pscore')   if `sgroup' == 0 & !mi(`pscore') & `_comsup'
    }
    else if `m' == 1 {
      qui replace `ipsweight' = (1 - `pscore') / `pscore'             if `sgroup' == 1 & !mi(`pscore') & `_comsup'
      qui replace `ipsweight' = 1                                      if `sgroup' == 0
    }
    else if `m' == 0 {
      qui replace `ipsweight' = `pscore' / (1 - `pscore')             if `sgroup' == 0 & !mi(`pscore') & `_comsup'
      qui replace `ipsweight' = 1                                      if `sgroup' == 1
    }
    else {
      di as error "m() must be 0, 1, or 2 (got `m')."
      exit 198
    }
    // Replace NAs in ipsweight with 0 so they drop out of the regression
    qui replace `ipsweight' = 0 if mi(`ipsweight')
  }

  // ── Estimate: long-form TWFE with unit FE absorbed via xtreg, fe.
  // Weights enter as pweights so SE machinery is consistent with sampling-
  // weight semantics.  When noipsw / no balance, ipsweight is identically 1.
  qui xtset `unit'
  local rhs `G0_Z' `G1_Z' `G0_post' `G1_post' `g0_covs' `g1_covs'
  xtreg `depvar' `rhs' [pw=`ipsweight'] if `touse' & `ipsweight' > 0, fe vce(cluster `unit')
  tempvar _did_esample
  gen byte `_did_esample' = e(sample)

  // ── Extract coefficients of interest
  scalar b_g0   = _b[`G0_Z']
  scalar b_g1   = _b[`G1_Z']
  scalar b_diff = b_g1 - b_g0
  scalar se_g0  = _se[`G0_Z']
  scalar se_g1  = _se[`G1_Z']
  // SE of the difference via the post-estimation vcov
  matrix V = e(V)
  local idx_g0 = colnumb(V, "`G0_Z'")
  local idx_g1 = colnumb(V, "`G1_Z'")
  scalar cov_01  = V[`idx_g0', `idx_g1']
  scalar se_diff = sqrt(se_g0^2 + se_g1^2 - 2*cov_01)

  scalar t_g0   = b_g0   / se_g0
  scalar t_g1   = b_g1   / se_g1
  scalar t_diff = b_diff / se_diff

  // df for analytical inference: from xtreg's degrees of freedom
  scalar df = e(df_r)
  scalar p_g0   = 2*ttail(df, abs(t_g0))
  scalar p_g1   = 2*ttail(df, abs(t_g1))
  scalar p_diff = 2*ttail(df, abs(t_diff))
  scalar ci_lb_g0   = b_g0   + invttail(df, 0.975)*se_g0
  scalar ci_ub_g0   = b_g0   + invttail(df, 0.025)*se_g0
  scalar ci_lb_g1   = b_g1   + invttail(df, 0.975)*se_g1
  scalar ci_ub_g1   = b_g1   + invttail(df, 0.025)*se_g1
  scalar ci_lb_diff = b_diff + invttail(df, 0.975)*se_diff
  scalar ci_ub_diff = b_diff + invttail(df, 0.025)*se_diff

  // ── Cluster bootstrap (refit PS per rep; Q1 + Q3 + Q4).
  // Whole units are resampled with replacement; fresh unit IDs are assigned
  // by bsample's idcluster() option so unit FE remain identified across
  // duplicate draws (the Cameron-Gelbach-Miller recipe).  The analytical
  // SE/CI/p scalars set above are overridden with bootstrap-based values.
  scalar B_ok    = 0
  scalar N_clust = .
  scalar use_bootstrap = ("`bootstrap'" != "nobootstrap")
  if use_bootstrap {
    qui levelsof `unit' if `touse', local(_clusters)
    scalar N_clust = `: word count `_clusters''
    if N_clust < 30 {
      di as text ///
"Note: cluster bootstrap with " as result %4.0f N_clust as text " clusters.  Pairs-cluster bootstrap can substantially understate SEs at < ~30 clusters.  Consider supplementing with a wild cluster bootstrap-t (Stata: boottest)."
    }
    tempfile _wsga_did_panel
    qui save `_wsga_did_panel'

    tempname _wsga_did_draws
    matrix `_wsga_did_draws' = J(`bsreps', 2, .)

    display as text "Cluster bootstrap (`bsreps' reps):"
    if "`seed'" != "" set seed `seed'
    forvalues _b = 1/`bsreps' {
      qui use `_wsga_did_panel', clear
      capture {
        if "`blockbootstrap'" != "" {
          qui bsample, cluster(`unit') idcluster(_wsga_new_unit) strata(`blockbootstrap')
        }
        else {
          qui bsample, cluster(`unit') idcluster(_wsga_new_unit)
        }

        tempvar _b_pscore _b_ipw
        qui gen double `_b_pscore' = .
        qui gen double `_b_ipw'    = 1
        if "`ipsw'" != "noipsw" & `: list sizeof balance' > 0 {
          if "`fixedps'" != "" {
            // fixedps: carry original-sample pscore; skip logit/probit refit
            qui replace `_b_pscore' = `pscore'
          }
          else {
            if "`probit'" != "" qui probit `sgroup' `balance'
            else                qui logit  `sgroup' `balance'
            qui predict double _wsga_ps_b, pr
            qui replace `_b_pscore' = _wsga_ps_b
            qui drop _wsga_ps_b
          }

          // Common support per bootstrap rep (mirrors main-sample logic)
          tempvar _b_comsup _b_uf
          if "`comsup'" != "" {
            qui by _wsga_new_unit (`sgroup'), sort: gen byte `_b_uf' = (_n == 1)
            qui sum `_b_pscore' if `sgroup' == 1 & `_b_uf' & !mi(`_b_pscore'), meanonly
            qui gen byte `_b_comsup' = (`_b_pscore' >= r(min) & `_b_pscore' <= r(max)) if !mi(`_b_pscore')
          }
          else {
            qui gen byte `_b_comsup' = 1 if !mi(`_b_pscore')
            qui by _wsga_new_unit (`sgroup'), sort: gen byte `_b_uf' = (_n == 1)
          }

          // IPW normalization: always recount units per bootstrap rep
          qui count if `_b_uf' & `sgroup' == 0 & !mi(`_b_pscore') & `_b_comsup'
          local _N0 = r(N)
          qui count if `_b_uf' & `sgroup' == 1 & !mi(`_b_pscore') & `_b_comsup'
          local _N1 = r(N)
          local _NT = `_N0' + `_N1'

          qui replace `_b_ipw' = .
          if `m' == 2 {
            qui replace `_b_ipw' = (`_N1' / `_NT') / `_b_pscore'         if `sgroup' == 1 & !mi(`_b_pscore') & `_b_comsup'
            qui replace `_b_ipw' = (`_N0' / `_NT') / (1 - `_b_pscore')   if `sgroup' == 0 & !mi(`_b_pscore') & `_b_comsup'
          }
          else if `m' == 1 {
            qui replace `_b_ipw' = (1 - `_b_pscore') / `_b_pscore'        if `sgroup' == 1 & !mi(`_b_pscore') & `_b_comsup'
            qui replace `_b_ipw' = 1                                       if `sgroup' == 0
          }
          else {
            qui replace `_b_ipw' = `_b_pscore' / (1 - `_b_pscore')        if `sgroup' == 0 & !mi(`_b_pscore') & `_b_comsup'
            qui replace `_b_ipw' = 1                                       if `sgroup' == 1
          }
          qui replace `_b_ipw' = 0 if mi(`_b_ipw')
        }

        qui xtset _wsga_new_unit
        qui xtreg `depvar' `rhs' [pw=`_b_ipw'] if `_b_ipw' > 0, ///
          fe vce(cluster _wsga_new_unit)
        matrix `_wsga_did_draws'[`_b', 1] = _b[`G0_Z']
        matrix `_wsga_did_draws'[`_b', 2] = _b[`G1_Z']
        scalar B_ok = B_ok + 1
      }
      if mod(`_b', 10) == 0 | `_b' == `bsreps' di "  " %4.0f `_b' "/`bsreps'"
    }
    qui use `_wsga_did_panel', clear

    if B_ok < 2 {
      di as error "Cluster bootstrap produced fewer than 2 successful reps (got " B_ok ").  Try increasing bsreps or relaxing balance() / m()."
      exit 498
    }

    // Aggregate: SEs, percentile CIs, (1+count)/(B+1) p-values
    preserve
    qui drop _all
    qui svmat double `_wsga_did_draws', names(_draw)
    qui keep if !mi(_draw1) & !mi(_draw2)
    qui gen double _drawdiff = _draw2 - _draw1

    qui summarize _draw1
    scalar se_g0   = r(sd)
    qui summarize _draw2
    scalar se_g1   = r(sd)
    qui summarize _drawdiff
    scalar se_diff = r(sd)
    qui correlate _draw1 _draw2, covariance
    matrix _wsga_cov = r(C)
    scalar cov_01 = _wsga_cov[1,2]
    matrix drop _wsga_cov
    scalar t_g0   = b_g0   / se_g0
    scalar t_g1   = b_g1   / se_g1
    scalar t_diff = b_diff / se_diff

    if "`normal'" != "" {
      // Normal approximation with bootstrap SE
      scalar p_g0      = 2*(1 - normal(abs(t_g0)))
      scalar p_g1      = 2*(1 - normal(abs(t_g1)))
      scalar p_diff    = 2*(1 - normal(abs(t_diff)))
      scalar ci_lb_g0   = b_g0   - invnormal(0.975)*se_g0
      scalar ci_ub_g0   = b_g0   + invnormal(0.975)*se_g0
      scalar ci_lb_g1   = b_g1   - invnormal(0.975)*se_g1
      scalar ci_ub_g1   = b_g1   + invnormal(0.975)*se_g1
      scalar ci_lb_diff = b_diff - invnormal(0.975)*se_diff
      scalar ci_ub_diff = b_diff + invnormal(0.975)*se_diff
    }
    else {
      // Empirical percentile CIs and (1+count)/(B+1) p-values
      qui _pctile _draw1,    percentiles(2.5 97.5)
      scalar ci_lb_g0    = r(r1)
      scalar ci_ub_g0    = r(r2)
      qui _pctile _draw2,    percentiles(2.5 97.5)
      scalar ci_lb_g1    = r(r1)
      scalar ci_ub_g1    = r(r2)
      qui _pctile _drawdiff, percentiles(2.5 97.5)
      scalar ci_lb_diff  = r(r1)
      scalar ci_ub_diff  = r(r2)

      qui count if abs(_draw1    - b_g0)   >= abs(b_g0)
      scalar p_g0    = (1 + r(N)) / (B_ok + 1)
      qui count if abs(_draw2    - b_g1)   >= abs(b_g1)
      scalar p_g1    = (1 + r(N)) / (B_ok + 1)
      qui count if abs(_drawdiff - b_diff) >= abs(b_diff)
      scalar p_diff  = (1 + r(N)) / (B_ok + 1)
    }
    restore
  }

  // ── Display
  local _stat_label = cond(use_bootstrap, "z", "t")
  local _p_label    = cond(use_bootstrap, "P>|z|", "P>|t|")
  local _ci_tag     = cond(use_bootstrap & "`normal'" == "", " (emp.)", "")
  di
  di as text "DiD subgroup analysis  |  unit: " as result "`unit'" ///
     as text "  |  time: " as result "`time'" ///
     as text "  |  post = " as result "`post_value'"
  di
  di as text "{hline 13}{c TT}{hline 64}"
  di as text %12s abbrev("`depvar'",12) " {c |}" ///
    _col(15) "{ralign 11:Coef.}" ///
    _col(26) "{ralign 12:Std. Err.}" ///
    _col(38) "{ralign 8:`_stat_label' }" ///
    _col(46) "{ralign 8:`_p_label'}" ///
    _col(54) "{ralign 25:[95% Conf. Interval`_ci_tag']}"
  di as text "{hline 13}{c +}{hline 64}"
  di as text "Subgroup" _col(14) "{c |}"
  forvalues g = 0/1 {
    display as text %12s abbrev("`g'",12) " {c |}" ///
      as result ///
      "  " %9.0g b_g`g' ///
      "  " %9.0g se_g`g' ///
      "    " %5.2f t_g`g' ///
      "   " %5.3f p_g`g' ///
      "    " %9.0g ci_lb_g`g' ///
      "   " %9.0g ci_ub_g`g'
  }
  di as text "{hline 13}{c +}{hline 64}"
  display as text "Difference   {c |}" ///
    as result ///
    "  " %9.0g b_diff ///
    "  " %9.0g se_diff ///
    "    " %5.2f t_diff ///
    "   " %5.3f p_diff ///
    "    " %9.0g ci_lb_diff ///
    "   " %9.0g ci_ub_diff
  di as text "{hline 13}{c BT}{hline 64}"

  qui count if `touse' & `sgroup' == 0
  scalar N_G0 = r(N)
  qui count if `touse' & `sgroup' == 1
  scalar N_G1 = r(N)
  di as text "N (G=0): " as result %4.0f N_G0 ///
     as text "   N (G=1): " as result %4.0f N_G1
  if use_bootstrap {
    di as text "Cluster bootstrap: " as result %4.0f B_ok ///
       as text " / " as result `bsreps' ///
       as text " reps, clustered on '" as result "`unit'" ///
       as text "' (" as result %4.0f N_clust as text " clusters)"
  }

  // ── Balance tables (aggregate + treated-only).  Per Q9a, in DiD mode
  // balance is reported both on the full active sample and conditional on
  // treat == 1 — the paper recommends checking the treated-only table since
  // the parallel-trends assumption is on the treated.
  if `: list sizeof balance' > 0 {
    forvalues _bidx = 1/2 {
      tempvar _bal_mask
      if `_bidx' == 1 {
        qui gen byte `_bal_mask' = `touse'
        local _bal_label "Aggregate"
      }
      else {
        qui gen byte `_bal_mask' = `touse' & (`treat' == 1)
        local _bal_label "Treated-only (D = 1)"
      }
      forvalues _widx = 1/2 {
        if `_widx' == 2 & "`ipsw'" == "noipsw" continue
        tempvar _bal_wt
        if `_widx' == 1 {
          qui gen byte `_bal_wt' = 1
          local _bal_wlabel "Unweighted"
        }
        else {
          qui gen double `_bal_wt' = `ipsweight'
          local _bal_wlabel "IPW-weighted"
        }
        di
        di as text "`_bal_label' balance ({txt}`_bal_wlabel'):"
        di as text "{hline 13}{c TT}{hline 44}"
        di as text "{ralign 12:Variable}  {c |}" ///
          _col(17) "{ralign 11:Mean G0}" ///
          _col(30) "{ralign 11:Mean G1}" ///
          _col(42) "{ralign 9:Std diff}" ///
          _col(54) "{ralign 8:p-value}"
        di as text "{hline 13}{c +}{hline 44}"
        foreach _v of local balance {
          qui sum `_v' [aw=`_bal_wt'] if `_bal_mask' & `sgroup' == 0, meanonly
          scalar _m0 = r(mean)
          qui sum `_v' [aw=`_bal_wt'] if `_bal_mask' & `sgroup' == 1, meanonly
          scalar _m1 = r(mean)
          qui sum `_v' if `_bal_mask'
          scalar _sd = r(sd)
          scalar _sdiff = cond(_sd > 0, (_m1 - _m0)/_sd, 0)
          qui reg `_v' `sgroup' [pw=`_bal_wt'] if `_bal_mask'
          scalar _pv = 2*ttail(e(df_r), abs(_b[`sgroup']/_se[`sgroup']))
          di as text "{ralign 12:`_v'}  {c |}" ///
            _col(17) as result %11.0g _m0 ///
            _col(30) %11.0g _m1 ///
            _col(42) %9.4f _sdiff ///
            _col(54) %8.4f _pv
        }
        di as text "{hline 13}{c BT}{hline 44}"
      }
    }
  }

  // ── Post clean e(b)/e(V) trimmed to the two treatment-effect columns.
  // Analytic case: cov_01 comes from xtreg's VCV.
  // Bootstrap case: cov_01 is overwritten above with the bootstrap covariance.
  matrix _b_did = (b_g0, b_g1)
  matrix colnames _b_did = G0_Z G1_Z
  matrix _V_did = (se_g0^2, cov_01 \ cov_01, se_g1^2)
  matrix rownames _V_did = G0_Z G1_Z
  matrix colnames _V_did = G0_Z G1_Z
  ereturn post _b_did _V_did, esample(`_did_esample')
  ereturn local cmd    "wsga"
  ereturn local subcmd "did"
  ereturn local depvar "`depvar'"
  ereturn scalar b_g0       = b_g0
  ereturn scalar b_g1       = b_g1
  ereturn scalar b_diff     = b_diff
  ereturn scalar se_g0      = se_g0
  ereturn scalar se_g1      = se_g1
  ereturn scalar se_diff    = se_diff
  ereturn scalar t_g0       = t_g0
  ereturn scalar t_g1       = t_g1
  ereturn scalar t_diff     = t_diff
  ereturn scalar p_g0       = p_g0
  ereturn scalar p_g1       = p_g1
  ereturn scalar p_diff     = p_diff
  ereturn scalar ci_lb_g0   = ci_lb_g0
  ereturn scalar ci_ub_g0   = ci_ub_g0
  ereturn scalar ci_lb_g1   = ci_lb_g1
  ereturn scalar ci_ub_g1   = ci_ub_g1
  ereturn scalar ci_lb_diff = ci_lb_diff
  ereturn scalar ci_ub_diff = ci_ub_diff
  ereturn scalar N_G0       = N_G0
  ereturn scalar N_G1       = N_G1
  ereturn scalar df         = df
  if use_bootstrap {
    ereturn scalar B_ok    = B_ok
    ereturn scalar N_clust = N_clust
  }
end

********************************************************************************

/*
CHANGE LOG
1.0.2
  - fix(Stata/did): wire up ipsweight()/pscore() as named output variables
  - fix(Stata/did): implement comsup restriction (G=1 pscore range, unit-level)
  - fix(Stata/did): implement blockbootstrap() as stratified cluster resample
  - fix(Stata/did): validate blockbootstrap() variable for unit-constancy
  - fix(Stata/did): apply comsup consistently in bootstrap reps
  - fix(Stata/did): remove stale "bootstrap not implemented" warning
1.0.1
  - bug(Stata/did): add ereturn post at end of _wsga_did to trim e(b)/e(V) to
    G0_Z and G1_Z columns; bootstrap path now also computes covariance of the
    two draws so V is fully bootstrap-derived when bootstrap is on
  - chore: align Stata *! version with unified 1.0.x versioning
1.3
  - Add design(did) path: sharp 2-period DiD-SGA via xtreg, fe with
    cluster-robust SEs (TWFE design matrix: subgroup-interacted unit FE
    is absorbed; G0_Z = G0*D*post and G1_Z = G1*D*post are the
    coefficients of interest).  IPW, balance tables, and the pairs
    cluster bootstrap will follow in subsequent patches; the
    corresponding options are accepted but trigger a one-line note.
1.2
  - Rename command: rddsga -> wsga (Weighted Subgroup Analysis)
  - rddsga retained as deprecated alias (separate .ado / .sthlp)
1.1.1
  - Fix: empirical/normal bootstrap CI display now actually differ
  - Fix: no-bootstrap diff p-value uses ttail(df, ...) for consistency
1.1
  - Re-estimate propensity score per bootstrap replicate (PS refit)
  - Empirical percentile CIs and (1+count)/(B+1) p-values for diff
1.0
  - Compute block bootstrapped variance-covariance matrix
  - Fix First Stage to compute bootstrap in IV
0.9
  - Compute bootstrapped variance-covariance matrix
  - Make program (and subprograms) e-class
  - Allow issuing no model
0.8
  - Add synthetic dataset for examples
0.7 
  - First alpha version ready for full usage
  - Implement nlcom hack to all models, detect diff coef position automatically
0.6
  - Implement nlcom hack to show difference as additional coefficient in ivreg
0.5
  - Fist working version with IVREG, reduced form and first stage equations
  - Implement output reporting with estimates table and estout
  - Default binarymodel is logit
0.4
  - First working version with IVREG equation
0.3
  - Standardize syntax to merge with original rddsga.ado
0.2
  - Implement _wsga_rdd_balancematrix as separate subroutine
  - Standardize _wsga_rdd_balancematrix output
0.1
	- First working version, independent of project
	- Remove any LaTeX output
	- Modify some option names and internal locals

KNOWN ISSUES/BUGS:
  - Should we use pweights or iweights? iw don't work with ivregress.

TODOS AND IDEAS:
  - Create subroutine of matlist formatting for display of _wsga_rdd_balancematrix output
  - Implement matrix manipulation in Mata
  - Get rid of sgroup0 hack
  - Allow that groupvar is not necessarily an indicator variable
  - Is it possible to allow for N subgroups?
*/
