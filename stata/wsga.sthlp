{smcl}
{* *! version 1.3.0 2026-05-11}{...}
{title:Title}

{pstd}
{hi:wsga} {hline 2} Weighted subgroup analysis (regression discontinuity and difference-in-differences)


{title:Syntax}

{pstd}
{ul:Regression discontinuity (RD) design} (default):

{p 8 16 2}
{cmd:wsga} {depvar} {it:assignvar} [{indepvars}] {ifin}
{cmd:,} {opt sg:roup(varname)} {opt bw:idth(real)} [{it:rdd_options} {it:common_options}]
{p_end}

{pstd}
{ul:Difference-in-differences (DiD) design}:

{p 8 16 2}
{cmd:wsga} {depvar} [{indepvars}] {ifin}
{cmd:,} {opt sg:roup(varname)} {opt des:ign(did)} {opt un:it(varname)} {opt ti:me(varname)} {opt tr:eat(varname)} [{it:did_options} {it:common_options}]
{p_end}

{phang}
For RD: {it:depvar} is the outcome; {it:assignvar} is the running variable with a known cutoff; {it:indepvars} are optional control/balance variables.{p_end}
{phang}
For DiD: {it:depvar} is the outcome; {it:indepvars} are optional covariates included in the outcome regression and used in the propensity score model.{p_end}

{synoptset 26 tabbed}{...}
{synopthdr}
{synoptline}
{syntab :Design (required)}
{p2coldent:* {opth sg:roup(varname)}}binary (0/1) subgroup indicator variable{p_end}
{synopt :{opt des:ign(string)}}{opt rdd} (default) or {opt did}; selects the estimation pipeline{p_end}

{syntab :RDD-specific options}
{p2coldent:* {opt bw:idth(real)}}symmetric bandwidth around the cutoff; required for {opt design(rdd)}{p_end}
{synopt :{opth f:uzzy(varname)}}actual treatment receipt for fuzzy RD; if omitted, sharp RD is assumed{p_end}
{synopt :{opt c:utoff(real)}}cutoff value in {it:assignvar}; default is {opt cutoff(0)}{p_end}
{synopt :{opt k:ernel(kernelfn)}}kernel function: {opt uniform} (default), {opt triangular}, or {opt epanechnikov}{p_end}
{synopt :{opt p:(int)}}polynomial order; default is {opt p(1)} (local linear){p_end}
{synopt :{opt first:stage}}estimate the first stage (discontinuity in treatment probability){p_end}
{synopt :{opt reduced:form}}estimate the reduced form RD effect{p_end}
{synopt :{opt iv:regress}}estimate the treatment effect via 2SLS; requires {opt fuzzy(varname)}{p_end}
{synopt :{opt fixed:bootstrap}}for fuzzy RD: bootstrap the reduced form only and divide by the fixed first-stage estimate{p_end}

{syntab :DiD-specific options}
{p2coldent:* {opt un:it(varname)}}panel unit identifier; required for {opt design(did)}{p_end}
{p2coldent:* {opt ti:me(varname)}}time variable (must have exactly 2 unique values); required for {opt design(did)}{p_end}
{p2coldent:* {opt tr:eat(varname)}}binary treatment indicator (unit-constant); required for {opt design(did)}{p_end}
{synopt :{opt post:_value(string)}}value of {opt time()} denoting the post period; default is the larger of the two values{p_end}

{syntab :Balance and weighting}
{p2coldent:+ {opth bal:ance(varlist)}}variables entering the propensity score model; default is {it:indepvars}{p_end}
{synopt :{opt probit}}estimate propensity score with {manhelp probit R:probit}; default is {manhelp logit R:logit}{p_end}
{synopt :{opt com:sup}}restrict sample to common propensity score support{p_end}
{synopt :{opt m:(int)}}weighting mode: {opt m(2)} (default) balances both subgroups toward the pooled distribution; {opt m(1)} weights G=1 toward G=0; {opt m(0)} weights G=0 toward G=1{p_end}
{synopt :{opt rbal:ance(int)}}conditional mean method for RD balance tables: {opt rbalance(0)} (default) extrapolates to the cutoff; {opt rbalance(1)} uses the within-bandwidth mean{p_end}
{synopt :{opt noipsw}}skip IPW entirely; subgroup estimates are computed without reweighting{p_end}
{synopt :{opt dibal:ance}}display unweighted and IPW-weighted balance tables{p_end}
{synopt :{opth ipsw:eight(newvar)}}save IPW weights to a new variable{p_end}
{synopt :{opth psc:ore(newvar)}}save propensity scores to a new variable{p_end}

{syntab :Inference}
{synopt :{opth vce(vcetype)}}analytical variance estimator; default is {opt robust}; see {help vce_option}{p_end}
{synopt :{opt weights(weightsvar)}}optional observation weights multiplied into the kernel (RDD) or used as frequency weights (DiD){p_end}
{synopt :{opt noboot:strap}}suppress bootstrap standard errors; use analytical SEs only{p_end}
{synopt :{opt bsr:eps(#)}}bootstrap replications; default is {opt bsreps(50)}{p_end}
{synopt :{opt norm:al}}report normal-approximation p-values and CIs instead of empirical (percentile){p_end}
{synopt :{opt block:bootstrap(varname)}}stratified bootstrap: resample within strata defined by {it:varname}{p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}* Required options.{p_end}
{p 4 6 2}+ Required if {it:indepvars} is empty and IPSW is used.{p_end}
{p 4 6 2}{it:indepvars} may contain factor variables; see {help fvvarlist}.{p_end}


{marker description}{...}
{title:Description}

{pstd}
{cmd:wsga} conducts binary subgroup analysis using inverse propensity score weighting (IPSW) in two research designs.

{pstd}
{ul:Regression discontinuity ({opt design(rdd)}, default).}
Observations near the cutoff in each subgroup are reweighted so that the two subgroups have the same distribution of observed moderators at the cutoff.
Any remaining difference in per-subgroup RD estimates can therefore be attributed to the subgroup characteristic rather than correlated moderators.
{cmd:wsga} supports sharp RD ({opt reducedform}), the first stage of a fuzzy RD ({opt firststage}), and fuzzy RD via 2SLS ({opt ivregress}).

{pstd}
{ul:Difference-in-differences ({opt design(did)}).}
For a balanced 2-period panel, {cmd:wsga} runs a long-form two-way fixed effects (TWFE) regression with subgroup × post-treatment interactions.
IPW reweighting ensures the two subgroups have the same distribution of observed moderators, so the difference in DiD estimates is attributable to the subgroup characteristic.
Treatment must be binary, sharp, and unit-constant; fuzzy DiD is not supported.
{opt design(did)} requires {opt unit()}, {opt time()}, and {opt treat()}.

{pstd}
{ul:Common features.}
The propensity score for subgroup membership is estimated via {manhelp logit R:logit} (default) or {manhelp probit R:probit} on {it:indepvars}, or a separate set specified with {opt balance(varlist)}.
Three weighting modes are available via {opt m()}: {opt m(2)} (default) targets the pooled distribution; {opt m(1)} and {opt m(0)} are ATT-style modes.
Standard errors are bootstrap-based by default (empirical p-values and percentile CIs); the {opt normal} option switches to normal-approximation inference; {opt nobootstrap} uses the analytical sandwich estimator only.
Balance tables report the standardized mean difference and joint F-statistic for both the unweighted and IPW-weighted samples.
For {opt design(did)}, a second balance table restricted to treated units (D=1) is also reported.


{marker options}{...}
{title:Options}

{dlgtab:Design (required)}

{phang}
{opt sgroup(varname)} binary (0/1) subgroup indicator. Required for all designs.

{phang}
{opt design(string)} selects the estimation pipeline. {opt rdd} (default) runs the regression discontinuity path; {opt did} runs the 2-period difference-in-differences path.


{dlgtab:RDD-specific options}

{phang}
{opt bwidth(real)} symmetric bandwidth around the cutoff. Required for {opt design(rdd)}.

{phang}
{opt fuzzy(varname)} actual treatment receipt indicator for fuzzy RD. Required when {opt ivregress} or {opt firststage} is specified.

{phang}
{opt cutoff(real)} cutoff value of {it:assignvar}; default is {opt cutoff(0)}.

{phang}
{opt kernel(kernelfn)} kernel for local-polynomial weighting: {opt uniform} (default), {opt triangular}, or {opt epanechnikov}.

{phang}
{opt p(int)} polynomial order; default is {opt p(1)} (local linear).

{phang}
{opt firststage} estimates the discontinuity in treatment probability (first stage of a fuzzy RD).

{phang}
{opt reducedform} estimates the reduced form RD effect on the outcome.

{phang}
{opt ivregress} estimates the treatment effect using 2SLS. Requires {opt fuzzy(varname)}.

{phang}
{opt fixedbootstrap} for fuzzy RD: bootstrap the reduced form only and divide by the fixed first-stage point estimate.


{dlgtab:DiD-specific options}

{phang}
{opt unit(varname)} panel unit identifier. Required for {opt design(did)}.

{phang}
{opt time(varname)} time variable. Must have exactly two unique non-missing values. Required for {opt design(did)}.

{phang}
{opt treat(varname)} binary treatment indicator. Must be unit-constant (same value in both periods for every unit). Required for {opt design(did)}.

{phang}
{opt post_value(string)} value of {opt time()} that denotes the post period. Defaults to the larger of the two time values.


{dlgtab:Balance and weighting}

{phang}
{opt balance(varlist)} variables entering the propensity score model.
Defaults to {it:indepvars}. Useful when the balancing set differs from the control set.
Required if {it:indepvars} is empty and IPSW is in use.

{phang}
{opt probit} use {manhelp probit R:probit} for the propensity score; default is {manhelp logit R:logit}.

{phang}
{opt comsup} restrict to common propensity score support: [{it:min}({it:pscore} | G=1), {it:max}({it:pscore} | G=1)].

{phang}
{opt m(int)} weighting mode.
{opt m(2)} (default): weights both groups toward the pooled distribution.
{opt m(1)}: weights group 1 toward the distribution of group 0.
{opt m(0)}: weights group 0 toward the distribution of group 1.

{phang}
{opt rbalance(int)} method for computing conditional means in RD balance tables.
{opt rbalance(0)} (default): extrapolates to the cutoff via local linear regression.
{opt rbalance(1)}: uses within-bandwidth sample means.

{phang}
{opt noipsw} skip IPW entirely; subgroup estimates are computed without reweighting.

{phang}
{opt dibalance} display unweighted and IPW-weighted balance tables.

{phang}
{opt ipsweight(newvar)} save IPW weights to a new variable.

{phang}
{opt pscore(newvar)} save propensity scores to a new variable.


{dlgtab:Inference}

{phang}
{opt vce(vcetype)} analytical variance estimator; default is {opt robust}. See {help vce_option}.
For {opt design(did)} with a {opt unit()} variable, the default upgrades to {opt vce(cluster unit)}.

{phang}
{opt weights(weightsvar)} optional observation-level weights.

{phang}
{opt nobootstrap} suppress bootstrap variance estimation; use the sandwich estimator from {opt vce()} only.

{phang}
{opt bsreps(#)} number of bootstrap replications; default is {opt bsreps(50)}.

{phang}
{opt normal} use normal-approximation p-values and CIs.
By default {cmd:wsga} uses the percentile method: the empirical p-value is the share of bootstrap draws with |coef| >= |estimate|, and CIs are the 2.5th/97.5th percentiles of the bootstrap distribution.

{phang}
{opt blockbootstrap(varname)} stratified bootstrap: resample within strata defined by {it:varname}.
For {opt design(did)}, the bootstrap is a pairs-cluster bootstrap over units (Cameron-Gelbach-Miller); {opt blockbootstrap()} adds an additional stratification layer within that.


{marker examples}{...}
{title:Examples}

{pstd}{ul:Regression discontinuity}

{pstd}Load RD synthetic dataset{p_end}
{phang2}{cmd:. use rddsga_synth}{p_end}

{pstd}Assess covariate balance and display balance tables{p_end}
{phang2}{cmd:. wsga Y Z, balance(X1) sgroup(G) bwidth(10) dibalance}{p_end}

{pstd}Reduced form with IPW balancing on X1 and X2, 200 bootstrap replications{p_end}
{phang2}{cmd:. wsga Y Z X1 X2, sgroup(G) bwidth(10) reducedform bsreps(200)}{p_end}

{pstd}Fuzzy RD via 2SLS, saving IPW weights{p_end}
{phang2}{cmd:. wsga Y Z X1 X2, sgroup(G) bwidth(6) ivregress fuzzy(T) ipsweight(wt)}{p_end}

{pstd}Sharp RD without IPW{p_end}
{phang2}{cmd:. wsga Y Z, sgroup(G) bwidth(10) reducedform noipsw nobootstrap}{p_end}

{pstd}{ul:Difference-in-differences}

{pstd}Load DiD synthetic panel dataset (500 units, 2 periods){p_end}
{phang2}{cmd:. use wsga_did_synth}{p_end}

{pstd}Unweighted DiD subgroup analysis (no IPW){p_end}
{phang2}{cmd:. wsga y, sgroup(sgroup) design(did) unit(unit) time(time) treat(D) noipsw}{p_end}

{pstd}IPW-weighted DiD balancing on moderator m, display balance tables{p_end}
{phang2}{cmd:. wsga y m, sgroup(sgroup) design(did) unit(unit) time(time) treat(D) dibalance}{p_end}

{pstd}DiD with 200 cluster bootstrap replications{p_end}
{phang2}{cmd:. wsga y m, sgroup(sgroup) design(did) unit(unit) time(time) treat(D) bsreps(200)}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:wsga} stores the following in {cmd:e()}:

{synoptset 28 tabbed}{...}
{p2col 5 28 32 2: Scalars}{p_end}
{synopt:{cmd:e(unw_N_G0)}}observations in subgroup 0 (unweighted){p_end}
{synopt:{cmd:e(unw_N_G1)}}observations in subgroup 1 (unweighted){p_end}
{synopt:{cmd:e(unw_Fstat)}}joint F-statistic (unweighted){p_end}
{synopt:{cmd:e(unw_pvalue)}}F-statistic p-value (unweighted){p_end}
{synopt:{cmd:e(unw_avgdiff)}}mean absolute standardized difference (unweighted){p_end}
{synopt:{cmd:e(ipsw_N_G0)}}observations in subgroup 0 (IPW){p_end}
{synopt:{cmd:e(ipsw_N_G1)}}observations in subgroup 1 (IPW){p_end}
{synopt:{cmd:e(ipsw_Fstat)}}joint F-statistic (IPW){p_end}
{synopt:{cmd:e(ipsw_pvalue)}}F-statistic p-value (IPW){p_end}
{synopt:{cmd:e(ipsw_avgdiff)}}mean absolute standardized difference (IPW){p_end}
{synopt:{cmd:e(N_reps)}}number of bootstrap replications{p_end}
{synopt:{cmd:e(pval0)}}empirical p-value, subgroup 0{p_end}
{synopt:{cmd:e(pval1)}}empirical p-value, subgroup 1{p_end}
{synopt:{cmd:e(lb_g0)}}bootstrap CI lower bound, subgroup 0{p_end}
{synopt:{cmd:e(ub_g0)}}bootstrap CI upper bound, subgroup 0{p_end}
{synopt:{cmd:e(lb_g1)}}bootstrap CI lower bound, subgroup 1{p_end}
{synopt:{cmd:e(ub_g1)}}bootstrap CI upper bound, subgroup 1{p_end}
{synopt:{cmd:e(design)}}design in use: {cmd:rdd} or {cmd:did}{p_end}
{synopt:{cmd:e(N_units)}}({opt design(did)} only) number of unique panel units{p_end}
{synopt:{cmd:e(N_clusters)}}number of clusters used in a cluster bootstrap{p_end}

{p2col 5 28 32 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}subgroup coefficient vector (G=0, G=1){p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix{p_end}
{synopt:{cmd:e(unw)}}aggregate balance table (unweighted){p_end}
{synopt:{cmd:e(ipsw)}}aggregate balance table (IPW){p_end}
{synopt:{cmd:e(unw_treated)}}({opt design(did)} only) treated-only balance table (unweighted){p_end}
{synopt:{cmd:e(ipsw_treated)}}({opt design(did)} only) treated-only balance table (IPW){p_end}
{p2colreset}{...}

{pstd}
Additionally, all scalars and macros from the underlying estimation command ({opt reducedform}/{opt firststage}/{opt ivregress} for RDD; {help xtreg:xtreg, fe} for DiD) are stored in {cmd:e()}.
See {help regress##results:Stored Results} and {help xtreg##results:xtreg Stored Results} for details.


{marker authors}{...}
{title:Authors}

{pstd}
Alvaro Carril (maintainer){break}
acarril@princeton.edu

{pstd}
Andre Cazor{break}

{pstd}
Maria Paula Gerardino{break}
Inter-American Development Bank{break}

{pstd}
Stephan Litschig{break}
National Graduate Institute for Policy Studies{break}

{pstd}
Dina Pomeranz{break}
University of Zurich{break}


{marker disclaimer}{...}
{title:Disclaimer}

{pstd}
This software is provided "as is", without warranty of any kind.
Please report issues at {browse "https://github.com/acarril/wsga/issues"}.


{marker references}{...}
{title:References}

{phang}
Carril, Alvaro, Andre Cazor, Maria Paula Gerardino, Stephan Litschig, and Dina Pomeranz. "Weighted Subgroup Analysis in Regression Discontinuity Designs." Working paper.{p_end}

{phang}
Cameron, A. Colin, Jonah B. Gelbach, and Douglas L. Miller. "Bootstrap-Based Improvements for Inference with Clustered Errors." {it:Review of Economics and Statistics} 90(3): 414-427, 2008.{p_end}
