{smcl}
{* *! version 1.2.0 2026-05-09}{...}
{title:Title}

{pstd}
{hi:wsga} {hline 2} Weighted subgroup analysis (regression discontinuity; sharp difference-in-differences planned)


{title:Syntax}

{p 8 16 2}
{cmd:wsga} {depvar} {it:assignvar} [{indepvars}] {ifin}
{cmd:,} {it:options}
{p_end}

{phang}
{it:depvar} is the outcome variable, {it:assignvar} is the assignment (running) variable for which there is a known cutoff, and {it:indepvars} are optional control variables included in the outcome regression.
{p_end}

{synoptset 24 tabbed}{...}
{synopthdr}
{synoptline}
{syntab :Subgroup and RDD}
{p2coldent:* {opth sg:roup(varname)}}binary subgroup indicator variable (must be 0/1){p_end}
{synopt :{opth f:uzzy(varname)}}actual treatment status indicator for fuzzy RD; if not specified, sharp RD is assumed{p_end}
{synopt :{opt c:utoff(real)}}cutoff value in {it:assignvar}; default is {opt cutoff(0)}{p_end}
{p2coldent:* {opt bw:idth(real)}}symmetric bandwidth around the cutoff{p_end}
{synopt :{opt k:ernel(kernelfn)}}kernel function for local-polynomial estimation; {opt uniform} (default), {opt triangular}, or {opt epanechnikov}{p_end}
{synopt :{opt p:(int)}}order of the local polynomial; default is {opt p(1)}{p_end}

{syntab :Balance and weighting}
{p2coldent:+ {opth bal:ance(varlist)}}variables entering the propensity score model; default is {it:indepvars}{p_end}
{synopt :{opt probit}}estimate propensity score with {manhelp probit R:probit}; default is {manhelp logit R:logit}{p_end}
{synopt :{opt com:sup}}restrict sample to common propensity score support{p_end}
{synopt :{opt m:(int)}}weighting mode: {opt m(2)} (default) balances both subgroups toward the pooled distribution; {opt m(1)} weights group 1 toward group 0; {opt m(0)} weights group 0 toward group 1{p_end}
{synopt :{opt rbal:ance:(int)}}how conditional means are computed for balance tables: {opt rbalance(0)} (default) evaluates at the cutoff; {opt rbalance(1)} uses the within-bandwidth sample mean{p_end}
{synopt :{opt noipsw}}skip inverse propensity score weighting entirely{p_end}
{synopt :{opt dibal:ance}}display unweighted and IPSW balance tables{p_end}
{synopt :{opth ipsw:eight(newvar)}}save IPW weights to a new variable{p_end}
{synopt :{opth psc:ore(newvar)}}save propensity scores to a new variable{p_end}

{syntab :Model}
{synopt :{opt first:stage}}estimate the first stage (discontinuity in treatment probability){p_end}
{synopt :{opt reduced:form}}estimate the reduced form effect{p_end}
{synopt :{opt iv:regress}}estimate the treatment effect via 2SLS; requires {opt fuzzy(varname)}{p_end}
{synopt :{opth vce(vcetype)}}{it:vcetype} may be {opt r:obust}, {opt cl:uster} {it:clustvar}, or any option accepted by {manhelp regress R:regress}; default is {opt robust}{p_end}
{synopt :{opt weights(weightsvar)}}optional observation weights multiplied into the kernel{p_end}

{syntab :Bootstrap}
{synopt :{opt noboot:strap}}suppress bootstrap standard errors; use sandwich SEs only{p_end}
{synopt :{opt bsr:eps(#)}}number of bootstrap replications; default is {opt bsreps(50)}{p_end}
{synopt :{opt norm:al}}report normal-approximation p-values and CIs instead of empirical (percentile){p_end}
{synopt :{opt fixed:bootstrap}}for fuzzy RD: bootstrap the reduced form only and divide by the fixed first-stage point estimate{p_end}
{synopt :{opt block:bootstrap(varname)}}stratified bootstrap: resample within strata defined by {it:varname}{p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}* Required options.{p_end}
{p 4 6 2}+ Required if {it:indepvars} is empty and IPSW is used.{p_end}
{p 4 6 2}{it:indepvars} may contain factor variables; see {help fvvarlist}.{p_end}


{marker description}{...}
{title:Description}

{pstd}
{cmd:wsga} conducts binary subgroup analysis in regression discontinuity (RD) settings using inverse propensity score weighting (IPSW).
Observations in each subgroup are reweighted by the inverse of their estimated probability of belonging to that subgroup, given a set of observed moderators.
This ensures that the two subgroups have the same distribution of observed moderators at the cutoff, so that any remaining difference in subgroup RD estimates can be attributed to the subgroup characteristic rather than to correlated moderators.

{pstd}
The propensity score is estimated via {manhelp logit R:logit} (default) or {manhelp probit R:probit} using {it:indepvars}, or a separate set specified with {opt balance(varlist)}.
Three weighting modes are available via {opt m()}: {opt m(2)} (default) targets the pooled distribution; {opt m(1)} and {opt m(0)} are ATT-style modes targeting group 0 and group 1 respectively.

{pstd}
Estimation supports sharp RD ({opt reducedform}), the first stage of a fuzzy RD ({opt firststage}), and fuzzy RD via 2SLS ({opt ivregress}).
All three estimators report per-subgroup RD gaps and their difference.
By default, standard errors are bootstrap-based with empirical p-values and percentile confidence intervals; the {opt normal} option switches to normal-approximation inference.

{pstd}
Balance tables report the standardized mean difference and p-value for each covariate, plus a joint F-statistic, for both the unweighted and IPW-weighted samples.


{marker options}{...}
{title:Options}

{dlgtab:Subgroup and RDD}

{phang}
{opt sgroup(varname)} binary (0/1) subgroup indicator. Required.

{phang}
{opt fuzzy(varname)} actual treatment receipt indicator for fuzzy RD. Required when {opt ivregress} or {opt firststage} is specified.

{phang}
{opt cutoff(real)} cutoff value of {it:assignvar}; default is {opt cutoff(0)}.

{phang}
{opt bwidth(real)} symmetric bandwidth around the cutoff. Required.

{phang}
{opt kernel(kernelfn)} kernel for local-polynomial weighting: {opt uniform} (default), {opt triangular}, or {opt epanechnikov}.

{phang}
{opt p(int)} polynomial order; default is {opt p(1)} (local linear).


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
{opt rbalance(int)} method for computing conditional means in balance tables.
{opt rbalance(0)} (default): extrapolates to the cutoff via local linear regression within each subgroup.
{opt rbalance(1)}: uses within-bandwidth sample means.

{phang}
{opt noipsw} skip IPW entirely; subgroup estimates are computed without reweighting.

{phang}
{opt dibalance} display unweighted and IPW-weighted balance tables.

{phang}
{opt ipsweight(newvar)} save IPW weights to a new variable.

{phang}
{opt pscore(newvar)} save propensity scores to a new variable.


{dlgtab:Model}

{phang}
{opt firststage} estimates the discontinuity in treatment probability (first stage).

{phang}
{opt reducedform} estimates the reduced form RD effect on the outcome.

{phang}
{opt ivregress} estimates the treatment effect using 2SLS. Requires {opt fuzzy(varname)}.

{phang}
{opt vce(vcetype)} variance estimator; default is {opt robust}. See {help vce_option}.

{phang}
{opt weights(weightsvar)} optional observation-level weights multiplied into the kernel weights.


{dlgtab:Bootstrap}

{phang}
{opt nobootstrap} suppress bootstrap variance estimation; use the sandwich estimator from {opt vce()} only.

{phang}
{opt bsreps(#)} number of bootstrap replications; default is {opt bsreps(50)}.

{phang}
{opt normal} use normal-approximation p-values and CIs.
By default {cmd:wsga} uses the percentile method: the empirical p-value is the share of bootstrap draws with |coef| >= |estimate|, and CIs are the 2.5th/97.5th percentiles of the bootstrap distribution.

{phang}
{opt fixedbootstrap} for fuzzy RD: bootstrap the reduced form only and divide by the fixed first-stage point estimate.

{phang}
{opt blockbootstrap(varname)} stratified bootstrap: resample within strata defined by {it:varname}.


{marker examples}{...}
{title:Examples}

{pstd}Load synthetic dataset{p_end}
{phang2}{cmd:. use rddsga_synth}{p_end}

{pstd}Assess covariate balance and display balance tables{p_end}
{phang2}{cmd:. wsga Y Z, balance(X1) sgroup(G) bwidth(10) dibalance}{p_end}

{pstd}Reduced form with IPW balancing on X1 and X2, 200 bootstrap replications{p_end}
{phang2}{cmd:. wsga Y Z X1 X2, sgroup(G) bwidth(10) reducedform bsreps(200)}{p_end}

{pstd}Fuzzy RD via 2SLS, saving IPW weights{p_end}
{phang2}{cmd:. wsga Y Z X1 X2, sgroup(G) bwidth(6) ivregress fuzzy(T) ipsweight(wt)}{p_end}

{pstd}Sharp RD without IPW{p_end}
{phang2}{cmd:. wsga Y Z, sgroup(G) bwidth(10) reducedform noipsw nobootstrap}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:wsga} stores the following in {cmd:e()}:

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Scalars}{p_end}
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

{p2col 5 24 28 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}subgroup RD coefficient vector (G=0, G=1){p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix{p_end}
{synopt:{cmd:e(unw)}}balance table matrix (unweighted){p_end}
{synopt:{cmd:e(ipsw)}}balance table matrix (IPW){p_end}
{p2colreset}{...}

{pstd}
Additionally, all scalars and macros from the underlying {opt reducedform}, {opt firststage}, or {opt ivregress} estimation are stored in {cmd:e()}.
See {help regress##results:Stored Results} for details.


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
