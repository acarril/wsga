{smcl}
{* *! version 1.2.0 2026-05-26}{...}
{viewerjumpto "Stored results" "wsga_rdd##results"}{...}
{title:Title}

{pstd}
{hi:wsga rdd} {hline 2} Weighted Subgroup Analysis for Regression Discontinuity Designs


{title:Syntax}

{p 8 16 2}
{cmd:wsga rdd} {it:depvar} [{it:covariates}] {ifin}{cmd:,}
  {opt sgroup(varname)}
  {opt running(varname)}
  {opt bwidth(#)}
  [{it:options}]

{synoptset 26 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{opt sgroup(varname)}}binary (0/1) subgroup indicator{p_end}
{synopt:{opt running(varname)}}running variable (centered at the cutoff internally){p_end}
{synopt:{opt bwidth(#)}}half-bandwidth (positive){p_end}
{syntab:RD model}
{synopt:{opt cutoff(#)}}cutoff value; default 0{p_end}
{synopt:{opt reducedform}}estimate reduced form (default){p_end}
{synopt:{opt firststage}}estimate first stage of fuzzy RD{p_end}
{synopt:{opt ivregress}}2SLS fuzzy RD{p_end}
{synopt:{opt fuzzy(varname)}}treatment take-up variable (required with {opt ivregress}/{opt firststage}){p_end}
{synopt:{opt p(#)}}polynomial order; default 1 (local linear){p_end}
{synopt:{opt kernel(string)}}{opt uni}form (default), {opt tri}angular, or {opt epa}nechnikov{p_end}
{syntab:IPW}
{synopt:{opt balance(varlist)}}moderators for propensity score; defaults to covariates{p_end}
{synopt:{opt noipsw}}skip IPW reweighting{p_end}
{synopt:{opt m(#)}}weighting mode: 2 = both groups (default), 1 = G1->G0, 0 = G0->G1{p_end}
{synopt:{opt probit}}use probit instead of logit for propensity score{p_end}
{synopt:{opt comsup}}restrict to common propensity score support{p_end}
{synopt:{opt dibalance}}display balance tables{p_end}
{synopt:{opt rbalance(#)}}balance table means: 0 = at cutoff (default), 1 = in sample{p_end}
{syntab:Inference}
{synopt:{opt nobootstrap}}report analytical SEs only{p_end}
{synopt:{opt bsreps(#)}}bootstrap replications; default 200{p_end}
{synopt:{opt normal}}normal-approximation CIs from bootstrap SE{p_end}
{synopt:{opt seed(#)}}RNG seed for reproducibility{p_end}
{synopt:{opt cluster(varname)}}clustering variable for pairs-cluster bootstrap{p_end}
{synopt:{opt wildcluster}}wild cluster bootstrap (WCB-U) with Rademacher signs; requires {opt cluster}{p_end}
{synopt:{opt fixedbootstrap}}fixed first-stage bootstrap for IV{p_end}
{synopt:{opt fixedps}}hold propensity score fixed across bootstrap reps (diagnostic){p_end}
{synopt:{opt blockbootstrap(varname)}}stratified block bootstrap{p_end}
{synopt:{opt vce(vcetype)}}SE type; default robust{p_end}
{synopt:{opt weights(varname)}}pre-existing observation weights{p_end}
{synoptline}


{title:Description}

{pstd}
{cmd:wsga rdd} estimates per-subgroup treatment effects and their difference
in a sharp or fuzzy regression discontinuity (RD) design. Inverse probability
weighting (IPW) balances observed moderators across subgroups, isolating the
subgroup-attributable component of the effect difference.


{title:Examples}

{pstd}Sharp RD, no IPW:{p_end}
{phang2}{cmd:. wsga rdd Y, sgroup(G) running(X) bwidth(0.5) noipsw nobootstrap}{p_end}

{pstd}Sharp RD with IPW on moderator M:{p_end}
{phang2}{cmd:. wsga rdd Y M, sgroup(G) running(X) bwidth(0.5) nobootstrap}{p_end}

{pstd}With bootstrap inference:{p_end}
{phang2}{cmd:. wsga rdd Y M, sgroup(G) running(X) bwidth(0.5) bsreps(200) seed(42)}{p_end}

{pstd}Fuzzy RD via 2SLS:{p_end}
{phang2}{cmd:. wsga rdd Y, sgroup(G) running(X) bwidth(0.5) ivregress fuzzy(D) nobootstrap}{p_end}

{pstd}Pairs-cluster bootstrap on school IDs:{p_end}
{phang2}{cmd:. wsga rdd Y M, sgroup(G) running(X) bwidth(0.5) bsreps(200) cluster(school) seed(42)}{p_end}

{pstd}Wild cluster bootstrap (recommended for small G):{p_end}
{phang2}{cmd:. wsga rdd Y M, sgroup(G) running(X) bwidth(0.5) bsreps(200) cluster(school) wildcluster seed(42)}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:wsga rdd} stores the following in {cmd:e()}:

{synoptset 23 tabbed}{...}
{p2col 5 23 26 2: Scalars}{p_end}
{synopt:{cmd:e(p_g0)}}p-value for subgroup-0{p_end}
{synopt:{cmd:e(p_g1)}}p-value for subgroup-1{p_end}
{synopt:{cmd:e(p_diff)}}p-value for the difference{p_end}
{synopt:{cmd:e(ci_lb_g0)}}95% confidence interval lower bound, subgroup-0{p_end}
{synopt:{cmd:e(ci_ub_g0)}}95% confidence interval upper bound, subgroup-0{p_end}
{synopt:{cmd:e(ci_lb_g1)}}95% confidence interval lower bound, subgroup-1{p_end}
{synopt:{cmd:e(ci_ub_g1)}}95% confidence interval upper bound, subgroup-1{p_end}
{synopt:{cmd:e(ci_lb_diff)}}95% confidence interval lower bound, difference{p_end}
{synopt:{cmd:e(ci_ub_diff)}}95% confidence interval upper bound, difference{p_end}
{synopt:{cmd:e(N_G0)}}number of observations in subgroup 0{p_end}
{synopt:{cmd:e(N_G1)}}number of observations in subgroup 1{p_end}
{synopt:{cmd:e(N_reps)}}(bootstrap) requested replications{p_end}
{synopt:{cmd:e(B_ok)}}(bootstrap) successful replications{p_end}
{synopt:{cmd:e(N_clust)}}(bootstrap, clustered) number of clusters{p_end}
{synopt:{cmd:e(level)}}confidence level (95){p_end}

{p2col 5 23 26 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:wsga}{p_end}
{synopt:{cmd:e(subcmd)}}{cmd:rdd}{p_end}
{synopt:{cmd:e(depvar)}}name of dependent variable{p_end}
{synopt:{cmd:e(boot_type)}}(bootstrap) {cmd:pairs} or {cmd:wild}{p_end}
{synopt:{cmd:e(clustvar)}}(bootstrap, clustered) clustering variable{p_end}
{synopt:{cmd:e(blockbootstrap)}}(bootstrap, stratified) stratification variable{p_end}

{p2col 5 23 26 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}1 x 2 coefficient row for ({cmd:G0_Z}, {cmd:G1_Z}){p_end}
{synopt:{cmd:e(V)}}2 x 2 variance-covariance matrix; bootstrap-derived when bootstrap is on{p_end}

{p2col 5 23 26 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}marks estimation sample{p_end}


{title:See also}

{pstd}
{help wsga}, {help wsga did}
