{smcl}
{* *! version 1.0.1 2026-05-11}{...}
{viewerjumpto "Stored results" "wsga_did##results"}{...}
{title:Title}

{pstd}
{hi:wsga did} {hline 2} Weighted Subgroup Analysis for Difference-in-Differences Designs


{title:Syntax}

{p 8 16 2}
{cmd:wsga did} {it:depvar} [{it:covariates}] {ifin}{cmd:,}
  {opt sgroup(varname)}
  {opt unit(varname)}
  {opt time(varname)}
  {opt treat(varname)}
  [{it:options}]

{synoptset 26 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{opt sgroup(varname)}}binary (0/1) subgroup indicator (unit-constant){p_end}
{synopt:{opt unit(varname)}}unit identifier (panel ID){p_end}
{synopt:{opt time(varname)}}time variable; must have exactly 2 unique values{p_end}
{synopt:{opt treat(varname)}}binary (0/1) treatment indicator (unit-constant){p_end}
{syntab:DiD options}
{synopt:{opt post_value(val)}}value of {it:time} denoting the post period; default max({it:time}){p_end}
{syntab:IPW}
{synopt:{opt balance(varlist)}}moderators for propensity score; defaults to covariates{p_end}
{synopt:{opt noipsw}}skip IPW reweighting{p_end}
{synopt:{opt m(#)}}weighting mode: 2 = both groups (default), 1 = G1→G0, 0 = G0→G1{p_end}
{synopt:{opt probit}}use probit instead of logit for propensity score{p_end}
{synopt:{opt comsup}}restrict to common propensity score support{p_end}
{synopt:{opt dibalance}}display balance tables{p_end}
{syntab:Inference}
{synopt:{opt nobootstrap}}report analytical cluster-robust SEs only{p_end}
{synopt:{opt bsreps(#)}}bootstrap replications; default 200{p_end}
{synopt:{opt normal}}normal-approximation CIs from bootstrap SE{p_end}
{synopt:{opt seed(#)}}RNG seed for reproducibility{p_end}
{synopt:{opt fixedps}}hold propensity score fixed across bootstrap reps (diagnostic){p_end}
{synopt:{opt blockbootstrap(varname)}}stratified block bootstrap{p_end}
{synopt:{opt weights(varname)}}pre-existing observation weights{p_end}
{synoptline}


{title:Description}

{pstd}
{cmd:wsga did} estimates per-subgroup treatment effects and their difference
in a sharp 2-period difference-in-differences design. The estimator is a
long-form TWFE regression (unit + time fixed effects) fully interacted by
subgroup. Inverse probability weighting (IPW) balances observed moderators
across subgroups.

{pstd}
The cluster bootstrap resamples whole units with replacement and assigns fresh
unit IDs per draw so fixed effects remain identified (Cameron–Gelbach–Miller).
Clustering is always on {it:unit}.


{title:Examples}

{pstd}Unweighted DiD subgroup analysis:{p_end}
{phang2}{cmd:. wsga did Y, sgroup(G) unit(id) time(t) treat(D) noipsw nobootstrap}{p_end}

{pstd}IPW-weighted DiD:{p_end}
{phang2}{cmd:. wsga did Y M, sgroup(G) unit(id) time(t) treat(D) nobootstrap}{p_end}

{pstd}With cluster bootstrap:{p_end}
{phang2}{cmd:. wsga did Y M, sgroup(G) unit(id) time(t) treat(D) bsreps(200) seed(1)}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:wsga did} stores the following in {cmd:e()}:

{synoptset 23 tabbed}{...}
{p2col 5 23 26 2: Scalars}{p_end}
{synopt:{cmd:e(b_g0)}}subgroup-0 treatment effect{p_end}
{synopt:{cmd:e(b_g1)}}subgroup-1 treatment effect{p_end}
{synopt:{cmd:e(b_diff)}}{cmd:e(b_g1)} minus {cmd:e(b_g0)}{p_end}
{synopt:{cmd:e(se_g0)}}standard error for subgroup-0{p_end}
{synopt:{cmd:e(se_g1)}}standard error for subgroup-1{p_end}
{synopt:{cmd:e(se_diff)}}standard error for the difference{p_end}
{synopt:{cmd:e(t_g0)}}t or z statistic for subgroup-0{p_end}
{synopt:{cmd:e(t_g1)}}t or z statistic for subgroup-1{p_end}
{synopt:{cmd:e(t_diff)}}t or z statistic for the difference{p_end}
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
{synopt:{cmd:e(df)}}degrees of freedom for analytical inference{p_end}
{synopt:{cmd:e(B_ok)}}(bootstrap) number of successful replications{p_end}
{synopt:{cmd:e(N_clust)}}(bootstrap) number of clusters{p_end}

{p2col 5 23 26 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:wsga}{p_end}
{synopt:{cmd:e(subcmd)}}{cmd:did}{p_end}
{synopt:{cmd:e(depvar)}}name of dependent variable{p_end}

{p2col 5 23 26 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}1{cmd:x}2 coefficient vector with columns named {cmd:G0_Z} and {cmd:G1_Z}{p_end}
{synopt:{cmd:e(V)}}2{cmd:x}2 variance-covariance matrix; bootstrap-derived when bootstrap is on, analytical (cluster-robust) otherwise{p_end}

{p2col 5 23 26 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}marks estimation sample{p_end}
{p2colreset}{...}

{pstd}
Named column access works after estimation:{p_end}
{phang2}{cmd:. matrix b = e(b)}{p_end}
{phang2}{cmd:. scalar b_g0 = b[1,"G0_Z"]}{p_end}
{phang2}{cmd:. scalar se_diff = sqrt(e(V)[1,1] + e(V)[2,2] - 2*e(V)[1,2])}{p_end}


{title:See also}

{pstd}
{help wsga}, {help wsga rdd}
