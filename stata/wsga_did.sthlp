{smcl}
{* *! version 1.5.0 2026-05-11}{...}
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


{title:See also}

{pstd}
{help wsga}, {help wsga rdd}
