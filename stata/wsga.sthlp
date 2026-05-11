{smcl}
{* *! version 1.5.0 2026-05-11}{...}
{title:Title}

{pstd}
{hi:wsga} {hline 2} Weighted Subgroup Analysis


{title:Description}

{pstd}
{cmd:wsga} implements weighted subgroup analysis for research designs that
require control variables for identification. It uses inverse probability
weighting (IPW) to balance observed moderators across subgroups, isolating
the subgroup-attributable component of the treatment effect difference.

{pstd}
{cmd:wsga} has two subcommands, one per research design:

{phang2}
{helpb wsga rdd} {hline 2} Regression discontinuity designs (sharp and fuzzy RD)

{phang2}
{helpb wsga did} {hline 2} Sharp 2-period difference-in-differences


{title:Quick syntax}

{pstd}
{cmd:wsga rdd} {it:depvar} [{it:covariates}]{cmd:,}
{opt sgroup(varname)} {opt running(varname)} {opt bwidth(#)} [...]

{pstd}
{cmd:wsga did} {it:depvar} [{it:covariates}]{cmd:,}
{opt sgroup(varname)} {opt unit(varname)} {opt time(varname)} {opt treat(varname)} [...]


{title:See also}

{pstd}
{help wsga rdd}, {help wsga did}
