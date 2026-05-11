{smcl}
{* *! version 1.0.2 2026-05-11}{...}
{title:Title}

{pstd}
{hi:rddsga} {hline 2} Removed alias (previously deprecated)


{title:Description}

{pstd}
{cmd:rddsga} has been removed. Use {cmd:wsga rdd} instead.

{pstd}
Note: the running variable is now a named option rather than the second
positional variable in the varlist:

{phang2}
{cmd:wsga rdd} {it:depvar} [{it:covariates}]{cmd:,}
{opt sgroup(varname)} {opt running(varname)} {opt bwidth(#)} [...]


{title:See also}

{pstd}
{help wsga}, {help wsga rdd}
