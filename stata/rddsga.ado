*! 1.0.3 Alvaro Carril 2026-05-12
program define rddsga
  version 11.1
  di as error "rddsga is removed. Use {cmd:wsga rdd} instead."
  di as text  "Note: the running variable is now a named option: {opt running(varname)}"
  di as text  "Example: wsga rdd `depvar' [covars], sgroup(`sgroup') running(`assignvar') bwidth(...)"
  di as text  "See {help wsga rdd} for full syntax."
  exit 199
end
