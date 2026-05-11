# wsga

**Weighted Subgroup Analysis**

Implements inverse probability weighted (IPW) subgroup analysis for research designs that require control variables for identification. When subgroups differ in observed moderators, a naive comparison of subgroup-specific treatment effects conflates the causal effect of the subgroup characteristic with the effect of correlated moderators. Reweighting observations via IPW balances observed moderators across subgroups, isolating the subgroup-attributable component of the treatment effect difference.

Both a **Stata** package and an **R** package are included in this repository. Both implement weighted subgroup analysis for **regression discontinuity (RD)** designs (sharp and fuzzy) and **sharp 2-period difference-in-differences (DiD)**.

---

## Stata

### Installation

Latest version from this repository:
```stata
net from https://raw.githubusercontent.com/acarril/wsga/main/stata/
net install wsga
net get wsga
```

### Quick start (RD)

```stata
use rddsga_synth
wsga Y Z X1 X2, sgroup(G) bwidth(10) reducedform bsreps(200)
```

### Quick start (DiD)

```stata
use wsga_did_synth
wsga y m, sgroup(sgroup) design(did) unit(unit) time(time) treat(D) bsreps(200)
```

See `help wsga` for full documentation. The previous command name `rddsga` is retained as a deprecated alias and is installed alongside `wsga`.

---

## R

### Installation

```r
devtools::install_github("acarril/wsga")
```

### Quick start (RD)

```r
library(wsga)
data(rddsga_synth)

fit <- wsga(
  y ~ m | sgroup,
  data    = rddsga_synth,
  running = ~ x,
  bwidth  = 0.5,
  bsreps  = 200,
  seed    = 42
)
print(fit)
summary(fit)   # also shows balance tables
```

### Quick start (DiD)

```r
library(wsga)
data(wsga_did_synth)

fit <- wsga(
  y ~ m | sgroup,
  data   = wsga_did_synth,
  design = "did",
  unit   = "unit",
  time   = "time",
  treat  = "D",
  bsreps = 200,
  seed   = 42
)
print(fit)
summary(fit)   # shows aggregate and treated-only balance tables
```

See `?wsga`, `vignette("wsga-intro")`, and `vignette("wsga-did-intro")` for full documentation. The previous function name `rddsga()` is retained as a deprecated alias.

---

## Reference

Carril, Alvaro, Andre Cazor, Maria Paula Gerardino, Stephan Litschig, and Dina Pomeranz. "Weighted Subgroup Analysis in Regression Discontinuity Designs." Working paper.

## Issues

Please report bugs at <https://github.com/acarril/wsga/issues>.
