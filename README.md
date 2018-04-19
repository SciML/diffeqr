# diffeqr

[![Build Status](https://travis-ci.org/JuliaDiffEq/diffeqr.svg?branch=master)](https://travis-ci.org/JuliaDiffEq/diffeqr)

diffeqr is a package for solving differential equations in R. It utilizes DifferentialEquations.jl
for its core routines to give high performance solving of ordinary differential equations (ODEs),
stochastic differential equations (SDEs), delay differential equations (DDEs), and
differential-algebraic equations (DAEs) directly in R.

## Installation

diffeqr is currently not registered into CRAN. Thus to use this package, use the following command:

```R
devtools::install_github('JuliaDiffEq/diffeqr', build_vignettes=T)
```
