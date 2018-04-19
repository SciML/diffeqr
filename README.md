# diffeqr

[![Build Status](https://travis-ci.org/JuliaDiffEq/diffeqr.svg?branch=master)](https://travis-ci.org/JuliaDiffEq/diffeqr)

diffeqr is a package for solving differential equations in R. It utilizes 
[DifferentialEquations.jl](http://docs.juliadiffeq.org/latest/) for its core routines 
to give high performance solving of ordinary differential equations (ODEs),
stochastic differential equations (SDEs), delay differential equations (DDEs), and
differential-algebraic equations (DAEs) directly in R.

## Installation

diffeqr is currently not registered into CRAN. Thus to use this package, use the following command:

```R
devtools::install_github('JuliaDiffEq/diffeqr', build_vignettes=T)
```

## Usage

diffeqr does not provide the full functionality of DifferentialEquations.jl. Instead, it supplies simplified
direct solving routines with an R interface. The most basic function is:

```R
ode.solve(f,u0,tspan,[p,abstol,reltol,saveat])
```

which solves the ODE `u' = f(u,p,t)` where `u(0)=u0` over the timespan `tspan`. 
The common interface arguments are documented 
[at the DifferentialEquations.jl page](http://docs.juliadiffeq.org/latest/basics/common_solver_opts.html).
Notice that not all options are allowed, but the most common arguments are supported.

## ODE Examples

### 1D Linear ODEs

Let's solv
```R
library(diffeqr)
```

f <- function(u,p,t) {
  return(1.01*u)
}
u0 = 1/2
tspan <- list(0.0,1.0)
sol = ode.solve(f,u0,tspan)

length(sol$u) == length(sol$t)
sol$u[length(sol$u)] - exp(1)/2 < 5e-2
exp(1)/2
plot(sol$t,sol$u,"l")
