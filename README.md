# diffeqr

[![Join the chat at https://gitter.im/JuliaDiffEq/Lobby](https://badges.gitter.im/JuliaDiffEq/Lobby.svg)](https://gitter.im/JuliaDiffEq/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
[![Build Status](https://travis-ci.com/SciML/diffeqr.svg?branch=master)](https://travis-ci.com/SciML/diffeqr)
[![Build status](https://ci.appveyor.com/api/projects/status/2pxp5kfu0uiddmpl?svg=true)](https://ci.appveyor.com/project/ChrisRackauckas/diffeqr)

diffeqr is a package for solving differential equations in R. It utilizes 
[DifferentialEquations.jl](http://diffeq.sciml.ai/dev/) for its core routines 
to give high performance solving of ordinary differential equations (ODEs),
stochastic differential equations (SDEs), delay differential equations (DDEs), and
differential-algebraic equations (DAEs) directly in R.

If you have any questions, or just want to chat about solvers/using the package,
please feel free to chat in the [Gitter channel](https://gitter.im/JuliaDiffEq/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge). For bug reports, feature requests, etc., please submit an issue. 

## Installation

[diffeqr is registered into CRAN](https://CRAN.R-project.org/package=diffeqr). Thus to add the package, use:

```R
install.packages("diffeqr")
```

To install the master branch of the package (for developers), use:

```R
devtools::install_github('SciML/diffeqr', build_vignettes=T)
```

You will need a working installation of Julia in your path. To install Julia, download a generic binary
from [the JuliaLang site](https://julialang.org/downloads/) and add it to your path. The download and
installation of DifferentialEquations.jl will happen on the first invocation of `diffeqr::diffeq_setup()`.

## Usage

diffeqr provides a direct wrapper over [DifferentialEquations.jl](diffeq.sciml.ai). 
The namespace is setup so that the standard syntax of Julia translates directly
over to the R environment. There are two things to keep in mind:

1. All DifferentialEquations.jl commands are prefaced by `de$`
2. All commands with a `!` are replaced with `_bang`, for example `solve!` becomes `solve_bang`.

## Ordinary Differential Equation (ODE) Examples

### 1D Linear ODEs

Let's solve the linear ODE `u'=1.01u`. First setup the package:

```R
de <- diffeqr::diffeq_setup()
```

Define the derivative function `f(u,p,t)`. 

```R
f <- function(u,p,t) {
  return(1.01*u)
}
```

Then we give it an initial condition and a time span to solve over:

```R
u0 <- 1/2
tspan <- c(0., 1.)
```

With those pieces we define the `ODEProblem` and `solve` the ODE:

```R
prob = de$ODEProblem(f, u0, tspan)
sol = de$solve(prob)
```

This gives back a solution object for which `sol$t` are the time points
and `sol$u` are the values. We can treat the solution as a continuous object
in time via 

```R
sol$.(0.2)
```

and a high order interpolation will compute the value at `t=0.2`. We can check 
the solution by plotting it:

```R 
plot(sol$t,sol$u,"l")
```

![linear_ode](https://user-images.githubusercontent.com/1814174/39011970-e04f1fe8-43c7-11e8-8da3-848362691783.png)

### Systems of ODEs

Now let's solve the Lorenz equations. In this case, our initial condition is a vector and our derivative functions
takes in the vector to return a vector (note: arbitrary dimensional arrays are allowed). We would define this as:

```R
f <- function(u,p,t) {
  du1 = p[1]*(u[2]-u[1])
  du2 = u[1]*(p[2]-u[3]) - u[2]
  du3 = u[1]*u[2] - p[3]*u[3]
  return(c(du1,du2,du3))
}
```

Here we utilized the parameter array `p`. Thus we use `diffeqr::ode.solve` like before, but also pass in parameters this time:

```R
u0 <- c(1.0,0.0,0.0)
tspan <- list(0.0,100.0)
p <- c(10.0,28.0,8/3)
prob <- de$ODEProblem(f, u0, tspan, p)
sol <- de$solve(prob)
```

The returned solution is like before except now `sol$u` is an array of arrays,
where `sol$u[i]` is the full system at time `sol$t[i]`. It can be convenient to
turn this into an R matrix through `sapply`:

```R
mat <- sapply(sol$u,identity)
```

This has each row as a time series. `t(mat)` makes each column a time series.
It is sometimes convenient to turn the output into a `data.frame` which is done
via:

```R
udf <- as.data.frame(t(mat))
```

Now we can use `matplot` to plot the timeseries together:

```R
matplot(sol$t,udf,"l",col=1:3)
```

![timeseries](https://user-images.githubusercontent.com/1814174/39012314-ef7a8fe2-43c8-11e8-9dde-1a8b87d3cfa4.png)

Now we can use the Plotly package to draw a phase plot:

```R
plotly::plot_ly(udf, x = ~V1, y = ~V2, z = ~V3, type = 'scatter3d', mode = 'lines')
```

![plotly_plot](https://user-images.githubusercontent.com/1814174/39012384-27ee7262-43c9-11e8-84d2-1edf937288ae.png)

Plotly is much prettier!

### Option Handling

If we want to have a more accurate solution, we can send `abstol` and `reltol`. Defaults are `1e-6` and `1e-3` respectively.
Generally you can think of the digits of accuracy as related to 1 plus the exponent of the relative tolerance, so the default is
two digits of accuracy. Absolute tolernace is the accuracy near 0. 

In addition, we may want to choose to save at more time points. We do this by giving an array of values to save at as `saveat`.
Together, this looks like:

```R
abstol <- 1e-8
reltol <- 1e-8
saveat <- 0:10000/100
sol <- de$solve(prob,abstol=abstol,reltol=reltol,saveat=saveat)
udf <- as.data.frame(t(sapply(sol$u,identity)))
plotly::plot_ly(udf, x = ~V1, y = ~V2, z = ~V3, type = 'scatter3d', mode = 'lines')
```

![precise_solution](https://user-images.githubusercontent.com/1814174/39012651-e03124e6-43c9-11e8-8496-bbee87987a37.png)

We can also choose to use a different algorithm. The choice is done using a string that matches the Julia syntax. See
[the ODE tutorial for details](http://diffeq.sciml.ai/dev/tutorials/ode_example.html#Choosing-a-Solver-Algorithm-1).
The list of choices for ODEs can be found at the [ODE Solvers page](http://diffeq.sciml.ai/dev/solvers/ode_solve.html).
For example, let's use a 9th order method due to Verner:

```R
sol <- de$solve(prob,de$Vern9(),abstol=abstol,reltol=reltol,saveat=saveat)
```

Note that each algorithm choice will cause a JIT compilation.

## Performance Enhancements

One way to enhance the performance of your code is to define the function in Julia 
so that way it is JIT compiled. diffeqr is built using 
[the JuliaCall package](https://github.com/Non-Contradiction/JuliaCall), and so 
you can utilize the Julia JIT compiler. We expose this automatically over ODE
functions via `jitoptimize_ode`, like in the following example:

```R
f <- function(u,p,t) {
  du1 = p[1]*(u[2]-u[1])
  du2 = u[1]*(p[2]-u[3]) - u[2]
  du3 = u[1]*u[2] - p[3]*u[3]
  return(c(du1,du2,du3))
}
u0 <- c(1.0,0.0,0.0)
tspan <- c(0.0,100.0)
p <- c(10.0,28.0,8/3)
prob <- de$ODEProblem(f, u0, tspan, p)
fastprob <- diffeqr::jitoptimize_ode(de,prob)
sol <- de$solve(fastprob,de$Tsit5())
```

Note that the first evaluation of the function will have an ~2 second lag since 
the compiler will run, and all subsequent runs will be orders of magnitude faster
than the pure R function. This means it's great for expensive functions (ex. large
PDEs) or functions called repeatedly, like during optimization of parameters.

We can also use the JuliaCall functions to directly define the function in Julia
to eliminate the R interpreter overhead and get full JIT compilation:

```R
julf <- JuliaCall::julia_eval("
function julf(du,u,p,t)
  du[1] = 10.0*(u[2]-u[1])
  du[2] = u[1]*(28.0-u[3]) - u[2]
  du[3] = u[1]*u[2] - (8/3)*u[3]
end")
JuliaCall::julia_assign("u0", u0)
JuliaCall::julia_assign("p", p)
JuliaCall::julia_assign("tspan", tspan)
prob3 = JuliaCall::julia_eval("ODEProblem(julf, u0, tspan, p)")
sol = de$solve(prob3,de$Tsit5())
```

To demonstrate the performance advantage, let's time them all:

```R
> system.time({ for (i in 1:100){ de$solve(prob    ,de$Tsit5()) }})
   user  system elapsed 
   6.69    0.06    6.78 
> system.time({ for (i in 1:100){ de$solve(fastprob,de$Tsit5()) }})
   user  system elapsed 
   0.11    0.03    0.14 
> system.time({ for (i in 1:100){ de$solve(prob3   ,de$Tsit5()) }})
   user  system elapsed 
   0.14    0.02    0.15 
```

This is about a 50x improvement!

## Stochastic Differential Equation (SDE) Examples

### 1D SDEs

Solving stochastic differential equations (SDEs) is the similar to ODEs. To solve an SDE, you use `diffeqr::sde.solve` and give
two functions: `f` and `g`, where `du = f(u,t)dt + g(u,t)dW_t`

```r
de <- diffeqr::diffeq_setup()
f <- function(u,p,t) {
  return(1.01*u)
}
g <- function(u,p,t) {
  return(0.87*u)
}
u0 <- 1/2
tspan <- list(0.0,1.0)
prob <- de$SDEProblem(f,g,u0,tspan)
sol <- de$solve(prob)
udf <- as.data.frame(t(sapply(sol$u,identity)))
plotly::plot_ly(udf, x = sol$t, y = sol$u, type = 'scatter', mode = 'lines')
```

![geometric_sdes](https://user-images.githubusercontent.com/1814174/39020683-9ea22b56-43e2-11e8-97f5-0c2a3ea69a2e.png)

### Systems of Diagonal Noise SDEs

Let's add diagonal multiplicative noise to the Lorenz attractor. diffeqr defaults to diagonal noise when a system of
equations is given. This is a unique noise term per system variable. Thus we generalize our previous functions as
follows:

```R
f <- function(u,p,t) {
  du1 = p[1]*(u[2]-u[1])
  du2 = u[1]*(p[2]-u[3]) - u[2]
  du3 = u[1]*u[2] - p[3]*u[3]
  return(c(du1,du2,du3))
}
g <- function(u,p,t) {
  return(c(0.3*u[1],0.3*u[2],0.3*u[3]))
}
u0 <- c(1.0,0.0,0.0)
tspan <- c(0.0,1.0)
p <- c(10.0,28.0,8/3)
prob <- de$SDEProblem(f,g,u0,tspan,p)
sol <- de$solve(prob,saveat=0.005)
udf <- as.data.frame(t(sapply(sol$u,identity)))
plotly::plot_ly(udf, x = ~V1, y = ~V2, z = ~V3, type = 'scatter3d', mode = 'lines')
```

Using a JIT compiled function for the drift and diffusion functions can greatly enhance the speed here.
With the speed increase we can comfortably solve over long time spans:

```R
tspan <- c(0.0,100.0)
prob <- de$SDEProblem(f,g,u0,tspan,p)
fastprob <- diffeqr::jitoptimize_sde(de,prob)
sol <- de$solve(fastprob,saveat=0.005)
udf <- as.data.frame(t(sapply(sol$u,identity)))
plotly::plot_ly(udf, x = ~V1, y = ~V2, z = ~V3, type = 'scatter3d', mode = 'lines')
```

![stochastic_lorenz](https://user-images.githubusercontent.com/1814174/39019723-216c3210-43df-11e8-82c0-2e676f53e235.png)

Let's see how much faster the JIT-compiled version was:

```R
> system.time({ for (i in 1:5){ de$solve(prob    ) }})
   user  system elapsed 
 146.40    0.75  147.22 
> system.time({ for (i in 1:5){ de$solve(fastprob) }})
   user  system elapsed 
   1.07    0.10    1.17
```

Holy Monster's Inc. that's about 145x faster.

### Systems of SDEs with Non-Diagonal Noise

In many cases you may want to share noise terms across the system. This is known as non-diagonal noise. The 
[DifferentialEquations.jl SDE Tutorial](http://diffeq.sciml.ai/dev/tutorials/sde_example.html#Example-4:-Systems-of-SDEs-with-Non-Diagonal-Noise-1)
explains how the matrix form of the diffusion term corresponds to the summation style of multiple Wiener processes. Essentially,
the row corresponds to which system the term is applied to, and the column is which noise term. So `du[i,j]` is the amount of
noise due to the `j`th Wiener process that's applied to `u[i]`. We solve the Lorenz system with correlated noise as follows:

```R
f <- JuliaCall::julia_eval("
function f(du,u,p,t)
  du[1] = 10.0*(u[2]-u[1])
  du[2] = u[1]*(28.0-u[3]) - u[2]
  du[3] = u[1]*u[2] - (8/3)*u[3]
end")
g <- JuliaCall::julia_eval("
function g(du,u,p,t)
  du[1,1] = 0.3u[1]
  du[2,1] = 0.6u[1]
  du[3,1] = 0.2u[1]
  du[1,2] = 1.2u[2]
  du[2,2] = 0.2u[2]
  du[3,2] = 0.3u[2]
end")
u0 <- c(1.0,0.0,0.0)
tspan <- c(0.0,100.0)
noise_rate_prototype <- matrix(c(0.0,0.0,0.0,0.0,0.0,0.0), nrow = 3, ncol = 2)

JuliaCall::julia_assign("u0", u0)
JuliaCall::julia_assign("tspan", tspan)
JuliaCall::julia_assign("noise_rate_prototype", noise_rate_prototype)
prob <- JuliaCall::julia_eval("SDEProblem(f, g, u0, tspan, p, noise_rate_prototype=noise_rate_prototype)")
sol <- de$solve(prob)
udf <- as.data.frame(t(sapply(sol$u,identity)))
plotly::plot_ly(udf, x = ~V1, y = ~V2, z = ~V3, type = 'scatter3d', mode = 'lines')
```

![noise_corr](https://user-images.githubusercontent.com/1814174/39022036-8958319a-43e8-11e8-849b-c21bcb2ec21e.png)

Here you can see that the warping effect of the noise correlations is quite visible!
Note that we applied JIT compilation since it's quite necessary for any difficult
stochastic example.

## Differential-Algebraic Equation (DAE) Examples

A differential-algebraic equation is defined by an implicit function `f(du,u,p,t)=0`. All of the controls are the
same as the other examples, except here you define a function which returns the residuals for each part of the equation
to define the DAE. The initial value `u0` and the initial derivative `du0` are required, though they do not necessarily
have to satisfy `f` (known as inconsistent initial conditions). The methods will automatically find consistent initial
conditions. In order for this to occur, `differential_vars` must be set. This vector states which of the variables are
differential (have a derivative term), with `false` meaning that the variable is purely algebraic.

This example shows how to solve the Robertson equation:

```R
f <- function (du,u,p,t) {
  resid1 = - 0.04*u[1]              + 1e4*u[2]*u[3] - du[1]
  resid2 = + 0.04*u[1] - 3e7*u[2]^2 - 1e4*u[2]*u[3] - du[2]
  resid3 = u[1] + u[2] + u[3] - 1.0
  c(resid1,resid2,resid3)
}
u0 <- c(1.0, 0, 0)
du0 <- c(-0.04, 0.04, 0.0)
tspan <- c(0.0,100000.0)
differential_vars <- c(TRUE,TRUE,FALSE)
prob <- de$DAEProblem(f,du0,u0,tspan,differential_vars=differential_vars)
sol <- de$solve(prob)
udf <- as.data.frame(t(sapply(sol$u,identity)))
plotly::plot_ly(udf, x = sol$t, y = ~V1, type = 'scatter', mode = 'lines') %>%
  plotly::add_trace(y = ~V2) %>%
  plotly::add_trace(y = ~V3)
```

Additionally, an in-place JIT compiled form for `f` can be used to enhance the speed:

```R
f = JuliaCall::julia_eval("function f(out,du,u,p,t)
  out[1] = - 0.04u[1]              + 1e4*u[2]*u[3] - du[1]
  out[2] = + 0.04u[1] - 3e7*u[2]^2 - 1e4*u[2]*u[3] - du[2]
  out[3] = u[1] + u[2] + u[3] - 1.0
end")
u0 <- c(1.0, 0, 0)
du0 <- c(-0.04, 0.04, 0.0)
tspan <- c(0.0,100000.0)
differential_vars <- c(TRUE,TRUE,FALSE)
JuliaCall::julia_assign("du0", du0)
JuliaCall::julia_assign("u0", u0)
JuliaCall::julia_assign("p", p)
JuliaCall::julia_assign("tspan", tspan)
JuliaCall::julia_assign("differential_vars", differential_vars)
prob = JuliaCall::julia_eval("DAEProblem(f, du0, u0, tspan, p, differential_vars=differential_vars)")
sol = de$solve(prob)
```

![daes](https://user-images.githubusercontent.com/1814174/39022955-d600814c-43ec-11e8-91bb-e096ff3d3fb7.png)

## Delay Differential Equation (DDE) Examples

A delay differential equation is an ODE which allows the use of previous values. In this case, the function
needs to be a JIT compiled Julia function. It looks just like the ODE, except in this case there is a function
`h(p,t)` which allows you to interpolate and grab previous values.

We must provide a history function `h(p,t)` that gives values for `u` before
`t0`. Here we assume that the solution was constant before the initial time point. Additionally, we pass
`constant_lags = c(20.0)` to tell the solver that only constant-time lags were used and what the lag length
was. This helps improve the solver accuracy by accurately stepping at the points of discontinuity. Together
this is:

```R
f <- JuliaCall::julia_eval("function f(du, u, h, p, t)
  du[1] = 1.1/(1 + sqrt(10)*(h(p, t-20)[1])^(5/4)) - 10*u[1]/(1 + 40*u[2])
  du[2] = 100*u[1]/(1 + 40*u[2]) - 2.43*u[2]
end")
h <- JuliaCall::julia_eval("function h(p, t)
  [1.05767027/3, 1.030713491/3]
end")
u0 <- c(1.05767027/3, 1.030713491/3)
tspan <- c(0.0, 100.0)
constant_lags <- c(20.0)
JuliaCall::julia_assign("u0", u0)
JuliaCall::julia_assign("tspan", tspan)
JuliaCall::julia_assign("constant_lags", tspan)
prob <- JuliaCall::julia_eval("DDEProblem(f, u0, h, tspan, constant_lags = constant_lags)")
sol <- de$solve(prob,de$MethodOfSteps(de$Tsit5()))
udf <- as.data.frame(t(sapply(sol$u,identity)))
plotly::plot_ly(udf, x = sol$t, y = ~V1, type = 'scatter', mode = 'lines') %>% plotly::add_trace(y = ~V2)
```

![delay](https://user-images.githubusercontent.com/1814174/39023532-10bdd750-43f0-11e8-837d-156d33ea2f99.png)

Notice that the solver accurately is able to simulate the kink (discontinuity) at `t=20` due to the discontinuity
of the derivative at the initial time point! This is why declaring discontinuities can enhance the solver accuracy.
