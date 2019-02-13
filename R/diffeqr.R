#' Setup diffeqr
#'
#' This function initializes Julia and the DifferentialEquations.jl package.
#' The first time will be long since it includes precompilation.
#'
#' @param ... Parameters are passed down to JuliaCall::julia_setup
#'
#' @examples
#'
#' \donttest{ ## diffeq_setup() is time-consuming and requires Julia+DifferentialEquations.jl
#'
#' diffeqr::diffeq_setup()
#'
#' }
#'
#' @export
diffeq_setup <- function (...){
  julia <- JuliaCall::julia_setup(...)
  JuliaCall::julia_install_package_if_needed("DifferentialEquations")
  JuliaCall::julia_library("DifferentialEquations")
}

#' Solve Ordinary Differential Equations (ODE)
#'
#' Solves an ODE with u'=f(u,p,t), for u(0)=u0 over the tspan
#' @param f the derivative function.
#' @param u0 the initial condition. Can be a number or (arbitrary dimension) array.
#' @param tspan the timespan to solve over. Should be a list of two values: (initial time, end time).
#' @param p the parameters. Defaults to no parameters. Can be a number or an array.
#' @param alg the algorithm used to solve the differential equation. Defaults to an adaptive choice.
#'        Algorithm choices are done through a string which matches the DifferentialEquations.jl form.
#' @param reltol the relative tolerance of the ODE solver. Defaults to 1e-3.
#' @param abstol the absolute tolerance of the ODE solver. Defaults to 1e-6.
#' @param maxiters the maximum number of iterations the adaptive solver is allowed to try before exiting.
#'        Defualt value is 1000000.
#' @param saveat the time points to save values at. Should be an array of times. Defaults to automatic.
#'
#' @return sol. Has the sol$t for the time points and sol$u for the values.
#'
#' @examples
#'
#' \donttest{ ## diffeq_setup() is time-consuming and requires Julia+DifferentialEquations.jl
#'
#' diffeqr::diffeq_setup()
#'
#' # Scalar ODEs
#'
#' f <- function(u,p,t) {
#' return(1.01*u)
#' }
#' u0 = 1/2
#' tspan <- list(0.0,1.0)
#' sol = diffeqr::ode.solve(f,u0,tspan)
#' plot(sol$t,sol$u,"l")
#'
#' saveat=1:10/10
#' sol2 = diffeqr::ode.solve(f,u0,tspan,saveat=saveat)
#' sol3 = diffeqr::ode.solve(f,u0,tspan,alg="Vern9()")
#' sol4 = diffeqr::ode.solve(f,u0,tspan,alg="Rosenbrock23()")
#'
#' # Systems of ODEs
#'
#' f <- function(u,p,t) {
#'   du1 = p[1]*(u[2]-u[1])
#'   du2 = u[1]*(p[2]-u[3]) - u[2]
#'   du3 = u[1]*u[2] - p[3]*u[3]
#' return(c(du1,du2,du3))
#' }
#'
#' u0 = c(1.0,0.0,0.0)
#' tspan <- list(0.0,100.0)
#' p = c(10.0,28.0,8/3)
#' sol = diffeqr::ode.solve(f,u0,tspan,p=p)
#' udf = as.data.frame(sol$u)
#' matplot(sol$t,udf,"l",col=1:3)
#' #plotly::plot_ly(udf, x = ~V1, y = ~V2, z = ~V3, type = 'scatter3d', mode = 'lines')
#'
#' f <- JuliaCall::julia_eval("
#' function f(du,u,p,t)
#'  du[1] = 10.0*(u[2]-u[1])
#'  du[2] = u[1]*(28.0-u[3]) - u[2]
#'  du[3] = u[1]*u[2] - (8/3)*u[3]
#' end")
#' sol = diffeqr::ode.solve('f',u0,tspan)
#'
#' }
#'
#' @export
ode.solve <- function(f,u0,tspan,p=NULL,alg="nothing",reltol=1e-3,abstol=1e-6,maxiters = 1000000,saveat=NULL){
  if (is.character(f)){
    fname = f
  } else {
    JuliaCall::julia_assign("___f", f)
    fname = "___f"
  }
  JuliaCall::julia_assign("u0", u0)
  tspan_tup = tspan
  class(tspan_tup) <- "JuliaTuple"
  JuliaCall::julia_assign("tspan", tspan_tup)
  if (is.null(p)){
    p_str = "nothing"
  } else {
    p_str = "p"
    JuliaCall::julia_assign("p", p)
  }
  if (is.null(saveat)){
    saveat_str = "eltype(prob.tspan)[]"
  } else {
    saveat_str = "saveat"
    JuliaCall::julia_assign("saveat", saveat)
  }
  jleval = stringr::str_interp("prob = ODEProblem(${fname},u0,tspan,${p_str})")
  JuliaCall::julia_eval(jleval)
  jleval = stringr::str_interp("sol = solve(prob,${alg},reltol=${reltol},abstol=${abstol}, maxiters = ${maxiters}, saveat=${saveat_str}); nothing")
  JuliaCall::julia_eval(jleval)
  u = JuliaCall::julia_eval("typeof(u0)<:Number ? Array(sol) : sol'")
  t = JuliaCall::julia_eval("sol.t")
  list(u=u,t=t)
}

#' Solve Stochastic Differential Equations (SDE)
#'
#' Solves an SDE with du=f(u,p,t)dt + g(u,p,t)dW_t, for u(0)=u0 over the tspan
#' @param f the drift function.
#' @param g the diffusion function.
#' @param u0 the initial condition. Can be a number or (arbitrary dimension) array.
#' @param tspan the timespan to solve over. Should be a list of two values: (initial time, end time).
#' @param p the parameters. Defaults to no parameters. Can be a number or an array.
#' @param alg the algorithm used to solve the differential equation. Defaults to an adaptive choice.
#'        Algorithm choices are done through a string which matches the DifferentialEquations.jl form.
#' @param reltol the relative tolerance of the ODE solver. Defaults to 1e-3.
#' @param abstol the absolute tolerance of the ODE solver. Defaults to 1e-6.
#' @param maxiters the maximum number of iterations the adaptive solver is allowed to try before exiting.
#'        Defualt value is 1000000.
#' @param saveat the time points to save values at. Should be an array of times. Defaults to automatic.
#' @param noise.dims list of the dimensions for the noise rate term. Defaults to NULL which gives diagonal noise.
#'
#' @return sol. Has the sol$t for the time points and sol$u for the values.
#'
#' @examples
#'
#' \donttest{ ## diffeq_setup() is time-consuming and requires Julia+DifferentialEquations.jl
#'
#' diffeqr::diffeq_setup()
#'
#' # Scalar SDEs
#'
#' f <- function(u,p,t) {
#'   return(1.01*u)
#' }
#' g <- function(u,p,t) {
#'   return(0.87*u)
#' }
#' u0 = 1/2
#' tspan <- list(0.0,1.0)
#' sol = diffeqr::sde.solve(f,g,u0,tspan)
#' #plotly::plot_ly(udf, x = sol$t, y = sol$u, type = 'scatter', mode = 'lines')
#'
#' # Diagonal Noise SDEs
#'
#' f <- JuliaCall::julia_eval("
#' function f(du,u,p,t)
#'   du[1] = 10.0*(u[2]-u[1])
#'   du[2] = u[1]*(28.0-u[3]) - u[2]
#'   du[3] = u[1]*u[2] - (8/3)*u[3]
#' end")
#'
#' g <- JuliaCall::julia_eval("
#' function g(du,u,p,t)
#'   du[1] = 0.3*u[1]
#'   du[2] = 0.3*u[2]
#'   du[3] = 0.3*u[3]
#' end")
#' tspan <- list(0.0,100.0)
#' sol = diffeqr::sde.solve('f','g',u0,tspan,p=p,saveat=0.05)
#' udf = as.data.frame(sol$u)
#' #plotly::plot_ly(udf, x = ~V1, y = ~V2, z = ~V3, type = 'scatter3d', mode = 'lines')
#'
#' # Non-Diagonal Noise SDEs
#'
#' f <- JuliaCall::julia_eval("
#' function f(du,u,p,t)
#'   du[1] = 10.0*(u[2]-u[1])
#'   du[2] = u[1]*(28.0-u[3]) - u[2]
#'   du[3] = u[1]*u[2] - (8/3)*u[3]
#' end")
#' g <- JuliaCall::julia_eval("
#' function g(du,u,p,t)
#'   du[1,1] = 0.3u[1]
#'   du[2,1] = 0.6u[1]
#'   du[3,1] = 0.2u[1]
#'   du[1,2] = 1.2u[2]
#'   du[2,2] = 0.2u[2]
#'   du[3,2] = 0.3u[2]
#' end")
#' u0 = c(1.0,0.0,0.0)
#' tspan <- list(0.0,100.0)
#' noise.dims = list(3,2)
#' sol = diffeqr::sde.solve('f','g',u0,tspan,saveat=0.005,noise.dims=noise.dims)
#' udf = as.data.frame(sol$u)
#' #plotly::plot_ly(udf, x = ~V1, y = ~V2, z = ~V3, type = 'scatter3d', mode = 'lines')
#'
#' }
#'
#' @export
sde.solve <- function(f,g,u0,tspan,p=NULL,alg="nothing",noise.dims=NULL,maxiters = 1000000,
                      reltol=1e-2,abstol=1e-2,saveat=NULL){
  if (is.character(f)){
    fname = f
  } else {
    JuliaCall::julia_assign("___f", f)
    fname = "___f"
  }
  if (is.character(g)){
    gname = g
  } else {
    JuliaCall::julia_assign("___g", g)
    gname = "___g"
  }
  JuliaCall::julia_assign("u0", u0)
  tspan_tup = tspan
  class(tspan_tup) <- "JuliaTuple"
  JuliaCall::julia_assign("tspan", tspan_tup)
  if (is.null(p)){
    p_str = "nothing"
  } else {
    p_str = "p"
    JuliaCall::julia_assign("p", p)
  }
  if (is.null(saveat)){
    saveat_str = "eltype(prob.tspan)[]"
  } else {
    saveat_str = "saveat"
    JuliaCall::julia_assign("saveat", saveat)
  }
  if (is.null(noise.dims)) {
    nrp_str = "nothing"
  } else {
    nrp_str = stringr::str_interp("zeros(${noise.dims[1]},${noise.dims[2]})")
  }
  jleval = stringr::str_interp("prob = SDEProblem(${fname},${gname},u0,tspan,${p_str},noise_rate_prototype=${nrp_str})")
  JuliaCall::julia_eval(jleval)
  jleval = stringr::str_interp("sol = solve(prob,${alg},reltol=${reltol},abstol=${abstol}, maxiters=${maxiters},saveat=${saveat_str}); nothing")
  JuliaCall::julia_eval(jleval)
  u = JuliaCall::julia_eval("typeof(u0)<:Number ? Array(sol) : sol'")
  t = JuliaCall::julia_eval("sol.t")
  list(u=u,t=t)
}

#' Solve Differential-Algebraic Equations (DAE)
#'
#' Solves a DAE with f(du,u,p,t)=0 for u(0)=u0 over the tspan
#' @param f the implicit ODE function.
#' @param du0 the initial derivative. Can be a number or (arbitrary dimension) array.
#' @param u0 the initial condition. Can be a number or (arbitrary dimension) array.
#' @param tspan the timespan to solve over. Should be a list of two values: (initial time, end time).
#' @param p the parameters. Defaults to no parameters. Can be a number or an array.
#' @param alg the algorithm used to solve the differential equation. Defaults to an adaptive choice.
#'        Algorithm choices are done through a string which matches the DifferentialEquations.jl form.
#' @param differential_vars boolean array declaring which variables are differential. All falses correspond to
#'        purely algebraic variables.
#' @param reltol the relative tolerance of the ODE solver. Defaults to 1e-3.
#' @param abstol the absolute tolerance of the ODE solver. Defaults to 1e-6.
#' @param maxiters the maximum number of iterations the adaptive solver is allowed to try before exiting.
#'        Defualt value is 1000000.
#' @param saveat the time points to save values at. Should be an array of times. Defaults to automatic.
#'
#' @return sol. Has the sol$t for the time points and sol$u for the values.
#'
#' @examples
#'
#' \donttest{ ## diffeq_setup() is time-consuming and requires Julia+DifferentialEquations.jl
#'
#' diffeqr::diffeq_setup()
#'
#' f <- function (du,u,p,t) {
#'   resid1 = - 0.04*u[1]              + 1e4*u[2]*u[3] - du[1]
#'   resid2 = + 0.04*u[1] - 3e7*u[2]^2 - 1e4*u[2]*u[3] - du[2]
#'   resid3 = u[1] + u[2] + u[3] - 1.0
#'   c(resid1,resid2,resid3)
#' }
#' u0 = c(1.0, 0, 0)
#' du0 = c(-0.04, 0.04, 0.0)
#' tspan = list(0.0,100000.0)
#' differential_vars = c(TRUE,TRUE,FALSE)
#' sol = diffeqr::dae.solve(f,du0,u0,tspan,differential_vars=differential_vars)
#' udf = as.data.frame(sol$u)
#' #plotly::plot_ly(udf, x = sol$t, y = ~V1, type = 'scatter', mode = 'lines') %>%
#' #plotly::add_trace(y = ~V2) %>%
#' #plotly::add_trace(y = ~V3)
#'
#' f = JuliaCall::julia_eval("function f(out,du,u,p,t)
#'   out[1] = - 0.04u[1]              + 1e4*u[2]*u[3] - du[1]
#'   out[2] = + 0.04u[1] - 3e7*u[2]^2 - 1e4*u[2]*u[3] - du[2]
#'   out[3] = u[1] + u[2] + u[3] - 1.0
#' end")
#' sol = diffeqr::dae.solve('f',du0,u0,tspan,differential_vars=differential_vars)
#'
#' }
#'
#' @export
dae.solve <- function(f,du0,u0,tspan,p=NULL,alg="nothing",reltol=1e-3,abstol=1e-6,maxiters = 1000000,saveat=NULL,differential_vars=NULL){
  if (is.character(f)){
    fname = f
  } else {
    JuliaCall::julia_assign("___f", f)
    fname = "___f"
  }
  JuliaCall::julia_assign("u0", u0)
  JuliaCall::julia_assign("du0", du0)
  tspan_tup = tspan
  class(tspan_tup) <- "JuliaTuple"
  JuliaCall::julia_assign("tspan", tspan_tup)
  if (is.null(p)){
    p_str = "nothing"
  } else {
    p_str = "p"
    JuliaCall::julia_assign("p", p)
  }
  if (is.null(saveat)){
    saveat_str = "eltype(prob.tspan)[]"
  } else {
    saveat_str = "saveat"
    JuliaCall::julia_assign("saveat", saveat)
  }
  if (is.null(differential_vars)){
    diffvar_str = "nothing"
  } else {
    diffvar_str = "diffvars"
    JuliaCall::julia_assign("diffvars",differential_vars)
  }
  jleval = stringr::str_interp("prob = DAEProblem(${fname},du0,u0,tspan,${p_str},differential_vars=${diffvar_str})")
  JuliaCall::julia_eval(jleval)
  jleval = stringr::str_interp("sol = solve(prob,${alg},reltol=${reltol},abstol=${abstol}, maxiters=${maxiters}, saveat=${saveat_str}); nothing")
  JuliaCall::julia_eval(jleval)
  u = JuliaCall::julia_eval("typeof(u0)<:Number ? Array(sol) : sol'")
  t = JuliaCall::julia_eval("sol.t")
  list(u=u,t=t)
}

#' Solve Delay Differential Equations (DDE)
#'
#' Solves a DDE with f(u,p,t)=0 for u(0)=u0 over the tspan
#' @param f the implicit ODE function.
#' @param u0 the initial condition. Can be a number or (arbitrary dimension) array.
#' @param h is the history function (p,t) which gives values of the solution before the initial time point.
#' @param tspan the timespan to solve over. Should be a list of two values: (initial time, end time).
#' @param p the parameters. Defaults to no parameters. Can be a number or an array.
#' @param alg the algorithm used to solve the differential equation. Defaults to an adaptive choice.
#'        Algorithm choices are done through a string which matches the DifferentialEquations.jl form.
#' @param reltol the relative tolerance of the ODE solver. Defaults to 1e-3.
#' @param abstol the absolute tolerance of the ODE solver. Defaults to 1e-6.
#' @param maxiters the maximum number of iterations the adaptive solver is allowed to try before exiting.
#'        Defualt value is 1000000.
#' @param saveat the time points to save values at. Should be an array of times. Defaults to automatic.
#' @param constant_lags a vector of floats for the constant-time lags. Defaults to NULL.
#'
#' @return sol. Has the sol$t for the time points and sol$u for the values.
#'
#' @examples
#'
#' \donttest{ ## diffeq_setup() is time-consuming and requires Julia+DifferentialEquations.jl
#'
#' diffeqr::diffeq_setup()
#'
#' f = JuliaCall::julia_eval("function f(du, u, h, p, t)
#'   du[1] = 1.1/(1 + sqrt(10)*(h(p, t-20)[1])^(5/4)) - 10*u[1]/(1 + 40*u[2])
#'   du[2] = 100*u[1]/(1 + 40*u[2]) - 2.43*u[2]
#' end")
#' u0 = c(1.05767027/3, 1.030713491/3)
#' h <- function (p,t){
#'   c(1.05767027/3, 1.030713491/3)
#' }
#' tspan = list(0.0, 100.0)
#' constant_lags = c(20.0)
#' sol = diffeqr::dde.solve('f',u0,h,tspan,constant_lags=constant_lags)
#' udf = as.data.frame(sol$u)
#' #plotly::plot_ly(udf, x = sol$t, y = ~V1, type = 'scatter', mode = 'lines') %>%
#' #plotly::add_trace(y = ~V2)
#' }
#'
#' @export
dde.solve <- function(f,u0,h,tspan,p=NULL,alg="nothing",reltol=1e-3,abstol=1e-6,maxiters = 1000000,saveat=NULL,constant_lags=NULL){
  if (is.character(f)){
    fname = f
  } else {
    JuliaCall::julia_assign("___f", f)
    fname = "___f"
  }
  if (is.character(h)){
    hname = h
  } else {
    JuliaCall::julia_assign("___h", h)
    hname = "___h"
  }

  JuliaCall::julia_assign("u0", u0)
  tspan_tup = tspan
  class(tspan_tup) <- "JuliaTuple"
  JuliaCall::julia_assign("tspan", tspan_tup)
  if (is.null(p)){
    p_str = "nothing"
  } else {
    p_str = "p"
    JuliaCall::julia_assign("p", p)
  }
  if (is.null(saveat)){
    saveat_str = "eltype(prob.tspan)[]"
  } else {
    saveat_str = "saveat"
    JuliaCall::julia_assign("saveat", saveat)
  }
  if (is.null(constant_lags)){
    cl_str = "nothing"
  } else {
    cl_str = "cl"
    JuliaCall::julia_assign("cl",constant_lags)
  }
  jleval = stringr::str_interp("prob = DDEProblem(${fname},u0,${hname},tspan,${p_str},constant_lags=${cl_str})")
  JuliaCall::julia_eval(jleval)
  jleval = stringr::str_interp("sol = solve(prob,${alg},reltol=${reltol},abstol=${abstol}, maxiters=${maxiters},saveat=${saveat_str}); nothing")
  JuliaCall::julia_eval(jleval)
  u = JuliaCall::julia_eval("typeof(u0)<:Number ? Array(sol) : sol'")
  t = JuliaCall::julia_eval("sol.t")
  list(u=u,t=t)
}
