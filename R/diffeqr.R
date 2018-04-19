library(JuliaCall)
library(stringr)

julia <- julia_setup()
julia_install_package_if_needed("DifferentialEquations")
julia_library("DifferentialEquations")

#' Solve Ordinary Differential Equations (ODE)
#'
#' Solves an ODE with u'=f(u,p,t), for u(0)=u0 over the tspan
#' @param f the derivative function.
#' @param u0 the initial condition. Can be a number or (arbitrary dimension) array.
#' @param tspan the timespan to solve over. Should be a list of two values: (initial time, end time).
#' @param p the parameters. Defaults to no parameters. Can be a number or an array.
#' @param reltol the relative tolerance of the ODE solver. Defaults to 1e-3.
#' @param abstol the absolute tolerance of the ODE solver. Defaults to 1e-6
#' @param saveat the time points to save values at. Should be an array of times. Defaults to automatic.
#' @param fname the name of the Julia function. Defaults to automatic, using the passed in f.
#'
#' @return sol. Has the sol$t for the time points and sol$u for the values.
#'
#' @export
ode.solve <- function(f,u0,tspan,p=NULL,alg="nothing",fname="___f",reltol=1e-3,abstol=1e-6,saveat=NULL){
  julia_assign("___f", f)
  julia_assign("u0", u0)
  tspan_tup = tspan
  class(tspan_tup) <- "JuliaTuple"
  julia_assign("tspan", tspan_tup)
  if (is.null(p)){
    p_str = "nothing"
  } else {
    p_str = "p"
    julia_assign("p", p)
  }
  if (is.null(saveat)){
    saveat_str = "eltype(prob.tspan)[]"
  } else {
    saveat_str = "saveat"
    julia_assign("saveat", saveat)
  }
  jleval = stringr::str_interp("prob = ODEProblem(${fname},u0,tspan,${p_str})")
  julia_eval(jleval)
  jleval = stringr::str_interp("sol = solve(prob,${alg},reltol=${reltol},abstol=${abstol}, saveat=${saveat_str}); nothing")
  julia_eval(jleval)
  u = julia_eval("typeof(u0)<:Number ? Array(sol) : sol'")
  t = julia_eval("sol.t")
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
#' @param reltol the relative tolerance of the ODE solver. Defaults to 1e-3.
#' @param abstol the absolute tolerance of the ODE solver. Defaults to 1e-6
#' @param saveat the time points to save values at. Should be an array of times. Defaults to automatic.
#' @param fname the name of the Julia function. Defaults to automatic, using the passed in f.
#'
#' @return sol. Has the sol$t for the time points and sol$u for the values.
#'
#' @export
sde.solve <- function(f,g,u0,tspan,p=NULL,alg="nothing",fname="___f",gname="___g",noise.dims=NULL,reltol=1e-3,abstol=1e-6,saveat=NULL){
  julia_assign("___f", f)
  julia_assign("___g", g)
  julia_assign("u0", u0)
  tspan_tup = tspan
  class(tspan_tup) <- "JuliaTuple"
  julia_assign("tspan", tspan_tup)
  if (is.null(p)){
    p_str = "nothing"
  } else {
    p_str = "p"
    julia_assign("p", p)
  }
  if (is.null(saveat)){
    saveat_str = "eltype(prob.tspan)[]"
  } else {
    saveat_str = "saveat"
    julia_assign("saveat", saveat)
  }
  if (is.null(noise.dims)) {
    nrp_str = "nothing"
  } else {
    nrp_str = stringr::str_interp("zeros(${noise.dims[1]},${noise.dims[2]})")
  }
  jleval = stringr::str_interp("prob = SDEProblem(${fname},${gname},u0,tspan,${p_str},noise_rate_prototype=${nrp_str})")
  julia_eval(jleval)
  jleval = stringr::str_interp("sol = solve(prob,${alg},reltol=${reltol},abstol=${abstol}, saveat=${saveat_str}); nothing")
  julia_eval(jleval)
  u = julia_eval("typeof(u0)<:Number ? Array(sol) : sol'")
  t = julia_eval("sol.t")
  list(u=u,t=t)
}
