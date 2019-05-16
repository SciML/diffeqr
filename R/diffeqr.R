#' Setup diffeqr
#'
#' This function initializes Julia and the DifferentialEquations.jl package.
#' The first time will be long since it includes precompilation.
#' Additonally, this will install Julia and the required packages
#' if they are missing.
#'
#' @param pkg.check logical, check for DifferentialEquations.jl package and install if necessary
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
diffeq_setup <- function (pkg_check=FALSE, ...){
  if is.null(julia_locate()){
    JuliaCall::julia_setup()
  }
  
  julia <- JuliaCall::julia_setup(...)
  if(pkg_check) JuliaCall::julia_install_package_if_needed("DifferentialEquations")
  JuliaCall::julia_library("DifferentialEquations")

  JuliaCall::julia_install_package_if_needed("ModelingToolkit")
  JuliaCall::julia_library("ModelingToolkit")

  functions <- JuliaCall::julia_eval("filter(isascii, replace.(string.(propertynames(DifferentialEquations)),\"!\"=>\"_bang\"))")
  de <- julia_pkg_import("DifferentialEquations",functions)
  JuliaCall::autowrap("DiffEqBase.AbstractODESolution", fields = c("t","u"))
  JuliaCall::autowrap("DiffEqBase.AbstractRODESolution", fields = c("t","u"))
  JuliaCall::autowrap("DiffEqBase.AbstractDDESolution", fields = c("t","u"))
  JuliaCall::autowrap("DiffEqBase.AbstractDAESolution", fields = c("t","u"))
  JuliaCall::autowrap("DiffEqBase.EnsembleSolution", fields = c("t","u"))
  de
}

julia_locate <- do.call(":::", list("JuliaCall", quote(julia_locate)))

#' Jit Optimize an ODEProblem
#'
#' This function JIT Optimizes and ODEProblem utilizing the Julia ModelingToolkit
#' and JIT compiler.
#'
#' @param de the current diffeqr environment
#' @param prob an ODEProblem
#'
#' @examples
#'
#' \donttest{ ## diffeq_setup() is time-consuming and requires Julia+DifferentialEquations.jl
#' de <- diffeqr::diffeq_setup()
#' f <- function(u,p,t) {
#'   du1 = p[1]*(u[2]-u[1])
#'   du2 = u[1]*(p[2]-u[3]) - u[2]
#'   du3 = u[1]*u[2] - p[3]*u[3]
#'   return(c(du1,du2,du3))
#' }
#' u0 <- c(1.0,0.0,0.0)
#' tspan <- c(0.0,100.0)
#' p <- c(10.0,28.0,8/3)
#' prob <- de$ODEProblem(f, u0, tspan, p)
#' fastprob <- diffeqr::jitoptimize_ode(de,prob)
#' sol <- de$solve(fastprob,de$Tsit5())
#' }
#'
#' @export
jitoptimize_ode <- function (de,prob){
  odesys = de$modelingtoolkitize(prob)
  JuliaCall::julia_assign("odesys", odesys)
  jul_f = JuliaCall::julia_eval("jitf = ODEFunction(odesys,jac=true)")
  JuliaCall::julia_assign("u0", prob$u0)
  JuliaCall::julia_assign("p", prob$p)
  JuliaCall::julia_assign("tspan", prob$tspan)
  new_prob <- JuliaCall::julia_eval("ODEProblem(jitf, u0, tspan, p)")
}

#' Jit Optimize an SDEProblem
#'
#' This function JIT Optimizes and SDEProblem utilizing the Julia ModelingToolkit
#' and JIT compiler.
#'
#' @param de the current diffeqr environment
#' @param prob an SDEProblem
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
jitoptimize_sde <- function (de,prob){
  sdesys = de$modelingtoolkitize(prob)
  JuliaCall::julia_assign("sdesys", sdesys)
  jul_f = JuliaCall::julia_eval("jitf = SDEFunction(sdesys,jac=true)")
  JuliaCall::julia_assign("u0", prob$u0)
  JuliaCall::julia_assign("p", prob$p)
  JuliaCall::julia_assign("tspan", prob$tspan)
  new_prob <- JuliaCall::julia_eval("SDEProblem(jitf, jitf.g, u0, tspan, p)")
}

#' Setup DiffEqGPU
#'
#' This function initializes the DiffEqGPU package for GPU-parallelized ensembles.
#' The first time will be long since it includes precompilation.
#'
#' @examples
#'
#' \donttest{ ## diffeq_setup() is time-consuming and requires Julia+DifferentialEquations.jl
#'
#' degpu <- diffeqr::diffeqgpu_setup()
#'
#' }
#'
#' @export
diffeqgpu_setup <- function (){
  JuliaCall::julia_install_package_if_needed("DiffEqGPU")
  JuliaCall::julia_library("DiffEqGPU")
  functions <- JuliaCall::julia_eval("filter(isascii, replace.(string.(propertynames(DiffEqGPU)),\"!\"=>\"_bang\"))")
  degpu <- julia_pkg_import("DiffEqGPU",functions)
}

julia_function <- function(func_name, pkg_name = "Main",
                           env = emptyenv()){
  fname <- paste0(pkg_name, ".", func_name)
  force(fname)
  f <- function(...,
                need_return = c("R", "Julia", "None"),
                show_value = FALSE){
    if (!isTRUE(env$initialized)) {
      env$setup()
    }
    JuliaCall::julia_do.call(func_name = fname, list(...),
                             need_return = match.arg(need_return),
                             show_value = show_value)
  }
  force(f)
  env[[func_name]] <- f
}

julia_pkg_import <- function(pkg_name, func_list){
  env <- new.env(parent = emptyenv())
  env$setup <- function(...){
    JuliaCall::julia_setup(...)
    JuliaCall::julia_library(pkg_name)
    env$initialized <- TRUE
  }
  for (fname in func_list) {
    julia_function(func_name = fname,
                   pkg_name = pkg_name,
                   env = env)
  }
  env
}
