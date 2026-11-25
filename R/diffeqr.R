#' Setup diffeqr
#'
#' This function initializes Julia and the DifferentialEquations.jl package.
#' The first time will be long since it includes precompilation.
#' Additionally, this will install Julia and the required packages
#' if they are missing. With DifferentialEquations.jl v8 and later, which no
#' longer bundles the SDE and DDE solvers, StochasticDiffEq.jl and
#' DelayDiffEq.jl are installed and loaded as well.
#'
#' With DifferentialEquations.jl v8 the returned object exposes fewer names than
#' with v7: for example \code{Rodas5}, \code{RK4}, \code{QNDF},
#' \code{JumpProblem} and \code{PresetTimeCallback} are not available through
#' it. To use one, install and load the package that owns it (for \code{Rodas5},
#' \code{JuliaCall::julia_install_package_if_needed("OrdinaryDiffEqRosenbrock")},
#' \code{JuliaCall::julia_library("OrdinaryDiffEqRosenbrock")} and then
#' \code{JuliaCall::julia_eval("Rodas5()")}), or stay on
#' DifferentialEquations.jl v7.
#'
#' @param pkg_check logical, check for the DifferentialEquations.jl package, and
#' with DifferentialEquations.jl v8 and later also StochasticDiffEq.jl and
#' DelayDiffEq.jl, and install them if necessary
#' @param ... Parameters are passed down to JuliaCall::julia_setup
#'
#' @examples
#'
#' \dontrun{ ## diffeq_setup() is time-consuming and requires Julia+DifferentialEquations.jl
#'
#' diffeqr::diffeq_setup()
#'
#' }
#'
#' @export
diffeq_setup <- function (pkg_check=TRUE,...){
  use_julia_libunwind(...)
  julia <- JuliaCall::julia_setup(installJulia=TRUE,...)
  if(pkg_check) JuliaCall::julia_install_package_if_needed("DifferentialEquations")
  JuliaCall::julia_library("DifferentialEquations")

  functions <- JuliaCall::julia_eval("filter(isascii, replace.(string.(propertynames(DifferentialEquations)),\"!\"=>\"_bang\"))")
  de <- julia_pkg_import("DifferentialEquations",functions)

  # DifferentialEquations.jl v8 no longer bundles the SDE and DDE solvers
  if(JuliaCall::julia_eval("pkgversion(DifferentialEquations) >= v\"8\"")){
    for(pkg in c("StochasticDiffEq","DelayDiffEq")){
      if(pkg_check) JuliaCall::julia_install_package_if_needed(pkg)
      JuliaCall::julia_library(pkg)
      pkg_functions <- JuliaCall::julia_eval(sprintf("filter(isascii, replace.(string.(propertynames(%s)),\"!\"=>\"_bang\"))",pkg))
      for(fname in setdiff(pkg_functions,functions)){
        julia_function(func_name = fname, pkg_name = pkg, env = de)
      }
      functions <- union(functions,pkg_functions)
    }
  }

  JuliaCall::autowrap("SciMLBase.AbstractODESolution", fields = c("t","u"))
  JuliaCall::autowrap("SciMLBase.AbstractRODESolution", fields = c("t","u"))
  JuliaCall::autowrap("SciMLBase.AbstractDDESolution", fields = c("t","u"))
  JuliaCall::autowrap("SciMLBase.AbstractDAESolution", fields = c("t","u"))
  JuliaCall::autowrap("SciMLBase.EnsembleSolution", fields = c("t","u"))
  de
}

julia_locate <- do.call(":::", list("JuliaCall", quote(julia_locate)))

# R puts the system library directory in LD_LIBRARY_PATH, which the loader searches
# before Julia's own, so Julia gets the system libunwind and segfaults.
# `...` matches julia_setup()'s arguments.
use_julia_libunwind <- function(JULIA_HOME = NULL, ...){
  if (Sys.info()[["sysname"]] != "Linux") return(invisible())
  JULIA_HOME <- julia_locate(JULIA_HOME)
  if (is.null(JULIA_HOME)) return(invisible())
  libunwind <- file.path(dirname(JULIA_HOME), "lib", "julia", "libunwind.so.8")
  if (!file.exists(libunwind)) return(invisible())
  # This process reads LD_LIBRARY_PATH only at startup, so load Julia's copy first.
  dyn.load(libunwind)
  # Julia processes started later (package installation, precompilation) must not
  # inherit the default system directories; the loader searches those anyway.
  dirs <- strsplit(Sys.getenv("LD_LIBRARY_PATH"), ":")[[1]]
  system_dir <- grepl("^(/usr)?/lib(64)?(/[^/]*-linux-gnu[^/]*)?/?$", dirs)
  Sys.setenv(LD_LIBRARY_PATH = paste(dirs[!system_dir], collapse = ":"))
  invisible()
}

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
#' \dontrun{ ## diffeq_setup() is time-consuming and requires Julia+DifferentialEquations.jl
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
  JuliaCall::julia_install_package_if_needed("ModelingToolkit")
  JuliaCall::julia_library("ModelingToolkit")
  functions <- JuliaCall::julia_eval("filter(isascii, replace.(string.(propertynames(ModelingToolkit)),\"!\"=>\"_bang\"))")

  # Can remove the de argument when breaking, but kept for backwards compat
  mtk <- julia_pkg_import("ModelingToolkit",functions)

  odesys = mtk$modelingtoolkitize(prob)
  
  JuliaCall::julia_assign("odesys", odesys)
  JuliaCall::julia_assign("tspan", prob$tspan)
  new_prob <- JuliaCall::julia_eval("ODEProblem(complete(odesys, split=false), [], tspan; jac=true)")
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
#' \dontrun{ ## diffeq_setup() is time-consuming and requires Julia+DifferentialEquations.jl
#'
#' diffeqr::diffeq_setup()
#'
#' }
#'
#' @export
jitoptimize_sde <- function (de,prob){
  JuliaCall::julia_install_package_if_needed("ModelingToolkit")
  JuliaCall::julia_library("ModelingToolkit")
  functions <- JuliaCall::julia_eval("filter(isascii, replace.(string.(propertynames(ModelingToolkit)),\"!\"=>\"_bang\"))")

  # Can remove the de argument when breaking, but kept for backwards compat
  mtk <- julia_pkg_import("ModelingToolkit",functions)

  sdesys = mtk$modelingtoolkitize(prob)
  JuliaCall::julia_assign("sdesys", sdesys)
  JuliaCall::julia_assign("tspan", prob$tspan)
  new_prob <- JuliaCall::julia_eval("SDEProblem(complete(sdesys, split=false), [], tspan; jac=true)")
}

#' Setup DiffEqGPU
#'
#' This function initializes the DiffEqGPU package for GPU-parallelized ensembles.
#' The first time will be long since it includes precompilation.
#'
#' @param backend the backend for the GPU computation. Choices are "CUDA", "AMDGPU", "Metal", or "oneAPI"
#'
#' @examples
#'
#' \dontrun{ ## diffeq_setup() is time-consuming and requires Julia+DifferentialEquations.jl
#'
#' degpu <- diffeqr::diffeqgpu_setup(backend="CUDA")
#'
#' }
#'
#' @export
diffeqgpu_setup <- function (backend){
  use_julia_libunwind()
  JuliaCall::julia_install_package_if_needed("DiffEqGPU")
  JuliaCall::julia_library("DiffEqGPU")
  functions <- JuliaCall::julia_eval("filter(isascii, replace.(string.(propertynames(DiffEqGPU)),\"!\"=>\"_bang\"))")
  degpu <- julia_pkg_import("DiffEqGPU",functions)

  if (backend == "CUDA") {
    JuliaCall::julia_install_package_if_needed("CUDA")
    JuliaCall::julia_library("CUDA")
    backend <- julia_pkg_import("CUDA",c("CUDABackend"))
    degpu$CUDABackend <- backend$CUDABackend
  } else if (backend == "AMDGPU") {
    JuliaCall::julia_install_package_if_needed("AMDGPU")
    JuliaCall::julia_library("AMDGPU")
    backend <- julia_pkg_import("AMDGPU",c("AMDGPUBackend"))
    degpu$AMDGPUBackend <- backend$AMDGPUBackend
  } else if (backend == "Metal") {
    JuliaCall::julia_install_package_if_needed("Metal")
    JuliaCall::julia_library("Metal")
    backend <- julia_pkg_import("Metal",c("MetalBackend"))
    degpu$MetalBackend <- backend$MetalBackend
  } else if (backend == "oneAPI") {
    JuliaCall::julia_install_package_if_needed("oneAPI")
    JuliaCall::julia_library("oneAPI")
    backend <- julia_pkg_import("oneAPI",c("oneAPIBackend"))
    degpu$oneAPIBackend <- backend$oneAPIBackend
  } else {
    stop(paste("Illegal backend choice found. Allowed choices: CUDA, AMDGPU, Metal, and oneAPI. Chosen backend: ", backend))
  }
  degpu
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
