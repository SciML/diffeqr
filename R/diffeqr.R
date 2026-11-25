#' Setup diffeqr
#'
#' This function initializes Julia and the DifferentialEquations.jl package.
#' The first time will be long since it includes precompilation.
#' Additionally, this will install Julia and the required packages
#' if they are missing. With DifferentialEquations.jl v8 and later, which no
#' longer bundles the SDE and DDE solvers, StochasticDiffEq.jl and
#' DelayDiffEq.jl are installed and loaded as well.
#'
#' On Linux, the \code{libunwind} that Julia ships with is loaded into R before
#' Julia starts, installing Julia first if none is found. A Julia already started
#' in the session, e.g. by \code{JuliaCall::julia_setup()}, keeps the unwinder it
#' started with.
#'
#' With DifferentialEquations.jl v8 the returned object exposes fewer names than
#' with v7: for example \code{Rodas5}, \code{RK4}, \code{QNDF},
#' \code{JumpProblem} and \code{PresetTimeCallback} are not available through
#' it. To use one, install and load the package that owns it (for \code{Rodas5},
#' \code{JuliaCall::julia_install_package_if_needed("OrdinaryDiffEqRosenbrock")},
#' \code{JuliaCall::julia_library("OrdinaryDiffEqRosenbrock")} and then
#' \code{JuliaCall::julia_eval("Rodas5()")}), or stay on
#' DifferentialEquations.jl v7. On Linux, install and load such a package with
#' \code{LD_LIBRARY_PATH} temporarily set to \code{""}, as in \code{NEWS.md};
#' otherwise the Julia processes this starts can segfault.
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
  restore_library_path <- filter_library_path()
  on.exit(restore_library_path())
  julia <- setup_julia(...)
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

# R puts the system library directory in LD_LIBRARY_PATH, which the loader searches
# before Julia's own, so Julia processes get the system libunwind and segfault.
# The processes started while Julia is set up or installs packages inherit this
# filtered path. The returned function puts the caller's value back.
filter_library_path <- function(){
  old <- Sys.getenv("LD_LIBRARY_PATH", unset = NA)
  if (Sys.info()[["sysname"]] != "Linux" || is.na(old)) return(function() invisible())
  dirs <- strsplit(old, ":")[[1]]
  system_dir <- grepl("^(/usr)?/lib(64)?(/[^/]*-linux-gnu[^/]*)?/?$", dirs)
  Sys.setenv(LD_LIBRARY_PATH = paste(dirs[!system_dir], collapse = ":"))
  function() Sys.setenv(LD_LIBRARY_PATH = old)
}

# JuliaCall::julia_setup() installs a missing Julia and embeds it in one call, so on Linux
# the install happens here, leaving room to load libunwind before Julia starts.
# `...` is bound to julia_setup()'s formals as JuliaCall would bind it.
setup_julia <- function(...){
  args <- as.list(match.call(JuliaCall::julia_setup,
                             as.call(c(quote(julia_setup), list(...), installJulia = TRUE))))[-1]
  if (Sys.info()[["sysname"]] == "Linux") {
    bindir <- julia_bindir(args$JULIA_HOME)
    # a Julia that is already running needs no install, even if it is not discoverable
    if (is.null(bindir) && !any(startsWith(names(getLoadedDLLs()), "libjulia"))) {
      if (!is.null(args$JULIA_HOME))
        stop("Julia is not found at JULIA_HOME = ", args$JULIA_HOME, ".")
      JuliaCall::install_julia(if ("version" %in% names(args)) args$version else formals(JuliaCall::julia_setup)$version)
      bindir <- julia_bindir()
      if (is.null(bindir)) stop("Julia is not found and automatic installation failed.")
    }
    if (!is.null(bindir)) args$JULIA_HOME <- bindir
    preload_julia_libunwind(bindir)
  }
  do.call(JuliaCall::julia_setup, args, quote = TRUE)
}

# The embedded Julia uses R's loader, which read LD_LIBRARY_PATH at R's startup,
# so Julia's libunwind has to be loaded first.
preload_julia_libunwind <- function(bindir){
  if (Sys.info()[["sysname"]] != "Linux" || is.null(bindir)) return(invisible())
  libunwind <- file.path(dirname(bindir), "lib", "julia", "libunwind.so.8")
  if (file.exists(libunwind)) dyn.load(libunwind)
  invisible()
}

# Sys.BINDIR of the Julia that JuliaCall::julia_setup(JULIA_HOME) selects on Linux,
# or NULL if there is none yet. The lookup order is JuliaCall's. A launcher may print
# other lines, so the path is read from the one tagged line and must hold a julia.
julia_bindir <- function(JULIA_HOME = NULL){
  if (is.null(JULIA_HOME)) JULIA_HOME <- getOption("JULIA_HOME")
  if (is.null(JULIA_HOME) && nzchar(Sys.getenv("JULIA_HOME"))) JULIA_HOME <- Sys.getenv("JULIA_HOME")
  if (is.null(JULIA_HOME)) {
    depot <- Sys.getenv("JULIA_DEPOT_PATH", unset = file.path(Sys.getenv("HOME"), ".julia"))
    prefs <- file.path(depot, "prefs", "JuliaCall")
    if (file.exists(prefs)) JULIA_HOME <- readChar(prefs, 256)
  }
  tag <- "DIFFEQR_JULIA_BINDIR="
  tryCatch({
    julia <- if (!is.null(JULIA_HOME)) {
      file.path(JULIA_HOME, "julia")
    } else if (nzchar(Sys.which("julia"))) {
      Sys.which("julia")
    } else {
      system2("bash", "-l -c 'which julia'", stdout = TRUE)[1]
    }
    out <- if (isTRUE(file.exists(julia)))
      system2(julia, c("--startup-file=no", "-e", shQuote(sprintf('println("%s", Sys.BINDIR)', tag))),
              stdout = TRUE)
    bindir <- substring(out[startsWith(out, tag)], nchar(tag) + 1)
    if (length(bindir) == 1 && file.access(file.path(bindir, "julia"), 1) == 0) bindir else NULL
  }, warning = function(w) NULL, error = function(e) NULL)
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
  restore_library_path <- filter_library_path()
  on.exit(restore_library_path())
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
  restore_library_path <- filter_library_path()
  on.exit(restore_library_path())
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
  restore_library_path <- filter_library_path()
  on.exit(restore_library_path())
  preload_julia_libunwind(julia_bindir())
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
