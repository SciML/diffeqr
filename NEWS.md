## Unreleased

Support DifferentialEquations.jl v8, which no longer bundles the SDE and DDE solvers or exposes
`DiffEqBase`. `diffeq_setup()` now also loads StochasticDiffEq.jl and DelayDiffEq.jl in that case.

With DifferentialEquations.jl v8, `de` exposes fewer names than with v7. For example `Rodas5`,
`RK4`, `Euler`, `QNDF`, `DFBDF`, `RadauIIA5`, `JumpProblem`, `MIRK4`, `Shooting`,
`KLUFactorization`, `PresetTimeCallback` and `TerminateSteadyState` are not on `de`. To use one,
install and load the package that owns it, or stay on DifferentialEquations.jl v7. `IDA` is not on
`de` with either version; it needs Sundials. For example:

```R
local({
  # Linux: keeps the Julia processes that install and precompile from segfaulting
  old <- Sys.getenv("LD_LIBRARY_PATH", unset = NA_character_)
  if (Sys.info()[["sysname"]] == "Linux") {
    Sys.setenv(LD_LIBRARY_PATH = "")
    on.exit(if (is.na(old)) Sys.unsetenv("LD_LIBRARY_PATH") else Sys.setenv(LD_LIBRARY_PATH = old))
  }
  JuliaCall::julia_install_package_if_needed("OrdinaryDiffEqRosenbrock")
  JuliaCall::julia_library("OrdinaryDiffEqRosenbrock")
})
sol <- de$solve(prob, JuliaCall::julia_eval("Rodas5()"))
```

The owning packages are OrdinaryDiffEqRosenbrock (`Rodas5`), OrdinaryDiffEqLowOrderRK (`RK4`,
`Euler`), OrdinaryDiffEqBDF (`QNDF`, `DFBDF`), OrdinaryDiffEqFIRK (`RadauIIA5`), Sundials (`IDA`),
JumpProcesses (`JumpProblem`), BoundaryValueDiffEq (`MIRK4`, `Shooting`), LinearSolve
(`KLUFactorization`) and DiffEqCallbacks (`PresetTimeCallback`, `TerminateSteadyState`).

Ensembles with SciMLBase.jl 2 and 3:

* With SciMLBase 3, pass `trajectories` as an integer (`10000L`, not `10000`).
* With SciMLBase 3, `prob_func` receives `(prob, ctx)` and the trajectory index is `ctx$sim_id`.
  With SciMLBase 2, it receives `(prob, i, rep)`, including when used with
  DifferentialEquations 8.0. `ctx` is an `EnsembleContext`; ``ctx$`repeat` `` is the repeat count.
  A `function (prob, ...)` that only uses `prob` works with both.

Fix segfaults on Linux systems with `libunwind` installed: `diffeq_setup()` and `diffeqgpu_setup()`
load the `libunwind` that Julia ships with into R. `diffeq_setup()` does so before Julia starts,
installing Julia first if none is found; `diffeqgpu_setup()` does so when Julia is already
installed. A Julia already started in the session keeps the unwinder it started with. While these
functions and `jitoptimize_ode()` / `jitoptimize_sde()` run, the default system library directories
are left out of `LD_LIBRARY_PATH`, so that Julia processes started by them use Julia's `libunwind`
too. The value of `LD_LIBRARY_PATH` is restored when they return, including on error. Julia
processes started at other times, such as package precompilation triggered by your own `JuliaCall`
calls, inherit your `LD_LIBRARY_PATH` and can segfault; the example above shows how to avoid this.

## Release v2.1.0

Better support for ModelingToolkit JIT tracing on SDEs.

## Release v2.0.1

Updated to support ModelingToolkit v9 from the Julia side with the JIT compilation.

## Release v2.0.0

Support new DiffEqGPU syntax. This requires passing a backend. Supports NVIDIA CUDA, Intel OneAPI,
AMD GPUs, and Apple Metal GPUs. Also much faster GPU compilation and runtime performance.

## Release v1.1.2

Bugfixes for newer Julia versions.

## Release v1.1.1

This package now ensures that the tests are not run on build so that 
installation of Julia will not occur unless the user specifically asks for it
via diffeqr::diffeq_setup()

## Release v1.0.0

Full recreation of the package. This provides a new simplified interface over
DifferentialEquations.jl that matches the Julia interface almost 1-1.

## Release v0.1.1

This is a quick patch to fix the vignettes of the v0.1.0 release.

## Release v0.1.0

This is the initial release of the package. It provides a simplified interface over DifferentialEquations.jl. Currently it's interfaced via 5 functions:

- diffeq_setup
- ode.solve
- sde.solve
- dae.solve
- dde.solve

The return is a list with sol$u and sol$t. In future updates this will be backwards compatibly updated to be a full solution object with the interpolation.
