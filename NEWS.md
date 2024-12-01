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
