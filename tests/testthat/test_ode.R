context("ODEs")

test_that('1D works',{
  skip_on_cran()
  de <- diffeqr::diffeq_setup()
  f <- function(u,p,t) {
    return(1.01*u)
  }
  u0 <- 1/2
  tspan <- c(0., 1.)
  prob = de$ODEProblem(f, u0, tspan)
  sol = de$solve(prob)
  sol$.(0.2)
  expect_true(length(sol$t)<200)
})

test_that('ODE system works',{

  skip_on_cran()
  de <- diffeqr::diffeq_setup()
  f <- function(u,p,t) {
    du1 = p[1]*(u[2]-u[1])
    du2 = u[1]*(p[2]-u[3]) - u[2]
    du3 = u[1]*u[2] - p[3]*u[3]
    return(c(du1,du2,du3))
  }
  u0 <- c(1.0,0.0,0.0)
  tspan <- list(0.0,100.0)
  p <- c(10.0,28.0,8/3)
  prob <- de$ODEProblem(f, u0, tspan, p)
  sol <- de$solve(prob)
  mat <- sapply(sol$u,identity)
  udf <- as.data.frame(t(mat))

  abstol <- 1e-8
  reltol <- 1e-8
  saveat <- 0:10000/100
  sol <- de$solve(prob,abstol=abstol,reltol=reltol,saveat=saveat)
  udf <- as.data.frame(t(sapply(sol$u,identity)))
  expect_true(length(sol$t)>200)
})

test_that('ODE JIT works',{

  skip_on_cran()
  de <- diffeqr::diffeq_setup()
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
  expect_true(length(sol$t)>200)
})
