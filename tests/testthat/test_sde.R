context("SDEs")

test_that('1D works',{

  skip_on_cran()

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

  expect_true(length(sol$t) > 10)
})

test_that('diagonal noise works',{

  skip_on_cran()
  de <- diffeqr::diffeq_setup()
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
  expect_true(length(sol$t)>10)
})

test_that('diagonal noise JIT works',{

  skip_on_cran()
  de <- diffeqr::diffeq_setup()
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
  fastprob <- diffeqr::jitoptimize_sde(de,prob)
  sol <- de$solve(fastprob,saveat=0.005)
  udf <- as.data.frame(t(sapply(sol$u,identity)))
  expect_true(length(sol$t)>10)
})

test_that('non-diagonal noise works',{

  skip_on_cran()
  de <- diffeqr::diffeq_setup()
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
  prob <- JuliaCall::julia_eval("SDEProblem(f, g, u0, tspan, nothing, noise_rate_prototype=noise_rate_prototype)")
  sol <- de$solve(prob)
  udf <- as.data.frame(t(sapply(sol$u,identity)))
  expect_true(length(sol$t)>10)
})
