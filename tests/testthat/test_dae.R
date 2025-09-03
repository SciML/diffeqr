context("DAEs")

test_that('DAEs work',{

  skip_on_cran()
  skip_if(Sys.getenv("CI") != "", "Skip on CI - Julia installation too time-consuming")

  de <- diffeqr::diffeq_setup()
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
  expect_true(length(sol$t)<200)
})
