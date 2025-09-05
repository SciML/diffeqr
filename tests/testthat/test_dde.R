context("DDEs")

test_that('DDEs work',{

  skip_on_cran()

  de <- diffeqr::diffeq_setup()
  f <- JuliaCall::julia_eval("function f(du, u, h, p, t)
  du[1] = 1.1/(1 + sqrt(10)*(h(p, t-20)[1])^(5/4)) - 10*u[1]/(1 + 40*u[2])
  du[2] = 100*u[1]/(1 + 40*u[2]) - 2.43*u[2]
end")
  h <- JuliaCall::julia_eval("function h(p, t)
  [1.05767027/3, 1.030713491/3]
end")
  u0 <- c(1.05767027/3, 1.030713491/3)
  tspan <- c(0.0, 100.0)
  constant_lags <- c(20.0)
  JuliaCall::julia_assign("u0", u0)
  JuliaCall::julia_assign("tspan", tspan)
  JuliaCall::julia_assign("constant_lags", constant_lags)
  prob <- JuliaCall::julia_eval("DDEProblem(f, u0, h, tspan, constant_lags = constant_lags)")
  sol <- de$solve(prob,de$MethodOfSteps(de$Tsit5()))
  udf <- as.data.frame(t(sapply(sol$u,identity)))
  expect_true(length(sol$t)<200)
  #plotly::plot_ly(udf, x = sol$t, y = ~V1, type = 'scatter', mode = 'lines') %>% plotly::add_trace(y = ~V2)
})
