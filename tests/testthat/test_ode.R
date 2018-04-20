context("ODEs")

f <- function(u,p,t) {
  return(1.01*u)
}
u0 = 1/2
tspan <- list(0.0,1.0)
sol = diffeqr::ode.solve(f,u0,tspan)
test_that('1D works',{
  expect_true(length(sol$t) < 10)
  expect_true(sol$u[length(sol$u)] - exp(1)/2 < 5e-2)
})


saveat=1:10/10
sol2 = diffeqr::ode.solve(f,u0,tspan,saveat=saveat)

test_that('saveat works',{
  expect_true(length(sol2$t) == 11)
  expect_equal(sol$u[length(sol$u)],sol2$u[length(sol2$u)])
})

sol3 = diffeqr::ode.solve(f,u0,tspan,alg="Vern9()")
test_that('alg works',{
  expect_true(length(sol3$t) < 5)
})

f <- function(u,p,t) {
  du1 = p[1]*(u[2]-u[1])
  du2 = u[1]*(p[2]-u[3]) - u[2]
  du3 = u[1]*u[2] - p[3]*u[3]
  return(c(du1,du2,du3))
}
u0 = c(1.0,0.0,0.0)
tspan <- list(0.0,100.0)
p = c(10.0,28.0,8/3)
sol = diffeqr::ode.solve(f,u0,tspan,p=p)
udf = as.data.frame(sol$u)
#matplot(sol$t,udf,"l",col=1:3)
#plotly::plot_ly(udf, x = ~V1, y = ~V2, z = ~V3, type = 'scatter3d', mode = 'lines')


f <- JuliaCall::julia_eval("
function f(du,u,p,t)
  du[1] = 10.0*(u[2]-u[1])
  du[2] = u[1]*(28.0-u[3]) - u[2]
  du[3] = u[1]*u[2] - (8/3)*u[3]
end")
u0 = c(1.0,0.0,0.0)
tspan <- list(0.0,100.0)
sol = diffeqr::ode.solve('f',u0,tspan)
udf = as.data.frame(sol$u)
#matplot(sol$t,udf,"l",col=1:3)
#plotly::plot_ly(udf, x = ~V1, y = ~V2, z = ~V3, type = 'scatter3d', mode = 'lines')
