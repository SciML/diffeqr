context("SDEs")

f <- function(u,p,t) {
  return(1.01*u)
}
g <- function(u,p,t) {
  return(0.87*u)
}
u0 = 1/2
tspan <- list(0.0,1.0)
sol = sde.solve(f,g,u0,tspan)
test_that('1D works',{
  expect_true(length(sol$t) > 10)
})
#plot(sol$t,sol$u,"l")

f <- function(u,p,t) {
  du1 = p[1]*(u[2]-u[1])
  du2 = u[1]*(p[2]-u[3]) - u[2]
  du3 = u[1]*u[2] - p[3]*u[3]
  return(c(du1,du2,du3))
}
g <- function(u,p,t) {
  return(c(0.3*u[1],0.3*u[2],0.3*u[3]))
}
u0 = c(1.0,0.0,0.0)
tspan <- list(0.0,1.0)
p = c(10.0,28.0,8/3)
sol = sde.solve(f,g,u0,tspan,p=p,saveat=0.05)
udf = as.data.frame(sol$u)
#plotly::plot_ly(udf, x = ~V1, y = ~V2, z = ~V3, type = 'scatter3d', mode = 'lines')

f <- julia_eval("
function f(du,u,p,t)
  du[1] = 10.0*(u[2]-u[1])
  du[2] = u[1]*(28.0-u[3]) - u[2]
  du[3] = u[1]*u[2] - (8/3)*u[3]
end")

g <- julia_eval("
function g(du,u,p,t)
  du[1] = 0.3*u[1]
  du[2] = 0.3*u[2]
  du[3] = 0.3*u[3]
end")
tspan <- list(0.0,100.0)
sol = sde.solve(f,g,u0,tspan,fname="f",gname="g",p=p,saveat=0.05)
udf = as.data.frame(sol$u)
#plotly::plot_ly(udf, x = ~V1, y = ~V2, z = ~V3, type = 'scatter3d', mode = 'lines')
