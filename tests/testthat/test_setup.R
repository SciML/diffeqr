context("Setup")

test_that('setup exposes the SDE and DDE solver interfaces',{
  skip_on_cran()
  de <- diffeqr::diffeq_setup()
  expect_true(is.function(de$SOSRI))
  expect_true(is.function(de$MethodOfSteps))
})
