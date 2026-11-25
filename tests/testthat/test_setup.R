context("Setup")

test_that('setup exposes the SDE and DDE solver interfaces',{
  skip_on_cran()
  de <- diffeqr::diffeq_setup()
  expect_true(is.function(de$MethodOfSteps))

  prob <- de$SDEProblem(function(u,p,t) 1.25, function(u,p,t) 0, 2, c(0,1))
  sol <- de$solve(prob, de$SOSRI(), saveat=0.25)
  t <- seq(0,1,0.25)
  expect_equal(sol$t, t)
  expect_equal(sol$u, 2 + 1.25*t)
})

test_that('the documented ensemble callbacks work',{
  skip_on_cran()
  de <- diffeqr::diffeq_setup()
  prob <- de$ODEProblem(function(u,p,t) 1, 0, c(0,1))
  endpoints <- function(prob_func){
    ensembleprob <- de$EnsembleProblem(prob, prob_func = prob_func, safetycopy=FALSE)
    sol <- de$solve(ensembleprob, de$Tsit5(), de$EnsembleSerial(), trajectories=3L)
    vapply(sol$u, function(s) s$u[[length(s$u)]], numeric(1))
  }

  # du/dt = 1 on (0,1), so each endpoint is u0 + 1
  calls <- 0
  index_free <- function(prob,...){
    calls <<- calls + 1
    de$remake(prob, u0=10*calls)
  }
  expect_equal(endpoints(index_free), c(11,21,31))

  by_index <- if (JuliaCall::julia_eval('pkgversion(SciMLBase) >= v"3"')) {
    function(prob,ctx) de$remake(prob, u0=2*ctx$sim_id)
  } else {
    function(prob,i,rep) de$remake(prob, u0=2*i)
  }
  expect_equal(endpoints(by_index), c(3,5,7))
})

test_that("setup leaves the caller's LD_LIBRARY_PATH alone",{
  skip_on_cran()
  skip_if_not(Sys.info()[["sysname"]] == "Linux")
  system_dir <- paste0("/usr/lib/", R.version$arch, "-linux-gnu")
  script <- tempfile(fileext = ".R")
  result <- tempfile()
  ld_log <- tempfile()
  dir.create(ld_log)
  # Julia can only be set up once per process, so this runs in a fresh one
  writeLines(deparse(quote({
    args <- commandArgs(TRUE)
    before <- Sys.getenv("LD_LIBRARY_PATH")
    # makes the loader of every child process log the libraries it loads
    Sys.setenv(LD_DEBUG = "libs", LD_DEBUG_OUTPUT = file.path(args[1], "ld"))
    error <- tryCatch({
      diffeqr::diffeq_setup(pkg_check = FALSE, sysimage_path = tempfile())
      ""
    }, error = conditionMessage)
    after_error <- Sys.getenv("LD_LIBRARY_PATH")
    de <- diffeqr::diffeq_setup()
    after_success <- Sys.getenv("LD_LIBRARY_PATH")
    Sys.unsetenv(c("LD_DEBUG", "LD_DEBUG_OUTPUT"))
    saveRDS(list(before = before, error = error, after_error = after_error,
                 after_success = after_success,
                 bindir = JuliaCall::julia_eval("Sys.BINDIR")), args[2])
  })), script)
  status <- system2(file.path(R.home("bin"), "Rscript"),
                    shQuote(c(script, ld_log, result)),
                    env = c(paste0("LD_LIBRARY_PATH=", system_dir), "R_TESTS=",
                            paste0("R_LIBS=", paste(.libPaths(), collapse = ":"))))
  expect_equal(status, 0)
  res <- readRDS(result)

  expect_true(system_dir %in% strsplit(res$before, ":")[[1]])
  expect_match(res$error, "sysimage")
  expect_identical(res$after_error, res$before)
  expect_identical(res$after_success, res$before)

  unwinders <- unlist(lapply(list.files(ld_log, full.names = TRUE), function(f){
    log <- readLines(f, warn = FALSE)
    if (!any(grepl("initialize program: (.*/)?julia$", log))) return(NULL)
    sub(".*calling init: ", "", grep("calling init: .*/libunwind\\.so\\.8$", log, value = TRUE))
  }))
  bundled <- file.path(dirname(res$bindir), "lib", "julia", "libunwind.so.8")
  expect_gt(length(unwinders), 0)
  expect_equal(unique(normalizePath(unwinders)), normalizePath(bundled))
})
