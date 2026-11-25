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

# Julia can only be set up once per process, so each case runs in a fresh R process
# whose startup search path has the system library directory.
run_fresh_setup <- function(case, julia_root, setup_args = list(), check_unwinder = TRUE){
  script <- tempfile(fileext = ".R")
  result <- tempfile()
  args_file <- tempfile()
  saveRDS(setup_args, args_file)
  ld_log <- tempfile()
  dir.create(ld_log)
  writeLines(deparse(quote({
    args <- commandArgs(TRUE)
    case <- args[1]
    julia_root <- args[2]
    before <- Sys.getenv("LD_LIBRARY_PATH")
    # makes the loader of every child process log the libraries it loads
    Sys.setenv(LD_DEBUG = "libs", LD_DEBUG_OUTPUT = file.path(args[3], "ld"))
    setup_args <- readRDS(args[5])
    installs <- 0
    installed <- NULL
    version_forced <- 0L
    julia_home_forced <- 0L
    warm_callable <- FALSE
    error <- ""
    installer <- function(version = "latest", ...){
      installs <<- installs + 1
      installed <<- version
    }
    if (case == "installer") {
      # No Julia is at JULIA_HOME until JuliaCall's installer puts the runtime there
      staging <- tempfile()
      options(JULIA_HOME = file.path(staging, "bin"))
      installer <- function(version = "latest", ...){
        installs <<- installs + 1
        installed <<- version
        file.symlink(julia_root, staging)
      }
    }
    utils::assignInNamespace("install_julia", installer, "JuliaCall")
    if (case == "blank_verbose") {
      error <- tryCatch({
        diffeqr::diffeq_setup(FALSE, file.path(julia_root, "bin"), ,
                              sysimage_path = setup_args$sysimage_path)
        ""
      }, error = conditionMessage)
    } else if (case == "lazy_version") {
      error <- tryCatch({
        diffeqr::diffeq_setup(FALSE,
                              version = {
                                version_forced <<- version_forced + 1L
                                stop("unused version evaluated")
                              },
                              sysimage_path = setup_args$sysimage_path)
        ""
      }, error = conditionMessage)
    } else if (case == "missing_image") {
      error <- tryCatch({
        do.call(diffeqr::diffeq_setup, setup_args)
        ""
      }, error = conditionMessage)
    } else {
      de <- do.call(diffeqr::diffeq_setup, setup_args)
      if (case == "installed") {
        # a running Julia that can no longer be discovered must not be installed again
        options(JULIA_HOME = tempfile())
        de <- diffeqr::diffeq_setup(FALSE, JULIA_HOME = {
          julia_home_forced <<- julia_home_forced + 1L
          stop("unused JULIA_HOME evaluated")
        })
        warm_callable <- is.function(de$ODEProblem)
      }
    }
    after <- Sys.getenv("LD_LIBRARY_PATH")
    Sys.unsetenv(c("LD_DEBUG", "LD_DEBUG_OUTPUT"))
    maps <- grep("/libunwind\\.so", readLines("/proc/self/maps"), value = TRUE)
    saveRDS(list(before = before, after = after, error = error, installs = installs,
                 installed = installed, version_forced = version_forced,
                 julia_home_forced = julia_home_forced,
                 warm_callable = warm_callable,
                 embedded = unique(sub("^[^/]*", "", maps))), args[4])
  })), script)
  system_dir <- paste0("/usr/lib/", R.version$arch, "-linux-gnu")
  status <- system2(file.path(R.home("bin"), "Rscript"),
                    shQuote(c(script, case, julia_root, ld_log, result, args_file)),
                    env = c(paste0("LD_LIBRARY_PATH=", system_dir), "R_TESTS=",
                            paste0("R_LIBS=", paste(.libPaths(), collapse = ":"))))
  expect_equal(status, 0)
  res <- readRDS(result)

  expect_true(system_dir %in% strsplit(res$before, ":")[[1]])
  expect_identical(res$after, res$before)

  unwinders <- unlist(lapply(list.files(ld_log, full.names = TRUE), function(f){
    log <- readLines(f, warn = FALSE)
    if (!any(grepl("initialize program: (.*/)?julia$", log))) return(NULL)
    sub(".*calling init: ", "", grep("calling init: .*/libunwind\\.so\\.8$", log, value = TRUE))
  }))
  bundled <- file.path(julia_root, "lib", "julia", "libunwind.so.8")
  if (check_unwinder) {
    expect_gt(length(unwinders), 0)
    expect_equal(unique(normalizePath(unwinders)), normalizePath(bundled))
    res$bundled <- normalizePath(bundled)
  }
  res
}

test_that("a failed setup leaves the caller's LD_LIBRARY_PATH alone",{
  skip_on_cran()
  skip_if_not(Sys.info()[["sysname"]] == "Linux")
  diffeqr::diffeq_setup()
  res <- run_fresh_setup("missing_image", dirname(JuliaCall::julia_eval("Sys.BINDIR")),
                         list(pkg_check = FALSE, sysimage_path = tempfile()))
  expect_match(res$error, "sysimage")
  expect_equal(res$installs, 0)
})

test_that("setup forwards version and JULIA_HOME the way julia_setup() binds them",{
  skip_on_cran()
  skip_if_not(Sys.info()[["sysname"]] == "Linux")
  diffeqr::diffeq_setup()
  julia_root <- dirname(JuliaCall::julia_eval("Sys.BINDIR"))
  image <- tempfile()
  # a launcher that prints a line before running Julia
  launcher <- tempfile()
  dir.create(launcher)
  writeLines(c("#!/bin/sh", "echo 'Julia launcher'",
               sprintf('exec "%s" "$@"', file.path(julia_root, "bin", "julia"))),
             file.path(launcher, "julia"))
  Sys.chmod(file.path(launcher, "julia"), "755")
  forms <- list(
    named = list(pkg_check = FALSE, sysimage_path = image, version = "1.10"),
    abbreviated = list(pkg_check = FALSE, sysimage_path = image, vers = "1.10"),
    positional = list(FALSE, NULL, TRUE, TRUE, FALSE, TRUE, FALSE, image, "1.10"),
    launcher = list(pkg_check = FALSE, JULIA_HOME = launcher, sysimage_path = image)
  )
  for (form in names(forms)) {
    res <- run_fresh_setup("missing_image", julia_root, forms[[form]])
    expect_match(res$error, paste("sysimage at path:", image), fixed = TRUE, info = form)
    expect_equal(res$installs, 0, info = form)
  }
})

test_that("setup preserves blank positional arguments and lazy version expressions",{
  skip_on_cran()
  skip_if_not(Sys.info()[["sysname"]] == "Linux")
  diffeqr::diffeq_setup()
  julia_root <- dirname(JuliaCall::julia_eval("Sys.BINDIR"))
  image <- tempfile()

  blank <- run_fresh_setup("blank_verbose", julia_root,
                           list(sysimage_path = image), check_unwinder = FALSE)
  expect_match(blank$error, paste("sysimage at path:", image), fixed = TRUE)
  expect_equal(blank$installs, 0)

  lazy <- run_fresh_setup("lazy_version", julia_root,
                          list(sysimage_path = image), check_unwinder = FALSE)
  expect_match(lazy$error, paste("sysimage at path:", image), fixed = TRUE)
  expect_equal(lazy$version_forced, 0L)
  expect_equal(lazy$installs, 0)
})

test_that("an explicit JULIA_HOME without Julia is reported instead of installing Julia",{
  skip_on_cran()
  skip_if_not(Sys.info()[["sysname"]] == "Linux")
  result <- tempfile()
  script <- tempfile(fileext = ".R")
  writeLines(deparse(quote({
    installs <- 0
    utils::assignInNamespace("install_julia", function(...) installs <<- installs + 1, "JuliaCall")
    error <- tryCatch({
      diffeqr::diffeq_setup(JULIA_HOME = tempfile())
      ""
    }, error = conditionMessage)
    saveRDS(list(error = error, installs = installs), commandArgs(TRUE)[1])
  })), script)
  status <- system2(file.path(R.home("bin"), "Rscript"), shQuote(c(script, result)),
                    env = c("R_TESTS=", paste0("R_LIBS=", paste(.libPaths(), collapse = ":"))))
  expect_equal(status, 0)
  res <- readRDS(result)
  expect_match(res$error, "Julia is not found at JULIA_HOME")
  expect_equal(res$installs, 0)
})

test_that("setup embeds Julia's libunwind and leaves the caller's LD_LIBRARY_PATH alone",{
  skip_on_cran()
  skip_if_not(Sys.info()[["sysname"]] == "Linux")
  diffeqr::diffeq_setup()
  julia_root <- dirname(JuliaCall::julia_eval("Sys.BINDIR"))
  for (case in c("installed", "installer")) {
    res <- run_fresh_setup(case, julia_root)
    expect_equal(res$installs, if (case == "installer") 1 else 0, info = case)
    expect_equal(res$embedded, res$bundled, info = case)
    if (case == "installed") {
      expect_true(res$warm_callable)
      expect_equal(res$julia_home_forced, 0L)
    }
  }
  res <- run_fresh_setup("installer", julia_root, list(vers = "1.10"))
  expect_equal(res$installs, 1)
  expect_identical(res$installed, "1.10")
  expect_equal(res$embedded, res$bundled)
})
