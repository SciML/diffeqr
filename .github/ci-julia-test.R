JuliaCall::julia_setup(installJulia = FALSE)

expected <- normalizePath(Sys.getenv("JULIA_HOME"), winslash = "/", mustWork = TRUE)
actual <- normalizePath(JuliaCall::julia_eval("Sys.BINDIR"), winslash = "/", mustWork = TRUE)
stopifnot(identical(actual, expected))
cat("JuliaCall uses the provisioned Julia", JuliaCall::julia_eval("string(VERSION)"), "\n")
