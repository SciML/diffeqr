workflow <- yaml::read_yaml(".github/workflows/cran-release.yml")
steps <- unlist(lapply(workflow$jobs, `[[`, "steps"), recursive = FALSE)
build_steps <- Filter(function(step) {
  isTRUE(step$name %in% c("Validate release metadata", "Build source tarball"))
}, steps)
stopifnot(length(build_steps) > 0L)

run_build <- function(event, tag = "") {
  work <- tempfile()
  log <- tempfile()
  dir.create(work)
  old_wd <- setwd(work)
  on.exit({
    setwd(old_wd)
    unlink(work, recursive = TRUE)
    unlink(log)
  })
  writeLines(c(
    "Package: releasefixture",
    "Version: 2.1.0",
    "Title: Release Workflow Test Fixture",
    "Description: A source tarball for testing the workflow.",
    'Authors@R: person("Test", "Author", email = "test@example.org", role = c("aut", "cre"))',
    "License: MIT"
  ), "DESCRIPTION")
  writeLines(character(), "NAMESPACE")
  env <- c(
    "LC_ALL=C.UTF-8", "GITHUB_OUTPUT=output",
    paste0("GITHUB_EVENT_NAME=", event), paste0("RELEASE_TAG=", tag)
  )
  output <- character()
  for (step in build_steps) {
    code <- shQuote(step$run)
    if (identical(step$shell, "Rscript --vanilla {0}")) {
      status <- system2("Rscript", c("--vanilla", "-e", code),
                        env = env, stdout = log, stderr = log)
    } else {
      tarball <- sub("^tarball=", "", output)
      status <- system2("bash", c("-e", "-o", "pipefail", "-c", code),
                        env = c(env, paste0("TARBALL=", tarball)),
                        stdout = log, stderr = log)
    }
    messages <- readLines(log)
    writeLines(messages)
    if (file.exists("output")) output <- readLines("output")
    if (status != 0L) break
  }
  list(status = status, output = output, messages = messages,
       tarball = file.exists("releasefixture_2.1.0.tar.gz"))
}

dispatch <- run_build("workflow_dispatch")
wrong_tag <- run_build("release", "v2.0.0")
cat("UTF-8 dispatch:", dispatch$status, dispatch$output, "tarball exists:", dispatch$tarball, "\n")
cat("Mismatched tag:", wrong_tag$status, "tarball exists:", wrong_tag$tarball, "\n")
stopifnot(
  dispatch$status == 0L,
  identical(dispatch$output, "tarball=releasefixture_2.1.0.tar.gz"),
  dispatch$tarball,
  wrong_tag$status != 0L,
  any(grepl("must match DESCRIPTION version v2.1.0", wrong_tag$messages, fixed = TRUE)),
  !wrong_tag$tarball,
  length(wrong_tag$output) == 0L
)
