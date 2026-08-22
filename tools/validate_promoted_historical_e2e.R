#!/usr/bin/env Rscript

# Source-tree gate for the historical methods that are currently promoted.
# It is intentionally split into independent runs: the broad Synopsis gate
# excludes the costly Gaussian LASSO focus, while Count has its own durable
# artifact contract.  Each run exercises K=2, K=3 and K=5 and verifies a
# plausible result plus Rock-backed replay.

args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) {
  cat(paste(
    "Usage: DSVERT_SERVER_SOURCE=/path/to/dsVert",
    "Rscript tools/validate_promoted_historical_e2e.R [--dry-run]\n",
    "Runs the promoted historical Synopsis, Gaussian LASSO and Count",
    "source-tree E2E gates at K=2,3,5.\n"
  ))
  quit(save = "no", status = 0L)
}

if (length(setdiff(args, "--dry-run"))) {
  stop("Unknown argument; use --help", call. = FALSE)
}
if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("validate_promoted_historical_e2e.R requires devtools", call. = FALSE)
}

.historical_e2e_script <- function() {
  value <- commandArgs(trailingOnly = FALSE)
  value <- value[startsWith(value, "--file=")]
  if (length(value) != 1L) {
    stop("Run this file with Rscript", call. = FALSE)
  }
  normalizePath(sub("^--file=", "", value), mustWork = TRUE)
}

.historical_e2e_client_root <- function() {
  root <- normalizePath(
    file.path(dirname(.historical_e2e_script()), ".."), mustWork = TRUE)
  if (!file.exists(file.path(root, "DESCRIPTION")) ||
      !dir.exists(file.path(root, "tests", "testthat"))) {
    stop("Cannot locate the dsVertClient source tree", call. = FALSE)
  }
  root
}

.historical_e2e_server_root <- function(client_root) {
  configured <- Sys.getenv("DSVERT_SERVER_SOURCE", unset = "")
  candidate <- if (nzchar(configured)) configured else
    file.path(dirname(client_root), "dsVert")
  root <- normalizePath(candidate, mustWork = FALSE)
  if (!dir.exists(root) || !file.exists(file.path(root, "DESCRIPTION"))) {
    stop(
      "Set DSVERT_SERVER_SOURCE to the sibling dsVert source tree",
      call. = FALSE)
  }
  root
}

client_root <- .historical_e2e_client_root()
server_root <- .historical_e2e_server_root(client_root)
gates <- list(
  list(name = "Synopsis", filter = "dp-synopsis-describe-rock-e2e",
       family = ""),
  list(name = "Gaussian LASSO", filter = "dp-synopsis-describe-rock-e2e",
       family = "gaussian_lasso_focal"),
  list(name = "Count", filter = "dp-count-synopsis-rock-e2e", family = "")
)

if (identical(args, "--dry-run")) {
  cat("dsVertClient: ", client_root, "\n", sep = "")
  cat("dsVert: ", server_root, "\n", sep = "")
  for (gate in gates) {
    cat("- ", gate$name, " (", gate$filter, ")\n", sep = "")
  }
  quit(save = "no", status = 0L)
}

Sys.setenv(
  DSVERT_SERVER_SOURCE = server_root,
  DSVERT_TEST_SYNOPSIS_E2E_K = "2,3,5"
)
for (gate in gates) {
  if (nzchar(gate$family)) {
    Sys.setenv(DSVERT_TEST_SYNOPSIS_E2E_FAMILY = gate$family)
  } else {
    Sys.unsetenv("DSVERT_TEST_SYNOPSIS_E2E_FAMILY")
  }
  started <- proc.time()[["elapsed"]]
  cat("\n== ", gate$name, " ==\n", sep = "")
  devtools::test(
    pkg = client_root, filter = gate$filter, reporter = "summary",
    stop_on_failure = TRUE)
  elapsed <- proc.time()[["elapsed"]] - started
  cat("== ", gate$name, " passed in ",
      format(round(elapsed, 1), nsmall = 1), " s ==\n", sep = "")
}
