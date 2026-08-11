.formal_dp_e2e_source_root <- function() {
  candidates <- unique(normalizePath(c(
    file.path(testthat::test_path(), "..", ".."),
    getwd(), file.path(getwd(), "..")),
    mustWork = FALSE))
  selected <- candidates[file.exists(file.path(
    candidates, "tools", "validate_formal_dp_e2e.R"))]
  if (!length(selected)) return(NULL)
  selected[[1L]]
}

.formal_dp_e2e_binary <- function(client_root) {
  server_root <- normalizePath(file.path(client_root, "..", "dsVert"),
                               mustWork = FALSE)
  system <- Sys.info()[["sysname"]]
  machine <- tolower(Sys.info()[["machine"]])
  platform <- if (identical(system, "Darwin") &&
      machine %in% c("arm64", "aarch64")) {
    "darwin-arm64"
  } else if (identical(system, "Darwin") &&
             machine %in% c("x86_64", "amd64")) {
    "darwin-amd64"
  } else if (identical(system, "Linux") &&
             machine %in% c("x86_64", "amd64")) {
    "linux-amd64"
  } else {
    return(NULL)
  }
  path <- file.path(server_root, "inst", "bin", platform, "dsvert-mpc")
  if (!file.exists(path)) NULL else normalizePath(path, mustWork = TRUE)
}

test_that("formal DP E2E harness is source-only and has a safe help path", {
  expect_identical(system.file(
    "scripts", "validate_formal_dp_e2e.R", package = "dsVertClient"), "")
  root <- .formal_dp_e2e_source_root()
  skip_if(is.null(root), "source-tree harness is not present after install")
  skip_if_not_installed("callr")
  script <- file.path(root, "tools", "validate_formal_dp_e2e.R")
  code <- paste(readLines(script, warn = FALSE), collapse = "\n")
  expect_false(grepl("server-private\\.log|DSVERT_E2E_PRIVATE", code))
  expect_false(grepl("Sys\\.chmod\\(path", code))
  expect_false(grepl(
    "dsvert\\.dp\\.total_delta[[:space:]]*=[[:space:]]*0", code))
  expect_false(grepl("dsVert:::\\.dsvert_dp_policy\\(\\)", code))
  expect_false(grepl("ds\\.vertDPStatus|ds\\.vertDPCapsulePlan", code))
  expect_false(grepl("noise_root|ledger_path|total_epsilon", code))
  expect_match(code, "dsvert.dp.implementation_delta = 1e-9", fixed = TRUE)
  expect_match(code, "first$artifact_key, replay$artifact_key", fixed = TRUE)
  expect_match(
    code, "sampler_recomputed_after_service_restart", fixed = TRUE)

  result <- callr::rscript(
    script, "--help", stdout = "|", stderr = "|", show = FALSE,
    wd = dirname(root))
  expect_identical(result$status, 0L)
  expect_match(result$stdout, "Source-tree", fixed = TRUE)
  expect_match(result$stdout, "--k=2,3,5", fixed = TRUE)
  expect_match(result$stdout, "--keep-state=PATH", fixed = TRUE)
})

test_that("formal DP E2E state cleans on success and failure unless retained", {
  root <- .formal_dp_e2e_source_root()
  skip_if(is.null(root), "source-tree harness is not present after install")
  skip_if_not_installed("callr")
  skip_if_not_installed("pkgload")
  script <- file.path(root, "tools", "validate_formal_dp_e2e.R")
  work <- withr::local_tempdir(.local_envir = environment())
  run <- function(mode, marker, keep = NULL, k = NULL) {
    cmdargs <- c(
      paste0("--cleanup-self-test=", mode),
      paste0("--cleanup-marker=", marker))
    if (!is.null(k)) cmdargs <- c(cmdargs, paste0("--k=", k))
    if (!is.null(keep)) cmdargs <- c(
      cmdargs, paste0("--keep-state=", keep))
    callr::rscript(
      script, cmdargs, stdout = "|", stderr = "|", show = FALSE,
      wd = dirname(root), fail_on_status = FALSE)
  }

  binary <- .formal_dp_e2e_binary(root)
  mode_before <- if (is.null(binary)) NULL else file.info(binary)$mode

  success_marker <- file.path(work, "success.marker")
  success <- run("success", success_marker)
  success_state <- readLines(success_marker, warn = FALSE)[[1L]]
  expect_identical(success$status, 0L)
  expect_false(file.exists(success_state))

  k5_marker <- file.path(work, "k5.marker")
  k5 <- run("success", k5_marker, k = 5L)
  k5_state <- readLines(k5_marker, warn = FALSE)[[1L]]
  expect_identical(k5$status, 0L)
  expect_false(file.exists(k5_state))

  failure_marker <- file.path(work, "failure.marker")
  failure <- run("failure", failure_marker)
  failure_state <- readLines(failure_marker, warn = FALSE)[[1L]]
  expect_false(identical(failure$status, 0L))
  expect_false(file.exists(failure_state))

  keep_marker <- file.path(work, "keep.marker")
  retained <- file.path(work, "explicitly-retained")
  kept <- run("success", keep_marker, retained)
  ephemeral <- readLines(keep_marker, warn = FALSE)[[1L]]
  expect_identical(kept$status, 0L)
  expect_false(file.exists(ephemeral))
  expect_true(dir.exists(retained))

  if (!is.null(binary)) {
    expect_identical(file.info(binary)$mode, mode_before)
  }
})
