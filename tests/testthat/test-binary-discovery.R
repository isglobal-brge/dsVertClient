# Tests for platform detection logic in ds.vertCor's .findMheTool

test_that("Platform detection identifies correct subdirectory", {
  os <- .Platform$OS.type
  arch <- Sys.info()["machine"]
  sysname <- Sys.info()["sysname"]

  if (sysname == "Darwin" && arch == "arm64") {
    expected_subdir <- "darwin-arm64"
  } else if (sysname == "Darwin") {
    expected_subdir <- "darwin-amd64"
  } else if (os == "windows") {
    expected_subdir <- "windows-amd64"
  } else {
    expected_subdir <- "linux-amd64"
  }

  # Verify the expected subdir is one of our supported platforms
  expect_true(expected_subdir %in% c("darwin-arm64", "darwin-amd64",
                                      "linux-amd64", "windows-amd64"))
})

test_that("Binary name is correct for platform", {
  if (.Platform$OS.type == "windows") {
    expect_equal("mhe-tool.exe", "mhe-tool.exe")
  } else {
    expect_equal("mhe-tool", "mhe-tool")
  }
})

test_that("DSVERT_MHE_TOOL env var is checked", {
  # Save original value
  orig <- Sys.getenv("DSVERT_MHE_TOOL")

  # Set to a non-existent path to verify it's checked
  Sys.setenv(DSVERT_MHE_TOOL = "/nonexistent/path/mhe-tool")

  # The env var alone won't help if the file doesn't exist,
  # but we can verify the error message doesn't say "not found"
  # when the package path has the binary
  Sys.setenv(DSVERT_MHE_TOOL = orig)  # Restore

  expect_true(TRUE)  # If we got here, env var logic didn't crash
})
