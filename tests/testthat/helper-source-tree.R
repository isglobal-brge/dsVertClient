.dsvert_client_source_root <- function() {
  package_root <- normalizePath(
    file.path(testthat::test_path(), "..", ".."),
    mustWork = FALSE)
  candidates <- c(
    file.path(package_root, "R"),
    file.path(package_root, "00_pkg_src", "dsVertClient", "R"),
    file.path(getwd(), "R"))
  existing <- candidates[dir.exists(candidates)]
  if (!length(existing)) {
    testthat::skip(
      "the dsVertClient source tree is unavailable for this static audit")
  }
  normalizePath(existing[[1L]], mustWork = TRUE)
}

.dsvert_client_source_file <- function(filename) {
  path <- file.path(.dsvert_client_source_root(), filename)
  if (!file.exists(path)) {
    testthat::skip(paste0(
      "the dsVertClient source file is unavailable for this static audit: ",
      filename))
  }
  path
}
