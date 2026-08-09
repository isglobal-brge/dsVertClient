.canonical_session_source <- function(filename) {
  paste(readLines(.dsvert_client_source_file(filename), warn = FALSE),
        collapse = "\n")
}

test_that("active typed model families use unprefixed canonical UUID sessions", {
  files <- c(
    multinomial = "ds.vertMultinomJointNewton.R",
    ordinal = "ds.vertOrdinalJointNewton.R",
    cox = "ds.vertCoxDiscreteNonDisclosive.R",
    negative_binomial = "ds.vertNBFullRegTheta.R")
  for (family in names(files)) {
    source <- .canonical_session_source(files[[family]])
    expect_match(
      source, "session_id\\s*<-\\s*\\.dsvert_uuid4\\(\\)",
      info = family)
    expect_false(
      grepl("session_id\\s*<-\\s*(?:paste0|sprintf)\\([^\\n]*uuid4",
            source, perl = TRUE),
      info = family)
  }
})

test_that("all package-generated MPC sessions satisfy the relay grammar", {
  session_re <- paste0(
    "^[0-9a-f]{8}-[0-9a-f]{4}-4[0-9a-f]{3}-",
    "[89ab][0-9a-f]{3}-[0-9a-f]{12}$")
  generated <- vapply(seq_len(32L), function(index) .dsvert_uuid4(),
                      character(1L))
  expect_true(all(grepl(session_re, generated)))
  expect_length(unique(generated), length(generated))
  expect_match(.mpc_session_id(), session_re)

  source_files <- list.files(.dsvert_client_source_root(),
                             pattern = "\\.[Rr]$", full.names = TRUE)
  source <- paste(unlist(lapply(
    source_files, readLines, warn = FALSE), use.names = FALSE), collapse = "\n")
  expect_false(grepl(
    "session_id\\s*<-\\s*(?:paste0|sprintf)\\([^\\n]*\\.dsvert_uuid4\\(\\)",
    source, perl = TRUE))
})
