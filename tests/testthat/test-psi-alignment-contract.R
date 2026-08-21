.client_padded_attestation <- function(contract_hash = strrep("4", 64L)) {
  list(
    attestation_version = 3L,
    alignment_attested = TRUE,
    alignment_protocol = "dsvert-pinned-padded-psi-v5",
    attestation_id = paste0("attest_", strrep("3", 64L)),
    contract_hash = contract_hash,
    policy_id = paste0("policy_", strrep("1", 64L)),
    alignment_purpose = "patient-record-alignment-v1",
    dataset_id = "test-logical-cohort",
    dataset_version = "v1",
    id_column = "patient_id",
    source_binding_id = paste0("source_", strrep("6", 64L)),
    pinset_id = paste0("pinset_", strrep("2", 64L)),
    capacity_bucket = 64L,
    relay_frame_bytes = 65536L,
    inline_max_bytes = 65536L,
    peer_count = 2L,
    reference_peer = "site_a",
    compute_peers = c("site_a", "site_b"))
}

test_that("alignment status accepts only identical padded attestations", {
  datasources <- list(site_a = list(), site_b = list())
  calls <- list()
  attestation <- .client_padded_attestation()
  fake_aggregate <- function(conns, expr, ...) {
    calls[[length(calls) + 1L]] <<- expr
    list(site_a = attestation, site_b = attestation)
  }

  result <- dsVertClient:::.psi_alignment_status(
    "DA", datasources, .aggregate = fake_aggregate)
  expect_true(result$aligned)
  expect_true(is.na(result$n_common))
  expect_length(calls, 1L)
  expect_identical(as.character(calls[[1L]][[1L]]),
                   "psiPaddedAttestationDS")
  expect_null(calls[[1L]]$session_id)
  expect_false(any(c("n", "hash", "token", "order_binding") %in%
                   names(result$manifests[[1L]])))
})

test_that("alignment status fails closed for mismatch, legacy or absence", {
  datasources <- list(site_a = list(), site_b = list())
  first <- .client_padded_attestation()
  second <- .client_padded_attestation(strrep("5", 64L))
  mismatch <- function(...) list(site_a = first, site_b = second)
  expect_false(dsVertClient:::.psi_alignment_status(
    "DA", datasources, .aggregate = mismatch)$aligned)

  legacy <- function(...) list(
    site_a = list(version = 2L, hash = strrep("a", 64L), n = 20L),
    site_b = list(version = 2L, hash = strrep("a", 64L), n = 20L))
  legacy_result <- dsVertClient:::.psi_alignment_status(
    "DA", datasources, .aggregate = legacy)
  expect_false(legacy_result$aligned)
  expect_true(is.na(legacy_result$n_common))

  absent <- function(...) stop("attestation absent")
  expect_false(dsVertClient:::.psi_alignment_status(
    "DA", datasources, .aggregate = absent)$aligned)
  expect_false(dsVertClient:::.psi_alignment_status(
    "DA", datasources[1L], .aggregate = mismatch)$aligned)
})

test_that("PSI resolves scalar and complete per-site routing aliases", {
  datasources <- list(site_a = list(), site_b = list(), site_c = list())
  expect_identical(
    dsVertClient:::.dsvert_site_character("D", datasources, "data_name"),
    c(site_a = "D", site_b = "D", site_c = "D"))
  expect_identical(
    dsVertClient:::.dsvert_site_character(
      c(site_c = "C", site_a = "A", site_b = "B"), datasources,
      "data_name"),
    c(site_a = "A", site_b = "B", site_c = "C"))
  expect_identical(
    dsVertClient:::.dsvert_site_character(
      list(site_a = "id_a", site_b = "id_b", site_c = "id_c"),
      datasources, "id_col"),
    c(site_a = "id_a", site_b = "id_b", site_c = "id_c"))
  expect_error(dsVertClient:::.dsvert_site_character(
    c(site_a = "A", site_b = "B"), datasources, "data_name"),
    "complete named")
  expect_error(dsVertClient:::.dsvert_site_character(
    c(site_a = "A", site_b = "", site_c = "C"), datasources,
    "data_name"), "non-empty")
})

test_that("public PSI front door delegates only to padded orchestration", {
  datasources <- list(site_a = list(), site_b = list())
  expected <- .client_padded_attestation()
  arguments <- NULL
  reset <- FALSE
  testthat::local_mocked_bindings(
    .dsvert_maybe_negotiate_dsi_chunk_size = function(conns) {
      expect_identical(conns, datasources)
      invisible(65536L)
    },
    .dsvert_reset_chunk_size = function() {
      reset <<- TRUE
      invisible(NULL)
    },
    .dsvert_psi_padded_align = function(...) {
      arguments <<- list(...)
      expected
    },
    .package = "dsVertClient")

  result <- ds.psiAlign(
    "D", "patient_id", "DA", verbose = FALSE,
    datasources = datasources, na.action = "none")
  expect_identical(result, expected)
  expect_identical(arguments$data_name, c(site_a = "D", site_b = "D"))
  expect_identical(arguments$id_col,
                   c(site_a = "patient_id", site_b = "patient_id"))
  expect_identical(arguments$newobj, "DA")
  expect_identical(arguments$datasources, datasources)
  expect_true(reset)
})

test_that("PSI client validates compatibility arguments and requires K >= 2", {
  expect_identical(eval(formals(ds.psiAlign)$na.action),
                   c("na.omit", "na.fail", "none"))
  expect_error(
    ds.psiAlign("D", "patient_id", datasources = list(),
                na.action = "silent"),
    "na.action must be one of")
  expect_error(ds.psiAlign("", "patient_id", datasources = list()),
               "data_name must be")
  expect_error(ds.psiAlign("D", NA_character_, datasources = list()),
               "id_col must be")
  expect_error(ds.psiAlign(
    "D", "patient_id", ref_server = "site_a",
    datasources = list(site_a = list(), site_b = list())),
    "ref_server must be NULL")
  expect_error(ds.psiAlign(
    "D", "patient_id", verbose = FALSE,
    datasources = list(site_a = list())),
    "at least two")
  expect_error(ds.psiAlign(
    "D", "patient_id", verbose = FALSE,
    datasources = list(list(), list())),
    "uniquely named")
})
