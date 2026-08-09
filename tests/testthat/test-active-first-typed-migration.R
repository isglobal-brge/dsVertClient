.client_r_text <- function(files) {
  paste(unlist(lapply(files, readLines, warn = FALSE), use.names = FALSE),
        collapse = "\n")
}

test_that("active first-wave callers no longer select generic blob slots", {
  root <- .dsvert_client_source_root()
  files <- list.files(root, pattern = "[.]R$", full.names = TRUE)
  quarantine <- file.path(root, "ds.vertLMM.closed_form.R")
  active <- setdiff(files, quarantine)
  text <- .client_r_text(active)
  retired_slots <- c(
    '"k2_peer_x_share"', '"k2_peer_y_share"',
    '"k2_grad_peer_r1"', '"k2_beaver_vecmul_peer_masked"',
    '"k2_peer_weight_share"', '"k2_peer_sqrt_weight_share"',
    '"cor_col_params"')
  for (slot in retired_slots) {
    expect_false(grepl(slot, text, fixed = TRUE), info = slot)
  }

  # Two old closed-form LMM relays remain explicitly quarantined; this test
  # makes any accidental spread of that generic slot visible.
  quarantine_text <- .client_r_text(quarantine)
  expect_identical(lengths(regmatches(
    quarantine_text,
    gregexpr('"k2_beaver_vecmul_peer_masked"', quarantine_text,
             fixed = TRUE))), 2L)
})

test_that("IKNP and correlation callers use typed/direct contracts", {
  root <- .dsvert_client_source_root()
  ot <- .client_r_text(file.path(root, "ot_beaver_helpers.R"))
  cor <- .client_r_text(file.path(root, "ds.vertCor.R"))
  expect_false(grepl(".dsvert_store_blob", ot, fixed = TRUE))
  expect_false(grepl("points_blob_key =", ot, fixed = TRUE))
  expect_false(grepl("ciphertexts_blob_key =", ot, fixed = TRUE))
  expect_false(grepl("u_matrix_blob_key =", ot, fixed = TRUE))
  expect_false(grepl("cor_col_params", cor, fixed = TRUE))
  expect_false(grepl("from_storage = TRUE", cor, fixed = TRUE))
})

test_that("typed relay wrappers retain access to their producer connection", {
  root <- .dsvert_client_source_root()
  files <- list.files(root, pattern = "[.]R$", full.names = TRUE)
  callers <- files[vapply(files, function(path) {
    grepl(".dsvert_store_transfer_or_legacy(",
          .client_r_text(path), fixed = TRUE)
  }, logical(1L))]
  callers <- setdiff(callers, file.path(root, "typed_blob_transport.R"))
  expect_gt(length(callers), 0L)
  for (path in callers) {
    expect_match(.client_r_text(path), "producer_conns = datasources",
                 fixed = TRUE, info = basename(path))
  }
})
