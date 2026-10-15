.irls_chain_fixtures <- function() {
  path <- system.file("dp-irls-fixtures", "synthetic-chains.json",
                      package = "dsVertClient")
  jsonlite::fromJSON(paste(readLines(path, warn = FALSE), collapse = ""),
                     simplifyVector = FALSE)
}

test_that("offline chain verification independently recomputes full and score-only paths", {
  for (f in .irls_chain_fixtures()) {
    verified <- .dsvert_dp_irls_verify_chain(f$signed_contract, f$pins, f$stages)
    expect_true(verified$complete)
    expect_identical(.dsvert_dp_irls_client_state_json(verified$beta),
      .dsvert_dp_irls_client_state_json(f$stages[[2L]]$body$update_result$beta))
    expect_identical(verified$accounting, f$accounting)
    expect_identical(verified$root, f$stages[[1L]]$body$root)
    partial <- .dsvert_dp_irls_verify_chain(f$signed_contract, f$pins, f$stages[1L])
    expect_false(partial$complete)
    expect_identical(partial$accounting$charged_epsilon, list(d = "1", n = "2"))
    expect_identical(partial$accounting$reserved_epsilon,
                     f$signed_contract$contract$accounting$epsilon_path)
    empty <- .dsvert_dp_irls_verify_chain(f$signed_contract, f$pins, list())
    expect_false(empty$complete)
    expect_identical(empty$beta, f$signed_contract$contract$path$start$beta)
    expect_error(.dsvert_dp_irls_verify_chain(f$signed_contract, f$pins,
      f$stages, expected_root = strrep("f", 64L)), "expected path")
  }
})

test_that("offline verifier rejects false arithmetic even with both valid signatures", {
  for (f in .irls_chain_fixtures()) {
    contract <- .dsvert_dp_irls_client_contract_validate(f$signed_contract, f$pins)
    expect_silent(.dsvert_dp_irls_client_stage_signatures(f$forged_update, contract))
    expect_error(.dsvert_dp_irls_verify_chain(f$signed_contract, f$pins,
      list(f$forged_update)), "incorrect update")
  }
})

test_that("offline chains reject signature, pin, predecessor, epoch and order tampering", {
  f <- .irls_chain_fixtures()[[1L]]
  bad <- f$stages
  bad[[1L]]$signatures[[1L]] <- strrep("A", 86L)
  expect_error(.dsvert_dp_irls_verify_chain(f$signed_contract, f$pins, bad), "signature")
  pins <- f$pins
  pins[[1L]] <- pins[[2L]]
  expect_error(.dsvert_dp_irls_verify_chain(f$signed_contract, pins, f$stages))
  bad <- f$stages
  bad[[2L]]$body$predecessor <- strrep("b", 64L)
  expect_error(.dsvert_dp_irls_verify_chain(f$signed_contract, f$pins, bad))
  bad <- f$signed_contract
  bad$contract$privacy_source$closure$dependencies[[1L]]$dependency$authorization_epoch <- strrep("a",64L)
  expect_error(.dsvert_dp_irls_verify_chain(bad, f$pins, f$stages))
  expect_error(.dsvert_dp_irls_verify_chain(f$signed_contract, f$pins, rev(f$stages)))
  expect_error(.dsvert_dp_irls_verify_chain(f$signed_contract, f$pins,
                                          c(f$stages, f$stages[3L])))
  bad <- f$stages
  bad[[2L]]$signatures <- bad[[2L]]$signatures[1L]
  expect_error(.dsvert_dp_irls_verify_chain(f$signed_contract, f$pins, bad), "authority")
  # No binary exists here; verification uses R/gmp and OpenSSL only.
  withr::with_envvar(c(DS_VERT_MPC = "/does/not/exist", DS_MPC_BIN = "/does/not/exist"),
    expect_true(.dsvert_dp_irls_verify_chain(f$signed_contract, f$pins, f$stages)$complete))
})
