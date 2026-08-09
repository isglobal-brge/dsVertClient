test_that("eta privacy auto-selection is deterministic by federation size", {
  select <- dsVertClient:::.dsvert_select_eta_privacy
  expect_identical(select("auto", 2L), "k2_beaver")
  expect_identical(select("auto", 3L), "secure_agg")
  expect_identical(select("auto", 20L), "secure_agg")
})

test_that("eta privacy rejects obsolete and incompatible routes", {
  select <- dsVertClient:::.dsvert_select_eta_privacy
  for (obsolete in c("transport", "he_link", "plaintext", "invalid")) {
    expect_error(select(obsolete, 2L), "eta_privacy")
  }
  expect_error(select("secure_agg", 2L), "at least 3 servers")
  expect_error(select("k2_beaver", 3L), "exactly 2 servers")
  expect_error(select("auto", 1L), "At least two")
  expect_error(select("auto", 2.5), "At least two")
})

test_that("only the two live eta protocols are accepted explicitly", {
  select <- dsVertClient:::.dsvert_select_eta_privacy
  expect_identical(select("k2_beaver", 2L), "k2_beaver")
  expect_identical(select("secure_agg", 3L), "secure_agg")
})
