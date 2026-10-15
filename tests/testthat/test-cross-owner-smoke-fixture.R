test_that("K3 smoke includes a same-owner pair and all three data owners", {
  source(test_path("..", "..", "tools", "cross_owner_smoke_fixture.R"), local = TRUE)
  fixture <- vertcor142_cor_fixture(3L, 64L)
  owner_counts <- table(sub("[$].*$", "", fixture$refs))
  expect_identical(as.integer(owner_counts), c(2L, 1L, 1L))
  expect_identical(names(owner_counts), c("site_a", "site_b", "site_c"))
  values <- vapply(fixture$refs, function(ref) {
    parts <- strsplit(ref, "$", fixed = TRUE)[[1L]]
    fixture$raw[[parts[[1L]]]][[parts[[2L]]]]
  }, numeric(64L))
  expect_equal(unname(colSums(values)), rep(32, 4L))
  expect_equal(unname(crossprod(sweep(values, 2L, colMeans(values)))), diag(16, 4L))
  expect_identical(fixture$registered, c("a", "b", "site_b$x", "site_c$x"))
  expect_true(all(vapply(fixture$raw, nrow, integer(1L)) == fixture$n))
  expect_identical(fixture$designated, c("site_a", "site_b"))
})

test_that("K5 keeps the requested block and includes two additional signed members", {
  source(test_path("..", "..", "tools", "cross_owner_smoke_fixture.R"), local = TRUE)
  fixture <- vertcor142_cor_fixture(5L, 64L)
  expect_identical(fixture$refs, vertcor142_cor_fixture(3L, 64L)$refs)
  expect_length(fixture$owners, 5L)
  expect_identical(names(fixture$raw$site_d), c("patient_id", "membership_marker"))
  expect_identical(names(fixture$raw$site_e), c("patient_id", "membership_marker"))
  expect_identical(fixture$raw$site_d$patient_id, fixture$raw$site_a$patient_id)
  expect_false(any(grepl("membership_marker", fixture$registered, fixed = TRUE)))
})

test_that("smoke fixtures reject incomplete balance and unsupported federations", {
  source(test_path("..", "..", "tools", "cross_owner_smoke_fixture.R"), local = TRUE)
  expect_error(vertcor142_cor_fixture(3L, 63L))
  expect_error(vertcor142_cor_fixture(4L, 64L))
  expect_error(vertcor142_cor_fixture(3L, 1025L))
  expect_identical(vertcor142_cor_fixture(2L, 8L)$registered, c("a", "x"))
})
