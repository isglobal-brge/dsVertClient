test_that("client DSI text framing uses the bounded twin implementation", {
  value <- paste0(
    "{\"payload\":\"", strrep("A_", 512L * 1024L), "\"}")
  original_safe <- .dsvert_dsi_text_safe_raw
  largest_scan <- 0L
  encoded <- testthat::with_mocked_bindings(
    .dsvert_dsi_text_encode(value),
    .dsvert_dsi_text_safe_raw = function(bytes) {
      largest_scan <<- max(largest_scan, length(bytes))
      original_safe(bytes)
    }, .package = "dsVertClient")
  expect_true(startsWith(encoded, "DSV1L_"))
  expect_lte(largest_scan, .DSVERT_DSI_TEXT_CHUNK_BYTES)

  largest_scan <- 0L
  expect_identical(testthat::with_mocked_bindings(
    .dsvert_dsi_text_decode(
      encoded, maximum_bytes = nchar(value, type = "bytes")),
    .dsvert_dsi_text_safe_raw = function(bytes) {
      largest_scan <<- max(largest_scan, length(bytes))
      original_safe(bytes)
    },
    .dsvert_dsi_text_encode = function(...) {
      stop("canonical decode rebuilt the complete frame", call. = FALSE)
    }, .package = "dsVertClient"), value)
  expect_lte(largest_scan, .DSVERT_DSI_TEXT_CHUNK_BYTES)
})

test_that("client and server codec twins stay byte-identical", {
  skip_if_not_installed("dsVert")
  server <- asNamespace("dsVert")
  required <- c(
    ".dsvert_dsi_text_encode", ".dsvert_dsi_text_decode",
    ".dsvert_dsi_text_encode_l", ".dsvert_dsi_text_decode_l_pass",
    ".DSVERT_DSI_TEXT_CHUNK_BYTES", ".DSVERT_DSI_TEXT_L_BATCH_RUNS",
    ".DSVERT_DSI_TEXT_HEX_VALUES")
  if (!all(vapply(required, exists, logical(1L), envir = server,
                  inherits = FALSE))) {
    skip("requires the paired current dsVert source namespace")
  }
  server_encode <- get(".dsvert_dsi_text_encode", server)
  server_decode <- get(".dsvert_dsi_text_decode", server)
  expect_identical(
    get(".DSVERT_DSI_TEXT_CHUNK_BYTES", server),
    .DSVERT_DSI_TEXT_CHUNK_BYTES)
  expect_identical(
    get(".DSVERT_DSI_TEXT_L_BATCH_RUNS", server),
    .DSVERT_DSI_TEXT_L_BATCH_RUNS)
  expect_identical(
    get(".DSVERT_DSI_TEXT_HEX_VALUES", server),
    .DSVERT_DSI_TEXT_HEX_VALUES)
  expect_identical(
    deparse(body(get(".dsvert_dsi_text_encode_l", server)),
            width.cutoff = 500L,
            control = c("keepNA", "keepInteger", "niceNames")),
    deparse(body(.dsvert_dsi_text_encode_l), width.cutoff = 500L,
            control = c("keepNA", "keepInteger", "niceNames")))
  expect_identical(
    deparse(body(get(".dsvert_dsi_text_decode_l_pass", server)),
            width.cutoff = 500L,
            control = c("keepNA", "keepInteger", "niceNames")),
    deparse(body(.dsvert_dsi_text_decode_l_pass), width.cutoff = 500L,
            control = c("keepNA", "keepInteger", "niceNames")))

  set.seed(20260809)
  alphabet <- c(
    letters, LETTERS, as.character(0:9), "_", "-", "{", "}",
    "\"", "\\", ":", ",", "[", "]", " ", "é")
  values <- c("", replicate(512L, paste0(sample(
    alphabet, sample.int(1024L, 1L) - 1L, replace = TRUE),
    collapse = "")))
  for (value in values) {
    client_wire <- .dsvert_dsi_text_encode(value)
    server_wire <- server_encode(value)
    expect_identical(client_wire, server_wire)
    expect_identical(.dsvert_dsi_text_decode(client_wire), enc2utf8(value))
    expect_identical(server_decode(client_wire), enc2utf8(value))
  }
})
