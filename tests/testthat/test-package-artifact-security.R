test_that("installed artifacts never contain serialized R sessions", {
  root <- testthat::test_path("..", "..")
  shipped <- list.files(
    file.path(root, "inst"), recursive = TRUE, full.names = TRUE,
    all.files = TRUE, include.dirs = FALSE)
  serialized <- shipped[grepl(
    paste0("\\.(rds|rda|rdata|rhistory|qs|fst|feather|parquet|",
           "sqlite|sqlite3|db|pem|key|p12|pfx)$"),
    shipped, ignore.case = TRUE)]

  expect_identical(serialized, character())

  forbidden_magic <- list(
    gzip = as.raw(c(0x1f, 0x8b)),
    bzip2 = charToRaw("BZh"),
    xz = as.raw(c(0xfd, 0x37, 0x7a, 0x58, 0x5a, 0x00)),
    r_serialized_xdr_v3 = as.raw(c(0x58, 0x0a, 0, 0, 0, 3)),
    r_serialized_ascii_v3 = charToRaw("A\n3\n"),
    r_workspace_v3 = charToRaw("RDX3\n"),
    sqlite = c(charToRaw("SQLite format 3"), as.raw(0)),
    parquet = charToRaw("PAR1"),
    arrow = charToRaw("ARROW1")
  )
  has_magic <- function(path, magic) {
    size <- file.info(path)$size
    if (is.na(size) || size < length(magic)) return(FALSE)
    con <- file(path, open = "rb")
    on.exit(close(con), add = TRUE)
    identical(readBin(con, what = "raw", n = length(magic)), magic)
  }
  disguised <- shipped[vapply(shipped, function(path) {
    any(vapply(forbidden_magic, has_magic, logical(1L), path = path))
  }, logical(1L))]

  expect_identical(disguised, character())
})
