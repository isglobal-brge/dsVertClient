test_that("installed scripts resolve the installed dsVertClient namespace", {
  scripts <- c(
    "benchmark_dsi_transport.R",
    "benchmark_typed_blob_memory.R",
    "benchmark_typed_source_streaming.R",
    "validate_dp_statistical_methods.R"
  )
  paths <- system.file("scripts", scripts, package = "dsVertClient")

  expect_true(all(nzchar(paths)))
  expect_true(all(file.exists(paths)))
  contents <- vapply(
    paths,
    function(path) paste(readLines(path, warn = FALSE), collapse = "\n"),
    character(1L)
  )
  expect_false(any(grepl("pkgload", contents, fixed = TRUE)))
  expect_false(any(grepl("package_root", contents, fixed = TRUE)))

  if (.Platform$OS.type != "windows") {
    expect_equal(unname(file.access(paths, mode = 1L)), rep(0L, length(paths)))
  }

  package_path <- normalizePath(
    find.package("dsVertClient"), mustWork = TRUE)
  skip_if_not(
    file.exists(file.path(package_path, "Meta", "package.rds")),
    "subprocess preflight requires an installed package tree")

  isolated_user_library <- tempfile("dsvert-client-user-library-")
  dir.create(isolated_user_library)
  on.exit(unlink(isolated_user_library, recursive = TRUE), add = TRUE)
  package_library <- dirname(package_path)
  process_libraries <- paste(
    unique(c(package_library, .libPaths())), collapse = .Platform$path.sep)
  process_env <- c(
    paste0("R_LIBS_USER=", shQuote(isolated_user_library)),
    paste0("R_LIBS=", shQuote(process_libraries))
  )
  rscript <- file.path(
    R.home("bin"),
    if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")

  for (path in paths) {
    output <- suppressWarnings(system2(
      rscript,
      c(shQuote(path), "--preflight"),
      stdout = TRUE,
      stderr = TRUE,
      env = process_env
    ))
    status <- attr(output, "status")
    expect_true(
      is.null(status) || identical(status, 0L),
      info = paste(output, collapse = "\n")
    )
    expect_match(
      paste(output, collapse = "\n"),
      "^dsVertClient preflight OK: ",
      info = basename(path)
    )
  }
})
