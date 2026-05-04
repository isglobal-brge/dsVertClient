
validation_data_path <- function(file) {
  candidates <- c(
    file.path("validation-data", file),
    file.path("..", "validation-data", file),
    file.path("vignettes", "validation-data", file),
    file.path("..", "vignettes", "validation-data", file)
  )
  hit <- candidates[file.exists(candidates)][1]
  if (is.na(hit)) stop("Validation data file not found: ", file, call. = FALSE)
  hit
}

load_validation_summary <- function() {
  x <- utils::read.csv(validation_data_path("validation_summary.csv"),
                       stringsAsFactors = FALSE)
  x$observed <- as.numeric(x$observed)
  x$tolerance <- as.numeric(x$tolerance)
  x$non_disclosive <- as.logical(x$non_disclosive)
  x
}

load_legacy_removed <- function() {
  utils::read.csv(validation_data_path("legacy_routes_removed.csv"),
                  stringsAsFactors = FALSE)
}

validation_rows <- function(method_id) {
  x <- load_validation_summary()
  out <- x[x$method_id == method_id, , drop = FALSE]
  if (!nrow(out)) stop("No validation rows for method_id=", method_id,
                       call. = FALSE)
  rownames(out) <- NULL
  out
}

assert_validation <- function(rows) {
  bad <- rows[
    rows$status != "PASS" |
      !isTRUE(all(rows$non_disclosive)) |
      rows$observed > rows$tolerance + sqrt(.Machine$double.eps),
    , drop = FALSE]
  if (nrow(bad)) {
    stop("Validation evidence outside accepted envelope:\n",
         paste(utils::capture.output(print(bad)), collapse = "\n"),
         call. = FALSE)
  }
  invisible(TRUE)
}

display_validation <- function(rows) {
  keep <- c("k_mode", "function_route", "dataset", "reference_target",
            "primary_metric", "observed", "tolerance", "tier", "status",
            "cache")
  out <- rows[, keep, drop = FALSE]
  out$observed <- signif(out$observed, 4)
  out$tolerance <- signif(out$tolerance, 4)
  knitr::kable(out)
}

assert_legacy_routes_removed <- function() {
  pkg_roots <- c(".", "..", file.path("..", ".."))
  pkg_roots <- pkg_roots[file.exists(file.path(pkg_roots, "DESCRIPTION"))]
  if (requireNamespace("pkgload", quietly = TRUE) && length(pkg_roots)) {
    try(pkgload::load_all(pkg_roots[[1]], quiet = TRUE), silent = TRUE)
  }
  ns <- asNamespace("dsVertClient")
  checks <- data.frame(
    route = c(
      "ds.vertCox.k3",
      "ds.vertCox(method=...)",
      "ds.vertNBFullRegTheta(variant='full_reg')",
      "ds.vertMultinom(method='warm')",
      "ds.vertOrdinal(method='warm')",
      "ds.vertGLMM(method=...)"
    ),
    pass = c(
      !exists("ds.vertCox.k3", envir = ns, inherits = FALSE),
      !"method" %in% names(formals(get("ds.vertCox", ns))),
      tryCatch({
        dsVertClient::ds.vertNBFullRegTheta(
          y ~ x, data = "D", variant = "full_reg", datasources = list())
        FALSE
      }, error = function(e) grepl("variant must", conditionMessage(e))),
      !"method" %in% names(formals(get("ds.vertMultinom", ns))),
      !"method" %in% names(formals(get("ds.vertOrdinal", ns))),
      !"method" %in% names(formals(get("ds.vertGLMM", ns)))
    )
  )
  if (!all(checks$pass)) {
    stop("A removed legacy route is still invokable.", call. = FALSE)
  }
  checks
}

