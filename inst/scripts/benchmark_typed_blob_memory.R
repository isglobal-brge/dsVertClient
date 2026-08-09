#!/usr/bin/env Rscript

# Isolated-process memory benchmark for the current typed producer-to-client
# pump.  The worker deliberately starts with one complete Base64url value,
# because that is the production contract being measured rather than the
# disk-backed relay substrate.

args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) {
  cat(paste(
    "Usage: benchmark_typed_blob_memory.R [--sizes-mib=N,...]",
    "[--output=PATH] [--preflight]\n"
  ))
  quit(save = "no", status = 0L)
}
if (!requireNamespace("dsVertClient", quietly = TRUE)) {
  stop("Install dsVertClient before running this script", call. = FALSE)
}
if ("--preflight" %in% args) {
  cat("dsVertClient preflight OK: ",
      as.character(utils::packageVersion("dsVertClient")), "\n", sep = "")
  quit(save = "no", status = 0L)
}

suppressPackageStartupMessages(library(jsonlite))

value_arg <- function(prefix, default = NULL) {
  hit <- args[startsWith(args, paste0(prefix, "="))]
  if (!length(hit)) return(default)
  sub(paste0("^", prefix, "="), "", hit[[1L]])
}
worker_mib <- suppressWarnings(as.numeric(value_arg("--worker-mib")))
sizes_mib <- suppressWarnings(as.numeric(strsplit(
  value_arg("--sizes-mib", "1,16,64"), ",", fixed = TRUE)[[1L]]))
output_path <- value_arg("--output", "")

script_arg <- commandArgs(FALSE)[startsWith(commandArgs(FALSE), "--file=")]
script_path <- if (length(script_arg)) {
  normalizePath(sub("^--file=", "", script_arg[[1L]]), mustWork = TRUE)
} else {
  system.file(
    "scripts", "benchmark_typed_blob_memory.R",
    package = "dsVertClient", mustWork = TRUE)
}

load_client <- function() {
  if (!requireNamespace("dsVertClient", quietly = TRUE)) {
    stop("Install dsVertClient before running this script", call. = FALSE)
  }
}

if (length(worker_mib) == 1L && !is.na(worker_mib)) {
  if (!is.finite(worker_mib) || worker_mib < 0 ||
      worker_mib != floor(worker_mib)) stop("Invalid worker size", call. = FALSE)
  load_client()
  if (worker_mib == 0) {
    cat("DSVERT_TYPED_RESULT ", jsonlite::toJSON(list(
      payload_mib = 0, payload_chars = 0, frame_chars = 0,
      frames = 0, dsi_calls = 0, wall_seconds = 0,
      integrity = "baseline"), auto_unbox = TRUE), "\n", sep = "")
    quit(save = "no", status = 0L)
  }

  payload_chars <- as.numeric(worker_mib * 1024^2)
  blob <- strrep("A", payload_chars)
  transfer <- list(
    ticket = strrep("A", 86L),
    transfer_id = "tb_00112233445566778899aabbccddeeff",
    capability_id = "blob.input.peer-x.v1",
    sender_name = "producer", recipient_name = "recipient",
    payload_chars = payload_chars,
    payload_sha256 = paste0(openssl::sha256(charToRaw(blob))))
  spool <- tempfile(pattern = "dsvert-typed-memory-")
  on.exit(unlink(spool), add = TRUE)
  file.create(spool)
  Sys.chmod(spool, mode = "0600")
  committed <- 0
  calls <- 0L
  aggregate <- function(conns, expr, async, error, errors.print) {
    calls <<- calls + 1L
    method <- as.character(expr[[1L]])
    if (identical(method, "mpcTypedBlobReceiptDS")) {
      return(list(producer = list(
        version = "dsvert-typed-blob-receipt-v1",
        transfer_id = transfer$transfer_id, confirmed = TRUE)))
    }
    if (!identical(method, "mpcTypedBlobStoreDS")) {
      stop("Unexpected typed benchmark method", call. = FALSE)
    }
    offset <- as.numeric(expr$offset)
    chunk <- as.character(expr$chunk)
    if (!identical(offset, committed)) stop("Benchmark offset gap", call. = FALSE)
    con <- file(spool, "ab")
    writeBin(charToRaw(chunk), con)
    flush(con)
    close(con)
    committed <<- committed + nchar(chunk, type = "bytes")
    sealed <- identical(committed, payload_chars)
    if (sealed && !identical(
        digest::digest(file = spool, algo = "sha256", serialize = FALSE),
        transfer$payload_sha256)) stop("Benchmark digest mismatch", call. = FALSE)
    list(recipient = list(
      version = "dsvert-typed-blob-v1",
      transfer_id = transfer$transfer_id,
      committed_chars = committed, total_chars = payload_chars,
      sealed = sealed, receipt = if (sealed) "cmVjZWlwdA" else NULL))
  }

  old <- options(
    dsvert.chunk_size = 640L * 1024L,
    dsvert.dsi_max_expression_bytes = 768L * 1024L - 1L)
  on.exit(options(old), add = TRUE)
  dsVertClient:::.dsvert_reset_chunk_size()
  started <- proc.time()[["elapsed"]]
  frames <- dsVertClient:::.dsvert_store_typed_blob(
    blob, transfer,
    conn = list(recipient = structure(list(), class = "benchmark")),
    session_id = "12345678-1234-4234-8234-123456789abc",
    producer_conn = list(
      producer = structure(list(), class = "benchmark")),
    .aggregate = aggregate)
  wall <- unname(proc.time()[["elapsed"]] - started)
  integrity <- if (identical(committed, payload_chars) && identical(
      digest::digest(file = spool, algo = "sha256", serialize = FALSE),
      transfer$payload_sha256)) "sha256_identical" else "mismatch"
  cat("DSVERT_TYPED_RESULT ", jsonlite::toJSON(list(
    payload_mib = worker_mib, payload_chars = payload_chars,
    frame_chars = 640L * 1024L, frames = as.integer(frames),
    dsi_calls = calls, wall_seconds = wall, integrity = integrity),
    auto_unbox = TRUE, digits = NA), "\n", sep = "")
  quit(save = "no", status = 0L)
}

if (!length(sizes_mib) || anyNA(sizes_mib) || any(!is.finite(sizes_mib)) ||
    any(sizes_mib < 1) || any(sizes_mib != floor(sizes_mib))) {
  stop("--sizes-mib must contain positive whole MiB values", call. = FALSE)
}
time_bin <- "/usr/bin/time"
if (!file.exists(time_bin)) stop("/usr/bin/time is required", call. = FALSE)
system_name <- Sys.info()[["sysname"]]
time_flag <- if (identical(system_name, "Darwin")) "-l" else "-v"
rscript <- file.path(R.home("bin"), "Rscript")

run_isolated <- function(size) {
  output <- suppressWarnings(system2(
    time_bin,
    c(time_flag, shQuote(rscript), shQuote(script_path),
      paste0("--worker-mib=", format(size, scientific = FALSE, trim = TRUE))),
    stdout = TRUE, stderr = TRUE))
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) {
    stop("Typed memory worker failed", call. = FALSE)
  }
  marker <- output[startsWith(output, "DSVERT_TYPED_RESULT ")]
  if (length(marker) != 1L) stop("Typed memory worker returned no result", call. = FALSE)
  result <- jsonlite::fromJSON(sub("^DSVERT_TYPED_RESULT ", "", marker))
  rss_line <- output[grepl("maximum resident set size", output,
                          ignore.case = TRUE)]
  if (length(rss_line) != 1L) stop("Could not parse maximum RSS", call. = FALSE)
  if (identical(system_name, "Darwin")) {
    rss <- suppressWarnings(as.numeric(sub(
      "^\\s*([0-9]+).*$", "\\1", rss_line)))
  } else {
    rss <- suppressWarnings(as.numeric(sub("^.*:\\s*([0-9]+).*$", "\\1",
                                          rss_line))) * 1024
  }
  if (!is.finite(rss) || rss <= 0) stop("Invalid maximum RSS", call. = FALSE)
  result$max_rss_bytes <- rss
  result
}

baseline <- run_isolated(0)
rows <- lapply(sizes_mib, function(size) {
  value <- run_isolated(size)
  data.frame(
    payload_mib = as.numeric(value$payload_mib),
    payload_chars = as.numeric(value$payload_chars),
    frame_chars = as.numeric(value$frame_chars),
    frames = as.integer(value$frames), dsi_calls = as.integer(value$dsi_calls),
    wall_seconds = as.numeric(value$wall_seconds),
    throughput_mib_s = as.numeric(value$payload_mib) /
      as.numeric(value$wall_seconds),
    max_rss_bytes = as.numeric(value$max_rss_bytes),
    baseline_max_rss_bytes = as.numeric(baseline$max_rss_bytes),
    rss_delta_bytes = as.numeric(value$max_rss_bytes - baseline$max_rss_bytes),
    source_memory_model = "O(payload)", integrity = value$integrity,
    stringsAsFactors = FALSE)
})
results <- do.call(rbind, rows)
rownames(results) <- NULL
print(results)
if (nzchar(output_path)) {
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(results, output_path, row.names = FALSE)
  message("Wrote ", normalizePath(output_path, mustWork = TRUE))
}
