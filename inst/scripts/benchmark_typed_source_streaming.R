#!/usr/bin/env Rscript

# Multiprocess benchmark for the data-free typed-source pilot. Two persistent
# TCP workers model the producer and recipient DataSHIELD processes; the
# controller runs the real dsVertClient source pump. Connections are closed and
# re-established mid-transfer, and selected committed responses are discarded
# so the exact immutable requests must be replayed.

args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) {
  cat(paste(
    "Usage: benchmark_typed_source_streaming.R --binary=PATH",
    "[--payload-mib=N] [--frame-kib=N] [--preflight]\n"
  ))
  quit(save = "no", status = 0L)
}
arg_value <- function(name, default = NULL) {
  values <- args[startsWith(args, paste0(name, "="))]
  if (!length(values)) return(default)
  sub(paste0("^", name, "="), "", values[[1L]])
}
worker <- arg_value("--worker", "")
needs_client <- "--preflight" %in% args ||
  !worker %in% c("producer", "recipient")
if (needs_client && !requireNamespace("dsVertClient", quietly = TRUE)) {
  stop("Install dsVertClient before running this script", call. = FALSE)
}
if ("--preflight" %in% args) {
  cat("dsVertClient preflight OK: ",
      as.character(utils::packageVersion("dsVertClient")), "\n", sep = "")
  quit(save = "no", status = 0L)
}

suppressPackageStartupMessages({
  library(digest)
  library(jsonlite)
  library(processx)
})

port <- suppressWarnings(as.integer(arg_value("--port", NA_character_)))
ready_path <- arg_value("--ready", "")
work_path <- arg_value("--work", "")
binary_path <- arg_value("--binary", "")
spool_cap <- suppressWarnings(as.numeric(arg_value("--spool-cap", NA_character_)))

script_arg <- commandArgs(FALSE)[startsWith(commandArgs(FALSE), "--file=")]
script_path <- if (length(script_arg)) {
  normalizePath(sub("^--file=", "", script_arg[[1L]]), mustWork = TRUE)
} else {
  system.file(
    "scripts", "benchmark_typed_source_streaming.R",
    package = "dsVertClient", mustWork = TRUE)
}

rss_bytes <- function() {
  value <- tryCatch(
    ps::ps_memory_info(ps::ps_handle(Sys.getpid()))[["rss"]],
    error = function(e) NA_real_)
  as.numeric(value)
}

read_chars_at <- function(path, offset, count) {
  con <- file(path, "rb")
  on.exit(close(con), add = TRUE)
  seek(con, where = offset, origin = "start")
  value <- readBin(con, "raw", n = count)
  if (length(value) != count) stop("source frame was truncated", call. = FALSE)
  rawToChar(value)
}

run_listener <- function(handler) {
  if (is.na(port) || port < 1024L || !nzchar(ready_path)) {
    stop("invalid worker listener arguments", call. = FALSE)
  }
  listener <- serverSocket(port)
  on.exit(close(listener), add = TRUE)
  if (!file.create(ready_path)) stop("could not publish worker readiness")
  Sys.chmod(ready_path, "0600")
  shutdown <- FALSE
  while (!shutdown) {
    connection <- socketAccept(listener, blocking = TRUE, open = "a+b")
    repeat {
      request <- tryCatch(unserialize(connection), error = function(e) NULL)
      if (is.null(request)) break
      response <- tryCatch(
        list(ok = TRUE, value = handler(request)),
        error = function(e) list(
          ok = FALSE, error = conditionMessage(e),
          call = paste(deparse(conditionCall(e)), collapse = " ")))
      tryCatch({
        serialize(response, connection)
        flush(connection)
      }, error = function(e) NULL)
      if (identical(request$command, "shutdown")) {
        shutdown <- TRUE
        break
      }
    }
    close(connection)
  }
  invisible(NULL)
}

if (identical(worker, "producer")) {
  if (!dir.exists(work_path) || !file.exists(binary_path) ||
      is.na(spool_cap) || spool_cap < 1) stop("invalid producer worker setup")
  state <- new.env(parent = emptyenv())
  state$source <- file.path(work_path, "producer-source.b64")
  state$transfer <- NULL
  state$frame_chars <- NULL
  state$high_water <- 0
  state$reads <- 0L
  state$duplicates <- 0L
  state$generation_count <- 0L
  state$backpressure <- 0L
  state$baseline_rss <- rss_bytes()
  state$peak_rss <- state$baseline_rss
  state$go_peak_rss <- 0
  state$source_removed <- FALSE
  update_peak <- function() {
    state$peak_rss <- max(state$peak_rss, rss_bytes(), na.rm = TRUE)
  }
  transfer_for <- function(metadata) list(
    ticket = strrep("A", 86L),
    transfer_id = "tb_00112233445566778899aabbccddeeff",
    capability_id = "blob.transport.source-probe.v1",
    sender_name = "producer", recipient_name = "recipient",
    payload_chars = as.numeric(metadata$payload_chars),
    payload_sha256 = metadata$payload_sha256)
  generate <- function(raw_bytes) {
    if (!is.null(state$transfer)) return(state$transfer)
    expected <- ceiling(4 * raw_bytes / 3)
    if (expected > spool_cap) stop("producer spool backpressure")
    input <- jsonlite::toJSON(list(
      output_path = file.path(normalizePath(work_path, mustWork = TRUE),
                              "producer-source.b64"),
      raw_bytes = as.numeric(raw_bytes)), auto_unbox = TRUE, digits = NA)
    process <- processx::process$new(
      binary_path, "typed-source-stream-probe",
      stdin = "|", stdout = "|", stderr = "|")
    process$write_input(input)
    processx::processx_conn_close(process$get_input_connection())
    while (process$is_alive()) {
      memory <- tryCatch(process$get_memory_info()[["rss"]],
                         error = function(e) 0)
      state$go_peak_rss <- max(state$go_peak_rss, as.numeric(memory))
      process$poll_io(5L)
      update_peak()
    }
    process$wait()
    output <- process$read_all_output()
    errors <- process$read_all_error()
    if (!identical(process$get_exit_status(), 0L)) {
      stop("typed-source runtime failed: ", errors, call. = FALSE)
    }
    metadata <- jsonlite::fromJSON(output, simplifyVector = TRUE)
    if (!identical(metadata$version, "dsvert-typed-source-stream-v1") ||
        !identical(as.numeric(file.size(state$source)),
                   as.numeric(metadata$payload_chars)) ||
        !identical(digest::digest(
          file = state$source, algo = "sha256", serialize = FALSE),
          metadata$payload_sha256)) stop("producer integrity mismatch")
    state$generation_count <- state$generation_count + 1L
    state$transfer <- transfer_for(metadata)
    update_peak()
    state$transfer
  }
  handler <- function(request) {
    update_peak()
    if (identical(request$command, "preflight")) {
      rejected <- as.numeric(request$total_chars) > spool_cap
      if (rejected) state$backpressure <- state$backpressure + 1L
      return(list(rejected = rejected, source_created = file.exists(state$source)))
    }
    if (identical(request$command, "stats")) return(list(
      baseline_rss = state$baseline_rss, peak_rss = state$peak_rss,
      go_peak_rss = state$go_peak_rss, reads = state$reads,
      duplicates = state$duplicates,
      generation_count = state$generation_count,
      backpressure_rejections = state$backpressure,
      high_water = state$high_water, source_removed = state$source_removed,
      payload_sha256 = state$transfer$payload_sha256))
    if (identical(request$command, "shutdown")) return(TRUE)
    expr <- request$expr
    method <- as.character(expr[[1L]])
    if (identical(method, "mpcTypedSourceProbeDS")) {
      return(list(source_transfer = generate(as.numeric(expr$payload_bytes))))
    }
    if (is.null(state$transfer)) stop("producer has no transfer")
    if (identical(method, "mpcTypedBlobReceiptDS")) {
      if (!identical(state$high_water, state$transfer$payload_chars)) {
        stop("receipt preceded complete source read")
      }
      unlink(state$source)
      state$source_removed <- !file.exists(state$source)
      return(list(
        version = "dsvert-typed-blob-receipt-v1",
        transfer_id = state$transfer$transfer_id, confirmed = TRUE))
    }
    if (!identical(expr$ticket, state$transfer$ticket)) {
      stop("unknown producer ticket")
    }
    if (identical(method, "mpcTypedBlobReadDS")) {
      offset <- as.numeric(expr$offset)
      frame_chars <- as.numeric(expr$max_chars)
      if (is.null(state$frame_chars)) state$frame_chars <- frame_chars
      if (!identical(frame_chars, state$frame_chars) ||
          offset %% frame_chars != 0 || offset > state$high_water) {
        stop("invalid source geometry")
      }
      count <- min(frame_chars, state$transfer$payload_chars - offset)
      chunk <- read_chars_at(state$source, offset, count)
      state$reads <- state$reads + 1L
      if (offset < state$high_water) state$duplicates <- state$duplicates + 1L
      if (identical(offset, state$high_water)) {
        state$high_water <- state$high_water + count
      }
      update_peak()
      return(list(
        version = "dsvert-typed-blob-source-v1",
        transfer_id = state$transfer$transfer_id, offset = offset,
        chunk = chunk, chunk_chars = as.numeric(count),
        chunk_sha256 = paste0(openssl::sha256(charToRaw(chunk))),
        total_chars = state$transfer$payload_chars,
        payload_sha256 = state$transfer$payload_sha256,
        final = identical(offset + count, state$transfer$payload_chars)))
    }
    stop("unexpected producer command")
  }
  run_listener(handler)
  quit(save = "no", status = 0L)
}

if (identical(worker, "recipient")) {
  if (!dir.exists(work_path) || is.na(spool_cap) || spool_cap < 1) {
    stop("invalid recipient worker setup")
  }
  state <- new.env(parent = emptyenv())
  state$path <- file.path(work_path, "recipient-spool.b64")
  state$transfer <- NULL
  state$committed <- 0
  state$stores <- 0L
  state$duplicates <- 0L
  state$backpressure <- 0L
  state$baseline_rss <- rss_bytes()
  state$peak_rss <- state$baseline_rss
  update_peak <- function() {
    state$peak_rss <- max(state$peak_rss, rss_bytes(), na.rm = TRUE)
  }
  handler <- function(request) {
    update_peak()
    if (identical(request$command, "preflight")) {
      rejected <- as.numeric(request$total_chars) > spool_cap
      if (rejected) state$backpressure <- state$backpressure + 1L
      return(list(rejected = rejected, spool_created = file.exists(state$path)))
    }
    if (identical(request$command, "stats")) return(list(
      baseline_rss = state$baseline_rss, peak_rss = state$peak_rss,
      stores = state$stores, duplicates = state$duplicates,
      backpressure_rejections = state$backpressure,
      committed = state$committed,
      payload_sha256 = if (file.exists(state$path)) digest::digest(
        file = state$path, algo = "sha256", serialize = FALSE) else NA_character_))
    if (identical(request$command, "shutdown")) return(TRUE)
    if (!identical(request$command, "store")) stop("unexpected recipient command")
    expr <- request$expr
    transfer <- request$transfer
    if (is.null(state$transfer)) {
      if (transfer$payload_chars > spool_cap) stop("recipient spool backpressure")
      if (!file.create(state$path)) stop("could not create recipient spool")
      Sys.chmod(state$path, "0600")
      state$transfer <- transfer
    }
    if (!identical(expr$ticket, state$transfer$ticket) ||
        !identical(transfer, state$transfer)) stop("recipient ticket conflict")
    offset <- as.numeric(expr$offset)
    chunk <- as.character(expr$chunk)
    count <- nchar(chunk, type = "bytes")
    if (offset > state$committed) stop("recipient offset gap")
    if (offset < state$committed) {
      if (!identical(read_chars_at(state$path, offset, count), chunk)) {
        stop("conflicting recipient replay")
      }
      state$duplicates <- state$duplicates + 1L
    } else {
      con <- file(state$path, "ab")
      writeBin(charToRaw(chunk), con, useBytes = TRUE)
      flush(con)
      close(con)
      state$committed <- state$committed + count
    }
    state$stores <- state$stores + 1L
    sealed <- identical(state$committed, state$transfer$payload_chars)
    if (sealed && !identical(digest::digest(
      file = state$path, algo = "sha256", serialize = FALSE),
      state$transfer$payload_sha256)) stop("recipient terminal hash mismatch")
    update_peak()
    list(
      version = "dsvert-typed-blob-v1",
      transfer_id = state$transfer$transfer_id,
      committed_chars = state$committed,
      total_chars = state$transfer$payload_chars, sealed = sealed,
      receipt = if (sealed) "cmVjZWlwdA" else NULL)
  }
  run_listener(handler)
  quit(save = "no", status = 0L)
}

payload_mib <- suppressWarnings(as.numeric(arg_value("--payload-mib", "128")))
frame_kib <- suppressWarnings(as.numeric(arg_value("--frame-kib", "512")))
if (!is.finite(payload_mib) || payload_mib < 1 ||
    payload_mib != floor(payload_mib) ||
    !is.finite(frame_kib) || frame_kib < 16 || frame_kib > 8192 ||
    frame_kib != floor(frame_kib)) stop("invalid benchmark geometry")
payload_chars <- as.numeric(payload_mib * 1024^2)
if (payload_chars %% 4 != 0) stop("payload size must be divisible by four")
raw_bytes <- payload_chars * 3 / 4
frame_chars <- as.numeric(frame_kib * 1024)
if (!file.exists(binary_path)) stop("--binary must name a built dsvert-mpc")

old_options <- options(
  dsvert.chunk_size = frame_chars,
  dsvert.dsi_max_expression_bytes = frame_chars + 4096,
  dsvert.dsi.retry_deadline_seconds = 30)
on.exit(options(old_options), add = TRUE)
dsVertClient:::.dsvert_reset_chunk_size()
on.exit(dsVertClient:::.dsvert_reset_chunk_size(), add = TRUE)

root <- tempfile(pattern = "dsvert-source-stream-benchmark-")
dir.create(root, mode = "0700")
producer_work <- file.path(root, "producer")
recipient_work <- file.path(root, "recipient")
dir.create(producer_work, mode = "0700")
dir.create(recipient_work, mode = "0700")
on.exit(unlink(root, recursive = TRUE), add = TRUE)
producer_ready <- file.path(root, "producer.ready")
recipient_ready <- file.path(root, "recipient.ready")
ports <- c(httpuv::randomPort(), httpuv::randomPort())
while (ports[[2L]] == ports[[1L]]) ports[[2L]] <- httpuv::randomPort()
rscript <- file.path(R.home("bin"), "Rscript")
worker_process <- function(role, selected_port, ready, work) {
  processx::process$new(rscript, c(
    script_path, paste0("--worker=", role),
    paste0("--port=", selected_port), paste0("--ready=", ready),
    paste0("--work=", work), paste0("--binary=", binary_path),
    paste0("--spool-cap=", format(payload_chars, scientific = FALSE))),
    stdout = "|", stderr = "|")
}
producer_process <- worker_process(
  "producer", ports[[1L]], producer_ready, producer_work)
recipient_process <- worker_process(
  "recipient", ports[[2L]], recipient_ready, recipient_work)
on.exit({
  if (producer_process$is_alive()) producer_process$kill()
  if (recipient_process$is_alive()) recipient_process$kill()
}, add = TRUE)
deadline <- Sys.time() + 15
while ((!file.exists(producer_ready) || !file.exists(recipient_ready)) &&
       Sys.time() < deadline) Sys.sleep(0.02)
if (!file.exists(producer_ready) || !file.exists(recipient_ready)) {
  stop("workers did not become ready: ",
       producer_process$read_all_error(), " ",
       recipient_process$read_all_error(), call. = FALSE)
}

connections <- new.env(parent = emptyenv())
connect_one <- function(role) {
  selected <- if (identical(role, "producer")) ports[[1L]] else ports[[2L]]
  deadline <- Sys.time() + 10
  repeat {
    value <- tryCatch(socketConnection(
      host = "127.0.0.1", port = selected, blocking = TRUE,
      open = "a+b", timeout = 30), error = function(e) NULL)
    if (!is.null(value)) return(value)
    if (Sys.time() >= deadline) stop("worker reconnect deadline exceeded")
    Sys.sleep(0.02)
  }
}
connections$producer <- connect_one("producer")
connections$recipient <- connect_one("recipient")
rpc <- function(role, request) {
  connection <- connections[[role]]
  serialize(request, connection)
  flush(connection)
  response <- unserialize(connection)
  if (!is.list(response) || !identical(response$ok, TRUE)) {
    detail <- if (is.list(response) && is.character(response$error)) {
      paste(response$error, if (is.character(response$call)) response$call else "")
    } else {
      "malformed worker response"
    }
    message("BENCHMARK_WORKER_REJECTION[", role, "]: ", detail)
    stop("worker rejected request", call. = FALSE)
  }
  response$value
}
on.exit({
  for (role in c("producer", "recipient")) {
    try(rpc(role, list(command = "shutdown")), silent = TRUE)
    try(close(connections[[role]]), silent = TRUE)
  }
}, add = TRUE)

producer_pressure <- rpc("producer", list(
  command = "preflight", total_chars = payload_chars + 1))
recipient_pressure <- rpc("recipient", list(
  command = "preflight", total_chars = payload_chars + 1))
if (!identical(producer_pressure, list(rejected = TRUE, source_created = FALSE)) ||
    !identical(recipient_pressure, list(rejected = TRUE, spool_created = FALSE))) {
  stop("backpressure preflight created state", call. = FALSE)
}

controller_baseline <- rss_bytes()
controller_peak <- controller_baseline
update_controller_peak <- function() {
  controller_peak <<- max(controller_peak, rss_bytes(), na.rm = TRUE)
}
generation_loss <- TRUE
source_loss <- TRUE
sink_loss <- TRUE
reconnected <- FALSE
reconnect_count <- 0L
transfer_seen <- NULL
aggregate <- function(conns, expr, async, error, errors.print) {
  update_controller_peak()
  expressions <- if (is.list(expr) && !is.call(expr)) {
    expr[names(conns)]
  } else {
    stats::setNames(rep(list(expr), length(conns)), names(conns))
  }
  methods <- vapply(expressions, function(value) {
    as.character(value[[1L]])[[1L]]
  }, character(1L))
  reconnect_ready <- vapply(seq_along(expressions), function(index) {
    method <- methods[[index]]
    value <- expressions[[index]]
    method %in% c("mpcTypedBlobReadDS", "mpcTypedBlobStoreDS") &&
      as.numeric(value$offset) >= payload_chars / 2
  }, logical(1L))
  if (!reconnected && any(reconnect_ready)) {
    close(connections$producer)
    close(connections$recipient)
    connections$producer <- connect_one("producer")
    connections$recipient <- connect_one("recipient")
    reconnected <<- TRUE
    reconnect_count <<- reconnect_count + 1L
  }
  values <- lapply(seq_along(expressions), function(index) {
    method <- methods[[index]]
    value_expr <- expressions[[index]]
    if (identical(method, "mpcTypedSourceProbeDS")) {
      value <- rpc("producer", list(command = "call", expr = value_expr))
      if (generation_loss) {
        generation_loss <<- FALSE
        return(NULL)
      }
      transfer_seen <<- value$source_transfer
      update_controller_peak()
      return(value)
    }
    if (identical(method, "mpcTypedBlobReadDS")) {
      value <- rpc("producer", list(command = "call", expr = value_expr))
      if (source_loss) {
        source_loss <<- FALSE
        return(NULL)
      }
      update_controller_peak()
      return(value)
    }
    if (identical(method, "mpcTypedBlobStoreDS")) {
      value <- rpc("recipient", list(
        command = "store", expr = value_expr, transfer = transfer_seen))
      if (sink_loss) {
        sink_loss <<- FALSE
        return(NULL)
      }
      update_controller_peak()
      return(value)
    }
    if (identical(method, "mpcTypedBlobReceiptDS")) {
      value <- rpc("producer", list(command = "call", expr = value_expr))
      update_controller_peak()
      return(value)
    }
    stop("unexpected benchmark aggregate method", call. = FALSE)
  })
  stats::setNames(values, names(conns))
}

started <- proc.time()[["elapsed"]]
result <- dsVertClient:::.dsvert_run_typed_source_probe(
  producer_conn = list(producer = structure(list(), class = "benchmark")),
  recipient_conn = list(recipient = structure(list(), class = "benchmark")),
  recipient_pk = "cGlubmVkLXJlY2lwaWVudC10cmFuc3BvcnQta2V5",
  payload_bytes = raw_bytes,
  session_id = "12345678-1234-4234-8234-123456789abc",
  .aggregate = aggregate)
wall_seconds <- unname(proc.time()[["elapsed"]] - started)
update_controller_peak()
producer_stats <- rpc("producer", list(command = "stats"))
recipient_stats <- rpc("recipient", list(command = "stats"))
expected_frames <- ceiling(payload_chars / frame_chars)
if (!identical(result$frames, as.integer(expected_frames)) ||
    !identical(producer_stats$payload_sha256,
               recipient_stats$payload_sha256) ||
    !identical(as.numeric(recipient_stats$committed), payload_chars) ||
    !identical(producer_stats$generation_count, 1L) ||
    producer_stats$duplicates < 1L || recipient_stats$duplicates < 1L ||
    !identical(reconnect_count, 1L) ||
    !identical(producer_stats$source_removed, TRUE) ||
    producer_stats$backpressure_rejections < 1L ||
    recipient_stats$backpressure_rejections < 1L) {
  stop("multiprocess source-stream invariants failed", call. = FALSE)
}

record <- list(
  payload_mib = payload_mib, payload_chars = payload_chars,
  raw_source_mib = raw_bytes / 1024^2,
  frame_kib = frame_kib, frames = as.integer(result$frames),
  wall_seconds = wall_seconds,
  throughput_mib_s = payload_mib / wall_seconds,
  controller_baseline_rss_bytes = controller_baseline,
  controller_peak_rss_bytes = controller_peak,
  controller_rss_delta_bytes = controller_peak - controller_baseline,
  controller_payload_spool_bytes = 0,
  redundant_integrity_io_bytes_avoided = 2 * payload_chars,
  producer_baseline_rss_bytes = producer_stats$baseline_rss,
  producer_peak_rss_bytes = producer_stats$peak_rss,
  producer_rss_delta_bytes = producer_stats$peak_rss - producer_stats$baseline_rss,
  go_producer_peak_rss_bytes = producer_stats$go_peak_rss,
  recipient_baseline_rss_bytes = recipient_stats$baseline_rss,
  recipient_peak_rss_bytes = recipient_stats$peak_rss,
  recipient_rss_delta_bytes = recipient_stats$peak_rss - recipient_stats$baseline_rss,
  producer_generation_count = producer_stats$generation_count,
  producer_frame_replays = producer_stats$duplicates,
  recipient_frame_replays = recipient_stats$duplicates,
  reconnect_count = reconnect_count,
  producer_backpressure_rejections = producer_stats$backpressure_rejections,
  recipient_backpressure_rejections = recipient_stats$backpressure_rejections,
  source_removed_after_receipt = producer_stats$source_removed,
  integrity = "sha256_identical",
  scope = "data-free-probe-only")
cat("DSVERT_TYPED_SOURCE_RESULT ", jsonlite::toJSON(
  record, auto_unbox = TRUE, digits = NA), "\n", sep = "")
