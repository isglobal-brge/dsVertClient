#!/usr/bin/env Rscript

# Reproducible in-process DSLite benchmark for dsVert's direct DSI fan-out.
# It registers a test-only method dynamically; no benchmark endpoint is added
# to dsVert's remote method surface.

args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) {
  cat(paste(
    "Usage: benchmark_dsi_transport.R [--quick] [--peers=N,...]",
    "[--sizes-mib=N,...] [--chunk-raw-kib=N,...]",
    "[--engines=direct,dsi-standard] [--repeats=N] [--balanced-ab]",
    "[--warmup] [--continue-on-error] [--output=PATH] [--preflight]\n"
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

suppressPackageStartupMessages({
  library(DSI)
  library(DSLite)
  library(digest)
  library(jsonlite)
})

quick <- "--quick" %in% args
continue_on_error <- "--continue-on-error" %in% args
balanced_ab <- "--balanced-ab" %in% args
warmup <- "--warmup" %in% args
value_arg <- function(prefix, default) {
  hit <- args[startsWith(args, paste0(prefix, "="))]
  if (!length(hit)) return(default)
  strsplit(sub(paste0("^", prefix, "="), "", hit[[1L]]), ",",
           fixed = TRUE)[[1L]]
}
peer_grid <- as.integer(value_arg("--peers", if (quick) "2" else c("2", "4", "8")))
size_grid_mib <- as.numeric(value_arg(
  "--sizes-mib", if (quick) "1" else c("1", "16", "64")))
chunk_grid_kib <- as.integer(value_arg(
  "--chunk-raw-kib", if (quick) "480" else "480"))
engine_grid <- value_arg("--engines", "direct")
repeats <- suppressWarnings(as.integer(value_arg("--repeats", "1")))
output_arg <- args[startsWith(args, "--output=")]
output_path <- if (length(output_arg)) sub("^--output=", "", output_arg[[1L]]) else ""
if (anyNA(peer_grid) || any(peer_grid < 1L) || anyNA(size_grid_mib) ||
    any(size_grid_mib <= 0) || anyNA(chunk_grid_kib) ||
    any(chunk_grid_kib < 16L) || any(chunk_grid_kib > 64L * 1024L) ||
    any(!engine_grid %in% c("direct", "dsi-standard")) ||
    length(repeats) != 1L || is.na(repeats) || repeats < 1L ||
    (balanced_ab &&
      (!setequal(engine_grid, c("direct", "dsi-standard")) ||
       length(engine_grid) != 2L || anyDuplicated(engine_grid)))) {
  stop("Invalid benchmark grid", call. = FALSE)
}

b64url_encode <- function(value) {
  encoded <- gsub("[\r\n]", "", jsonlite::base64_enc(value))
  chartr("+/", "-_", sub("=+$", "", encoded))
}

b64url_decode <- function(value) {
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !grepl("^[A-Za-z0-9_-]+$", value) || nchar(value) %% 4L == 1L) {
    stop("invalid Base64url", call. = FALSE)
  }
  standard <- chartr("-_", "+/", value)
  padding <- nchar(standard) %% 4L
  if (padding == 2L) standard <- paste0(standard, "==")
  if (padding == 3L) standard <- paste0(standard, "=")
  decoded <- jsonlite::base64_dec(standard)
  if (!identical(b64url_encode(decoded), value)) {
    stop("non-canonical Base64url", call. = FALSE)
  }
  decoded
}

rss_mib <- function() {
  output <- suppressWarnings(system2(
    "ps", c("-o", "rss=", "-p", as.character(Sys.getpid())),
    stdout = TRUE, stderr = FALSE))
  value <- suppressWarnings(as.numeric(trimws(output[[1L]])))
  if (!length(value) || !is.finite(value)) return(NA_real_)
  value / 1024
}

write_source <- function(path, total_bytes, block_bytes = 1024^2) {
  con <- file(path, "wb")
  on.exit(close(con), add = TRUE)
  offset <- 0
  while (offset < total_bytes) {
    n <- min(block_bytes, total_bytes - offset)
    value <- as.raw(((offset + seq_len(n) - 1) %% 251))
    writeBin(value, con)
    offset <- offset + n
  }
  flush(con)
  invisible(path)
}

read_at <- function(path, offset, n) {
  con <- file(path, "rb")
  on.exit(close(con), add = TRUE)
  seek(con, where = offset, origin = "start")
  value <- readBin(con, "raw", n = n)
  if (length(value) != n) stop("benchmark source read truncated", call. = FALSE)
  value
}

new_peer_method <- function(spool_dir, spool_cap_bytes, state_id) {
  state <- new.env(parent = emptyenv())
  state$operations <- list()
  state$spool_dir <- spool_dir
  state$spool_cap_bytes <- spool_cap_bytes
  registry_name <- ".dsvert_dsi_benchmark_registry"
  if (!exists(registry_name, envir = .GlobalEnv, inherits = FALSE)) {
    assign(registry_name, new.env(parent = emptyenv()), envir = .GlobalEnv)
  }
  registry <- get(registry_name, envir = .GlobalEnv, inherits = FALSE)
  registry[[state_id]] <- state
  method <- function(state_id, operation_id, payload, offset, total_bytes,
                     payload_hash, final) {
    registry <- get(".dsvert_dsi_benchmark_registry", envir = .GlobalEnv,
                    inherits = FALSE)
    state <- registry[[state_id]]
    if (!is.environment(state)) stop("unknown benchmark state", call. = FALSE)
    spool_dir <- state$spool_dir
    spool_cap_bytes <- state$spool_cap_bytes
    if (!is.character(operation_id) || length(operation_id) != 1L ||
        !grepl("^bench_[0-9a-f]{32}$", operation_id)) {
      stop("invalid operation", call. = FALSE)
    }
    offset <- as.numeric(offset)
    total_bytes <- as.numeric(total_bytes)
    if (!is.finite(offset) || !is.finite(total_bytes) || offset < 0 ||
        total_bytes < 1 || offset != floor(offset) ||
        total_bytes != floor(total_bytes)) stop("invalid geometry", call. = FALSE)
    if (!is.character(payload_hash) || length(payload_hash) != 1L ||
        !grepl("^[0-9a-f]{64}$", payload_hash)) stop("invalid hash", call. = FALSE)
    if (!is.logical(final) || length(final) != 1L || is.na(final)) {
      stop("invalid final flag", call. = FALSE)
    }
    encoded_bytes <- nchar(payload, type = "bytes")
    chunk <- b64url_decode(payload)
    chunk_bytes <- length(chunk)
    if (encoded_bytes != ceiling(4 * chunk_bytes / 3) || chunk_bytes < 1L ||
        offset + chunk_bytes > total_bytes ||
        !identical(isTRUE(final), offset + chunk_bytes == total_bytes)) {
      stop("invalid chunk shape", call. = FALSE)
    }

    operation <- state$operations[[operation_id]]
    if (is.null(operation)) {
      retained <- sum(vapply(state$operations, function(value)
        value$total_bytes, numeric(1L)))
      if (retained + total_bytes > spool_cap_bytes) {
        stop("benchmark spool backpressure", call. = FALSE)
      }
      path <- file.path(spool_dir, paste0(operation_id, ".bin"))
      file.create(path)
      Sys.chmod(path, mode = "0600")
      operation <- list(
        path = path, total_bytes = total_bytes, payload_hash = payload_hash,
        next_offset = 0, complete = FALSE)
    } else if (!identical(operation$total_bytes, total_bytes) ||
               !identical(operation$payload_hash, payload_hash)) {
      stop("conflicting operation replay", call. = FALSE)
    }

    if (offset < operation$next_offset) {
      existing <- read_at(operation$path, offset, chunk_bytes)
      if (!identical(existing, chunk)) stop("conflicting chunk replay", call. = FALSE)
    } else {
      if (offset != operation$next_offset) stop("offset gap", call. = FALSE)
      con <- file(operation$path, "ab")
      writeBin(chunk, con)
      flush(con)
      close(con)
      operation$next_offset <- offset + chunk_bytes
    }
    if (isTRUE(final)) {
      if (file.size(operation$path) != total_bytes ||
          !identical(digest::digest(
            file = operation$path, algo = "sha256", serialize = FALSE),
            payload_hash)) stop("terminal hash mismatch", call. = FALSE)
      operation$complete <- TRUE
    }
    state$operations[[operation_id]] <- operation
    list(
      version = "dsvert-dsi-benchmark-v1", operation_id = operation_id,
      ack_offset = operation$next_offset, complete = operation$complete,
      payload_hash = operation$payload_hash)
  }
  list(method = method, state = state)
}

connect_peers <- function(server_names) {
  builder <- DSI::newDSLoginBuilder()
  for (site in names(server_names)) {
    builder$append(
      server = site, url = server_names[[site]], table = "t",
      driver = "DSLiteDriver")
  }
  DSI::datashield.login(builder$build(), assign = FALSE)
}

validate_ack <- function(value, operation_id, expected_offset,
                         expected_complete, payload_hash) {
  required <- c("version", "operation_id", "ack_offset", "complete",
                "payload_hash")
  if (!is.list(value) || !identical(sort(names(value)), sort(required)) ||
      !identical(value$version, "dsvert-dsi-benchmark-v1") ||
      !identical(value$operation_id, operation_id) ||
      !identical(as.numeric(value$ack_offset), as.numeric(expected_offset)) ||
      !identical(value$complete, expected_complete) ||
      !identical(value$payload_hash, payload_hash)) {
    stop("malformed benchmark acknowledgement", call. = FALSE)
  }
  TRUE
}

run_case <- function(n_peers, size_mib, chunk_bytes = 480 * 1024,
                     dispatch_engine = "direct") {
  total_bytes <- as.numeric(size_mib * 1024^2)
  if (total_bytes != floor(total_bytes)) stop("size must resolve to bytes")
  # The benchmark deliberately evaluates candidate geometries above the
  # package's portable floor.  It raises only the local pre-parse guard, to a
  # fixed value known before transmission; the real connector/parser still has
  # to accept the data-free frame.  Production code must negotiate this bound.
  expression_cap <- ceiling(4 * chunk_bytes / 3) + 64 * 1024
  old_options <- options(dsvert.dsi_max_expression_bytes = expression_cap)
  on.exit(options(old_options), add = TRUE)
  case_dir <- tempfile(pattern = "dsvert-dsi-benchmark-")
  dir.create(case_dir, mode = "0700")
  on.exit(unlink(case_dir, recursive = TRUE), add = TRUE)
  source_path <- file.path(case_dir, "source.bin")
  write_source(source_path, total_bytes)
  payload_hash <- digest::digest(
    file = source_path, algo = "sha256", serialize = FALSE)

  peers <- sprintf("site%02d", seq_len(n_peers))
  server_objects <- list()
  server_names <- character()
  states <- list()
  state_ids <- character()
  chunk_kib <- as.integer(chunk_bytes / 1024)
  method_name <- paste0("dsvertTransportBenchmarkDS", Sys.getpid(),
                        sprintf("%02d%03d%05d", n_peers,
                                as.integer(size_mib), chunk_kib))
  for (index in seq_along(peers)) {
    site <- peers[[index]]
    spool <- file.path(case_dir, paste0("spool-", site))
    dir.create(spool, mode = "0700")
    state_id <- paste(
      Sys.getpid(), n_peers, as.integer(size_mib), chunk_kib, index,
      sep = "_")
    implementation <- new_peer_method(
      spool, total_bytes + chunk_bytes, state_id)
    server <- DSLite::newDSLiteServer(tables = list(t = data.frame(x = 1)))
    server$aggregateMethod(method_name, implementation$method)
    object_name <- paste0("dsvert_dsi_bench_", Sys.getpid(), "_", index)
    assign(object_name, server, envir = .GlobalEnv)
    server_objects[[site]] <- object_name
    server_names[[site]] <- object_name
    states[[site]] <- implementation$state
    state_ids[[site]] <- state_id
  }
  on.exit(rm(list = unname(unlist(server_objects, use.names = FALSE)),
             envir = .GlobalEnv), add = TRUE)
  on.exit(rm(
    list = unname(state_ids),
    envir = get(".dsvert_dsi_benchmark_registry", envir = .GlobalEnv,
                inherits = FALSE)), add = TRUE)
  conns <- connect_peers(server_names)
  on.exit(try(DSI::datashield.logout(conns), silent = TRUE), add = TRUE)

  operation_id <- paste0("bench_", substr(digest::digest(
    paste(n_peers, size_mib, chunk_kib, payload_hash, sep = "|"),
    algo = "sha256", serialize = FALSE), 1L, 32L))
  offsets <- stats::setNames(rep(0, n_peers), peers)
  calls <- 0L
  replay_cycles <- 0L
  raw_submitted <- 0
  encoded_submitted <- 0
  expression_submitted <- 0
  max_expression <- 0
  lost_injected <- FALSE
  reconnected <- FALSE
  rss_start <- rss_mib()
  rss_peak <- rss_start
  spool_peak <- 0
  started <- proc.time()[[3L]]
  aggregate_engine <- if (identical(dispatch_engine, "direct")) {
    DSI::datashield.aggregate
  } else {
    # A non-identical wrapper deliberately exercises DSI's standard aggregate
    # orchestration for an A/B comparison with dsVert's immutable-frame path.
    function(...) DSI::datashield.aggregate(...)
  }

  while (any(offsets < total_bytes)) {
    expressions <- stats::setNames(lapply(peers, function(site) {
      offset <- offsets[[site]]
      n <- min(chunk_bytes, total_bytes - offset)
      chunk <- read_at(source_path, offset, n)
      encoded <- b64url_encode(chunk)
      call(
        name = method_name, state_id = state_ids[[site]],
        operation_id = operation_id, payload = encoded,
        offset = offset, total_bytes = total_bytes,
        payload_hash = payload_hash, final = offset + n == total_bytes)
    }), peers)
    sizes <- dsVertClient:::.dsvert_validate_dsi_expression_sizes(expressions)
    cycle <- dsVertClient:::.dsvert_fanout_cycle(
      conns, expressions, operation = "DSI transport benchmark",
      .aggregate = aggregate_engine)
    calls <- calls + 1L
    per_site_raw <- vapply(peers, function(site)
      min(chunk_bytes, total_bytes - offsets[[site]]), numeric(1L))
    per_site_encoded <- vapply(expressions, function(expr)
      nchar(expr[["payload"]], type = "bytes"), numeric(1L))
    raw_submitted <- raw_submitted + sum(per_site_raw)
    encoded_submitted <- encoded_submitted + sum(per_site_encoded)
    expression_submitted <- expression_submitted + sum(sizes)
    max_expression <- max(max_expression, sizes)
    if (!identical(cycle$state, "ok")) {
      stop("DSLite benchmark unexpectedly unavailable: ",
           paste(DSI::datashield.errors(), collapse = " | "), call. = FALSE)
    }

    # Simulate one lost complete response. The next loop replays byte-identical
    # expressions and absolute offsets; no acknowledgement is advanced.
    if (!lost_injected) {
      lost_injected <- TRUE
      replay_cycles <- replay_cycles + 1L
    } else {
      for (site in peers) {
        next_offset <- offsets[[site]] + per_site_raw[[site]]
        validate_ack(
          cycle$result[[site]], operation_id, next_offset,
          next_offset == total_bytes, payload_hash)
        offsets[[site]] <- next_offset
      }
    }

    if (!reconnected && all(offsets >= total_bytes / 2)) {
      DSI::datashield.logout(conns)
      conns <- connect_peers(server_names)
      reconnected <- TRUE
    }
    current_spool <- sum(vapply(states, function(state) {
      if (!length(state$operations)) return(0)
      sum(vapply(state$operations, function(operation)
        as.numeric(file.size(operation$path)), numeric(1L)))
    }, numeric(1L)))
    spool_peak <- max(spool_peak, current_spool)
    rss_peak <- max(rss_peak, rss_mib(), na.rm = TRUE)
  }
  wall <- proc.time()[[3L]] - started
  complete <- vapply(states, function(state) {
    operation <- state$operations[[operation_id]]
    isTRUE(operation$complete) && identical(digest::digest(
      file = operation$path, algo = "sha256", serialize = FALSE), payload_hash)
  }, logical(1L))
  if (!all(complete)) stop("benchmark byte-integrity failure", call. = FALSE)

  data.frame(
    dispatch_engine = dispatch_engine,
    peers = n_peers, payload_mib_per_peer = size_mib,
    unique_raw_mib = n_peers * size_mib,
    chunk_raw_kib = chunk_bytes / 1024,
    fanout_cycles = calls, replay_cycles = replay_cycles,
    reconnects = as.integer(reconnected), wall_seconds = wall,
    throughput_unique_mib_s = n_peers * size_mib / wall,
    raw_submitted_mib = raw_submitted / 1024^2,
    base64_submitted_mib = encoded_submitted / 1024^2,
    expression_submitted_mib = expression_submitted / 1024^2,
    base64_expansion = encoded_submitted / raw_submitted,
    max_expression_kib = max_expression / 1024,
    rss_start_mib = rss_start, rss_peak_mib = rss_peak,
    rss_delta_mib = rss_peak - rss_start,
    spool_peak_mib = spool_peak / 1024^2,
    source_disk_mib = as.numeric(file.size(source_path)) / 1024^2,
    integrity = "sha256_identical", stringsAsFactors = FALSE)
}

failed_case <- function(n_peers, size_mib, chunk_kib, dispatch_engine,
                        wall_seconds) {
  data.frame(
    dispatch_engine = dispatch_engine,
    peers = n_peers, payload_mib_per_peer = size_mib,
    unique_raw_mib = n_peers * size_mib,
    chunk_raw_kib = chunk_kib,
    fanout_cycles = NA_integer_, replay_cycles = NA_integer_,
    reconnects = NA_integer_, wall_seconds = wall_seconds,
    throughput_unique_mib_s = NA_real_, raw_submitted_mib = NA_real_,
    base64_submitted_mib = NA_real_, expression_submitted_mib = NA_real_,
    base64_expansion = NA_real_, max_expression_kib =
      ceiling(4 * chunk_kib / 3) + 64,
    rss_start_mib = NA_real_, rss_peak_mib = NA_real_,
    rss_delta_mib = NA_real_, spool_peak_mib = NA_real_,
    source_disk_mib = size_mib, integrity = "transport_rejected",
    stringsAsFactors = FALSE)
}

results <- list()
for (size in size_grid_mib) {
  for (peers in peer_grid) {
    for (chunk_kib in chunk_grid_kib) {
      if (isTRUE(warmup)) {
        for (engine in unique(engine_grid)) {
          message(sprintf(
            "Warming %s: %d peers x %g MiB, %d KiB raw frames ...",
            engine, peers, size, chunk_kib))
          invisible(run_case(
            peers, size, chunk_bytes = chunk_kib * 1024,
            dispatch_engine = engine))
          gc()
        }
      }
      for (replicate_index in seq_len(repeats)) {
        engines <- if (isTRUE(balanced_ab) && replicate_index %% 2L == 0L) {
          rev(engine_grid)
        } else {
          engine_grid
        }
        for (order_position in seq_along(engines)) {
          engine <- engines[[order_position]]
          message(sprintf(
            paste0(
              "Benchmarking %s: %d peers x %g MiB, %d KiB raw frames ",
              "(replicate %d/%d, order %d) ..."),
            engine, peers, size, chunk_kib, replicate_index, repeats,
            order_position))
          case_started <- proc.time()[[3L]]
          result <- tryCatch(
            run_case(
              peers, size, chunk_bytes = chunk_kib * 1024,
              dispatch_engine = engine),
            error = function(e) {
              if (!isTRUE(continue_on_error)) stop(e)
              message("Candidate geometry was rejected; continuing the sweep.")
              failed_case(
                peers, size, chunk_kib, engine,
                proc.time()[[3L]] - case_started)
            })
          if (isTRUE(balanced_ab) || repeats > 1L) {
            result <- cbind(
              result["dispatch_engine"],
              replicate = as.integer(replicate_index),
              order_position = as.integer(order_position),
              result[setdiff(names(result), "dispatch_engine")],
              stringsAsFactors = FALSE)
          }
          results[[length(results) + 1L]] <- result
          print(results[[length(results)]])
          if (nzchar(output_path)) {
            dir.create(
              dirname(output_path), recursive = TRUE, showWarnings = FALSE)
            checkpoint <- do.call(rbind, results)
            rownames(checkpoint) <- NULL
            utils::write.csv(checkpoint, output_path, row.names = FALSE)
          }
          gc()
        }
      }
    }
  }
}
results <- do.call(rbind, results)
rownames(results) <- NULL
print(results)
if (nzchar(output_path)) {
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(results, output_path, row.names = FALSE)
  message("Wrote ", normalizePath(output_path, mustWork = TRUE))
}
