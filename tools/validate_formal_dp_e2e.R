#!/usr/bin/env Rscript

# Source-tree-only release validation of dsVert's canonical sticky Count.
# Validation dependencies: callr, digest, DSI, DSLite, jsonlite, methods,
# openssl and pkgload. This script is not an installed runtime entry point.
#
# This harness deliberately uses only exported analyst APIs and documented
# custodian boundaries.  Every peer is a separate persistent R process with an
# exact DSLite allowlist. Persistent identities bootstrap automatically. There
# is no client-supplied/test bypass, analyst secret, plaintext exception or raw
# remote route: the harness acts as custodian administrator and derives the
# server-side surface token only after verifying the exact effective inventory.

args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) {
  cat(paste(
    "Source-tree release harness; run from a dsVert/dsVertClient checkout.\n",
    "Usage: validate_formal_dp_e2e.R [--k=2,3,5] [--output=PATH]",
    "[--keep-state=PATH]\n"
  ))
  quit(save = "no", status = 0L)
}

suppressPackageStartupMessages({
  library(callr)
  library(digest)
  library(DSI)
  library(DSLite)
  library(jsonlite)
  library(methods)
  library(pkgload)
})

.e2e_arg <- function(prefix, default = NULL) {
  value <- args[startsWith(args, paste0(prefix, "="))]
  if (!length(value)) return(default)
  sub(paste0("^", prefix, "="), "", value[[1L]])
}

.e2e_k <- as.integer(strsplit(.e2e_arg("--k", "2,3"), ",", fixed = TRUE)[[1L]])
.e2e_output <- .e2e_arg("--output", "")
.e2e_keep_path_arg <- .e2e_arg("--keep-state", "")
.e2e_cleanup_test <- .e2e_arg("--cleanup-self-test", "")
.e2e_cleanup_marker <- .e2e_arg("--cleanup-marker", "")
if ("--keep-state" %in% args) {
  stop("--keep-state requires an explicit =PATH destination",
       call. = FALSE)
}
if (anyNA(.e2e_k) || any(!.e2e_k %in% c(2L, 3L, 5L)) ||
    anyDuplicated(.e2e_k)) {
  stop("--k must be a comma-separated subset of 2,3,5", call. = FALSE)
}
if (!.e2e_cleanup_test %in% c("", "success", "failure")) {
  stop("Invalid cleanup self-test mode", call. = FALSE)
}

.e2e_script_path <- function() {
  file_arg <- commandArgs(trailingOnly = FALSE)
  file_arg <- file_arg[startsWith(file_arg, "--file=")]
  if (!length(file_arg)) return(NULL)
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
}

.e2e_repo_root <- function() {
  starts <- c(getwd(), if (!is.null(.e2e_script_path())) {
    dirname(.e2e_script_path())
  } else character())
  for (start in starts) {
    current <- normalizePath(start, mustWork = TRUE)
    repeat {
      if (file.exists(file.path(current, "dsVert", "DESCRIPTION")) &&
          file.exists(file.path(current, "dsVertClient", "DESCRIPTION"))) {
        return(current)
      }
      parent <- dirname(current)
      if (identical(parent, current)) break
      current <- parent
    }
  }
  stop("Cannot locate sibling dsVert and dsVertClient source trees",
       call. = FALSE)
}

.e2e_root <- .e2e_repo_root()
.e2e_server_dir <- normalizePath(file.path(.e2e_root, "dsVert"),
                                 mustWork = TRUE)
.e2e_client_dir <- normalizePath(file.path(.e2e_root, "dsVertClient"),
                                 mustWork = TRUE)

.e2e_mpc_binary <- function() {
  system <- Sys.info()[["sysname"]]
  machine <- tolower(Sys.info()[["machine"]])
  platform <- if (identical(system, "Darwin") &&
      machine %in% c("arm64", "aarch64")) {
    "darwin-arm64"
  } else if (identical(system, "Darwin") &&
             machine %in% c("x86_64", "amd64")) {
    "darwin-amd64"
  } else if (identical(system, "Linux") &&
             machine %in% c("x86_64", "amd64")) {
    "linux-amd64"
  } else if (identical(.Platform$OS.type, "windows") &&
             machine %in% c("x86_64", "amd64")) {
    "windows-amd64"
  } else {
    stop("No packaged dsvert-mpc binary for this validation platform",
         call. = FALSE)
  }
  executable <- if (identical(.Platform$OS.type, "windows")) {
    "dsvert-mpc.exe"
  } else {
    "dsvert-mpc"
  }
  path <- normalizePath(
    file.path(.e2e_server_dir, "inst", "bin", platform, executable),
    mustWork = TRUE)
  if (.Platform$OS.type != "windows") {
    if (file.access(path, mode = 1L) != 0L) {
      stop("The packaged dsvert-mpc binary is not executable", call. = FALSE)
    }
  }
  path
}

.e2e_real_mpc <- .e2e_mpc_binary()
.e2e_real_mpc_sha256 <- digest::digest(
  file = .e2e_real_mpc, algo = "sha256", serialize = FALSE)
.e2e_real_mpc_mode <- if (.Platform$OS.type == "windows") NA_integer_ else
  as.integer(file.info(.e2e_real_mpc)$mode)

.e2e_keep_path <- if (!nzchar(.e2e_keep_path_arg)) NULL else {
  requested <- path.expand(.e2e_keep_path_arg)
  parent <- normalizePath(dirname(requested), mustWork = TRUE)
  base <- basename(requested)
  requested_link <- Sys.readlink(requested)
  if (!nzchar(base) || base %in% c(".", "..") || file.exists(requested) ||
      (!is.na(requested_link) && nzchar(requested_link))) {
    stop("--keep-state destination must not already exist", call. = FALSE)
  }
  file.path(parent, base)
}

.e2e_state_parent <- normalizePath(tempdir(), mustWork = TRUE)
.e2e_state_root <- file.path(
  .e2e_state_parent,
  paste0(".dsvert-formal-dp-e2e-", Sys.getpid(), "-",
         paste(sprintf("%02x", as.integer(openssl::rand_bytes(8L))),
               collapse = "")))
dir.create(.e2e_state_root, recursive = FALSE, mode = "0700")
Sys.chmod(.e2e_state_root, mode = "0700")

.e2e_state_guard <- new.env(parent = emptyenv())
.e2e_state_guard$runtime <- .e2e_state_root
.e2e_state_guard$runtime_parent <- .e2e_state_parent
.e2e_state_guard$retain <- .e2e_keep_path
.e2e_state_guard$finalized <- FALSE
.e2e_state_guard$retained <- NULL

.e2e_finish_state <- function(guard, quiet = FALSE) {
  if (isTRUE(guard$finalized)) return(guard$retained)
  runtime <- guard$runtime
  valid_name <- grepl(
    "^\\.dsvert-formal-dp-e2e-[0-9]+-[0-9a-f]{16}$",
    basename(runtime))
  valid_parent <- identical(
    normalizePath(dirname(runtime), mustWork = TRUE),
    guard$runtime_parent)
  if (!isTRUE(valid_name) || !isTRUE(valid_parent)) {
    if (isTRUE(quiet)) return(NULL)
    stop("Refusing an unsafe harness state cleanup target", call. = FALSE)
  }
  retained <- NULL
  if (dir.exists(runtime)) {
    if (is.null(guard$retain)) {
      status <- unlink(runtime, recursive = TRUE, force = FALSE)
      if (!identical(as.integer(status), 0L) || dir.exists(runtime)) {
        if (isTRUE(quiet)) return(NULL)
        stop("Could not remove the ephemeral harness state", call. = FALSE)
      }
    } else {
      retain_link <- Sys.readlink(guard$retain)
      if (file.exists(guard$retain) ||
          (!is.na(retain_link) && nzchar(retain_link)) ||
          !isTRUE(file.rename(runtime, guard$retain))) {
        # Never leave an untracked secret-bearing runtime directory merely
        # because an explicitly requested diagnostic move failed.
        unlink(runtime, recursive = TRUE, force = FALSE)
        if (isTRUE(quiet)) return(NULL)
        stop("Could not retain harness state at the requested destination",
             call. = FALSE)
      }
      Sys.chmod(guard$retain, mode = "0700")
      retained <- guard$retain
    }
  }
  guard$retained <- retained
  guard$finalized <- TRUE
  retained
}

invisible(reg.finalizer(.e2e_state_guard, function(guard) {
  try(.e2e_finish_state(guard, quiet = TRUE), silent = TRUE)
}, onexit = TRUE))

if (nzchar(.e2e_cleanup_test)) {
  if (!nzchar(.e2e_cleanup_marker)) {
    stop("Cleanup self-test requires an explicit marker path", call. = FALSE)
  }
  marker_parent <- normalizePath(
    dirname(.e2e_cleanup_marker), mustWork = TRUE)
  marker <- file.path(marker_parent, basename(.e2e_cleanup_marker))
  writeLines(.e2e_state_root, marker, useBytes = TRUE)
  if (identical(.e2e_cleanup_test, "failure")) {
    stop("Forced harness cleanup self-test failure", call. = FALSE)
  }
  retained <- .e2e_finish_state(.e2e_state_guard)
  cat(if (is.null(retained)) "cleanup-self-test: removed\n" else
    paste0("cleanup-self-test: retained ", retained, "\n"))
  quit(save = "no", status = 0L)
}

.e2e_sampler_commands <- c(
  "exact-gc-worker",
  "joint-dp-vector-convolution-share-v3",
  "joint-dp-vector-convolution-finalize-v3",
  "joint-dp-vector-gaussian-share-v2",
  "joint-dp-vector-gaussian-finalize-v2",
  "joint-dp-convolution-share-v1",
  "dp-noise-int64",
  "dp-gaussian-int64")

.e2e_create_wrapper <- function(peer_dir) {
  wrapper <- file.path(peer_dir, "dsvert-mpc-observer")
  log_path <- file.path(peer_dir, "mpc-command.log")
  writeLines(c(
    "#!/bin/sh",
    "cmd=\"${1-}\"",
    paste0("printf '%s\\n' \"$cmd\" >> ", shQuote(log_path)),
    paste0("exec ", shQuote(.e2e_real_mpc), " \"$@\"")), wrapper,
    useBytes = TRUE)
  Sys.chmod(wrapper, mode = "0700")
  file.create(log_path)
  Sys.chmod(log_path, mode = "0600")
  list(wrapper = normalizePath(wrapper, mustWork = TRUE),
       log = normalizePath(log_path, mustWork = TRUE))
}

# The connector inherits DSLiteConnection so dsVertClient applies its normal
# in-process DSLite transport contract.  Only execution is moved into a
# separate persistent R process, which gives each peer independent options,
# secrets, filesystem state and concurrent MPC execution.
setClass(
  "DSVertE2EProcessConnection", contains = "DSLiteConnection",
  slots = c(worker = "ANY"))
setClass(
  "DSVertE2EProcessResult", contains = "DSResult",
  slots = c(conn = "DSVertE2EProcessConnection", consumed = "environment"))

.e2e_worker_submit <- function(worker, operation, args, async) {
  if (!identical(worker$get_state(), "idle")) {
    stop("The isolated DSLite peer already has a pending command",
         call. = FALSE)
  }
  job <- switch(
    operation,
    aggregate = function(expr) {
      get(".dsvert_e2e_server", .GlobalEnv)$aggregate(
        get(".dsvert_e2e_sid", .GlobalEnv), expr)
    },
    assign = function(symbol, expr) {
      get(".dsvert_e2e_server", .GlobalEnv)$assignExpr(
        get(".dsvert_e2e_sid", .GlobalEnv), symbol, expr)
      NULL
    })
  if (isTRUE(async)) {
    worker$call(job, args = args)
    return(NULL)
  }
  worker$run(job, args = args)
}

setMethod("dsIsAsync", "DSVertE2EProcessConnection", function(conn) {
  list(session = NULL, aggregate = TRUE, assignTable = FALSE,
       assignExpr = TRUE, assignResource = FALSE)
})
setMethod(
  "dsAggregate", "DSVertE2EProcessConnection",
  function(conn, expr, async = TRUE) {
    payload <- .e2e_worker_submit(
      conn@worker, "aggregate", list(expr = expr), async)
    state <- new.env(parent = emptyenv())
    if (!isTRUE(async)) state$sync_payload <- list(value = payload)
    new("DSVertE2EProcessResult", conn = conn,
        consumed = state)
  })
setMethod(
  "dsAssignExpr", "DSVertE2EProcessConnection",
  function(conn, symbol, expr, async = TRUE) {
    payload <- .e2e_worker_submit(
      conn@worker, "assign", list(symbol = symbol, expr = expr), async)
    state <- new.env(parent = emptyenv())
    if (!isTRUE(async)) state$sync_payload <- list(value = payload)
    new("DSVertE2EProcessResult", conn = conn,
        consumed = state)
  })
setMethod("dsIsCompleted", "DSVertE2EProcessResult", function(res) {
  identical(res@conn@worker$poll_process(0), "ready")
})
setMethod("dsFetch", "DSVertE2EProcessResult", function(res) {
  if (isTRUE(res@consumed$yes)) {
    stop("The isolated DSLite result was already consumed", call. = FALSE)
  }
  if (exists("sync_payload", envir = res@consumed, inherits = FALSE)) {
    payload <- res@consumed$sync_payload$value
    res@consumed$yes <- TRUE
    return(payload)
  }
  while (!identical(res@conn@worker$poll_process(1000), "ready")) {
    invisible(NULL)
  }
  payload <- res@conn@worker$read()
  res@consumed$yes <- TRUE
  if (!is.null(payload$error)) {
    stop(conditionMessage(payload$error), call. = FALSE)
  }
  payload$result
})
setMethod("dsKeepAlive", "DSVertE2EProcessConnection", function(conn) {
  invisible(TRUE)
})
setMethod(
  "dsDisconnect",
  signature(conn = "DSVertE2EProcessConnection", save = "ANY"),
  function(conn, save = NULL) invisible(TRUE))

.e2e_dummy_server <- DSLite::newDSLiteServer(
  config = list(AggregateMethods = data.frame(),
                AssignMethods = data.frame(), Options = list()),
  strict = TRUE, home = file.path(.e2e_state_root, "connector-dummy"))

.e2e_boot_peer <- function(spec, pins = NULL, raw = NULL,
                           configure_count = FALSE) {
  worker <- callr::r_session$new(wait = TRUE)
  boot <- tryCatch(worker$run(
    function(server_dir, spec, pins, raw, configure_count) {
      options(
        dsvert.state_dir = spec$state_dir,
        dsvert.peer_name = spec$name,
        dsvert.mpc_binary = spec$wrapper,
        dsvert.psi.max_input_ids = 64L)
      pkgload::load_all(server_dir, quiet = TRUE)

      description <- read.dcf(file.path(server_dir, "DESCRIPTION"))
      parse_field <- function(field, type) {
        methods <- trimws(strsplit(
          description[1L, field], ",", fixed = TRUE)[[1L]])
        methods <- methods[nzchar(methods)]
        if (anyDuplicated(methods) ||
            any(grepl("=", methods, fixed = TRUE)) ||
            any(!grepl("^[A-Za-z][A-Za-z0-9.]*DS$", methods))) {
          stop("The source DESCRIPTION has an invalid remote allowlist",
               call. = FALSE)
        }
        data.frame(
          name = methods, value = paste0("dsVert::", methods),
          package = "dsVert", version = unname(description[1L, "Version"]),
          type = type, class = "function", stringsAsFactors = FALSE)
      }
      config <- list(
        AggregateMethods = parse_field("AggregateMethods", "aggregate"),
        AssignMethods = parse_field("AssignMethods", "assign"),
        Options = list())

      expected <- dsVert:::.dsvert_remote_surface_contract(
        file.path(server_dir, "DESCRIPTION"))
      actual <- rbind(
        config$AggregateMethods[, c("name", "type", "value")],
        config$AssignMethods[, c("name", "type", "value")])
      order_contract <- function(value) value[
        order(value$type, value$name, method = "radix"), , drop = FALSE]
      expected <- order_contract(expected)
      actual <- order_contract(actual)
      rownames(expected) <- rownames(actual) <- NULL
      if (!identical(actual, expected) ||
          nrow(actual) != length(unique(actual$name)) ||
          !all(vapply(actual$name, function(name) {
            exists(name, envir = asNamespace("dsVert"), mode = "function",
                   inherits = FALSE) &&
              name %in% getNamespaceExports("dsVert")
          }, logical(1L)))) {
        stop("The effective DSLite surface does not match the source contract",
             call. = FALSE)
      }
      canonical <- paste(
        paste(actual$type, actual$name, actual$value, sep = "\x1f"),
        collapse = "\x1e")
      # Documented Rock administrator boundary: verify the effective inventory
      # above, then derive the deterministic public token from that exact
      # installed server contract.
      attestation <- dsVert:::.dsvert_remote_surface_attestation_expected(
        file.path(server_dir, "DESCRIPTION"))
      options(dsvert.remote_surface_attestation = attestation)

      if (!is.null(pins)) {
        options(dsvert.trusted_peers = pins)
      }
      if (!is.null(raw)) {
        authorization <- dsVert::dsvertPSISourceDescriptor(
          raw, id_col = "patient_id", id = spec$dataset_id,
          version = spec$dataset_version,
          purpose = "patient-record-alignment-v1")
        options(dsvert.psi.authorized_sources = list(D = authorization))
      }
      if (isTRUE(configure_count)) {
        options(
          dsvert.dp.epsilon = 1,
          dsvert.dp.delta = 1e-6,
          dsvert.dp.implementation_delta = 1e-9,
          dsvert.dp.domain = spec$domain,
          dsvert.dp.cohort_id = spec$cohort_id,
          dsvert.dp.adjacency = "add_remove_patient",
          dsvert.dp.unit_capacity = 64L,
          dsvert.dp.fixed_cohort_size = NULL)
      }

      tables <- if (is.null(raw)) list() else list(D = raw)
      server <- DSLite::newDSLiteServer(
        tables = tables, config = config, strict = TRUE,
        home = spec$dslite_home)
      sid <- server$newSession(profile = "default")
      if (!is.null(raw)) server$assignTable(sid, "D", "D")
      assign(".dsvert_e2e_server", server, envir = .GlobalEnv)
      assign(".dsvert_e2e_sid", sid, envir = .GlobalEnv)
      list(
        surface_count = nrow(actual),
        surface_sha256 = digest::digest(
          canonical, algo = "sha256", serialize = FALSE),
        source_present = !is.null(raw),
        count_configured = isTRUE(configure_count))
    }, args = list(
      server_dir = .e2e_server_dir, spec = spec, pins = pins, raw = raw,
      configure_count = configure_count)),
    error = identity)
  if (inherits(boot, "error")) {
    try(worker$close(), silent = TRUE)
    stop(boot)
  }
  connection <- new(
    "DSVertE2EProcessConnection", name = spec$name,
    sid = paste0("process:", spec$name, ":", worker$get_pid()),
    server = .e2e_dummy_server, worker = worker)
  list(worker = worker, connection = connection, boot = boot)
}

.e2e_close_peers <- function(peers) {
  for (peer in peers) try(peer$worker$close(), silent = TRUE)
  invisible(TRUE)
}

.e2e_connections <- function(peers) {
  stats::setNames(lapply(peers, `[[`, "connection"), names(peers))
}

.e2e_load_client <- function() {
  if ("dsVertClient" %in% loadedNamespaces()) {
    pkgload::unload("dsVertClient")
  }
  pkgload::load_all(.e2e_client_dir, quiet = TRUE)
  invisible(TRUE)
}

.e2e_raw_data <- function(index) {
  common <- sprintf("patient-%03d", 1:12)
  extra <- sprintf("site-%02d-extra-%02d", index, 1:2)
  ids <- c(common, extra)
  order <- if (index %% 2L) rev(seq_along(ids)) else {
    c(seq(2L, length(ids), by = 2L), seq(1L, length(ids), by = 2L))
  }
  data.frame(
    patient_id = ids[order],
    value = as.numeric((seq_along(ids) * (index + 3L)) %% 97L)[order],
    stringsAsFactors = FALSE, check.names = FALSE)
}

.e2e_make_specs <- function(k, scenario_dir) {
  names <- paste0("site_", letters[seq_len(k)])
  stats::setNames(lapply(seq_along(names), function(index) {
    peer_dir <- file.path(scenario_dir, names[[index]])
    dir.create(peer_dir, recursive = TRUE, mode = "0700")
    Sys.chmod(peer_dir, mode = "0700")
    observer <- .e2e_create_wrapper(peer_dir)
    list(
      name = names[[index]], state_dir = peer_dir,
      dslite_home = file.path(peer_dir, "dslite"),
      wrapper = observer$wrapper, command_log = observer$log,
      dataset_id = paste0("formal-dp-e2e-k", k),
      dataset_version = "v1", domain = paste0("formal-dp-e2e-k", k),
      cohort_id = paste0("formal-dp-e2e-k", k, "-cohort"))
  }), names)
}

.e2e_truncate_logs <- function(specs) {
  for (spec in specs) {
    connection <- file(spec$command_log, open = "wb")
    close(connection)
    Sys.chmod(spec$command_log, mode = "0600")
  }
  invisible(TRUE)
}

.e2e_command_counts <- function(specs) {
  stats::setNames(lapply(specs, function(spec) {
    commands <- readLines(spec$command_log, warn = FALSE)
    commands <- commands[nzchar(commands)]
    sampler <- commands[commands %in% .e2e_sampler_commands]
    list(sampler = length(sampler),
         exact_gc_workers = sum(sampler == "exact-gc-worker"),
         sampler_histogram = as.list(table(factor(
           sampler, levels = .e2e_sampler_commands))))
  }), names(specs))
}

.e2e_assert <- function(value, message) {
  if (!isTRUE(value)) stop(message, call. = FALSE)
  invisible(TRUE)
}

.e2e_run_scenario <- function(k) {
  message("[formal-dp-e2e] K=", k, ": bootstrapping isolated peers")
  scenario_dir <- file.path(.e2e_state_root, paste0("k", k))
  dir.create(scenario_dir, recursive = TRUE, mode = "0700")
  Sys.chmod(scenario_dir, mode = "0700")
  specs <- .e2e_make_specs(k, scenario_dir)
  peer_names <- names(specs)
  raw <- stats::setNames(lapply(seq_len(k), function(index) {
    .e2e_raw_data(index)
  }), peer_names)

  bootstrap <- stats::setNames(lapply(peer_names, function(peer) {
    .e2e_boot_peer(specs[[peer]], raw = raw[[peer]])
  }), peer_names)
  on.exit(.e2e_close_peers(bootstrap), add = TRUE)
  .e2e_load_client()
  bootstrap_conns <- .e2e_connections(bootstrap)
  identities <- dsVertClient::ds.getIdentityPks(bootstrap_conns)
  .e2e_assert(identical(names(identities), peer_names),
              "Identity discovery returned the wrong peer set")
  .e2e_assert(length(unique(unlist(identities, use.names = FALSE))) == k,
              "Every peer must have a unique identity key")
  .e2e_close_peers(bootstrap)

  pins <- stats::setNames(lapply(peer_names, function(peer) {
    values <- unlist(identities[setdiff(peer_names, peer)], use.names = FALSE)
    stats::setNames(values, setdiff(peer_names, peer))
  }), peer_names)
  peers <- stats::setNames(lapply(peer_names, function(peer) {
    .e2e_boot_peer(specs[[peer]], pins = pins[[peer]], raw = raw[[peer]])
  }), peer_names)
  on.exit(.e2e_close_peers(peers), add = TRUE)
  conns <- .e2e_connections(peers)
  identities_after_pin <- dsVertClient::ds.getIdentityPks(conns)
  .e2e_assert(identical(identities_after_pin, identities),
              "Identity keys changed across the provisioning restart")

  message("[formal-dp-e2e] K=", k, ": running public pinned PSI")
  dsVertClient::ds.psiAlign(
    "D", "patient_id", "DA", datasources = conns, verbose = FALSE)
  aligned <- dsVertClient::ds.isPsiAligned("DA", datasources = conns)
  .e2e_assert(isTRUE(aligned$aligned) && is.na(aligned$n_common),
              "Pinned PSI did not return its non-disclosive attestation")

  # Restart once into the exact custodian-owned Count configuration. The raw
  # data are retained because a stateless current-snapshot artifact is
  # deliberately recomputed rather than served from durable release storage.
  .e2e_close_peers(peers)
  peers <- stats::setNames(lapply(peer_names, function(peer) {
    .e2e_boot_peer(
      specs[[peer]], pins = pins[[peer]], raw = raw[[peer]],
      configure_count = TRUE)
  }), peer_names)
  conns <- .e2e_connections(peers)

  # Recreate the aligned object through the public protocol after the service
  # restart. PSI run identifiers may change; Count semantics must not.
  dsVertClient::ds.psiAlign(
    "D", "patient_id", "DA", datasources = conns, verbose = FALSE)
  .e2e_assert(isTRUE(dsVertClient::ds.isPsiAligned(
    "DA", datasources = conns)$aligned),
    "Pinned PSI did not survive Count configuration bootstrap")
  .e2e_truncate_logs(specs)

  message("[formal-dp-e2e] K=", k, ": creating one canonical Count release")
  first <- dsVertClient::ds.vertDPCount("DA", datasources = conns)
  legacy_fields <- c(
    "capsule_id", "final_vector_root", "history_gate", "request_limit",
    "operation_limit")
  .e2e_assert(
    identical(first$released, TRUE) &&
      identical(first$adjacency, "add_remove_patient") &&
      identical(first$sticky_replay, TRUE) &&
      identical(first$intermediate_values_exposed, FALSE) &&
      identical(first$privacy$finite_global_composition_claim, FALSE) &&
      identical(first$privacy$distinct_artifacts_compose, TRUE) &&
      identical(as.numeric(first$privacy$public_openings), 1) &&
      !length(intersect(names(first), legacy_fields)),
    "Count did not expose the canonical per-artifact release contract")
  identity_values <- unlist(identities, use.names = FALSE)
  authority_pks <- c(
    first$signed_release$source_identity_pk,
    first$signed_release$finalizer_identity_pk)
  authorities <- names(identities)[match(authority_pks, identity_values)]
  .e2e_assert(
    length(authorities) == 2L && !anyNA(authorities) &&
      !anyDuplicated(authorities),
    "The signed Count authorities do not map to two pinned peers")
  first_counts <- .e2e_command_counts(specs)
  .e2e_assert(all(vapply(authorities, function(peer) {
    identical(first_counts[[peer]]$sampler, 1L) &&
      identical(first_counts[[peer]]$exact_gc_workers, 1L)
  }, logical(1L))),
  "Each Count authority must invoke exactly one exact-GC sampler worker")
  non_authorities <- setdiff(peer_names, authorities)
  .e2e_assert(all(vapply(non_authorities, function(peer) {
    identical(first_counts[[peer]]$sampler, 0L)
  }, logical(1L))),
  "A non-authority peer invoked a Count sampler")

  # A real client-process restart creates a new ephemeral MPC session. The
  # sampler runs again, but persistent identity-derived subseeds must reproduce
  # exactly the same public artifact and release bytes.
  .e2e_load_client()
  replay <- dsVertClient::ds.vertDPCount("DA", datasources = conns)
  replay_counts <- .e2e_command_counts(specs)
  replay_incremented <- all(vapply(peer_names, function(peer) {
    increment <- if (peer %in% authorities) 1L else 0L
    identical(replay_counts[[peer]]$sampler,
              first_counts[[peer]]$sampler + increment) &&
      identical(replay_counts[[peer]]$exact_gc_workers,
                first_counts[[peer]]$exact_gc_workers + increment)
  }, logical(1L)))
  .e2e_assert(
    identical(first$value, replay$value) &&
      identical(first$artifact_key, replay$artifact_key) &&
      identical(first$release_sha256, replay$release_sha256) &&
      identical(first$signed_release, replay$signed_release),
    "A new client session changed the sticky Count artifact or release")
  .e2e_assert(replay_incremented,
              "A new Count session did not recompute on exactly two peers")

  message("[formal-dp-e2e] K=", k,
          ": restarting every server with the current protected snapshot")
  .e2e_close_peers(peers)
  peers <- stats::setNames(lapply(peer_names, function(peer) {
    .e2e_boot_peer(
      specs[[peer]], pins = pins[[peer]], raw = raw[[peer]],
      configure_count = TRUE)
  }), peer_names)
  conns <- .e2e_connections(peers)
  source_present <- vapply(peers, function(peer) {
    identical(peer$boot$source_present, TRUE)
  }, logical(1L))
  .e2e_assert(all(source_present),
              "The current protected source was not restored for recomputation")
  restart_identities <- dsVertClient::ds.getIdentityPks(conns)
  .e2e_assert(identical(restart_identities, identities),
              "Identity keys changed across the service restart")
  dsVertClient::ds.psiAlign(
    "D", "patient_id", "DA", datasources = conns, verbose = FALSE)
  .e2e_assert(isTRUE(dsVertClient::ds.isPsiAligned(
    "DA", datasources = conns)$aligned),
    "Pinned PSI did not reconstruct the current Count snapshot")
  .e2e_load_client()
  restarted <- dsVertClient::ds.vertDPCount("DA", datasources = conns)
  restart_counts <- .e2e_command_counts(specs)
  restart_incremented <- all(vapply(peer_names, function(peer) {
    increment <- if (peer %in% authorities) 1L else 0L
    identical(restart_counts[[peer]]$sampler,
              replay_counts[[peer]]$sampler + increment) &&
      identical(restart_counts[[peer]]$exact_gc_workers,
                replay_counts[[peer]]$exact_gc_workers + increment)
  }, logical(1L)))
  .e2e_assert(
    identical(first$value, restarted$value) &&
      identical(first$artifact_key, restarted$artifact_key) &&
      identical(first$release_sha256, restarted$release_sha256) &&
      identical(first$signed_release, restarted$signed_release),
    "A service restart changed the sticky Count artifact or release")
  .e2e_assert(restart_incremented,
              "Restarted Count did not recompute on exactly two peers")

  .e2e_close_peers(peers)
  list(
    k = k, passed = TRUE,
    surface_method_count = bootstrap[[1L]]$boot$surface_count,
    surface_sha256 = bootstrap[[1L]]$boot$surface_sha256,
    mpc_binary_sha256 = .e2e_real_mpc_sha256,
    psi_aligned_without_cardinality = TRUE,
    protected_source_present_on_restart = all(source_present),
    identity_persistent = identical(restart_identities, identities),
    dp_value = first$value,
    artifact_key = first$artifact_key,
    release_sha256 = first$release_sha256,
    sticky_same_output = identical(first$value, replay$value),
    restart_same_output = identical(first$value, restarted$value),
    finite_global_composition_claim =
      first$privacy$finite_global_composition_claim,
    sampler_invocations = stats::setNames(vapply(
      first_counts, `[[`, integer(1L), "sampler"), peer_names),
    sampler_recomputed_after_client_restart = replay_incremented,
    sampler_recomputed_after_service_restart = restart_incremented)
}

.e2e_load_client()
.e2e_started <- Sys.time()
.e2e_results <- lapply(.e2e_k, .e2e_run_scenario)
.e2e_real_mpc_mode_after <- if (.Platform$OS.type == "windows") NA_integer_ else
  as.integer(file.info(.e2e_real_mpc)$mode)
.e2e_assert(identical(.e2e_real_mpc_mode_after, .e2e_real_mpc_mode),
            "The validation harness changed the packaged MPC binary mode")
.e2e_retained_state <- .e2e_finish_state(.e2e_state_guard)
.e2e_report <- list(
  version = "dsvert-canonical-count-e2e-v1",
  passed = all(vapply(.e2e_results, `[[`, logical(1L), "passed")),
  started_utc = format(.e2e_started, tz = "UTC", usetz = TRUE),
  completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  scenarios = .e2e_results,
  mpc_binary_mode_unchanged = TRUE,
  state_retained = !is.null(.e2e_retained_state),
  state_path = .e2e_retained_state)

encoded <- jsonlite::toJSON(
  .e2e_report, auto_unbox = TRUE, pretty = TRUE, null = "null", digits = 17)
cat(encoded, "\n")
if (nzchar(.e2e_output)) {
  output <- normalizePath(dirname(.e2e_output), mustWork = TRUE)
  path <- file.path(output, basename(.e2e_output))
  writeLines(encoded, path, useBytes = TRUE)
}
if (!isTRUE(.e2e_report$passed)) quit(save = "no", status = 1L)
