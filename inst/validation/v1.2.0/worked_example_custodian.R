# CUSTODIAN-SIDE SCAFFOLD for the chapter-4 worked example (dsVert v1.2.0).
#
# Adapted from dsVertClient v1.2.0 tools/validate_formal_dp_e2e.R (the
# released source-tree validation harness).  Each peer is a separate
# persistent R process with an exact DSLite allowlist parsed from the dsVert
# DESCRIPTION, because dsVert server state (peer name, identity keys, policy
# options, key store) is process-scoped by design.  There is no analyst
# bypass: the analyst surface used by the worked example is only the
# exported dsVertClient API over these DSI connections.
#
# The scaffold exposes:
#   we_boot_peer(spec, server_dir, pins, raw, policy)  -> peer handle
#   we_connections(peers)                              -> named DSI conns
#   we_close_peers(peers)
# and the S4 connection classes routing DSI aggregate/assign calls into the
# per-peer R process (verbatim from the release harness).

suppressPackageStartupMessages({
  library(callr)
  library(digest)
  library(DSI)
  library(DSLite)
  library(jsonlite)
  library(methods)
  library(pkgload)
})

setClass(
  "DSVertE2EProcessConnection", contains = "DSLiteConnection",
  slots = c(worker = "ANY"))
setClass(
  "DSVertE2EProcessResult", contains = "DSResult",
  slots = c(conn = "DSVertE2EProcessConnection", consumed = "environment"))

.we_worker_submit <- function(worker, operation, args, async) {
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
    payload <- .we_worker_submit(
      conn@worker, "aggregate", list(expr = expr), async)
    state <- new.env(parent = emptyenv())
    if (!isTRUE(async)) state$sync_payload <- list(value = payload)
    new("DSVertE2EProcessResult", conn = conn, consumed = state)
  })
setMethod(
  "dsAssignExpr", "DSVertE2EProcessConnection",
  function(conn, symbol, expr, async = TRUE) {
    payload <- .we_worker_submit(
      conn@worker, "assign", list(symbol = symbol, expr = expr), async)
    state <- new.env(parent = emptyenv())
    if (!isTRUE(async)) state$sync_payload <- list(value = payload)
    new("DSVertE2EProcessResult", conn = conn, consumed = state)
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

.we_dummy_server_home <- file.path(tempdir(), "we-connector-dummy")
.we_dummy_server <- DSLite::newDSLiteServer(
  config = list(AggregateMethods = data.frame(),
                AssignMethods = data.frame(), Options = list()),
  strict = TRUE, home = .we_dummy_server_home)

# Boot one custodian peer process.  `policy` is the per-peer custodian DP
# Synopsis configuration (a named list of dsvert.dp.* option values); `raw`
# is that custodian's private table, published to the DSLite session as "D".
we_boot_peer <- function(spec, server_dir, pins = NULL, raw = NULL,
                         policy = NULL) {
  worker <- callr::r_session$new(wait = TRUE)
  boot <- tryCatch(worker$run(
    function(server_dir, spec, pins, raw, policy) {
      options(
        dsvert.state_dir = spec$state_dir,
        dsvert.peer_name = spec$name,
        dsvert.mpc_binary = spec$mpc_binary)
      pkgload::load_all(server_dir, quiet = TRUE)

      description <- read.dcf(file.path(server_dir, "DESCRIPTION"))
      parse_field <- function(field, type) {
        methods <- trimws(strsplit(
          description[1L, field], ",", fixed = TRUE)[[1L]])
        methods <- methods[nzchar(methods)]
        data.frame(
          name = methods, value = paste0("dsVert::", methods),
          package = "dsVert", version = unname(description[1L, "Version"]),
          type = type, class = "function", stringsAsFactors = FALSE)
      }
      config <- list(
        AggregateMethods = parse_field("AggregateMethods", "aggregate"),
        AssignMethods = parse_field("AssignMethods", "assign"),
        Options = list())

      # Documented Rock administrator boundary: verify the effective
      # inventory, then derive the deterministic public surface token from
      # the exact installed server contract (as in the release harness).
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
      if (!identical(actual, expected)) {
        stop("The effective DSLite surface does not match the source contract",
             call. = FALSE)
      }
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
      if (!is.null(policy)) {
        names(policy) <- paste0("dsvert.dp.", names(policy))
        do.call(options, policy)
      }

      tables <- if (is.null(raw)) list() else list(D = raw)
      server <- DSLite::newDSLiteServer(
        tables = tables, config = config, strict = TRUE,
        home = spec$dslite_home)
      sid <- server$newSession(profile = "default")
      if (!is.null(raw)) server$assignTable(sid, "D", "D")
      assign(".dsvert_e2e_server", server, envir = .GlobalEnv)
      assign(".dsvert_e2e_sid", sid, envir = .GlobalEnv)
      list(surface_count = nrow(actual), source_present = !is.null(raw))
    }, args = list(server_dir = server_dir, spec = spec, pins = pins,
                   raw = raw, policy = policy)),
    error = identity)
  if (inherits(boot, "error")) {
    try(worker$close(), silent = TRUE)
    stop(boot)
  }
  connection <- new(
    "DSVertE2EProcessConnection", name = spec$name,
    sid = paste0("process:", spec$name, ":", worker$get_pid()),
    server = .we_dummy_server, worker = worker)
  list(worker = worker, connection = connection, boot = boot)
}

we_connections <- function(peers) {
  stats::setNames(lapply(peers, `[[`, "connection"), names(peers))
}

we_close_peers <- function(peers) {
  for (peer in peers) try(peer$worker$close(), silent = TRUE)
  invisible(TRUE)
}
