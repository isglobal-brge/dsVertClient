# Fail-closed wrappers around DSI 1.8 aggregate fan-out. DSI deliberately
# returns a named list containing NULL for a site that failed; callers must not
# mistake that partial list for a successful protocol phase. Connector error
# text is intentionally discarded client-side because it may contain server
# implementation details; privileged logging belongs on the data node.

.dsvert_client_parse_release_instance_retry <- function(message) {
  if (!is.character(message) || length(message) != 1L || is.na(message)) {
    return(NULL)
  }
  token <- regmatches(message, regexpr(
    paste0(
      "\\[dsvert_retry_current_release_instance:",
      "(retry_unpublished_instance|new_release_instance)\\]"),
    message, perl = TRUE))
  if (!length(token) || !nzchar(token)) return(NULL)
  action <- sub(
    "^\\[dsvert_retry_current_release_instance:([^]]+)\\]$", "\\1",
    token, perl = TRUE)
  structure(list(
    message = paste0(
      token,
      " A pinned server requested one refresh of the current release roots."),
    call = NULL, retry_action = action),
    class = c("dsvert_release_instance_retry", "error", "condition"))
}

.dsvert_client_parse_resource_backpressure <- function(message) {
  if (!is.character(message) || length(message) != 1L || is.na(message) ||
      !grepl("[dsvert_resource_backpressure:v1]", message, fixed = TRUE)) {
    return(NULL)
  }
  structure(list(
    message = paste0(
      "[dsvert_resource_backpressure:v1] A pinned server reported transient ",
      "byte-capacity backpressure for an idempotent transport operation."),
    call = NULL, code = "resource_backpressure", retryable = TRUE),
    class = c("dsvert_resource_backpressure", "error", "condition"))
}

.dsvert_client_parse_phase_not_ready <- function(message) {
  if (!is.character(message) || length(message) != 1L || is.na(message) ||
      !grepl("[dsvert_phase_not_ready:v1]", message, fixed = TRUE)) {
    return(NULL)
  }
  structure(list(
    message = paste0(
      "[dsvert_phase_not_ready:v1] A pinned server reported that this ",
      "idempotent protocol phase is not ready."),
    call = NULL, code = "phase_not_ready", retryable = FALSE),
    class = c("dsvert_phase_not_ready", "error", "condition"))
}

.dsvert_client_parse_dp_lifetime_budget_exhausted <- function(message) {
  token <- "[dsvert_dp_lifetime_budget_exhausted:v1]"
  if (!is.character(message) || length(message) != 1L || is.na(message) ||
      !grepl(token, message, fixed = TRUE)) {
    return(NULL)
  }
  structure(list(
    message = token, call = NULL,
    code = "dp_lifetime_budget_exhausted", retryable = FALSE),
    class = c(
      "dsvert_dp_lifetime_budget_exhausted", "error", "condition"))
}

.dsvert_client_has_other_bracketed_typed_token <- function(message, token) {
  if (!is.character(message) || length(message) != 1L || is.na(message)) {
    return(FALSE)
  }
  remainder <- gsub(token, "", message, fixed = TRUE)
  grepl("\\[dsvert_[^]\\r\\n]+\\]", remainder, perl = TRUE)
}

.dsvert_client_resource_oversize <- function(
    message = NULL, requested_bytes = NULL, capacity_bytes = NULL,
    scope = "DataSHIELD transport") {
  if (is.null(message)) {
    message <- paste0(
      "[dsvert_resource_oversize:v1] resource_oversize: ", scope,
      " cannot fit within its fixed byte/resource policy.")
  }
  fields <- list(
    message = message, call = NULL, code = "resource_oversize",
    retryable = FALSE, scope = scope)
  if (!is.null(requested_bytes)) fields$requested_bytes <- requested_bytes
  if (!is.null(capacity_bytes)) fields$capacity_bytes <- capacity_bytes
  structure(fields,
    class = c("dsvert_resource_oversize", "error", "condition"))
}

.dsvert_client_parse_resource_oversize <- function(message) {
  if (!is.character(message) || length(message) != 1L || is.na(message) ||
      !grepl("[dsvert_resource_oversize:v1]", message, fixed = TRUE)) {
    return(NULL)
  }
  .dsvert_client_resource_oversize(
    message = paste0(
      "[dsvert_resource_oversize:v1] A pinned server rejected a payload ",
      "that cannot fit its fixed public resource policy."),
    scope = "pinned server transport")
}

.dsvert_client_remote_contract_failure <- function(operation) {
  structure(list(
    message = paste0(
      "A pinned server returned a terminal negative outcome during '",
      operation, "'; no transport retry was attempted."),
    call = NULL, code = "remote_contract_failure", retryable = FALSE,
    operation = operation),
    class = c("dsvert_remote_contract_failure", "error", "condition"))
}

#' Portable upper bound for one deparsed DSI expression (internal)
#'
#' DSI 1.8 with DSLite 1.4.1 accepts a 512 KiB raw frame as about 683 KiB of
#' Base64url expression text. Production therefore probes 688 KiB of request
#' text while retaining a complete-expression bound below 768 KiB. No payload
#' uses that geometry until every target accepts the stateless probe, and the
#' resulting common geometry is frozen before its first payload transmission.
#' @keywords internal
.dsvert_dsi_max_expression_bytes <- function() {
  negotiated <- if (exists(".dsvert_chunk_env", inherits = TRUE) &&
      is.environment(.dsvert_chunk_env)) {
    .dsvert_chunk_env$effective_expression_bytes
  } else {
    NULL
  }
  value <- if (is.null(negotiated)) {
    getOption("dsvert.dsi_max_expression_bytes", 768L * 1024L - 1L)
  } else {
    negotiated
  }
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      !is.finite(value) || value != floor(value) || value < 16 * 1024L ||
      value > 256 * 1024^2) {
    stop("dsvert.dsi_max_expression_bytes must be one integer between 16 KiB and 256 MiB",
         call. = FALSE)
  }
  as.numeric(value)
}

#' Fail before DSI parsing when a request exceeds negotiated geometry.
#' @param expr One DSI call or a named list of calls.
#' @param capacity_bytes Optional scalar or per-expression named byte ceiling.
#'   \code{NULL} uses the current negotiated ceiling.
#' @keywords internal
.dsvert_validate_dsi_expression_sizes <- function(
    expr, capacity_bytes = NULL) {
  expressions <- if (is.list(expr) && !is.call(expr)) expr else list(expr)
  if (!length(expressions) || any(!vapply(expressions, is.call, logical(1L)))) {
    stop("DSI expressions must contain R calls", call. = FALSE)
  }
  sizes <- vapply(expressions, function(value) {
    text <- paste(deparse(value, width.cutoff = 500L), collapse = "\n")
    nchar(text, type = "bytes")
  }, numeric(1L))
  if (is.null(capacity_bytes)) {
    capacities <- rep(.dsvert_dsi_max_expression_bytes(), length(sizes))
  } else {
    capacities <- suppressWarnings(as.numeric(capacity_bytes))
    if (length(capacities) == 1L) {
      capacities <- rep(capacities, length(sizes))
    } else {
      expression_names <- names(expressions)
      capacity_names <- names(capacity_bytes)
      if (is.null(expression_names) || is.null(capacity_names) ||
          anyNA(expression_names) || anyNA(capacity_names) ||
          any(!nzchar(expression_names)) || any(!nzchar(capacity_names)) ||
          anyDuplicated(expression_names) || anyDuplicated(capacity_names) ||
          !setequal(expression_names, capacity_names)) {
        stop("Per-site DSI expression capacities must map exactly to expressions",
             call. = FALSE)
      }
      capacities <- capacities[match(expression_names, capacity_names)]
    }
    if (length(capacities) != length(sizes) || anyNA(capacities) ||
        any(!is.finite(capacities)) || any(capacities != floor(capacities)) ||
        any(capacities < 16 * 1024L) || any(capacities > 256 * 1024^2)) {
      stop("Invalid per-site DSI expression capacities", call. = FALSE)
    }
  }
  oversized <- sizes > capacities
  if (any(oversized)) {
    stop(.dsvert_client_resource_oversize(
      requested_bytes = max(sizes[oversized]),
      capacity_bytes = min(capacities[oversized]),
      scope = "DataSHIELD negotiated expression"))
  }
  invisible(sizes)
}

# DSI::datashield.aggregate() deparses every expression again while polling
# asynchronous jobs and while formatting progress messages.  For an immutable
# 640 KiB transport frame that repeated formatting is both avoidable CPU work
# and an avoidable full-size allocation.  The transport paths below already
# validate and freeze their calls, so use DSI's public job primitives directly:
# launch every asynchronous peer first, execute synchronous connectors second,
# then poll/fetch by logical site name.  Test harnesses keep receiving the
# injected aggregate function unchanged.
.dsvert_dsi_poll_seconds <- function(checks) {
  checks <- suppressWarnings(as.numeric(checks))
  if (length(checks) != 1L || is.na(checks) || !is.finite(checks) ||
      checks < 1 || checks != floor(checks)) {
    stop("Invalid DSI polling state", call. = FALSE)
  }
  t0 <- getOption("datashield.polling.sleep.0", 0.05)
  t1 <- getOption("datashield.polling.sleep.1", 1)
  t10 <- getOption("datashield.polling.sleep.10", t1 * 2)
  t60 <- getOption("datashield.polling.sleep.60", t1 * 10)
  t600 <- getOption("datashield.polling.sleep.600", t1 * 60)
  t3600 <- getOption("datashield.polling.sleep.3600", t1 * 600)
  timings <- suppressWarnings(as.numeric(c(t0, t1, t10, t60, t600, t3600)))
  if (length(timings) != 6L || anyNA(timings) ||
      any(!is.finite(timings)) || any(timings <= 0)) {
    stop("Invalid DataSHIELD polling interval", call. = FALSE)
  }
  n1 <- 1 / timings[[1L]]
  n10 <- n1 + 9 / timings[[2L]]
  n60 <- n10 + 50 / timings[[3L]]
  n600 <- n60 + 540 / timings[[4L]]
  n3600 <- n600 + 3000 / timings[[5L]]
  if (checks >= n3600) return(timings[[6L]])
  if (checks >= n600) return(timings[[5L]])
  if (checks >= n60) return(timings[[4L]])
  if (checks >= n10) return(timings[[3L]])
  if (checks >= n1) return(timings[[2L]])
  timings[[1L]]
}

# A remote statistical job is not an idempotent transport retry. Reusing the
# five-minute frame-availability deadline here caused valid long-running jobs
# to be abandoned and, for Armadillo, correctly poisoned their singleton
# last-result lifecycle. By default DSI jobs therefore run until completion or
# explicit user interruption. Deployments that need an operational wall-clock
# ceiling may set one without changing request or retry accounting.
.dsvert_dsi_job_timeout_seconds <- function() {
  value <- suppressWarnings(as.numeric(getOption(
    "dsvert.dsi.job_timeout_seconds", Inf)))
  if (length(value) != 1L || is.na(value) || value <= 0) {
    stop("dsvert.dsi.job_timeout_seconds must be positive or Inf",
         call. = FALSE)
  }
  value
}

# DSI has no portable cancellation primitive for DSResult jobs. In
# particular, Armadillo 3.0.x exposes one last-command/last-result lifecycle
# per authenticated handle. Reusing a handle after an unresolved async job can
# therefore associate a later fetch with the wrong submission. Remember such
# handles for this client process and require a fresh login/session instead of
# silently replaying on an ambiguous lifecycle.
.dsvert_dsi_job_state <- new.env(parent = emptyenv())
.dsvert_dsi_job_state$poisoned <- new.env(hash = TRUE, parent = emptyenv())

.dsvert_dsi_job_session_key <- function(connection) {
  details <- tryCatch(
    .dsvert_dsi_inspect_connection(connection), error = function(e) NULL)
  kind <- if (is.list(details) && is.character(details$kind) &&
      length(details$kind) == 1L) details$kind else "unknown"
  endpoint <- if (is.list(details) && is.list(details$endpoint) &&
      is.character(details$endpoint$scope) &&
      length(details$endpoint$scope) == 1L) {
    details$endpoint$scope
  } else {
    "in-process-or-unavailable"
  }
  session <- if (identical(kind, "armadillo")) {
    # Bearer-token refresh does not create a new Armadillo backend command
    # lifecycle. Bind poison state to the stable libcurl handle instead: a
    # fresh dsConnect/login creates a new handle, whereas token refresh and S4
    # copies retain this pointer.
    .dsvert_dsi_armadillo_handle_identity(connection)
  } else {
    .dsvert_dsi_connection_session(connection)
  }
  if (is.null(session)) return(NULL)
  fields <- c(
    "dsVert/dsi-pending-job/v1",
    paste(class(connection), collapse = "/"), kind, endpoint, session)
  digest::digest(
    paste0(nchar(fields, type = "bytes"), ":", fields, collapse = ""),
    algo = "sha256", serialize = FALSE)
}

.dsvert_dsi_poison_sessions <- function(keys) {
  keys <- unique(keys[is.character(keys) & !is.na(keys) & nzchar(keys)])
  for (key in keys) assign(key, TRUE, envir = .dsvert_dsi_job_state$poisoned)
  invisible(keys)
}

.dsvert_dsi_session_is_poisoned <- function(key) {
  is.character(key) && length(key) == 1L && !is.na(key) && nzchar(key) &&
    exists(key, envir = .dsvert_dsi_job_state$poisoned, inherits = FALSE)
}

.dsvert_dsi_clear_poisoned_sessions <- function() {
  rm(list = ls(.dsvert_dsi_job_state$poisoned, all.names = TRUE),
     envir = .dsvert_dsi_job_state$poisoned)
  invisible(TRUE)
}

.dsvert_dsi_poisoned_session_condition <- function(sites) {
  structure(list(
    message = paste0(
      "A previous DSI submission has an unresolved result lifecycle for ",
      paste(sites, collapse = ", "),
      ". DSI cannot cancel or safely disambiguate that job. Create fresh ",
      "DSI login connections before retrying this phase."),
    call = NULL, sites = sites),
    class = c("dsvert_dsi_poisoned_session", "error", "condition"))
}

.dsvert_dsi_poison_ambiguous_sessions <- function(
    conns, sites = names(conns)) {
  sites <- intersect(as.character(sites), names(conns))
  keys <- unlist(lapply(conns[sites], function(connection) {
    tryCatch(.dsvert_dsi_job_session_key(connection), error = function(e) NULL)
  }), use.names = FALSE)
  .dsvert_dsi_poison_sessions(keys)
  .dsvert_dsi_poisoned_session_condition(sites)
}

.dsvert_dsi_pre_job_rejection_is_settled <- function(condition) {
  message <- substr(conditionMessage(condition), 1L, 8192L)
  grepl(
    "^\\[Client error: \\(413\\) Request Entity Too Large\\](?:\\s|$)",
    message, perl = TRUE)
}

.dsvert_dsi_probe_bound_failure_is_settled <- function(condition) {
  message <- conditionMessage(condition)
  if (!startsWith(message, "Command 'dsvertTransportProbeDS(")) return(FALSE)
  tail_start <- max(1L, nchar(message, type = "chars") - 8191L)
  tail <- substr(message, tail_start, nchar(message, type = "chars"))
  any(vapply(c(
    "Transport-probe padding exceeds the server byte bound.",
    "Transport-probe padding exceeds the server byte bound",
    "Transport-response-probe padding exceeds the server byte bound.",
    "Transport-response-probe padding exceeds the server byte bound"),
    function(suffix) endsWith(tail, suffix), logical(1L)))
}

.dsvert_dsi_sync_failure_is_settled <- function(connection, condition) {
  if (inherits(condition, c(
      "dsvert_peer_not_recognized", "dsvert_resource_oversize",
      "dsvert_phase_not_ready", "dsvert_dp_lifetime_budget_exhausted"))) {
    return(TRUE)
  }
  message <- substr(conditionMessage(condition), 1L, 8192L)
  if (!is.null(.dsvert_client_parse_peer_not_recognized(message)) ||
      !is.null(.dsvert_client_parse_resource_oversize(message)) ||
      !is.null(.dsvert_client_parse_phase_not_ready(message)) ||
      !is.null(.dsvert_client_parse_dp_lifetime_budget_exhausted(message))) {
    return(TRUE)
  }
  details <- tryCatch(
    .dsvert_dsi_inspect_connection(connection), error = function(e) NULL)
  kind <- if (is.list(details)) details$kind else NULL
  if (identical(kind, "dslite")) return(TRUE)
  unsupported_argument <- grepl(
    paste0(
      "unused argument[^\\n]*(transport_contract|response_padding_chars)|",
      "(transport_contract|response_padding_chars)[^\\n]*unused argument"),
    message, ignore.case = TRUE, perl = TRUE)
  if (isTRUE(unsupported_argument)) return(TRUE)
  if (.dsvert_dsi_probe_bound_failure_is_settled(condition)) return(TRUE)
  if (identical(kind, "armadillo")) {
    return(grepl(
      paste0(
        "^(Command .* failed on .*: Error whilst evaluating |",
        "Bad request: |Unauthorized$)"),
      message, perl = TRUE))
  }
  if (identical(kind, "opal")) {
    # A reverse proxy can report 502/503/504 after Opal accepted a command but
    # before its result identifier reached the client. Only an unequivocal
    # request-size rejection is known to precede job creation. Expected remote
    # application failures were matched by their exact public contracts above.
    return(.dsvert_dsi_pre_job_rejection_is_settled(condition))
  }
  FALSE
}

.dsvert_dsi_direct_aggregate <- function(
    conns, expr, async = TRUE, error = NULL, errors.print = FALSE,
    timeout_seconds = .dsvert_dsi_job_timeout_seconds(),
    .clock = .dsvert_monotonic_seconds, .sleep = Sys.sleep,
    require_settled_sync_failure = FALSE) {
  targets <- names(conns)
  if (!is.list(conns) || !length(conns) || is.null(targets) ||
      anyNA(targets) || any(!nzchar(targets)) || anyDuplicated(targets)) {
    stop("DSI direct fan-out requires uniquely named connections",
         call. = FALSE)
  }
  if (!is.logical(async) || length(async) != 1L || is.na(async) ||
      !is.logical(errors.print) || length(errors.print) != 1L ||
      is.na(errors.print) || (!is.null(error) && !is.function(error)) ||
      !is.function(.clock) || !is.function(.sleep) ||
      !is.logical(require_settled_sync_failure) ||
      length(require_settled_sync_failure) != 1L ||
      is.na(require_settled_sync_failure)) {
    stop("Invalid DSI direct fan-out contract", call. = FALSE)
  }
  timeout_seconds <- suppressWarnings(as.numeric(timeout_seconds))
  started <- .clock()
  if (length(timeout_seconds) != 1L || is.na(timeout_seconds) ||
      timeout_seconds <= 0 || !is.numeric(started) ||
      length(started) != 1L || is.na(started) || !is.finite(started)) {
    stop("Invalid DSI direct fan-out deadline", call. = FALSE)
  }
  deadline <- if (is.finite(timeout_seconds)) {
    started + timeout_seconds
  } else {
    Inf
  }
  expressions <- if (is.list(expr) && !is.call(expr)) expr else
    stats::setNames(rep(list(expr), length(targets)), targets)
  if (!is.list(expressions) || is.null(names(expressions)) ||
      anyNA(names(expressions)) || any(!nzchar(names(expressions))) ||
      anyDuplicated(names(expressions)) ||
      !setequal(names(expressions), targets) ||
      any(!vapply(expressions, is.call, logical(1L)))) {
    stop("DSI direct fan-out expressions must map exactly to connections",
         call. = FALSE)
  }
  expressions <- expressions[targets]

  failed <- stats::setNames(rep(FALSE, length(targets)), targets)
  completed <- stats::setNames(rep(FALSE, length(targets)), targets)
  is_async <- stats::setNames(rep(FALSE, length(targets)), targets)
  jobs <- stats::setNames(vector("list", length(targets)), targets)
  values <- stats::setNames(vector("list", length(targets)), targets)
  session_keys <- stats::setNames(rep(NA_character_, length(targets)), targets)
  launched <- stats::setNames(rep(FALSE, length(targets)), targets)
  settled <- stats::setNames(rep(FALSE, length(targets)), targets)
  on.exit({
    unresolved <- targets[is_async & launched & !settled]
    if (length(unresolved)) {
      .dsvert_dsi_poison_sessions(session_keys[unresolved])
    }
  }, add = TRUE)
  report <- function(site, condition) {
    failed[[site]] <<- TRUE
    completed[[site]] <<- TRUE
    if (is.function(error)) error(site, conditionMessage(condition))
    invisible(NULL)
  }

  # Preserve DSI's lazy-session contract. A failed session is subsequently
  # surfaced by dsAggregate for that site and cannot be mistaken for success.
  DSI::datashield.sessions(
    conns, async = async,
    error = function(site, message) {
      if (length(site) && site[[1L]] %in% targets && is.function(error)) {
        error(as.character(site[[1L]]), message)
      }
      invisible(NULL)
    }, errors.print = FALSE)

  for (site in targets) {
    capability <- tryCatch(
      DSI::dsIsAsync(conns[[site]])$aggregate,
      interrupt = function(e) stop(e),
      error = function(e) {
        report(site, e)
        NULL
      })
    if (isTRUE(failed[[site]])) next
    if (!is.logical(capability) || length(capability) != 1L ||
        is.na(capability)) {
      report(site, simpleError("Connector returned an invalid async capability"))
      next
    }
    is_async[[site]] <- isTRUE(async) && isTRUE(capability)
  }

  for (site in targets[is_async & !failed]) {
    session_key <- .dsvert_dsi_job_session_key(conns[[site]])
    if (is.null(session_key)) {
      stop(
        "An asynchronous DSI connector does not expose a stable authenticated session identity; its pending jobs cannot be retried safely",
        call. = FALSE)
    }
    session_keys[[site]] <- session_key
  }
  poisoned <- targets[is_async & !failed & vapply(
    session_keys, .dsvert_dsi_session_is_poisoned, logical(1L))]
  if (length(poisoned)) {
    stop(.dsvert_dsi_poisoned_session_condition(poisoned))
  }

  # Launch all genuinely asynchronous connectors before any synchronous work.
  for (site in targets[is_async & !failed]) {
    jobs[[site]] <- tryCatch(
      DSI::dsAggregate(conns[[site]], expressions[[site]], async = TRUE),
      interrupt = function(e) {
        .dsvert_dsi_poison_sessions(session_keys[[site]])
        stop(e)
      },
      error = function(e) {
        # A connector can reject an oversized request before creating a job.
        # Only an explicitly settled provider error permits the public probe
        # ladder to descend on the same login; every ambiguous failure still
        # poisons that exact authenticated session.
        settled_rejection <- isTRUE(require_settled_sync_failure) &&
          .dsvert_dsi_pre_job_rejection_is_settled(e)
        if (!settled_rejection) {
          .dsvert_dsi_poison_sessions(session_keys[[site]])
        }
        report(site, e)
        NULL
      })
    if (!isTRUE(failed[[site]])) launched[[site]] <- TRUE
  }
  for (site in targets[!is_async & !failed]) {
    jobs[[site]] <- tryCatch(
      DSI::dsAggregate(conns[[site]], expressions[[site]], async = FALSE),
      interrupt = function(e) {
        if (isTRUE(require_settled_sync_failure)) {
          .dsvert_dsi_poison_ambiguous_sessions(conns, site)
        }
        stop(e)
      },
      error = function(e) {
        if (isTRUE(require_settled_sync_failure) &&
            !.dsvert_dsi_sync_failure_is_settled(conns[[site]], e)) {
          stop(.dsvert_dsi_poison_ambiguous_sessions(conns, site))
        }
        report(site, e)
        NULL
      })
    if (isTRUE(failed[[site]])) next
    fetched <- tryCatch(
      DSI::dsFetch(jobs[[site]]),
      interrupt = function(e) {
        if (isTRUE(require_settled_sync_failure)) {
          .dsvert_dsi_poison_ambiguous_sessions(conns, site)
        }
        stop(e)
      },
      error = function(e) {
        if (isTRUE(require_settled_sync_failure) &&
            !.dsvert_dsi_sync_failure_is_settled(conns[[site]], e)) {
          stop(.dsvert_dsi_poison_ambiguous_sessions(conns, site))
        }
        report(site, e)
        NULL
      })
    values[site] <- list(fetched)
    if (!isTRUE(failed[[site]])) completed[[site]] <- TRUE
  }

  checks <- 1L
  while (any(!completed)) {
    for (site in targets[!completed]) {
      ready <- tryCatch(
        DSI::dsIsCompleted(jobs[[site]]),
        interrupt = function(e) stop(e),
        error = function(e) {
          report(site, e)
          FALSE
        })
      if (isTRUE(failed[[site]])) next
      if (!is.logical(ready) || length(ready) != 1L || is.na(ready)) {
        report(site, simpleError("Connector returned an invalid job state"))
        next
      }
      if (isTRUE(ready)) {
        fetched <- tryCatch(
          DSI::dsFetch(jobs[[site]]),
          interrupt = function(e) stop(e),
          error = function(e) {
            if (.dsvert_dsi_sync_failure_is_settled(conns[[site]], e)) {
              settled[[site]] <<- TRUE
            }
            report(site, e)
            NULL
          })
        values[site] <- list(fetched)
        if (!isTRUE(failed[[site]])) {
          completed[[site]] <- TRUE
          settled[[site]] <- TRUE
        }
      }
    }
    if (any(!completed)) {
      now <- .clock()
      if (!is.numeric(now) || length(now) != 1L || is.na(now) ||
          !is.finite(now) || now < started) {
        stop("Monotonic DSI fan-out clock is unavailable", call. = FALSE)
      }
      if (now >= deadline) {
        .dsvert_dsi_poison_sessions(session_keys[targets[!completed]])
        for (site in targets[!completed]) {
          report(site, simpleError(
            "Connector aggregate job exceeded its monotonic deadline"))
        }
        next
      }
      # Match DSI's session-liveness behaviour for peers that completed early.
      for (site in targets[completed & !failed]) {
        tryCatch(
          DSI::dsKeepAlive(conns[[site]]),
          interrupt = function(e) stop(e),
          error = function(e) {
            if (is.function(error)) error(site, conditionMessage(e))
            invisible(NULL)
          })
      }
      .sleep(min(.dsvert_dsi_poll_seconds(checks), deadline - now))
      checks <- checks + 1L
    }
  }
  values
}

.dsvert_transport_aggregate <- function(
    .aggregate, conns, expr, async = TRUE, error = NULL,
    errors.print = FALSE, require_settled_sync_failure = FALSE) {
  if (!is.function(.aggregate)) stop(".aggregate must be a function", call. = FALSE)
  .dsvert_dsi_text_require_framed(expr)
  if (identical(.aggregate, DSI::datashield.aggregate)) {
    return(.dsvert_dsi_direct_aggregate(
      conns = conns, expr = expr, async = async, error = error,
      errors.print = errors.print,
      require_settled_sync_failure = require_settled_sync_failure))
  }
  .aggregate(
    conns = conns, expr = expr, async = async, error = error,
    errors.print = errors.print)
}

.DSVERT_IDEMPOTENT_TYPED_PRODUCERS <- c(
  "k2ShareInputDS", "glmRing63ShareExtraInputDS",
  "glmRing63ExportOwnShareDS", "k2GradientR1DS",
  "k2BeaverShareVectorDS", "k2BeaverVecmulR1DS",
  "k2ShareWeightsDS",
  "k2IknpBaseSenderChoicesDS", "k2IknpBaseReceiverEncryptDS",
  "k2IknpReceiverExtendDS", "k2IknpSenderEncryptDS",
  "mpcTypedSourceProbeDS",
  "psiPaddedRelayExchangeDS",
  "exactGCTransportInitDS", "exactGCBindPeersDS",
  "exactGCChisqProductPrepareDS", "exactGCVecmulClaimInputsDS",
  "exactGCVecmulStartDS", "exactGCVecmulValidityDS",
  "exactGCVecmulValidityReceiveDS", "exactGCVecmulCommitDS",
  "k2ChisqCrossAccumulateCountDS",
  "dsvertDPCountCompileDS",
  "dsvertDPCountAuthorizeDS",
  "dsvertDPCountStartDS",
  "dsvertDPCountFinalShareDS",
  "dsvertDPCountReleaseDS",
  "dsvertDPFrequencyClaimDS", "dsvertDPFrequencyCompileDS",
  "dsvertDPFrequencyAuthorizeDS", "dsvertDPFrequencyCleanupDS",
  "dsvertDPFrequencySourceWindowDS", "dsvertDPFrequencyFinalizeWindowDS",
  "dsvertDPFrequencyReplayDS",
  "dsvertDPCapsuleManifestDraftDS",
  "dsvertDPCapsuleManifestSignDS",
  "dsvertDPCapsuleManifestBuildDS",
  "dsvertDPCapsuleSourceTicketDS",
  "dsvertDPCapsuleSourcePrepareDS",
  "dsvertDPCapsuleSourceChunkDS",
  "dsvertDPCapsuleSourceAcceptDS",
  "dsvertDPAlignmentMaskStartDS",
  "dsvertDPAlignmentMaskStoreDS",
  "dsvertDPAlignmentMaskSealDS",
  "dsvertDPAlignmentMaskReceiveDS",
  "dsvertDPCategoricalCrossBindDS",
  "dsvertDPCategoricalCrossPrepareDS",
  "dsvertDPCategoricalCrossFinalizeDS",
  "dsvertDPGaussianCrossBindDS",
  "dsvertDPGaussianCrossPrepareDS",
  "dsvertDPGaussianCrossFinalizeDS",
  "dsvertJointDPVectorPrepareDS",
  "dsvertJointDPVectorStartDS",
  "dsvertJointDPVectorResultDS",
  "dsvertJointDPVectorFinalShareDS",
  "dsvertJointDPVectorReleaseDS",
  "dsvertJointDPVectorReplayDS",
  "dsvertJointDPVectorFinalizeAckDS",
  "dsvertDPSynopsisClaimDS", "dsvertDPSynopsisCompileDS",
  "dsvertDPSynopsisPrepareDS", "dsvertDPSynopsisStartDS",
  "dsvertDPSynopsisResultDS", "dsvertDPSynopsisFinalShareDS",
  "dsvertDPSynopsisReleaseDS",
  "dsvertDPSynopsisSourceTicketDS", "dsvertDPSynopsisSourcePrepareDS",
  "dsvertDPSynopsisSourceChunkDS", "dsvertDPSynopsisSourceAcceptDS",
  "dsvertDPSynopsisCategoricalCrossBindDS",
  "dsvertDPSynopsisCategoricalCrossFinalizeDS",
  "dsvertDPSynopsisGaussianCrossBindDS",
  "dsvertDPSynopsisGaussianCrossFinalizeDS",
  "dsvertDPSynopsisGaussianCrossEvidenceDS",
  "dsvertDPSynopsisAlignmentMaskStartDS",
  "dsvertDPSynopsisAlignmentMaskStoreDS",
  "dsvertDPSynopsisAlignmentMaskSealDS",
  "dsvertDPSynopsisAlignmentMaskReceiveDS",
  "dsvertDPSynopsisBootstrapDS", "dsvertDPSynopsisBindDS",
  "dsvertDPSynopsisPublicationDS", "dsvertDPSynopsisPublishedReplayDS",
  "dsvertDPSynopsisFinalizeAckDS")

.dsvert_dsi_call_names <- function(expr) {
  expressions <- if (is.list(expr) && !is.call(expr)) expr else list(expr)
  vapply(expressions, function(value) {
    if (!is.call(value) || length(value) < 1L) return(NA_character_)
    name <- as.character(value[[1L]])
    if (length(name) != 1L || is.na(name) || !nzchar(name)) NA_character_
    else name
  }, character(1L))
}

.dsvert_dsi_named_argument <- function(expr, name) {
  if (!is.call(expr) || !is.character(name) || length(name) != 1L ||
      is.na(name) || !nzchar(name)) return(list(present = FALSE, value = NULL))
  arguments <- as.list(expr)[-1L]
  argument_names <- names(arguments)
  if (is.null(argument_names)) return(list(present = FALSE, value = NULL))
  matches <- which(argument_names == name)
  if (length(matches) != 1L) return(list(present = FALSE, value = NULL))
  list(present = TRUE, value = arguments[[matches]])
}

.dsvert_dsi_arguments_uniquely_named <- function(expr) {
  if (!is.call(expr)) return(FALSE)
  arguments <- as.list(expr)[-1L]
  if (!length(arguments)) return(TRUE)
  argument_names <- names(arguments)
  !is.null(argument_names) && !anyNA(argument_names) &&
    all(nzchar(argument_names)) && !anyDuplicated(argument_names)
}

.dsvert_dsi_absent_null_or_empty <- function(expr, name) {
  argument <- .dsvert_dsi_named_argument(expr, name)
  if (!isTRUE(argument$present)) return(TRUE)
  identical(argument$value, quote(NULL)) ||
    (is.character(argument$value) && length(argument$value) == 1L &&
       !is.na(argument$value) && !nzchar(argument$value))
}

.dsvert_is_idempotent_typed_producer <- function(expr) {
  expressions <- if (is.list(expr) && !is.call(expr)) expr else list(expr)
  names <- .dsvert_dsi_call_names(expressions)
  if (!length(names) || anyNA(names) ||
      any(!names %in% .DSVERT_IDEMPOTENT_TYPED_PRODUCERS)) return(FALSE)

  # These producers retain transitional output-to-legacy-blob branches.  Only
  # their active-first form reaches the server-side operation replay cache;
  # classifying by function name alone could replay a mutating legacy call.
  typed_form <- Map(function(call, producer) switch(
    producer,
    k2IknpBaseReceiverEncryptDS =
      .dsvert_dsi_arguments_uniquely_named(call) &&
      .dsvert_dsi_absent_null_or_empty(call, "points_blob_key") &&
      .dsvert_dsi_absent_null_or_empty(call, "ciphertexts_blob_key"),
    k2IknpReceiverExtendDS =
      .dsvert_dsi_arguments_uniquely_named(call) &&
      .dsvert_dsi_absent_null_or_empty(call, "u_matrix_blob_key"),
    k2IknpSenderEncryptDS =
      .dsvert_dsi_arguments_uniquely_named(call) &&
      .dsvert_dsi_absent_null_or_empty(call, "u_matrix_blob_key") &&
      .dsvert_dsi_absent_null_or_empty(call, "ciphertexts_blob_key"),
    TRUE), expressions, names)
  all(unlist(typed_form, use.names = FALSE))
}

.dsvert_has_typed_remote_rejection <- function(result) {
  if (!is.list(result)) return(FALSE)
  required <- c("version", "operation", "rejected")
  if (!is.null(names(result)) && all(required %in% names(result)) &&
      identical(result$version, "dsvert-typed-blob-rejection-v1") &&
      identical(result$rejected, TRUE)) return(TRUE)
  manifest_required <- c("version", "phase", "reason_code", "rejected")
  if (!is.null(names(result)) &&
      all(manifest_required %in% names(result)) &&
      identical(
        result$version,
        "dsvert-biomedical-capsule-manifest-rejection-v1") &&
      identical(result$rejected, TRUE)) return(TRUE)
  any(vapply(result, .dsvert_has_typed_remote_rejection, logical(1L)))
}

#' Strict DataSHIELD aggregate result contract (internal)
#'
#' @param conns Named DataSHIELD connection list.
#' @param expr One aggregate call or a named list of calls.
#' @param operation Public, non-sensitive phase label used in errors.
#' @param allow_null Named sites whose documented protocol contract permits a
#'   NULL result. Empty by default.
#' @param result_contract Either \code{"non_null"} or
#'   \code{"logical_true"}.
#' @param .aggregate Aggregate implementation, injectable for tests.
#' @param .expression_capacity_bytes Optional scalar or per-site named
#'   expression ceilings in bytes, applied after DSV1 framing. \code{NULL}
#'   uses the current negotiated ceiling.
#' @return A result list ordered exactly like \code{conns}.
#' @keywords internal
.dsvert_aggregate_strict <- function(
    conns, expr, operation = "DataSHIELD aggregate",
    allow_null = character(),
    result_contract = c("non_null", "logical_true"),
    .aggregate = DSI::datashield.aggregate,
    .expression_capacity_bytes = NULL) {
  result_contract <- match.arg(result_contract)
  targets <- names(conns)
  if (!is.list(conns) || length(conns) < 1L || is.null(targets) ||
      anyNA(targets) || any(!nzchar(targets)) || anyDuplicated(targets)) {
    stop("conns must be a non-empty list with unique logical site names",
         call. = FALSE)
  }
  if (!is.character(operation) || length(operation) != 1L ||
      is.na(operation) || !nzchar(operation)) {
    stop("operation must be one non-empty public phase label", call. = FALSE)
  }
  if (!is.character(allow_null) || anyNA(allow_null) ||
      anyDuplicated(allow_null) || !all(allow_null %in% targets)) {
    stop("allow_null must contain unique names from conns", call. = FALSE)
  }
  if (!is.function(.aggregate)) {
    stop(".aggregate must be a function", call. = FALSE)
  }
  if (is.list(expr) && !is.call(expr)) {
    expr_names <- names(expr)
    if (is.null(expr_names) || anyNA(expr_names) ||
        any(!nzchar(expr_names)) || anyDuplicated(expr_names) ||
        !setequal(expr_names, targets)) {
      stop("per-site expressions must be named exactly once for every connection",
           call. = FALSE)
    }
    expr <- expr[targets]
  }
  expr <- .dsvert_dsi_text_frame_expressions(expr)
  if (is.null(.expression_capacity_bytes)) {
    .dsvert_validate_dsi_expression_sizes(expr)
  } else {
    .dsvert_validate_dsi_expression_sizes(
      expr, capacity_bytes = .expression_capacity_bytes)
  }
  .dsvert_validate_real_dsi_transport(conns, .aggregate)

  idempotent_producer <- .dsvert_is_idempotent_typed_producer(expr)
  attempt <- function() {
    callback_failed <- character()
    callback_not_phase <- character()
    callback_not_lifetime <- character()
    peer_rejections <- list()
    instance_retries <- list()
    resource_backpressure <- list()
    resource_oversize <- list()
    phase_not_ready <- list()
    lifetime_exhausted <- list()
    lifetime_token_conflict <- FALSE
    error_callback <- function(site, message) {
      site <- if (length(site)) as.character(site[[1L]]) else "unknown"
      callback_failed <<- union(callback_failed, site)
      phase <- .dsvert_client_parse_phase_not_ready(message)
      if (is.null(phase)) {
        callback_not_phase <<- union(callback_not_phase, site)
      } else {
        phase_not_ready[[site]] <<- phase
      }
      lifetime <-
        .dsvert_client_parse_dp_lifetime_budget_exhausted(message)
      if (is.null(lifetime)) {
        callback_not_lifetime <<- union(callback_not_lifetime, site)
      } else {
        lifetime_exhausted[[site]] <<- lifetime
        lifetime_token_conflict <<-
          lifetime_token_conflict ||
          .dsvert_client_has_other_bracketed_typed_token(
            message, conditionMessage(lifetime))
      }
      rejection <- .dsvert_client_parse_peer_not_recognized(message)
      if (!is.null(rejection)) peer_rejections[[site]] <<- rejection
      retry <- .dsvert_client_parse_release_instance_retry(message)
      if (!is.null(retry)) instance_retries[[site]] <<- retry
      backpressure <- .dsvert_client_parse_resource_backpressure(message)
      if (!is.null(backpressure)) {
        resource_backpressure[[site]] <<- backpressure
      }
      oversize <- .dsvert_client_parse_resource_oversize(message)
      if (!is.null(oversize)) resource_oversize[[site]] <<- oversize
      invisible(NULL)
    }
    top_error <- FALSE
    top_rejection <- NULL
    top_instance_retry <- NULL
    top_resource_backpressure <- NULL
    top_resource_oversize <- NULL
    top_phase_not_ready <- NULL
    top_lifetime_exhausted <- NULL
    top_bracketed_typed_token <- FALSE
    result <- tryCatch(
      .dsvert_transport_aggregate(
        .aggregate = .aggregate, conns = conns, expr = expr,
        error = error_callback, async = TRUE, errors.print = FALSE),
      interrupt = function(e) stop(e),
      error = function(e) {
        if (inherits(e, "dsvert_dsi_poisoned_session")) stop(e)
        top_rejection <<- if (inherits(e, "dsvert_peer_not_recognized")) {
          e
        } else {
          .dsvert_client_parse_peer_not_recognized(conditionMessage(e))
        }
        top_instance_retry <<-
          .dsvert_client_parse_release_instance_retry(conditionMessage(e))
        top_resource_backpressure <<-
          .dsvert_client_parse_resource_backpressure(conditionMessage(e))
        top_resource_oversize <<-
          .dsvert_client_parse_resource_oversize(conditionMessage(e))
        top_phase_not_ready <<-
          .dsvert_client_parse_phase_not_ready(conditionMessage(e))
        top_lifetime_exhausted <<-
          .dsvert_client_parse_dp_lifetime_budget_exhausted(
            conditionMessage(e))
        top_bracketed_typed_token <<- grepl(
          "\\[dsvert_[^]\\r\\n]+\\]", conditionMessage(e), perl = TRUE)
        if (inherits(
              top_lifetime_exhausted,
              "dsvert_dp_lifetime_budget_exhausted")) {
          lifetime_token_conflict <<-
            lifetime_token_conflict ||
            .dsvert_client_has_other_bracketed_typed_token(
              conditionMessage(e),
              conditionMessage(top_lifetime_exhausted))
        }
        top_error <<- TRUE
        NULL
      })
    list(
      result = result, callback_failed = callback_failed,
      callback_not_phase = callback_not_phase,
      callback_not_lifetime = callback_not_lifetime,
      top_error = top_error,
      phase_not_ready = phase_not_ready,
      top_phase_not_ready = top_phase_not_ready,
      lifetime_exhausted = lifetime_exhausted,
      top_lifetime_exhausted = top_lifetime_exhausted,
      top_bracketed_typed_token = top_bracketed_typed_token,
      lifetime_token_conflict = lifetime_token_conflict,
      peer_rejections = c(peer_rejections,
        if (is.null(top_rejection)) list() else list(top = top_rejection)),
      instance_retries = c(instance_retries,
        if (is.null(top_instance_retry)) list() else
          list(top = top_instance_retry)),
      resource_backpressure = c(resource_backpressure,
        if (is.null(top_resource_backpressure)) list() else
          list(top = top_resource_backpressure)),
      resource_oversize = c(resource_oversize,
        if (is.null(top_resource_oversize)) list() else
          list(top = top_resource_oversize)))
  }
  classify <- function(value) {
    phase_seen <- length(value$phase_not_ready) > 0L ||
      inherits(value$top_phase_not_ready, "dsvert_phase_not_ready")
    lifetime_seen <- length(value$lifetime_exhausted) > 0L ||
      inherits(
        value$top_lifetime_exhausted,
        "dsvert_dp_lifetime_budget_exhausted")
    other_typed <- length(value$peer_rejections) > 0L ||
      length(value$instance_retries) > 0L ||
      length(value$resource_backpressure) > 0L ||
      length(value$resource_oversize) > 0L
    if ((isTRUE(phase_seen) &&
         (isTRUE(lifetime_seen) || isTRUE(other_typed))) ||
        (isTRUE(lifetime_seen) &&
         (isTRUE(other_typed) || isTRUE(value$lifetime_token_conflict)))) {
      return(list(state = "fatal", result = value$result))
    }
    if (length(value$peer_rejections)) {
      return(list(
        state = "peer_not_recognized", result = value$result,
        rejection = value$peer_rejections[[order(
          names(value$peer_rejections), method = "radix")[[1L]]]]))
    }
    if (length(value$resource_oversize)) {
      terminal <- value$resource_oversize[
        order(names(value$resource_oversize), method = "radix")]
      return(list(
        state = "resource_oversize", result = value$result,
        rejection = terminal[[1L]]))
    }
    callback_retry <- length(value$callback_failed) > 0L &&
      setequal(value$callback_failed, names(value$instance_retries))
    top_retry <- isTRUE(value$top_error) &&
      "top" %in% names(value$instance_retries) &&
      (length(value$callback_failed) == 0L || isTRUE(callback_retry))
    if (isTRUE(callback_retry) || isTRUE(top_retry)) {
      retries <- value$instance_retries[
        order(names(value$instance_retries), method = "radix")]
      return(list(
        state = "release_instance_retry", result = value$result,
        rejection = retries[[1L]]))
    }
    callback_phase <- length(value$callback_failed) > 0L &&
      !length(value$callback_not_phase) &&
      all(value$callback_failed %in% targets) &&
      setequal(value$callback_failed, names(value$phase_not_ready))
    top_phase <- isTRUE(value$top_error) &&
      inherits(value$top_phase_not_ready, "dsvert_phase_not_ready")
    named_result <- is.list(value$result) &&
      length(value$result) == length(targets) &&
      !is.null(names(value$result)) && !anyDuplicated(names(value$result)) &&
      setequal(names(value$result), targets)
    unreported_null <- if (isTRUE(named_result)) {
      null_sites <- targets[vapply(
        value$result[targets], is.null, logical(1L))]
      setdiff(null_sites, union(value$callback_failed, allow_null))
    } else {
      targets
    }
    phase_only <- (isTRUE(callback_phase) || isTRUE(top_phase)) &&
      !isTRUE(other_typed) &&
      (!length(value$callback_failed) || isTRUE(callback_phase)) &&
      (!isTRUE(value$top_error) || isTRUE(top_phase)) &&
      (isTRUE(value$top_error) ||
       (isTRUE(named_result) && !length(unreported_null)))
    if (isTRUE(phase_only)) {
      terminal <- if (length(value$phase_not_ready)) {
        value$phase_not_ready[
          order(names(value$phase_not_ready), method = "radix")][[1L]]
      } else {
        value$top_phase_not_ready
      }
      return(list(
        state = "phase_not_ready", result = value$result,
        rejection = terminal))
    }
    callback_lifetime <- length(value$callback_failed) > 0L &&
      !length(value$callback_not_lifetime) &&
      all(value$callback_failed %in% targets) &&
      setequal(value$callback_failed, names(value$lifetime_exhausted))
    top_lifetime <- isTRUE(value$top_error) &&
      inherits(
        value$top_lifetime_exhausted,
        "dsvert_dp_lifetime_budget_exhausted")
    lifetime_only <- (isTRUE(callback_lifetime) || isTRUE(top_lifetime)) &&
      !isTRUE(phase_seen) && !isTRUE(other_typed) &&
      (!length(value$callback_failed) || isTRUE(callback_lifetime)) &&
      (!isTRUE(value$top_error) || isTRUE(top_lifetime)) &&
      (isTRUE(value$top_error) ||
       (isTRUE(named_result) && !length(unreported_null)))
    if (isTRUE(lifetime_only)) {
      terminal <- if (length(value$lifetime_exhausted)) {
        value$lifetime_exhausted[
          order(names(value$lifetime_exhausted), method = "radix")][[1L]]
      } else {
        value$top_lifetime_exhausted
      }
      return(list(
        state = "dp_lifetime_budget_exhausted", result = value$result,
        rejection = terminal))
    }
    callback_backpressure <- length(value$callback_failed) > 0L &&
      setequal(value$callback_failed,
               setdiff(names(value$resource_backpressure), "top"))
    top_backpressure <- isTRUE(value$top_error) &&
      "top" %in% names(value$resource_backpressure) &&
      (length(value$callback_failed) == 0L ||
       isTRUE(callback_backpressure))
    if (isTRUE(idempotent_producer) &&
        (isTRUE(callback_backpressure) || isTRUE(top_backpressure))) {
      return(list(state = "missing", result = value$result))
    }
    # A per-site callback is a present negative outcome. It can represent a
    # server contract/signature error, so it is never converted into retry.
    if (length(value$callback_failed) > 0L) {
      return(list(state = "fatal", result = value$result))
    }
    if (isTRUE(value$top_error) &&
        isTRUE(value$top_bracketed_typed_token)) {
      return(list(state = "fatal", result = value$result))
    }
    if (isTRUE(value$top_error) || is.null(value$result)) {
      return(list(
        state = if (idempotent_producer) "missing" else "fatal",
        result = value$result))
    }
    if (.dsvert_has_typed_remote_rejection(value$result)) {
      return(list(state = "rejected", result = value$result))
    }
    if (isTRUE(idempotent_producer) && named_result &&
        any(vapply(value$result, is.null, logical(1L)))) {
      return(list(state = "missing", result = value$result))
    }
    list(state = "ack", result = value$result)
  }
  outcome <- if (isTRUE(idempotent_producer)) {
    .dsvert_retry_idempotent(
      attempt = attempt, classify = classify,
      operation = paste0(
        "idempotent typed producer ",
        paste(unique(.dsvert_dsi_call_names(expr)), collapse = "+")))
  } else {
    classify(attempt())
  }
  if (!identical(outcome$state, "ack")) {
    if (identical(outcome$state, "peer_not_recognized") &&
        inherits(outcome$rejection, "dsvert_peer_not_recognized")) {
      stop(outcome$rejection)
    }
    if (identical(outcome$state, "release_instance_retry") &&
        inherits(outcome$rejection, "dsvert_release_instance_retry")) {
      stop(outcome$rejection)
    }
    if (identical(outcome$state, "resource_oversize") &&
        inherits(outcome$rejection, "dsvert_resource_oversize")) {
      stop(outcome$rejection)
    }
    if (identical(outcome$state, "phase_not_ready") &&
        inherits(outcome$rejection, "dsvert_phase_not_ready")) {
      stop(outcome$rejection)
    }
    if (identical(outcome$state, "dp_lifetime_budget_exhausted") &&
        inherits(
          outcome$rejection, "dsvert_dp_lifetime_budget_exhausted")) {
      stop(outcome$rejection)
    }
    stop("DataSHIELD transport failed during '", operation,
         "'; partial or invalid site results were rejected.",
         call. = FALSE)
  }
  result <- outcome$result

  names_valid <- is.list(result) && length(result) == length(targets) &&
    !is.null(names(result)) && !anyDuplicated(names(result)) &&
    setequal(names(result), targets)
  if (names_valid) result <- result[targets]

  null_sites <- if (names_valid) {
    targets[vapply(result, is.null, logical(1L))]
  } else {
    targets
  }
  unexpected_null <- setdiff(null_sites, allow_null)
  contract_failed <- FALSE
  if (names_valid && identical(result_contract, "logical_true")) {
    required <- targets[!vapply(result, is.null, logical(1L))]
    contract_failed <- any(!vapply(
      result[required], identical, logical(1L), TRUE))
  }

  if (!names_valid || length(unexpected_null) > 0L || contract_failed) {
    stop("DataSHIELD transport failed during '", operation,
         "'; partial or invalid site results were rejected.",
         call. = FALSE)
  }

  result
}

#' One idempotent DSI fan-out cycle with an explicit availability state
#'
#' This helper does not retry and never advances protocol state.  It is for
#' offset-addressed, idempotent exchanges whose caller can replay the identical
#' request.  DSI 1.8 represents a missing site result as a named NULL; that and
#' a top-level transport exception are the only non-fatal unavailable states.
#' A present but misassociated outer result is a protocol failure.
#'
#' @keywords internal
.dsvert_fanout_cycle <- function(
    conns, expressions, operation = "DataSHIELD exchange",
    .aggregate = DSI::datashield.aggregate) {
  targets <- names(conns)
  if (!is.list(conns) || length(conns) < 1L || is.null(targets) ||
      anyNA(targets) || any(!nzchar(targets)) || anyDuplicated(targets)) {
    stop("conns must be a non-empty list with unique logical site names",
         call. = FALSE)
  }
  if (!is.list(expressions) || is.null(names(expressions)) ||
      anyNA(names(expressions)) || any(!nzchar(names(expressions))) ||
      anyDuplicated(names(expressions)) ||
      !setequal(names(expressions), targets)) {
    stop("expressions must be named exactly once for every connection",
         call. = FALSE)
  }
  expressions <- expressions[targets]
  if (any(!vapply(expressions, is.call, logical(1L)))) {
    stop("every per-site expression must be an R call", call. = FALSE)
  }
  if (!is.character(operation) || length(operation) != 1L ||
      is.na(operation) || !nzchar(operation)) {
    stop("operation must be one non-empty public phase label", call. = FALSE)
  }
  if (!is.function(.aggregate)) stop(".aggregate must be a function", call. = FALSE)
  .dsvert_validate_dsi_expression_sizes(expressions)
  .dsvert_validate_real_dsi_transport(conns, .aggregate)

  callback_failed <- FALSE
  callback_backpressure <- FALSE
  callback_contract_failure <- FALSE
  peer_rejection <- NULL
  resource_oversize <- NULL
  error_callback <- function(site, message) {
    callback_failed <<- TRUE
    parsed <- .dsvert_client_parse_peer_not_recognized(message)
    if (is.null(peer_rejection) && !is.null(parsed)) {
      peer_rejection <<- parsed
    }
    resource <- .dsvert_client_parse_resource_backpressure(message)
    oversized <- .dsvert_client_parse_resource_oversize(message)
    if (is.null(resource_oversize) && !is.null(oversized)) {
      resource_oversize <<- oversized
    }
    if (!is.null(resource)) {
      callback_backpressure <<- TRUE
    } else if (is.null(parsed) && is.null(oversized)) {
      callback_contract_failure <<- TRUE
    }
    invisible(NULL)
  }
  top_error <- FALSE
  top_resource_oversize <- NULL
  result <- tryCatch(
    .dsvert_transport_aggregate(
      .aggregate = .aggregate,
      conns = conns, expr = expressions, async = TRUE,
      error = error_callback, errors.print = FALSE),
    interrupt = function(e) stop(e),
    error = function(e) {
      if (inherits(e, "dsvert_dsi_poisoned_session")) stop(e)
      candidate <- if (inherits(e, "dsvert_peer_not_recognized")) {
        e
      } else {
        .dsvert_client_parse_peer_not_recognized(conditionMessage(e))
      }
      if (!inherits(peer_rejection, "dsvert_peer_not_recognized") &&
          inherits(candidate, "dsvert_peer_not_recognized")) {
        peer_rejection <<- candidate
      }
      top_resource_oversize <<-
        .dsvert_client_parse_resource_oversize(conditionMessage(e))
      top_error <<- TRUE
      NULL
    })
  if (inherits(peer_rejection, "dsvert_peer_not_recognized")) {
    stop(peer_rejection)
  }
  if (inherits(resource_oversize, "dsvert_resource_oversize")) {
    stop(resource_oversize)
  }
  if (inherits(top_resource_oversize, "dsvert_resource_oversize")) {
    stop(top_resource_oversize)
  }
  if (isTRUE(callback_contract_failure)) {
    stop("DataSHIELD exchange returned a server contract failure during '",
         operation, "'; it was not retried.", call. = FALSE)
  }
  if (isTRUE(top_error) || is.null(result)) {
    return(list(state = "unavailable", result = NULL))
  }
  names_valid <- is.list(result) && length(result) == length(targets) &&
    !is.null(names(result)) && !anyNA(names(result)) &&
    !anyDuplicated(names(result)) && setequal(names(result), targets)
  if (!names_valid) {
    stop("DataSHIELD transport returned a malformed or misassociated result during '",
         operation, "'.", call. = FALSE)
  }
  result <- result[targets]
  if (isTRUE(callback_backpressure) || isTRUE(callback_failed) ||
      any(vapply(result, is.null, logical(1L)))) {
    return(list(state = "unavailable", result = NULL))
  }
  list(state = "ok", result = result)
}

#' One-request named per-site DSI fan-out (internal)
#'
#' DSI 1.8 maps a named expression list to connections by logical site name
#' and submits async-capable connectors before fetching their results. This
#' helper validates the map before invoking that path and then applies the
#' strict result contract above.
#'
#' @param conns Named DataSHIELD connection list.
#' @param expressions Calls named exactly once for every connection.
#' @param operation Public, non-sensitive phase label used in errors.
#' @param allow_null Named sites whose protocol permits a NULL result.
#' @param result_contract Either \code{"non_null"} or
#'   \code{"logical_true"}.
#' @param .aggregate Aggregate implementation, injectable for tests.
#' @param .expression_capacity_bytes Optional scalar or per-site named
#'   expression ceilings in bytes, applied after DSV1 framing. \code{NULL}
#'   uses the current negotiated ceiling.
#' @keywords internal
.dsvert_fanout_by_site <- function(
    conns, expressions, operation = "DataSHIELD fan-out",
    allow_null = character(),
    result_contract = c("non_null", "logical_true"),
    .aggregate = DSI::datashield.aggregate,
    .expression_capacity_bytes = NULL) {
  targets <- names(conns)
  if (!is.list(expressions) || is.null(names(expressions)) ||
      anyNA(names(expressions)) || any(!nzchar(names(expressions))) ||
      anyDuplicated(names(expressions)) ||
      !setequal(names(expressions), targets)) {
    stop("expressions must be named exactly once for every connection",
         call. = FALSE)
  }
  expressions <- expressions[targets]
  if (any(!vapply(expressions, is.call, logical(1L)))) {
    stop("every per-site expression must be an R call", call. = FALSE)
  }
  .dsvert_aggregate_strict(
    conns = conns, expr = expressions, operation = operation,
    allow_null = allow_null, result_contract = result_contract,
    .aggregate = .aggregate,
    .expression_capacity_bytes = .expression_capacity_bytes)
}

#' Fail-closed DataSHIELD assignment fan-out (internal)
#'
#' DSI assignments do not return a typed acknowledgement.  This wrapper uses
#' the connector error callback as the completion contract, suppresses remote
#' error text, and rejects any partial fan-out.  A named expression list is
#' mapped to sites in one asynchronous DSI phase.
#' @keywords internal
.dsvert_assign_strict <- function(
    conns, symbol, values, operation = "DataSHIELD assignment",
    .assign = DSI::datashield.assign.expr) {
  targets <- names(conns)
  if (!is.list(conns) || length(conns) < 1L || is.null(targets) ||
      anyNA(targets) || any(!nzchar(targets)) || anyDuplicated(targets)) {
    stop("conns must be a non-empty list with unique logical site names",
         call. = FALSE)
  }
  if (!is.character(symbol) || length(symbol) != 1L || is.na(symbol) ||
      !nzchar(symbol)) {
    stop("symbol must be one non-empty string", call. = FALSE)
  }
  if (!is.character(operation) || length(operation) != 1L ||
      is.na(operation) || !nzchar(operation)) {
    stop("operation must be one non-empty public phase label", call. = FALSE)
  }
  expressions <- if (is.list(values) && !is.call(values)) values else {
    stats::setNames(rep(list(values), length(targets)), targets)
  }
  if (is.null(names(expressions)) || anyNA(names(expressions)) ||
      any(!nzchar(names(expressions))) || anyDuplicated(names(expressions)) ||
      !setequal(names(expressions), targets) ||
      any(!vapply(expressions, is.call, logical(1L)))) {
    stop("values must contain exactly one R call per connection", call. = FALSE)
  }
  expressions <- expressions[targets]
  if (!is.function(.assign)) stop(".assign must be a function", call. = FALSE)
  .dsvert_validate_dsi_expression_sizes(expressions)
  real_dsi <- identical(.assign, DSI::datashield.assign.expr)
  async_sites <- character()
  if (isTRUE(real_dsi)) {
    .dsvert_require_dsi_fanout_version()
    .dsvert_validate_dsi_transport_security(conns)
    async_capability <- vapply(targets, function(site) {
      value <- tryCatch(
        DSI::dsIsAsync(conns[[site]])$assignExpr,
        error = function(error) NULL)
      if (!is.logical(value) || length(value) != 1L || is.na(value)) {
        stop("A DSI connector returned an invalid assignment capability",
             call. = FALSE)
      }
      isTRUE(value)
    }, logical(1L))
    async_sites <- targets[async_capability]
    session_keys <- stats::setNames(lapply(async_sites, function(site) {
      .dsvert_dsi_job_session_key(conns[[site]])
    }), async_sites)
    stable <- vapply(session_keys, function(key) {
      is.character(key) && length(key) == 1L && !is.na(key) && nzchar(key)
    }, logical(1L))
    if (any(!stable)) {
      stop(
        "An asynchronous DSI connector does not expose a stable authenticated session identity; its pending assignments cannot be retried safely",
        call. = FALSE)
    }
    poisoned <- async_sites[vapply(
      session_keys, .dsvert_dsi_session_is_poisoned, logical(1L))]
    if (length(poisoned)) {
      stop(.dsvert_dsi_poisoned_session_condition(poisoned))
    }
  }

  callback_failed <- character()
  callback_succeeded <- character()
  callback_contract_failed <- FALSE
  peer_rejection <- NULL
  resource_oversize <- NULL
  success_callback <- function(site) {
    site <- if (length(site) == 1L && !is.na(site)) {
      as.character(site[[1L]])
    } else {
      ""
    }
    if (!nzchar(site) || !site %in% targets ||
        site %in% callback_succeeded) {
      callback_contract_failed <<- TRUE
    } else {
      callback_succeeded <<- c(callback_succeeded, site)
    }
    invisible(NULL)
  }
  error_callback <- function(site, message) {
    site <- if (length(site) == 1L && !is.na(site)) {
      as.character(site[[1L]])
    } else {
      ""
    }
    if (!nzchar(site) || !site %in% targets) {
      callback_contract_failed <<- TRUE
    } else {
      callback_failed <<- union(callback_failed, site)
    }
    parsed <- .dsvert_client_parse_peer_not_recognized(message)
    if (is.null(peer_rejection) && !is.null(parsed)) {
      peer_rejection <<- parsed
    }
    oversized <- .dsvert_client_parse_resource_oversize(message)
    if (is.null(resource_oversize) && !is.null(oversized)) {
      resource_oversize <<- oversized
    }
    invisible(NULL)
  }
  old_print <- getOption("datashield.errors.print", NULL)
  options(datashield.errors.print = FALSE)
  on.exit(options(datashield.errors.print = old_print), add = TRUE)
  top_error <- FALSE
  top_resource_oversize <- NULL
  tryCatch(
    .assign(
      conns = conns, symbol = symbol, expr = expressions,
      async = TRUE, success = success_callback,
      error = error_callback, errors.print = FALSE),
    interrupt = function(e) {
      if (isTRUE(real_dsi) && length(async_sites)) {
        .dsvert_dsi_poison_ambiguous_sessions(conns, async_sites)
      }
      stop(e)
    },
    error = function(e) {
      candidate <- if (inherits(e, "dsvert_peer_not_recognized")) {
        e
      } else {
        .dsvert_client_parse_peer_not_recognized(conditionMessage(e))
      }
      if (!inherits(peer_rejection, "dsvert_peer_not_recognized") &&
          inherits(candidate, "dsvert_peer_not_recognized")) {
        peer_rejection <<- candidate
      }
      top_resource_oversize <<-
        .dsvert_client_parse_resource_oversize(conditionMessage(e))
      top_error <<- TRUE
      invisible(NULL)
    })
  poison_condition <- NULL
  if (isTRUE(real_dsi)) {
    unconfirmed <- setdiff(targets, callback_succeeded)
    ambiguous <- if (isTRUE(callback_contract_failed) || isTRUE(top_error)) {
      async_sites
    } else {
      intersect(async_sites, union(callback_failed, unconfirmed))
    }
    if (length(ambiguous)) {
      poison_condition <- .dsvert_dsi_poison_ambiguous_sessions(
        conns, ambiguous)
    }
  }
  if (inherits(peer_rejection, "dsvert_peer_not_recognized")) {
    stop(peer_rejection)
  }
  if (inherits(resource_oversize, "dsvert_resource_oversize")) {
    stop(resource_oversize)
  }
  if (inherits(top_resource_oversize, "dsvert_resource_oversize")) {
    stop(top_resource_oversize)
  }
  if (!is.null(poison_condition)) stop(poison_condition)
  if (isTRUE(top_error) || length(callback_failed) ||
      isTRUE(callback_contract_failed) ||
      !setequal(callback_succeeded, targets)) {
    stop("DataSHIELD transport failed during '", operation,
         "'; partial assignment was rejected.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Final best-effort cleanup without exposing connector details (internal)
#'
#' Cleanup is the sole phase allowed to ignore missing replies: it runs only
#' from an on-exit handler and cannot advance a statistical protocol.
#' @keywords internal
.dsvert_cleanup_best_effort <- function(
    conns, expressions, .aggregate = DSI::datashield.aggregate) {
  targets <- names(conns)
  if (!is.list(conns) || !length(conns) || is.null(targets) ||
      anyNA(targets) || any(!nzchar(targets)) || anyDuplicated(targets) ||
      !is.function(.aggregate)) return(invisible(FALSE))
  expressions <- if (is.list(expressions) && !is.call(expressions)) {
    expressions
  } else {
    stats::setNames(rep(list(expressions), length(targets)), targets)
  }
  if (is.null(names(expressions)) || !setequal(names(expressions), targets) ||
      any(!vapply(expressions, is.call, logical(1L)))) {
    return(invisible(FALSE))
  }
  expressions <- expressions[targets]
  expressions <- .dsvert_dsi_text_frame_expressions(expressions)
  if (inherits(try(.dsvert_validate_dsi_expression_sizes(expressions),
                   silent = TRUE), "try-error")) return(invisible(FALSE))
  if (inherits(try(.dsvert_validate_real_dsi_transport(conns, .aggregate),
                   silent = TRUE), "try-error")) return(invisible(FALSE))
  tryCatch({
    .dsvert_transport_aggregate(
      .aggregate = .aggregate, conns = conns, expr = expressions, async = TRUE,
      error = function(site, message) invisible(NULL), errors.print = FALSE)
    invisible(TRUE)
  }, error = function(e) invisible(FALSE))
}
