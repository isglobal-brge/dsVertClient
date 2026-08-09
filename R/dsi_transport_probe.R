# Data-free DSI request-size negotiation.
#
# Candidate probes contain only public ASCII padding. They run before the first
# opaque payload and may safely fall back after a missing/transport response
# because the server endpoint is stateless. Once a payload begins, geometry is
# locked and only byte-identical replay is permitted.

.DSVERT_DSI_PROBE_VERSION <- "dsvert-transport-probe-v1"
.DSVERT_DSI_RESPONSE_PROBE_VERSION <-
  "dsvert-transport-response-probe-v1"
.DSVERT_DSI_PROBE_CANDIDATES <- as.numeric(c(
  8L * 1024L^2, 4L * 1024L^2, 2L * 1024L^2,
  1L * 1024L^2, 640L * 1024L, 320L * 1024L,
  160L * 1024L, 80L * 1024L, 32L * 1024L, 16L * 1024L))
.DSVERT_DSI_PROBE_ABSOLUTE_MAX <- 8L * 1024L^2
.DSVERT_DSI_PROBE_MIN <- 16L * 1024L
.DSVERT_DSI_PROBE_EXPRESSION_RESERVE <- 64L * 1024L
.DSVERT_DSI_RESPONSE_PROBE_HEADROOM <- 128L * 1024L
.DSVERT_DSI_RESPONSE_PROBE_USABLE_NUMERATOR <- 7L
.DSVERT_DSI_RESPONSE_PROBE_USABLE_DENOMINATOR <- 8L
.DSVERT_DSI_PORTABLE_CHUNK_CHARS <- 640L * 1024L
.DSVERT_DSI_PORTABLE_EXPRESSION_BYTES <- 768L * 1024L - 1L
.DSVERT_DSI_PROBE_CACHE_ENTRIES <- 64L
.DSVERT_DSI_KNOWN_DSLITE_VERSION <- "1.4.1"

.dsvert_dsi_probe_nonce <- function() {
  bytes <- openssl::rand_bytes(16L)
  paste0("tp_", paste(sprintf("%02x", as.integer(bytes)), collapse = ""))
}

.dsvert_dsi_probe_hash <- function(value) {
  # `as.character.hash()` retains openssl's class attribute; tickets and DSI
  # responses are plain JSON/R strings, so normalise to an unclassed scalar.
  paste0(openssl::sha256(charToRaw(value)))
}

.dsvert_dsi_cache_scalar <- function(value) {
  if (is.character(value) && length(value) == 1L && !is.na(value) &&
      nzchar(value)) return(value)
  if (is.numeric(value) && length(value) == 1L && !is.na(value) &&
      is.finite(value)) return(format(value, scientific = FALSE, trim = TRUE))
  NULL
}

.dsvert_dsi_connection_session <- function(connection) {
  if (inherits(connection, "ArmadilloConnection")) {
    # A bearer-token refresh is not a new connector/login lifecycle, while a
    # fresh dsConnect() creates a new libcurl handle (and may select a different
    # Armadillo profile).  The authoritative probe cache must therefore follow
    # the exact in-process handle, not the reusable/refreshable credential.
    return(.dsvert_dsi_armadillo_handle_identity(connection))
  }
  direct <- .dsvert_dsi_cache_scalar(
    .dsvert_dsi_member(connection, "sid"))
  if (!is.null(direct)) return(paste0("sid:", direct))
  opal <- .dsvert_dsi_member(connection, "opal")
  rid <- .dsvert_dsi_cache_scalar(.dsvert_dsi_member(opal, "rid"))
  if (!is.null(rid)) return(paste0("rid:", rid))
  handle <- .dsvert_dsi_member(connection, "handle")
  for (field in c("session_id", "session", "sid", "rid", "id")) {
    value <- .dsvert_dsi_cache_scalar(.dsvert_dsi_member(handle, field))
    if (!is.null(value)) return(paste0(field, ":", value))
  }
  NULL
}

.dsvert_dsi_probe_cache_key <- function(conns) {
  sites <- names(conns)
  if (!is.list(conns) || !length(conns) || is.null(sites) ||
      anyNA(sites) || any(!nzchar(sites)) || anyDuplicated(sites)) {
    stop("Probe connections must have unique logical site names",
         call. = FALSE)
  }
  sites <- sort(sites)
  components <- vapply(sites, function(site) {
    connection <- conns[[site]]
    details <- .dsvert_dsi_inspect_connection(connection)
    session <- .dsvert_dsi_connection_session(connection)
    if (is.null(session)) return(NA_character_)
    endpoint <- if (is.null(details$endpoint)) {
      "in-process"
    } else {
      details$endpoint$scope
    }
    fields <- c(site, paste(class(connection), collapse = "/"),
                details$kind, endpoint, session)
    paste0(nchar(fields, type = "bytes"), ":", fields, collapse = "")
  }, character(1L), USE.NAMES = FALSE)
  if (anyNA(components)) return(NULL)
  .dsvert_dsi_probe_hash(paste0(
    "dsVert/dsi-probe-cache/v1|", paste0(components, collapse = "")))
}

.dsvert_clear_transport_probe_cache <- function() {
  .dsvert_chunk_env$probe_cache <- list()
  .dsvert_chunk_env$probe_hint_cache <- list()
  .dsvert_chunk_env$probe_clock <- 0
  .dsvert_reset_chunk_size()
  invisible(TRUE)
}

# A session cache is authoritative only for that exact authenticated handle.
# This coarser profile cache is merely a probe-order hint: a new session still
# has to accept one public-padding request before the geometry is used. A stale
# or colliding hint can therefore cost one failed stateless probe, but can never
# authorize an untested payload size.
.dsvert_dsi_probe_profile_key <- function(conns) {
  sites <- names(conns)
  if (!is.list(conns) || !length(conns) || is.null(sites) ||
      anyNA(sites) || any(!nzchar(sites)) || anyDuplicated(sites)) {
    return(NULL)
  }
  sites <- sort(sites)
  components <- vapply(sites, function(site) {
    connection <- conns[[site]]
    details <- .dsvert_dsi_inspect_connection(connection)
    if (!details$kind %in% c("dslite", "opal", "armadillo")) {
      return(NA_character_)
    }
    endpoint <- if (identical(details$kind, "dslite")) {
      version <- tryCatch(
        as.character(utils::packageVersion("DSLite")),
        error = function(e) NA_character_)
      if (length(version) != 1L || is.na(version)) return(NA_character_)
      paste0("in-process/", version)
    } else {
      if (is.null(details$endpoint) ||
          !identical(details$endpoint$scheme, "https")) {
        return(NA_character_)
      }
      details$endpoint$scope
    }
    fields <- c(site, paste(class(connection), collapse = "/"),
                details$kind, endpoint)
    paste0(nchar(fields, type = "bytes"), ":", fields, collapse = "")
  }, character(1L), USE.NAMES = FALSE)
  if (anyNA(components)) return(NULL)
  .dsvert_dsi_probe_hash(paste0(
    "dsVert/dsi-probe-profile/v1|", paste0(components, collapse = "")))
}

.dsvert_dsi_probe_hint_get <- function(key) {
  if (is.null(key)) return(NULL)
  entry <- .dsvert_chunk_env$probe_hint_cache[[key]]
  if (is.null(entry)) return(NULL)
  .dsvert_chunk_env$probe_clock <- .dsvert_chunk_env$probe_clock + 1
  entry$touched <- .dsvert_chunk_env$probe_clock
  .dsvert_chunk_env$probe_hint_cache[[key]] <- entry
  entry
}

.dsvert_dsi_probe_hint_put <- function(key, entry) {
  if (is.null(key)) return(invisible(FALSE))
  cache <- .dsvert_chunk_env$probe_hint_cache
  if (is.null(cache[[key]]) &&
      length(cache) >= .DSVERT_DSI_PROBE_CACHE_ENTRIES) {
    touched <- vapply(cache, `[[`, numeric(1L), "touched")
    cache[[names(cache)[[which.min(touched)]]]] <- NULL
  }
  .dsvert_chunk_env$probe_clock <- .dsvert_chunk_env$probe_clock + 1
  entry$touched <- .dsvert_chunk_env$probe_clock
  cache[[key]] <- entry
  .dsvert_chunk_env$probe_hint_cache <- cache
  invisible(TRUE)
}

.dsvert_dsi_known_probe_hint <- function(
    conns, candidates, .package_version = utils::packageVersion) {
  dslite <- vapply(conns, inherits, logical(1L), "DSLiteConnection")
  if (!any(dslite)) return(NULL)
  version <- tryCatch(
    as.character(.package_version("DSLite")), error = function(e) NA_character_)
  # The recorded parser sweep is specific to this exact connector release.
  # Unknown/newer versions are probed normally instead of inheriting a brittle
  # class-wide ceiling.
  if (length(version) != 1L || is.na(version) ||
      !identical(version, .DSVERT_DSI_KNOWN_DSLITE_VERSION) ||
      !.DSVERT_DSI_PORTABLE_CHUNK_CHARS %in% candidates) return(NULL)
  as.numeric(.DSVERT_DSI_PORTABLE_CHUNK_CHARS)
}

.dsvert_dsi_prioritized_probe_candidates <- function(
    conns, candidates, profile_hint = NULL) {
  hint <- if (is.list(profile_hint)) profile_hint$chunk_chars else NULL
  if (is.null(hint)) hint <- .dsvert_dsi_known_probe_hint(conns, candidates)
  hint <- suppressWarnings(as.numeric(hint))
  if (length(hint) != 1L || is.na(hint) || !is.finite(hint) ||
      !hint %in% candidates) return(candidates)
  c(hint, candidates[candidates < hint])
}

.dsvert_dsi_probe_cache_get <- function(key) {
  if (is.null(key)) return(NULL)
  entry <- .dsvert_chunk_env$probe_cache[[key]]
  if (is.null(entry)) return(NULL)
  .dsvert_chunk_env$probe_clock <- .dsvert_chunk_env$probe_clock + 1
  entry$touched <- .dsvert_chunk_env$probe_clock
  .dsvert_chunk_env$probe_cache[[key]] <- entry
  entry
}

.dsvert_dsi_probe_cache_put <- function(key, entry) {
  if (is.null(key)) return(invisible(FALSE))
  cache <- .dsvert_chunk_env$probe_cache
  if (is.null(cache[[key]]) &&
      length(cache) >= .DSVERT_DSI_PROBE_CACHE_ENTRIES) {
    touched <- vapply(cache, `[[`, numeric(1L), "touched")
    cache[[names(cache)[[which.min(touched)]]]] <- NULL
  }
  .dsvert_chunk_env$probe_clock <- .dsvert_chunk_env$probe_clock + 1
  entry$touched <- .dsvert_chunk_env$probe_clock
  cache[[key]] <- entry
  .dsvert_chunk_env$probe_cache <- cache
  invisible(TRUE)
}

.dsvert_apply_dsi_probe_geometry <- function(entry, cache_key) {
  .dsvert_chunk_env$effective_chunk_size <- as.integer(entry$chunk_chars)
  .dsvert_chunk_env$effective_expression_bytes <-
    as.numeric(entry$max_expression_bytes)
  .dsvert_chunk_env$probed <- identical(entry$source, "probe")
  .dsvert_chunk_env$active_probe_cache_key <- cache_key
  invisible(entry)
}

.dsvert_validate_dsi_probe_ack <- function(value, nonce, padding_chars,
                                            padding_sha256) {
  required <- c(
    "version", "nonce", "padding_chars", "padding_sha256",
    "server_max_padding_chars")
  if (!is.list(value) || is.null(names(value)) || anyDuplicated(names(value)) ||
      !identical(sort(names(value)), sort(required)) ||
      !identical(value$version, .DSVERT_DSI_PROBE_VERSION) ||
      !identical(value$nonce, nonce) ||
      !identical(value$padding_sha256, padding_sha256)) return(FALSE)
  observed <- suppressWarnings(as.numeric(value$padding_chars))
  server_max <- suppressWarnings(as.numeric(value$server_max_padding_chars))
  length(observed) == 1L && !is.na(observed) && is.finite(observed) &&
    observed == floor(observed) && observed == padding_chars &&
    length(server_max) == 1L && !is.na(server_max) && is.finite(server_max) &&
    server_max == floor(server_max) && server_max >= padding_chars &&
    server_max <= .DSVERT_DSI_PROBE_ABSOLUTE_MAX
}

.dsvert_dsi_response_probe_usable_bytes <- function(padding_chars) {
  padding_chars <- suppressWarnings(as.numeric(padding_chars))
  if (length(padding_chars) != 1L || is.na(padding_chars) ||
      !is.finite(padding_chars) || padding_chars != floor(padding_chars) ||
      padding_chars < .DSVERT_DSI_PROBE_MIN ||
      padding_chars > .DSVERT_DSI_PROBE_ABSOLUTE_MAX) {
    return(NA_real_)
  }
  max(0, floor(
    padding_chars * .DSVERT_DSI_RESPONSE_PROBE_USABLE_NUMERATOR /
      .DSVERT_DSI_RESPONSE_PROBE_USABLE_DENOMINATOR) -
      .DSVERT_DSI_RESPONSE_PROBE_HEADROOM)
}

.dsvert_validate_dsi_response_probe_ack <- function(
    value, nonce, request_padding, request_sha256, response_padding_chars) {
  required <- c(
    "version", "nonce", "padding_chars", "padding_sha256",
    "server_max_padding_chars", "response_padding_chars",
    "response_padding_sha256", "server_max_response_padding_chars",
    "response_padding")
  if (!is.list(value) || is.null(names(value)) || anyNA(names(value)) ||
      anyDuplicated(names(value)) ||
      !identical(sort(names(value)), sort(required)) ||
      !identical(value$version, .DSVERT_DSI_RESPONSE_PROBE_VERSION) ||
      !identical(value$nonce, nonce) ||
      !identical(value$padding_sha256, request_sha256) ||
      !is.character(value$response_padding) ||
      length(value$response_padding) != 1L ||
      is.na(value$response_padding)) return(FALSE)
  request_observed <- suppressWarnings(as.numeric(value$padding_chars))
  request_max <- suppressWarnings(as.numeric(value$server_max_padding_chars))
  response_observed <- suppressWarnings(
    as.numeric(value$response_padding_chars))
  response_max <- suppressWarnings(
    as.numeric(value$server_max_response_padding_chars))
  valid_numbers <- all(vapply(
    list(request_observed, request_max, response_observed, response_max),
    function(item) length(item) == 1L && !is.na(item) && is.finite(item) &&
      item == floor(item), logical(1L)))
  if (!isTRUE(valid_numbers) || request_observed != nchar(
        request_padding, type = "bytes") || request_max < request_observed ||
      request_max > .DSVERT_DSI_PROBE_ABSOLUTE_MAX ||
      response_observed != response_padding_chars ||
      response_max < response_observed ||
      response_max > .DSVERT_DSI_PROBE_ABSOLUTE_MAX ||
      nchar(value$response_padding, type = "bytes") != response_observed ||
      !grepl("^R+\\z", value$response_padding, perl = TRUE,
             useBytes = TRUE)) return(FALSE)
  identical(
    .dsvert_dsi_probe_hash(value$response_padding),
    value$response_padding_sha256)
}

.dsvert_dsi_probe_once <- function(conns, padding_chars, .aggregate) {
  padding <- strrep("A", padding_chars)
  padding_sha256 <- .dsvert_dsi_probe_hash(padding)
  nonce <- .dsvert_dsi_probe_nonce()
  expression <- call(
    name = "dsvertTransportProbeDS", nonce = nonce, padding = padding,
    padding_sha256 = padding_sha256)
  expressions <- stats::setNames(
    rep(list(expression), length(conns)), names(conns))
  expression_sizes <- vapply(expressions, function(value) {
    nchar(paste(deparse(value, width.cutoff = 500L), collapse = "\n"),
          type = "bytes")
  }, numeric(1L))
  if (any(expression_sizes >
          padding_chars + .DSVERT_DSI_PROBE_EXPRESSION_RESERVE)) {
    stop("Transport-probe expression overhead exceeds its fixed reserve",
         call. = FALSE)
  }

  callback_failed <- character()
  top_failed <- FALSE
  response <- tryCatch(
    .dsvert_transport_aggregate(
      .aggregate = .aggregate,
      # A deliberately oversized public probe may be rejected after an async
      # connector has accepted a job.  DSI has no portable cancellation or
      # fetch-resume primitive, so that ambiguity would poison the login and
      # prevent the safe candidate ladder from descending.  Run only this
      # stateless, data-free negotiation call synchronously; all statistical
      # and MPC phases retain named asynchronous fan-out.
      conns = conns, expr = expressions, async = FALSE,
      error = function(site, message) {
        callback_failed <<- union(
          callback_failed,
          if (length(site)) as.character(site[[1L]]) else "unknown")
        invisible(NULL)
      }, errors.print = FALSE,
      require_settled_sync_failure = TRUE),
    error = function(e) {
      if (inherits(e, "dsvert_dsi_poisoned_session")) stop(e)
      top_failed <<- TRUE
      NULL
    })
  if (is.null(response) ||
      (is.list(response) && length(response) == 0L)) {
    return(list(
      success = FALSE, successful_sites = character(),
      expression_bytes = expression_sizes))
  }
  targets <- names(conns)
  response_names <- names(response)
  if (!is.list(response) || is.null(response_names) ||
      anyNA(response_names) || any(!nzchar(response_names)) ||
      anyDuplicated(response_names) ||
      any(!response_names %in% targets)) {
    stop("Transport probe returned malformed target routing", call. = FALSE)
  }
  present <- response_names[!vapply(response, is.null, logical(1L))]
  for (site in present) {
    if (!.dsvert_validate_dsi_probe_ack(
          response[[site]], nonce, padding_chars, padding_sha256)) {
      stop("Transport probe returned a malformed acknowledgement for target '",
           site, "'", call. = FALSE)
    }
  }
  success <- !isTRUE(top_failed) && !length(callback_failed) &&
    length(response) == length(targets) && setequal(response_names, targets) &&
    length(present) == length(targets)
  list(
    success = success, successful_sites = present,
    expression_bytes = expression_sizes)
}

.dsvert_dsi_response_probe_failure_kind <- function(message) {
  if (!is.character(message) || !length(message) || is.na(message[[1L]])) {
    return("other")
  }
  message <- substr(as.character(message[[1L]]), 1L, 8192L)
  if (grepl(
      paste0(
        "unused argument[^\\n]*response_padding_chars|",
        "response_padding_chars[^\\n]*unused argument"),
      message, ignore.case = TRUE, perl = TRUE)) return("unsupported")
  if (!is.null(.dsvert_client_parse_resource_oversize(message)) ||
      grepl(
        paste0(
          "(^|[^0-9])413([^0-9]|$)|payload too large|",
          "request entity too large|response[^\\n]{0,64}too large|",
          "transport-response-probe padding exceeds"),
        message, ignore.case = TRUE, perl = TRUE)) return("oversize")
  "other"
}

.dsvert_dsi_response_probe_once <- function(
    conns, response_padding_chars, .aggregate) {
  request_padding <- "A"
  request_sha256 <- .dsvert_dsi_probe_hash(request_padding)
  nonce <- .dsvert_dsi_probe_nonce()
  expression <- call(
    name = "dsvertTransportProbeDS", nonce = nonce,
    padding = request_padding, padding_sha256 = request_sha256,
    response_padding_chars = as.numeric(response_padding_chars))
  expressions <- stats::setNames(
    rep(list(expression), length(conns)), names(conns))
  .dsvert_validate_dsi_expression_sizes(expressions)

  targets <- names(conns)
  failure_kind <- stats::setNames(rep("", length(targets)), targets)
  callback <- function(site, message) {
    site <- if (length(site)) as.character(site[[1L]]) else ""
    if (site %in% targets) {
      failure_kind[[site]] <<-
        .dsvert_dsi_response_probe_failure_kind(message)
    }
    invisible(NULL)
  }
  result <- tryCatch(
    .dsvert_transport_aggregate(
      .aggregate = .aggregate, conns = conns, expr = expressions,
      async = FALSE, error = callback, errors.print = FALSE,
      require_settled_sync_failure = TRUE),
    interrupt = function(error) {
      .dsvert_dsi_poison_ambiguous_sessions(conns)
      stop(error)
    },
    error = function(error) {
      if (inherits(error, "dsvert_dsi_poisoned_session")) stop(error)
      kind <- .dsvert_dsi_response_probe_failure_kind(
        conditionMessage(error))
      if (kind %in% c("unsupported", "oversize")) {
        failure_kind[] <<- kind
        return(NULL)
      }
      stop(.dsvert_dsi_poison_ambiguous_sessions(conns))
    })
  if (!is.null(result) && (!is.list(result) || is.null(names(result)) ||
      anyNA(names(result)) || any(!nzchar(names(result))) ||
      anyDuplicated(names(result)) || any(!names(result) %in% targets))) {
    stop("Transport response probe returned malformed target routing",
         call. = FALSE)
  }
  successful <- character()
  if (is.list(result)) {
    present <- names(result)[!vapply(result, is.null, logical(1L))]
    for (site in present) {
      if (!.dsvert_validate_dsi_response_probe_ack(
            result[[site]], nonce, request_padding, request_sha256,
            response_padding_chars)) {
        stop("Transport response probe returned a malformed acknowledgement for target '",
             site, "'", call. = FALSE)
      }
      successful <- c(successful, site)
    }
  }
  unresolved <- setdiff(targets, successful)
  if (length(unresolved) && any(!failure_kind[unresolved] %in%
                               c("unsupported", "oversize"))) {
    # A synchronous remote failure may be settled yet unrelated to size. It
    # cannot authorize a response geometry, and silently reusing the handle
    # would be unsafe for singleton-result connectors.
    stop(.dsvert_dsi_poison_ambiguous_sessions(conns, unresolved))
  }
  list(
    successful_sites = successful,
    unsupported_sites = names(failure_kind)[failure_kind == "unsupported"],
    oversized_sites = names(failure_kind)[failure_kind == "oversize"])
}

.dsvert_negotiate_dsi_response_sizes <- function(
    conns, .aggregate, candidates) {
  targets <- names(conns)
  selected <- stats::setNames(rep(NA_real_, length(targets)), targets)
  supported <- stats::setNames(rep(TRUE, length(targets)), targets)
  remaining <- targets
  for (candidate in candidates) {
    if (!length(remaining)) break
    attempt <- .dsvert_dsi_response_probe_once(
      conns[remaining], candidate, .aggregate)
    accepted <- intersect(remaining, attempt$successful_sites)
    selected[accepted] <- as.numeric(candidate)
    unsupported <- intersect(remaining, attempt$unsupported_sites)
    supported[unsupported] <- FALSE
    remaining <- setdiff(remaining, c(accepted, unsupported))
  }
  list(
    response_padding_chars = selected,
    response_usable_bytes = stats::setNames(vapply(
      selected, .dsvert_dsi_response_probe_usable_bytes, numeric(1L)),
      targets),
    response_probe_supported = supported)
}

#' Negotiate one common request geometry before opaque payload transfer
#' @keywords internal
.dsvert_negotiate_dsi_chunk_size <- function(
    conns, .aggregate = DSI::datashield.aggregate,
    .candidates = .DSVERT_DSI_PROBE_CANDIDATES) {
  if (isTRUE(.dsvert_chunk_env$geometry_locked)) {
    stop("DSI payload geometry cannot be renegotiated after transfer begins",
         call. = FALSE)
  }
  if (!is.numeric(.candidates) || !length(.candidates) ||
      anyNA(.candidates) || any(!is.finite(.candidates)) ||
      any(.candidates != floor(.candidates)) || any(.candidates < 1) ||
      any(.candidates > .DSVERT_DSI_PROBE_ABSOLUTE_MAX) ||
      is.unsorted(-.candidates, strictly = TRUE)) {
    stop("Invalid descending DSI transport-probe candidates",
         call. = FALSE)
  }
  .dsvert_validate_real_dsi_transport(conns, .aggregate)
  response_candidates <- .candidates[
    .candidates >= .DSVERT_DSI_PROBE_MIN]
  cache_key <- .dsvert_dsi_probe_cache_key(conns)
  cached <- .dsvert_dsi_probe_cache_get(cache_key)
  if (!is.null(cached)) {
    return(.dsvert_apply_dsi_probe_geometry(cached, cache_key))
  }

  profile_key <- .dsvert_dsi_probe_profile_key(conns)
  profile_hint <- .dsvert_dsi_probe_hint_get(profile_key)
  candidates <- .dsvert_dsi_prioritized_probe_candidates(
    conns, .candidates, profile_hint)
  selected <- NULL
  site_request_padding <- stats::setNames(
    rep(NA_real_, length(conns)), names(conns))
  site_expression_bytes <- site_request_padding
  for (candidate in candidates) {
    attempt <- .dsvert_dsi_probe_once(conns, candidate, .aggregate)
    first_success <- attempt$successful_sites[
      is.na(site_request_padding[attempt$successful_sites])]
    if (length(first_success)) {
      site_request_padding[first_success] <- as.numeric(candidate)
      site_expression_bytes[first_success] <- as.numeric(
        candidate + .DSVERT_DSI_PROBE_EXPRESSION_RESERVE)
    }
    if (isTRUE(attempt$success)) {
      selected <- list(
        chunk_chars = as.numeric(candidate),
        max_expression_bytes = as.numeric(
          candidate + .DSVERT_DSI_PROBE_EXPRESSION_RESERVE),
        site_request_padding_chars = site_request_padding,
        site_max_expression_bytes = site_expression_bytes,
        source = "probe")
      break
    }
  }
  if (is.null(selected)) {
    stop(
      "The outer DSI channel rejected every data-free transport probe; ",
      "refusing to send an opaque payload at an unverified size.",
      call. = FALSE)
  }
  response_limits <- if (length(response_candidates)) {
    .dsvert_negotiate_dsi_response_sizes(
      conns, .aggregate, response_candidates)
  } else {
    list(
      response_padding_chars = stats::setNames(
        rep(NA_real_, length(conns)), names(conns)),
      response_usable_bytes = stats::setNames(
        rep(NA_real_, length(conns)), names(conns)),
      response_probe_supported = stats::setNames(
        rep(FALSE, length(conns)), names(conns)))
  }
  selected <- c(selected, response_limits)
  # Opal may create its lazy DSI session during the first probe. Recompute the
  # key so a cached geometry is never attached to a pre-session descriptor.
  final_key <- .dsvert_dsi_probe_cache_key(conns)
  .dsvert_dsi_probe_cache_put(final_key, selected)
  if (identical(selected$source, "probe")) {
    .dsvert_dsi_probe_hint_put(profile_key, selected)
  }
  .dsvert_apply_dsi_probe_geometry(selected, final_key)
}

.dsvert_maybe_negotiate_dsi_chunk_size <- function(
    conns, .aggregate = DSI::datashield.aggregate) {
  cache_key <- .dsvert_dsi_probe_cache_key(conns)
  if (!is.null(.dsvert_chunk_env$effective_chunk_size) &&
      identical(.dsvert_chunk_env$active_probe_cache_key, cache_key)) {
    return(invisible(NULL))
  }
  if (isTRUE(.dsvert_chunk_env$geometry_locked)) {
    stop("DSI payload geometry cannot change during an active transfer",
         call. = FALSE)
  }
  if (!is.null(.dsvert_chunk_env$effective_chunk_size) ||
      !identical(.dsvert_chunk_env$active_probe_cache_key, cache_key)) {
    .dsvert_reset_chunk_size()
  }
  # Dependency-injected aggregators are test/protocol harnesses. They opt into
  # probing explicitly; production DSI takes the automatic data-free path.
  if (!identical(.aggregate, DSI::datashield.aggregate)) {
    # Even when a harness deliberately skips probing, bind its configured
    # geometry to the exact connection/session set.  A subsequent cohort must
    # not inherit this value merely because it runs in the same R process.
    .dsvert_chunk_env$active_probe_cache_key <- cache_key
    return(invisible(NULL))
  }
  .dsvert_negotiate_dsi_chunk_size(conns, .aggregate)
  invisible(NULL)
}

.dsvert_dsi_transport_site_limits <- function(conns, .aggregate) {
  targets <- names(conns)
  if (!is.list(conns) || !length(conns) || is.null(targets) ||
      anyNA(targets) || any(!nzchar(targets)) || anyDuplicated(targets) ||
      !is.function(.aggregate)) {
    stop("Invalid connections for per-site DSI transport limits",
         call. = FALSE)
  }
  if (!identical(.aggregate, DSI::datashield.aggregate)) {
    # Dependency-injected protocol harnesses do not traverse an outer DSI
    # path. Give them the fixed application maximum; their aggregate callback
    # remains authoritative and all production connections take the verified
    # cache branch below.
    maximum <- as.numeric(.DSVERT_DSI_PROBE_ABSOLUTE_MAX)
    return(list(
      request_payload_bytes = stats::setNames(
        rep(maximum, length(targets)), targets),
      expression_bytes = stats::setNames(
        rep(maximum + .DSVERT_DSI_PROBE_EXPRESSION_RESERVE,
            length(targets)), targets),
      response_bytes = stats::setNames(rep(
        .dsvert_dsi_response_probe_usable_bytes(maximum), length(targets)),
        targets),
      response_probe_supported = stats::setNames(
        rep(TRUE, length(targets)), targets)))
  }
  key <- .dsvert_dsi_probe_cache_key(conns)
  entry <- if (is.null(key)) NULL else
    .dsvert_chunk_env$probe_cache[[key]]
  valid_map <- function(value, allow_na = FALSE) {
    is.numeric(value) && !is.null(names(value)) &&
      !anyNA(names(value)) && !anyDuplicated(names(value)) &&
      setequal(names(value), targets) &&
      (allow_na || !anyNA(value)) &&
      all(is.na(value) | (is.finite(value) & value >= 0 &
                            value == floor(value)))
  }
  if (is.null(entry) ||
      !valid_map(entry$site_request_padding_chars) ||
      !valid_map(entry$site_max_expression_bytes) ||
      !valid_map(entry$response_usable_bytes, allow_na = TRUE) ||
      !is.logical(entry$response_probe_supported) ||
      is.null(names(entry$response_probe_supported)) ||
      !setequal(names(entry$response_probe_supported), targets) ||
      anyNA(entry$response_probe_supported)) {
    return(list(
      request_payload_bytes = stats::setNames(
        rep(NA_real_, length(targets)), targets),
      expression_bytes = stats::setNames(
        rep(NA_real_, length(targets)), targets),
      response_bytes = stats::setNames(
        rep(NA_real_, length(targets)), targets),
      response_probe_supported = stats::setNames(
        rep(FALSE, length(targets)), targets)))
  }
  list(
    request_payload_bytes = entry$site_request_padding_chars[targets],
    expression_bytes = entry$site_max_expression_bytes[targets],
    response_bytes = entry$response_usable_bytes[targets],
    response_probe_supported = entry$response_probe_supported[targets])
}
