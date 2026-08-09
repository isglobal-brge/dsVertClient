probe_connection <- function(sid = "session-1", name = "site") {
  structure(list(sid = sid, name = name),
            class = c("DSLiteConnection", "list"))
}

probe_ack <- function(expr, server_max = 8 * 1024^2) {
  fields <- as.list(expr)
  list(
    version = "dsvert-transport-probe-v1",
    nonce = fields$nonce,
    padding_chars = as.numeric(nchar(fields$padding, type = "bytes")),
    padding_sha256 = fields$padding_sha256,
    server_max_padding_chars = as.numeric(server_max))
}

response_probe_ack <- function(expr, server_max = 8 * 1024^2) {
  fields <- as.list(expr)
  response_padding <- strrep("R", fields$response_padding_chars)
  list(
    version = "dsvert-transport-response-probe-v1",
    nonce = fields$nonce,
    padding_chars = as.numeric(nchar(fields$padding, type = "bytes")),
    padding_sha256 = fields$padding_sha256,
    server_max_padding_chars = as.numeric(server_max),
    response_padding_chars = as.numeric(fields$response_padding_chars),
    response_padding_sha256 = .dsvert_dsi_probe_hash(response_padding),
    server_max_response_padding_chars = as.numeric(server_max),
    response_padding = response_padding)
}

test_that("probe descends on absent transport responses and selects common min", {
  .dsvert_clear_transport_probe_cache()
  on.exit(.dsvert_clear_transport_probe_cache(), add = TRUE)
  conns <- list(
    alpha = probe_connection("session-a", "alpha"),
    beta = probe_connection("session-b", "beta"))
  attempted <- numeric()
  aggregate <- function(conns, expr, async, error, errors.print) {
    size <- nchar(as.list(expr[[1L]])$padding, type = "bytes")
    attempted <<- c(attempted, size)
    answer <- lapply(expr, probe_ack)
    if (size > 1024L) answer$beta <- NULL
    answer
  }

  selected <- .dsvert_negotiate_dsi_chunk_size(
    conns, aggregate, .candidates = c(4096, 2048, 1024))
  expect_identical(attempted, c(4096, 2048, 1024))
  expect_identical(selected$source, "probe")
  expect_identical(selected$chunk_chars, 1024)
  expect_identical(.dsvert_get_chunk_size(), 1024L)
  expect_identical(
    .dsvert_chunk_env$effective_expression_bytes,
    as.numeric(1024L + .DSVERT_DSI_PROBE_EXPRESSION_RESERVE))
})

test_that("completed remote probe rejection descends without poisoning login", {
  .dsvert_clear_transport_probe_cache()
  .dsvert_dsi_clear_poisoned_sessions()
  on.exit({
    .dsvert_clear_transport_probe_cache()
    .dsvert_dsi_clear_poisoned_sessions()
  }, add = TRUE)
  conns <- list(site = probe_connection("remote-probe-session"))
  attempts <- numeric()
  async_flags <- logical()
  observed_fields <- list()

  selected <- testthat::with_mocked_bindings(
    .dsvert_negotiate_dsi_chunk_size(
      conns, DSI::datashield.aggregate,
      .candidates = c(4096, 2048, 1024)),
    datashield.sessions = function(...) invisible(NULL),
    dsIsAsync = function(...) list(aggregate = TRUE),
    dsAggregate = function(conn, expr, async) {
      fields <- as.list(expr)
      size <- nchar(fields$padding, type = "bytes")
      attempts <<- c(attempts, size)
      async_flags <<- c(async_flags, async)
      observed_fields[[length(observed_fields) + 1L]] <<- names(fields)
      structure(list(expr = expr, size = size), class = "probe_job")
    },
    dsFetch = function(job) {
      if (job$size > 1024L) {
        stop("HTTP 413 / completed parser rejection")
      }
      probe_ack(job$expr)
    },
    dsIsCompleted = function(...) stop("synchronous probe must not poll"),
    dsKeepAlive = function(...) stop("synchronous probe must not keep alive"),
    .package = "DSI")

  expect_identical(attempts, c(4096, 2048, 1024))
  expect_identical(async_flags, rep(FALSE, 3L))
  expect_true(all(vapply(
    observed_fields, identical, logical(1L),
    c("", "nonce", "padding", "padding_sha256"))))
  expect_identical(selected$chunk_chars, 1024)
  expect_false(.dsvert_dsi_session_is_poisoned(
    .dsvert_dsi_job_session_key(conns[[1L]])))
})

test_that("probe never downgrades a poisoned session to a smaller candidate", {
  .dsvert_clear_transport_probe_cache()
  on.exit(.dsvert_clear_transport_probe_cache(), add = TRUE)
  poisoned <- .dsvert_dsi_poisoned_session_condition("site")
  calls <- 0L
  condition <- tryCatch(
    .dsvert_negotiate_dsi_chunk_size(
      list(site = probe_connection()),
      .aggregate = function(...) {
        calls <<- calls + 1L
        stop(poisoned)
      },
      .candidates = c(4096, 2048, 1024)),
    dsvert_dsi_poisoned_session = identity)
  expect_identical(condition, poisoned)
  expect_identical(calls, 1L)
})

test_that("probe cache is bound to the exact connection session set", {
  .dsvert_clear_transport_probe_cache()
  on.exit(.dsvert_clear_transport_probe_cache(), add = TRUE)
  calls <- 0L
  aggregate <- function(conns, expr, async, error, errors.print) {
    calls <<- calls + 1L
    lapply(expr, probe_ack)
  }
  conns <- list(site = probe_connection("session-1"))
  first <- .dsvert_negotiate_dsi_chunk_size(
    conns, aggregate, .candidates = c(4096, 2048))
  expect_identical(first$chunk_chars, 4096)
  expect_identical(calls, 1L)

  .dsvert_reset_chunk_size()
  cached <- .dsvert_negotiate_dsi_chunk_size(
    conns, aggregate, .candidates = c(4096, 2048))
  expect_identical(cached$chunk_chars, 4096)
  expect_identical(calls, 1L)

  .dsvert_reset_chunk_size()
  changed <- list(site = probe_connection("session-2"))
  .dsvert_negotiate_dsi_chunk_size(
    changed, aggregate, .candidates = c(4096, 2048))
  expect_identical(calls, 2L)
})

test_that("response ceilings are negotiated per site with 413 fallback", {
  .dsvert_clear_transport_probe_cache()
  on.exit(.dsvert_clear_transport_probe_cache(), add = TRUE)
  conns <- list(
    alpha = probe_connection("response-a", "alpha"),
    beta = probe_connection("response-b", "beta"),
    gamma = probe_connection("response-c", "gamma"))
  ceilings <- c(alpha = 64L, beta = 32L, gamma = 16L) * 1024L
  response_attempts <- list()
  aggregate <- function(conns, expr, async, error, errors.print) {
    fields <- as.list(expr[[1L]])[-1L]
    if (is.null(fields$response_padding_chars)) {
      return(lapply(expr, probe_ack))
    }
    candidate <- as.numeric(fields$response_padding_chars)
    response_attempts[[length(response_attempts) + 1L]] <<-
      list(candidate = candidate, sites = names(conns))
    stats::setNames(lapply(names(conns), function(site) {
      if (candidate > ceilings[[site]]) {
        error(site, "HTTP 413 Payload Too Large")
        return(NULL)
      }
      response_probe_ack(expr[[site]], server_max = ceilings[[site]])
    }), names(conns))
  }

  selected <- .dsvert_negotiate_dsi_chunk_size(
    conns, aggregate, .candidates = c(64, 32, 16) * 1024L)
  expect_identical(
    unname(selected$response_padding_chars), as.numeric(ceilings))
  expect_identical(names(selected$response_padding_chars), names(ceilings))
  expect_true(all(selected$response_probe_supported))
  expect_identical(
    lapply(response_attempts, `[[`, "sites"),
    list(c("alpha", "beta", "gamma"), c("beta", "gamma"), "gamma"))
  expect_identical(
    vapply(response_attempts, `[[`, numeric(1L), "candidate"),
    as.numeric(c(64, 32, 16) * 1024L))

  attempt_count <- length(response_attempts)
  .dsvert_reset_chunk_size()
  cached <- .dsvert_negotiate_dsi_chunk_size(
    conns, aggregate, .candidates = c(64, 32, 16) * 1024L)
  expect_identical(cached$response_padding_chars,
                   selected$response_padding_chars)
  expect_identical(length(response_attempts), attempt_count)

  .dsvert_reset_chunk_size()
  ceilings[["beta"]] <- 16L * 1024L
  reconnected <- lapply(conns, function(value) {
    value$sid <- paste0(value$sid, "-fresh")
    value
  })
  fresh <- .dsvert_negotiate_dsi_chunk_size(
    reconnected, aggregate, .candidates = c(64, 32, 16) * 1024L)
  expect_identical(
    unname(fresh$response_padding_chars), as.numeric(ceilings))
  expect_gt(length(response_attempts), attempt_count)
})

test_that("request ceilings retain heterogeneous per-site evidence", {
  .dsvert_clear_transport_probe_cache()
  on.exit(.dsvert_clear_transport_probe_cache(), add = TRUE)
  conns <- list(
    alpha = probe_connection("request-a", "alpha"),
    beta = probe_connection("request-b", "beta"),
    gamma = probe_connection("request-c", "gamma"))
  ceilings <- c(alpha = 64L, beta = 32L, gamma = 16L) * 1024L
  aggregate <- function(conns, expr, async, error, errors.print) {
    fields <- as.list(expr[[1L]])[-1L]
    if (!is.null(fields$response_padding_chars)) {
      for (site in names(conns)) {
        error(site, "unused argument (response_padding_chars = 65536)")
      }
      return(stats::setNames(vector("list", length(conns)), names(conns)))
    }
    size <- nchar(fields$padding, type = "bytes")
    stats::setNames(lapply(names(conns), function(site) {
      if (size > ceilings[[site]]) {
        error(site, "HTTP 413 Payload Too Large")
        return(NULL)
      }
      probe_ack(expr[[site]], server_max = ceilings[[site]])
    }), names(conns))
  }

  selected <- .dsvert_negotiate_dsi_chunk_size(
    conns, aggregate, .candidates = c(64, 32, 16) * 1024L)
  expect_identical(selected$chunk_chars, 16 * 1024)
  expect_identical(
    unname(selected$site_request_padding_chars), as.numeric(ceilings))
  expect_identical(
    unname(selected$site_max_expression_bytes),
    as.numeric(ceilings + .DSVERT_DSI_PROBE_EXPRESSION_RESERVE))
})

test_that("request hints do not truncate the independent response ladder", {
  .dsvert_clear_transport_probe_cache()
  on.exit(.dsvert_clear_transport_probe_cache(), add = TRUE)
  conns <- list(site = probe_connection("asymmetric-channel"))
  request_attempts <- numeric()
  response_attempts <- numeric()
  aggregate <- function(conns, expr, async, error, errors.print) {
    fields <- as.list(expr[[1L]])[-1L]
    if (is.null(fields$response_padding_chars)) {
      size <- nchar(fields$padding, type = "bytes")
      request_attempts <<- c(request_attempts, size)
      if (size > 640L * 1024L) {
        error("site", "HTTP 413 Payload Too Large")
        return(list(site = NULL))
      }
      return(list(site = probe_ack(
        expr[[1L]], server_max = 640L * 1024L)))
    }
    size <- as.numeric(fields$response_padding_chars)
    response_attempts <<- c(response_attempts, size)
    acknowledgement <- response_probe_ack(
      expr[[1L]], server_max = 8L * 1024L^2)
    acknowledgement$server_max_padding_chars <- 640L * 1024L
    list(site = acknowledgement)
  }

  selected <- testthat::with_mocked_bindings(
    .dsvert_negotiate_dsi_chunk_size(
      conns, aggregate,
      .candidates = as.numeric(c(
        4L * 1024L^2, 640L * 1024L, 320L * 1024L, 16L * 1024L))),
    .dsvert_dsi_known_probe_hint = function(...) 640L * 1024L,
    .package = "dsVertClient")

  expect_identical(request_attempts, 640 * 1024)
  expect_identical(response_attempts, 4 * 1024^2)
  expect_identical(selected$chunk_chars, 640 * 1024)
  expect_identical(
    selected$response_padding_chars[["site"]], 4 * 1024^2)
})

test_that("response acknowledgements reject a trailing newline", {
  response_padding_chars <- 16L * 1024L
  request_padding <- "A"
  request_sha256 <- .dsvert_dsi_probe_hash(request_padding)
  nonce <- paste0("tp_", strrep("1", 32L))
  expression <- call(
    name = "dsvertTransportProbeDS", nonce = nonce,
    padding = request_padding, padding_sha256 = request_sha256,
    response_padding_chars = response_padding_chars)
  acknowledgement <- response_probe_ack(expression)
  acknowledgement$response_padding <- paste0(
    strrep("R", response_padding_chars - 1L), "\n")
  acknowledgement$response_padding_sha256 <- .dsvert_dsi_probe_hash(
    acknowledgement$response_padding)

  expect_false(.dsvert_validate_dsi_response_probe_ack(
    acknowledgement, nonce, request_padding, request_sha256,
    response_padding_chars))
})

test_that("old response probes fall back to scalar without poisoning", {
  .dsvert_clear_transport_probe_cache()
  .dsvert_dsi_clear_poisoned_sessions()
  on.exit({
    .dsvert_clear_transport_probe_cache()
    .dsvert_dsi_clear_poisoned_sessions()
  }, add = TRUE)
  conns <- list(site = probe_connection("old-response-probe"))
  aggregate <- function(conns, expr, async, error, errors.print) {
    fields <- as.list(expr[[1L]])[-1L]
    if (is.null(fields$response_padding_chars)) return(lapply(expr, probe_ack))
    error("site", "unused argument (response_padding_chars = 65536)")
    list(site = NULL)
  }
  selected <- .dsvert_negotiate_dsi_chunk_size(
    conns, aggregate, .candidates = c(64, 32, 16) * 1024L)
  expect_false(selected$response_probe_supported[["site"]])
  expect_true(is.na(selected$response_padding_chars[["site"]]))
  expect_identical(.dsvert_dp_capsule_source_effective_window(
    "site", "site", .dsvert_dp_capsule_source_window_capability(),
    list(
      request_payload_bytes = selected$site_request_padding_chars,
      response_bytes = selected$response_usable_bytes,
      response_probe_supported = selected$response_probe_supported)), 1L)
  expect_false(.dsvert_dsi_session_is_poisoned(
    .dsvert_dsi_job_session_key(conns[[1L]])))
})

test_that("ambiguous response-probe failure poisons the exact handle", {
  .dsvert_dsi_clear_poisoned_sessions()
  on.exit(.dsvert_dsi_clear_poisoned_sessions(), add = TRUE)
  local_mocked_bindings(
    .dsvert_dsi_job_session_key = function(connection) "response-handle",
    .package = "dsVertClient")
  condition <- tryCatch(
    .dsvert_dsi_response_probe_once(
      list(site = probe_connection("ambiguous-response")), 16L * 1024L,
      function(...) stop("connection reset after execute")),
    dsvert_dsi_poisoned_session = identity)
  expect_s3_class(condition, "dsvert_dsi_poisoned_session")
  expect_true(.dsvert_dsi_session_is_poisoned("response-handle"))
})

test_that("request and response probes cross a real DSLite connector", {
  skip_if_not_installed("DSLite")
  .dsvert_clear_transport_probe_cache()
  on.exit(.dsvert_clear_transport_probe_cache(), add = TRUE)
  server_name <- paste0("dsvert_duplex_probe_", Sys.getpid())
  server <- DSLite::newDSLiteServer(tables = list(t = data.frame(x = 1)))
  server$aggregateMethod("dsvertTransportProbeDS", function(
      nonce, padding, padding_sha256, response_padding_chars = NULL) {
    request <- list(
      version = "dsvert-transport-probe-v1", nonce = nonce,
      padding_chars = as.numeric(nchar(padding, type = "bytes")),
      padding_sha256 = padding_sha256,
      server_max_padding_chars = 8 * 1024^2)
    if (is.null(response_padding_chars)) return(request)
    response_padding <- strrep("R", as.integer(response_padding_chars))
    list(
      version = "dsvert-transport-response-probe-v1", nonce = nonce,
      padding_chars = request$padding_chars,
      padding_sha256 = padding_sha256,
      server_max_padding_chars = 8 * 1024^2,
      response_padding_chars = as.numeric(response_padding_chars),
      response_padding_sha256 = digest::digest(
        response_padding, "sha256", serialize = FALSE),
      server_max_response_padding_chars = 8 * 1024^2,
      response_padding = response_padding)
  })
  assign(server_name, server, envir = .GlobalEnv)
  withr::defer(rm(list = server_name, envir = .GlobalEnv))
  builder <- DSI::newDSLoginBuilder()
  builder$append(
    server = "site", url = server_name, table = "t",
    driver = "DSLiteDriver")
  conns <- DSI::datashield.login(builder$build(), assign = FALSE)
  withr::defer(DSI::datashield.logout(conns))

  selected <- .dsvert_negotiate_dsi_chunk_size(
    conns, DSI::datashield.aggregate,
    .candidates = c(64, 32, 16) * 1024L)
  expect_identical(selected$chunk_chars, 64 * 1024)
  expect_identical(
    selected$response_padding_chars[["site"]], 64 * 1024)
  expect_true(selected$response_probe_supported[["site"]])
})

test_that("a connector profile hint is re-probed after reconnect", {
  .dsvert_clear_transport_probe_cache()
  on.exit(.dsvert_clear_transport_probe_cache(), add = TRUE)
  armadillo <- function(token) structure(
    list(handle = list(url = "https://armadillo.example.org"), token = token),
    class = c("ArmadilloConnection", "list"))
  accepted <- 2048
  attempted <- numeric()
  aggregate <- function(conns, expr, async, error, errors.print) {
    size <- nchar(as.list(expr[[1L]])$padding, type = "bytes")
    attempted <<- c(attempted, size)
    if (size > accepted) {
      return(stats::setNames(rep(list(NULL), length(conns)), names(conns)))
    }
    lapply(expr, probe_ack)
  }
  candidates <- c(4096, 2048, 1024)
  first <- .dsvert_negotiate_dsi_chunk_size(
    list(site = armadillo("token-one")), aggregate,
    .candidates = candidates)
  expect_identical(first$chunk_chars, 2048)
  expect_identical(attempted, c(4096, 2048))

  .dsvert_reset_chunk_size()
  attempted <- numeric()
  reconnected <- .dsvert_negotiate_dsi_chunk_size(
    list(site = armadillo("token-two")), aggregate,
    .candidates = candidates)
  expect_identical(reconnected$chunk_chars, 2048)
  expect_identical(attempted, 2048)

  # A stale hint is never trusted: a smaller proxy ceiling is discovered with
  # public padding before any opaque payload uses the new geometry.
  .dsvert_reset_chunk_size()
  attempted <- numeric()
  accepted <- 1024
  smaller <- .dsvert_negotiate_dsi_chunk_size(
    list(site = armadillo("token-three")), aggregate,
    .candidates = candidates)
  expect_identical(smaller$chunk_chars, 1024)
  expect_identical(attempted, c(2048, 1024))
})

test_that("probe caches distinguish port and base path on one host", {
  skip_if_not_installed("httr")
  .dsvert_clear_transport_probe_cache()
  on.exit(.dsvert_clear_transport_probe_cache(), add = TRUE)
  armadillo <- function(url) structure(
    list(handle = httr::handle(url), token = "same-token"),
    class = c("ArmadilloConnection", "list"))
  one <- list(site = armadillo("https://example.org:8443/one"))
  two <- list(site = armadillo("https://example.org:9443/two"))
  expect_false(identical(
    .dsvert_dsi_probe_cache_key(one), .dsvert_dsi_probe_cache_key(two)))
  expect_false(identical(
    .dsvert_dsi_probe_profile_key(one),
    .dsvert_dsi_probe_profile_key(two)))
})

test_that("Armadillo authoritative probe cache follows the curl handle", {
  skip_if_not_installed("httr")
  armadillo <- function(handle, token) structure(
    list(handle = handle, token = token),
    class = c("ArmadilloConnection", "list"))
  first_handle <- httr::handle("https://armadillo.example.org")
  first <- armadillo(first_handle, "token-one")
  refreshed <- armadillo(first_handle, "token-two")
  fresh_login <- armadillo(
    httr::handle("https://armadillo.example.org"), "token-two")

  expect_identical(
    .dsvert_dsi_probe_cache_key(list(site = first)),
    .dsvert_dsi_probe_cache_key(list(site = refreshed)))
  expect_false(identical(
    .dsvert_dsi_probe_cache_key(list(site = first)),
    .dsvert_dsi_probe_cache_key(list(site = fresh_login))))
})

test_that("only the measured DSLite release receives a class-derived hint", {
  conns <- list(site = probe_connection())
  candidates <- .DSVERT_DSI_PROBE_CANDIDATES
  expect_identical(.dsvert_dsi_known_probe_hint(
    conns, candidates, .package_version = function(package) "1.4.1"),
    as.numeric(640 * 1024))
  expect_null(.dsvert_dsi_known_probe_hint(
    conns, candidates, .package_version = function(package) "1.5.0"))
  expect_null(.dsvert_dsi_known_probe_hint(
    conns, candidates, .package_version = function(package) stop("missing")))
})

test_that("automatic geometry cannot leak across connection session sets", {
  .dsvert_clear_transport_probe_cache()
  on.exit(.dsvert_clear_transport_probe_cache(), add = TRUE)
  first <- list(site = probe_connection("cohort-one"))
  second <- list(site = probe_connection("cohort-two"))
  calls <- character()
  aggregate <- function(conns, expr, async, error, errors.print) {
    calls <<- c(calls, .dsvert_dsi_connection_session(conns[[1L]]))
    lapply(expr, probe_ack)
  }

  .dsvert_negotiate_dsi_chunk_size(
    first, aggregate, .candidates = c(4096, 2048))
  first_key <- .dsvert_chunk_env$active_probe_cache_key
  .dsvert_chunk_env$geometry_locked <- FALSE
  .dsvert_maybe_negotiate_dsi_chunk_size(second, aggregate)

  # Injected aggregators intentionally skip an automatic probe, but the old
  # cohort's negotiated geometry and cache identity must still be discarded.
  expect_null(.dsvert_chunk_env$effective_chunk_size)
  expect_false(identical(
    .dsvert_chunk_env$active_probe_cache_key, first_key))
  expect_identical(
    .dsvert_chunk_env$active_probe_cache_key,
    .dsvert_dsi_probe_cache_key(second))
  expect_identical(calls, "sid:cohort-one")
})

test_that("geometry unlocks only after a payload transfer returns", {
  .dsvert_clear_transport_probe_cache()
  on.exit(.dsvert_clear_transport_probe_cache(), add = TRUE)
  conns <- list(site = probe_connection("transfer-owner"))
  .dsvert_chunk_env$active_probe_cache_key <-
    .dsvert_dsi_probe_cache_key(conns)
  observed_locked <- logical()
  response <- .dsvert_adaptive_send(
    strrep("A", 20L),
    function(chunk_str, chunk_index, n_chunks) {
      observed_locked <<- c(
        observed_locked, isTRUE(.dsvert_chunk_env$geometry_locked))
      list(site = TRUE)
    }, target = "site", idempotent = TRUE)
  expect_identical(response, 1L)
  expect_true(all(observed_locked))
  expect_false(.dsvert_chunk_env$geometry_locked)
})

test_that("present malformed probe acknowledgements are fatal", {
  .dsvert_clear_transport_probe_cache()
  on.exit(.dsvert_clear_transport_probe_cache(), add = TRUE)
  calls <- 0L
  malformed <- function(conns, expr, async, error, errors.print) {
    calls <<- calls + 1L
    list(site = list(version = "wrong"))
  }
  expect_error(
    .dsvert_negotiate_dsi_chunk_size(
      list(site = probe_connection()), malformed,
      .candidates = c(4096, 2048)),
    "malformed acknowledgement")
  expect_identical(calls, 1L)
})

test_that("unavailable probe never authorizes an untested geometry", {
  .dsvert_clear_transport_probe_cache()
  on.exit(.dsvert_clear_transport_probe_cache(), add = TRUE)
  attempted <- numeric()
  unavailable <- function(conns, expr, async, error, errors.print) {
    attempted <<- c(
      attempted, nchar(as.list(expr[[1L]])$padding, type = "bytes"))
    stats::setNames(rep(list(NULL), length(conns)), names(conns))
  }
  expect_error(.dsvert_negotiate_dsi_chunk_size(
    list(site = probe_connection()), unavailable,
    .candidates = c(4096, 2048, 1024)),
    "rejected every data-free transport probe")
  expect_identical(attempted, c(4096, 2048, 1024))
  expect_null(.dsvert_chunk_env$effective_chunk_size)
})

test_that("payload geometry cannot downshift after its first attempt", {
  .dsvert_clear_transport_probe_cache()
  on.exit(.dsvert_clear_transport_probe_cache(), add = TRUE)
  withr::local_options(list(dsvert.dsi.retry_deadline_seconds = 0.001))
  conns <- list(site = probe_connection())
  aggregate <- function(conns, expr, async, error, errors.print) {
    lapply(expr, probe_ack)
  }
  .dsvert_negotiate_dsi_chunk_size(
    conns, aggregate, .candidates = c(1024, 512))
  observed <- list()
  expect_error(.dsvert_adaptive_send(
    strrep("A", 2048L),
    function(chunk_str, chunk_index, n_chunks) {
      observed[[length(observed) + 1L]] <<- list(
        chunk = chunk_str, index = chunk_index, count = n_chunks)
      NULL
    }, target = "site", idempotent = TRUE),
    class = "retry_deadline_exceeded")
  expect_gte(length(observed), 1L)
  expect_true(all(vapply(observed, identical, logical(1L), observed[[1L]])))
  expect_identical(nchar(observed[[1L]]$chunk), 1024L)
  expect_false(.dsvert_chunk_env$geometry_locked)
  # The failed transfer is over from the client's perspective. Reusing the
  # same connection/session retains the immutable cached geometry; a different
  # connection is negotiated only between transfer attempts.
  cached <- .dsvert_negotiate_dsi_chunk_size(
    conns, aggregate, .candidates = c(512, 256))
  expect_identical(cached$chunk_chars, 1024)
})

test_that("production candidate ladder and expression reserve are bounded", {
  expect_identical(
    .DSVERT_DSI_PROBE_CANDIDATES,
    as.numeric(c(
      8 * 1024^2, 4 * 1024^2, 2 * 1024^2, 1024^2,
      640 * 1024, 320 * 1024, 160 * 1024, 80 * 1024,
      32 * 1024, 16 * 1024)))
  expect_true(all(diff(.DSVERT_DSI_PROBE_CANDIDATES) < 0))
  expect_lte(max(.DSVERT_DSI_PROBE_CANDIDATES), 8 * 1024^2)
  expect_gte(.DSVERT_DSI_PROBE_EXPRESSION_RESERVE, 32 * 1024L)
  expect_identical(class(.dsvert_dsi_probe_hash("public")), "character")
  withr::local_options(list(dsvert.chunk_size = 8 * 1024^2 + 1))
  .dsvert_reset_chunk_size()
  on.exit(.dsvert_reset_chunk_size(), add = TRUE)
  expect_error(.dsvert_get_chunk_size(), "no larger than 8 MiB")
})
