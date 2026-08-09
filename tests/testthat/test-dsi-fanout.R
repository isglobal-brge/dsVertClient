test_that("named per-site fan-out uses one aggregate call and preserves mapping", {
  conns <- list(s1 = list(), s2 = list(), s3 = list())
  expressions <- list(
    s1 = call("phaseDS", value = 1L),
    s2 = call("phaseDS", value = 2L),
    s3 = call("phaseDS", value = 3L))
  calls <- list()
  aggregate <- function(conns, expr, error, async, errors.print, ...) {
    calls[[length(calls) + 1L]] <<- list(
      conns = conns, expr = expr, async = async,
      errors.print = errors.print)
    stats::setNames(lapply(expr, function(x) x[["value"]]), names(expr))
  }

  result <- dsVertClient:::.dsvert_fanout_by_site(
    conns, expressions, operation = "test phase", .aggregate = aggregate)

  expect_length(calls, 1L)
  expect_identical(names(calls[[1L]]$conns), names(conns))
  expect_identical(calls[[1L]]$expr, expressions)
  expect_true(calls[[1L]]$async)
  expect_identical(calls[[1L]]$errors.print, FALSE)
  expect_identical(result, list(s1 = 1L, s2 = 2L, s3 = 3L))
})

test_that("immutable transport fan-out maps real DSLite results by site", {
  skip_if_not_installed("DSLite")
  method <- paste0("dsvertDirectTransport", Sys.getpid())
  refs <- character()
  for (index in 1:2) {
    server <- DSLite::newDSLiteServer(tables = list(t = data.frame(x = 1)))
    server$aggregateMethod(method, function(payload, ordinal) list(
      payload_chars = nchar(payload, type = "bytes"), ordinal = ordinal))
    object <- paste0("dsvert_direct_transport_", Sys.getpid(), "_", index)
    assign(object, server, envir = .GlobalEnv)
    refs[[paste0("site", index)]] <- object
  }
  on.exit(rm(list = unname(refs), envir = .GlobalEnv), add = TRUE)
  builder <- DSI::newDSLoginBuilder()
  for (site in names(refs)) {
    builder$append(
      server = site, url = refs[[site]], table = "t",
      driver = "DSLiteDriver")
  }
  conns <- DSI::datashield.login(builder$build(), assign = FALSE)
  on.exit(DSI::datashield.logout(conns), add = TRUE)
  expressions <- stats::setNames(lapply(seq_along(conns), function(index) {
    call(name = method, payload = strrep("A", 16 * 1024L), ordinal = index)
  }), names(conns))

  result <- .dsvert_transport_aggregate(
    DSI::datashield.aggregate, conns, expressions,
    async = TRUE, errors.print = FALSE)
  expect_named(result, names(conns))
  expect_identical(unname(vapply(
    result, `[[`, numeric(1L), "payload_chars")), rep(16 * 1024, 2L))
  expect_identical(unname(vapply(
    result, `[[`, integer(1L), "ordinal")), 1:2)

  failures <- character()
  partial_expr <- expressions
  partial_expr$site2 <- call(name = paste0(method, "Missing"))
  partial <- .dsvert_transport_aggregate(
    DSI::datashield.aggregate, conns, partial_expr, async = TRUE,
    error = function(site, message) failures <<- c(failures, site),
    errors.print = FALSE)
  expect_named(partial, names(conns))
  expect_identical(partial$site1$ordinal, 1L)
  expect_null(partial$site2)
  expect_identical(failures, "site2")

  reused <- .dsvert_transport_aggregate(
    DSI::datashield.aggregate, conns, expressions,
    async = TRUE, errors.print = FALSE)
  expect_identical(unname(vapply(
    reused, `[[`, integer(1L), "ordinal")), 1:2)
})

test_that("direct async fan-out stops polling at a monotonic deadline", {
  .dsvert_dsi_clear_poisoned_sessions()
  on.exit(.dsvert_dsi_clear_poisoned_sessions(), add = TRUE)
  times <- c(0, 0.5, 1.1)
  clock <- function() {
    value <- times[[1L]]
    times <<- times[-1L]
    value
  }
  failures <- character()
  polls <- 0L
  result <- testthat::with_mocked_bindings(
    .dsvert_dsi_direct_aggregate(
      conns = list(site = list(sid = "deadline-session")),
      expr = call("exchangeDS"),
      error = function(site, message) failures <<- c(failures, site),
      timeout_seconds = 1, .clock = clock, .sleep = function(...) NULL),
    datashield.sessions = function(...) invisible(NULL),
    dsIsAsync = function(...) list(aggregate = TRUE),
    dsAggregate = function(...) structure(list(), class = "test_job"),
    dsIsCompleted = function(...) {
      polls <<- polls + 1L
      FALSE
    },
    dsFetch = function(...) stop("a pending job must not be fetched"),
    dsKeepAlive = function(...) invisible(NULL),
    .package = "DSI")

  expect_identical(result, list(site = NULL))
  expect_identical(failures, "site")
  expect_identical(polls, 2L)
})

test_that("remote jobs have no automatic wall-clock cutoff by default", {
  withr::local_options(list(dsvert.dsi.job_timeout_seconds = NULL))
  expect_identical(.dsvert_dsi_job_timeout_seconds(), Inf)

  withr::local_options(list(dsvert.dsi.job_timeout_seconds = 7200))
  expect_identical(.dsvert_dsi_job_timeout_seconds(), 7200)

  for (invalid in list(0, -1, NA_real_, NaN)) {
    withr::local_options(list(dsvert.dsi.job_timeout_seconds = invalid))
    expect_error(.dsvert_dsi_job_timeout_seconds(), "positive or Inf")
  }
})

test_that("an unresolved async job poisons only that authenticated session", {
  .dsvert_dsi_clear_poisoned_sessions()
  on.exit(.dsvert_dsi_clear_poisoned_sessions(), add = TRUE)
  launches <- 0L
  poll <- FALSE
  aggregate <- function(conns, clock) {
    testthat::with_mocked_bindings(
      .dsvert_dsi_direct_aggregate(
        conns = conns, expr = call("exchangeDS"),
        timeout_seconds = 1, .clock = clock,
        .sleep = function(...) NULL),
      datashield.sessions = function(...) invisible(NULL),
      dsIsAsync = function(...) list(aggregate = TRUE),
      dsAggregate = function(...) {
        launches <<- launches + 1L
        structure(list(), class = "test_job")
      },
      dsIsCompleted = function(...) poll,
      dsFetch = function(...) "settled",
      dsKeepAlive = function(...) invisible(NULL),
      .package = "DSI")
  }

  times <- c(0, 0.5, 1.1)
  old <- list(site = list(sid = "old-session"))
  expect_identical(aggregate(old, function() {
    value <- times[[1L]]
    times <<- times[-1L]
    value
  }), list(site = NULL))
  expect_identical(launches, 1L)

  expect_error(
    aggregate(old, local({ value <- -0.1; function() {
      value <<- value + 0.1
      value
    }})),
    "fresh DSI login connections")
  expect_identical(launches, 1L)

  expect_error(
    aggregate(list(renamed = old[[1L]]), local({
      value <- -0.1
      function() {
        value <<- value + 0.1
        value
      }
    })),
    "fresh DSI login connections")
  expect_identical(launches, 1L)

  poll <- TRUE
  fresh <- list(site = list(sid = "fresh-session"))
  expect_identical(
    aggregate(fresh, local({ value <- -0.1; function() {
      value <<- value + 0.1
      value
    }})),
    list(site = "settled"))
  expect_identical(launches, 2L)
})

test_that("direct async fan-out launches, polls, fetches, and keeps peers alive", {
  events <- character()
  polls <- c(site1 = 0L, site2 = 0L)
  result <- testthat::with_mocked_bindings(
    .dsvert_dsi_direct_aggregate(
      conns = list(
        site1 = list(sid = "lifecycle-site1"),
        site2 = list(sid = "lifecycle-site2")),
      expr = list(
        site1 = call("exchangeDS", value = 1L),
        site2 = call("exchangeDS", value = 2L)),
      timeout_seconds = 2,
      .clock = local({ value <- -0.1; function() {
        value <<- value + 0.1
        value
      }}),
      .sleep = function(...) events <<- c(events, "sleep")),
    datashield.sessions = function(...) {
      events <<- c(events, "sessions")
      invisible(NULL)
    },
    dsIsAsync = function(conn) list(aggregate = TRUE),
    dsAggregate = function(conn, expr, async) {
      site <- if (identical(expr[["value"]], 1L)) "site1" else "site2"
      events <<- c(events, paste0("launch:", site))
      structure(list(site = site), class = "test_job")
    },
    dsIsCompleted = function(job) {
      site <- job$site
      polls[[site]] <<- polls[[site]] + 1L
      events <<- c(events, paste0("poll:", site))
      identical(site, "site1") || polls[[site]] >= 2L
    },
    dsFetch = function(job) {
      events <<- c(events, paste0("fetch:", job$site))
      job$site
    },
    dsKeepAlive = function(conn) {
      events <<- c(events, "keepalive:site1")
      invisible(NULL)
    },
    .package = "DSI")

  expect_identical(result, list(site1 = "site1", site2 = "site2"))
  expect_lt(match("launch:site2", events), match("poll:site1", events))
  expect_lt(match("fetch:site1", events), match("keepalive:site1", events))
  expect_lt(match("keepalive:site1", events), match("fetch:site2", events))
  expect_identical(unname(polls), c(1L, 2L))
})

test_that("a consumed terminal async error does not poison Armadillo", {
  .dsvert_dsi_clear_poisoned_sessions()
  withr::defer(.dsvert_dsi_clear_poisoned_sessions())
  key <- "terminal-async-armadillo"
  launches <- 0L
  fetches <- 0L
  failures <- character()
  local_mocked_bindings(
    .dsvert_dsi_inspect_connection = function(connection) {
      list(kind = "armadillo", endpoint = list(scope = "https://example.test"))
    },
    .dsvert_dsi_job_session_key = function(connection) key,
    .package = "dsVertClient")

  run <- function() testthat::with_mocked_bindings(
    .dsvert_dsi_direct_aggregate(
      list(site = list()), call("publicCapabilityDS"),
      error = function(site, message) failures <<- c(failures, site)),
    datashield.sessions = function(...) invisible(NULL),
    dsIsAsync = function(...) list(aggregate = TRUE),
    dsAggregate = function(...) {
      launches <<- launches + 1L
      structure(list(), class = "test_job")
    },
    dsIsCompleted = function(...) TRUE,
    dsFetch = function(...) {
      fetches <<- fetches + 1L
      if (fetches == 1L) stop(simpleError(paste0(
        "Command 'publicCapabilityDS' failed on armadillo: Error whilst ",
        "evaluating request -> expected remote miss")))
      "reused"
    },
    dsKeepAlive = function(...) invisible(NULL),
    .package = "DSI")

  expect_identical(run(), list(site = NULL))
  expect_identical(failures, "site")
  expect_false(.dsvert_dsi_session_is_poisoned(key))
  expect_identical(run(), list(site = "reused"))
  expect_identical(c(launches, fetches), c(2L, 2L))
})

test_that("a fetched phase-not-ready terminal is settled and does not poison", {
  .dsvert_dsi_clear_poisoned_sessions()
  withr::defer(.dsvert_dsi_clear_poisoned_sessions())
  key <- "terminal-async-phase-not-ready"
  launches <- 0L
  fetches <- 0L
  failures <- character()
  local_mocked_bindings(
    .dsvert_dsi_job_session_key = function(connection) key,
    .package = "dsVertClient")

  run <- function() testthat::with_mocked_bindings(
    .dsvert_dsi_direct_aggregate(
      list(site = list()), call("dsvertJointDPVectorResultDS"),
      error = function(site, message) failures <<- c(failures, site)),
    datashield.sessions = function(...) invisible(NULL),
    dsIsAsync = function(...) list(aggregate = TRUE),
    dsAggregate = function(...) {
      launches <<- launches + 1L
      structure(list(), class = "test_job")
    },
    dsIsCompleted = function(...) TRUE,
    dsFetch = function(...) {
      fetches <<- fetches + 1L
      if (fetches == 1L) {
        stop(simpleError(paste0(
          "[connector (400)] [dsvert_phase_not_ready:v1] ",
          "SECRET remote state")))
      }
      "reused"
    },
    dsKeepAlive = function(...) invisible(NULL),
    .package = "DSI")

  expect_identical(run(), list(site = NULL))
  expect_identical(failures, "site")
  expect_false(.dsvert_dsi_session_is_poisoned(key))
  expect_identical(run(), list(site = "reused"))
  expect_false(.dsvert_dsi_session_is_poisoned(key))
  expect_identical(c(launches, fetches), c(2L, 2L))
})

test_that("a fetched lifetime exhaustion terminal is settled and does not poison", {
  .dsvert_dsi_clear_poisoned_sessions()
  withr::defer(.dsvert_dsi_clear_poisoned_sessions())
  key <- "terminal-async-lifetime-exhausted"
  launches <- 0L
  fetches <- 0L
  failures <- character()
  local_mocked_bindings(
    .dsvert_dsi_job_session_key = function(connection) key,
    .package = "dsVertClient")

  run <- function() testthat::with_mocked_bindings(
    .dsvert_dsi_direct_aggregate(
      list(site = list()), call("dsvertJointDPVectorResultDS"),
      error = function(site, message) failures <<- c(failures, site)),
    datashield.sessions = function(...) invisible(NULL),
    dsIsAsync = function(...) list(aggregate = TRUE),
    dsAggregate = function(...) {
      launches <<- launches + 1L
      structure(list(), class = "test_job")
    },
    dsIsCompleted = function(...) TRUE,
    dsFetch = function(...) {
      fetches <<- fetches + 1L
      if (fetches == 1L) {
        stop(simpleError(paste0(
          "[connector (400)] ",
          "[dsvert_dp_lifetime_budget_exhausted:v1] ",
          "SECRET remote state")))
      }
      "reused"
    },
    dsKeepAlive = function(...) invisible(NULL),
    .package = "DSI")

  expect_identical(run(), list(site = NULL))
  expect_identical(failures, "site")
  expect_false(.dsvert_dsi_session_is_poisoned(key))
  expect_identical(run(), list(site = "reused"))
  expect_false(.dsvert_dsi_session_is_poisoned(key))
  expect_identical(c(launches, fetches), c(2L, 2L))
})

test_that("ambiguous async fetch reset and 500 poison Armadillo", {
  local_mocked_bindings(
    .dsvert_dsi_inspect_connection = function(connection) {
      list(kind = "armadillo", endpoint = list(scope = "https://example.test"))
    },
    .dsvert_dsi_job_session_key = function(connection) connection$sid,
    .package = "dsVertClient")

  for (message in c(
      "connection reset during dsFetch",
      "Internal server error: HTTP 500 after dsFetch")) {
    .dsvert_dsi_clear_poisoned_sessions()
    withr::defer(.dsvert_dsi_clear_poisoned_sessions())
    key <- paste0("ambiguous-async-", digest::digest(
      message, algo = "sha256", serialize = FALSE))
    launches <- 0L
    outcome <- testthat::with_mocked_bindings({
      first <- .dsvert_dsi_direct_aggregate(
        list(site = list(sid = key)), call("publicCapabilityDS"))
      second <- tryCatch(
        .dsvert_dsi_direct_aggregate(
          list(site = list(sid = key)), call("publicCapabilityDS")),
        dsvert_dsi_poisoned_session = identity)
      list(first = first, second = second)
    },
    datashield.sessions = function(...) invisible(NULL),
    dsIsAsync = function(...) list(aggregate = TRUE),
    dsAggregate = function(...) {
      launches <<- launches + 1L
      structure(list(), class = "test_job")
    },
    dsIsCompleted = function(...) TRUE,
    dsFetch = function(...) stop(simpleError(message)),
    dsKeepAlive = function(...) invisible(NULL),
    .package = "DSI")

    expect_identical(outcome$first, list(site = NULL), info = message)
    expect_s3_class(outcome$second, "dsvert_dsi_poisoned_session")
    expect_true(.dsvert_dsi_session_is_poisoned(key), info = message)
    expect_identical(launches, 1L, info = message)
  }
})

test_that("an async fetch interrupt poisons the launched session", {
  .dsvert_dsi_clear_poisoned_sessions()
  withr::defer(.dsvert_dsi_clear_poisoned_sessions())
  key <- "interrupted-async-fetch"
  local_mocked_bindings(
    .dsvert_dsi_job_session_key = function(connection) key,
    .package = "dsVertClient")
  interruption <- structure(
    list(message = "user interrupt during dsFetch", call = NULL),
    class = c("interrupt", "condition"))

  outcome <- testthat::with_mocked_bindings(
    tryCatch(
      .dsvert_dsi_direct_aggregate(
        list(site = list()), call("publicCapabilityDS")),
      interrupt = identity),
    datashield.sessions = function(...) invisible(NULL),
    dsIsAsync = function(...) list(aggregate = TRUE),
    dsAggregate = function(...) structure(list(), class = "test_job"),
    dsIsCompleted = function(...) TRUE,
    dsFetch = function(...) stop(interruption),
    dsKeepAlive = function(...) invisible(NULL),
    .package = "DSI")

  expect_s3_class(outcome, "interrupt")
  expect_true(.dsvert_dsi_session_is_poisoned(key))
})

test_that("named fan-out rejects DSI versions without the required mapping API", {
  expect_invisible(.dsvert_require_dsi_fanout_version(
    .package_version = function(package) "1.8.0"))
  expect_error(
    .dsvert_require_dsi_fanout_version(
      .package_version = function(package) "1.7.9"),
    "DSI >= 1.8.0")
  expect_error(
    .dsvert_require_dsi_fanout_version(
      .package_version = function(package) stop("not installed")),
    "DSI >= 1.8.0")
})

test_that("strict aggregate rejects partial and total NULL results", {
  conns <- list(s1 = list(), s2 = list())
  expr <- call("phaseDS")

  partial <- function(conns, expr, error, ...) {
    error("s2", "SECRET backend detail")
    list(s1 = TRUE, s2 = NULL)
  }
  condition <- tryCatch(
    dsVertClient:::.dsvert_aggregate_strict(
      conns, expr, operation = "safe phase", .aggregate = partial),
    error = identity)
  expect_s3_class(condition, "error")
  expect_match(conditionMessage(condition), "partial or invalid site results")
  expect_false(grepl("SECRET", conditionMessage(condition), fixed = TRUE))

  total <- function(conns, expr, error, ...) NULL
  expect_error(
    dsVertClient:::.dsvert_aggregate_strict(
      conns, expr, operation = "safe phase", .aggregate = total),
    "partial or invalid site results")
})

test_that("strict aggregate routes real DSI through the direct dispatcher", {
  called <- FALSE
  conns <- list(site = structure(
    list(sid = "strict-direct-session"),
    class = c("DSLiteConnection", "list")))
  result <- testthat::with_mocked_bindings(
    .dsvert_aggregate_strict(
      conns, call("phaseDS"), result_contract = "logical_true",
      .aggregate = DSI::datashield.aggregate),
    .dsvert_dsi_direct_aggregate = function(
        conns, expr, async, error, errors.print, ...) {
      called <<- TRUE
      list(site = TRUE)
    })
  expect_true(called)
  expect_identical(result, list(site = TRUE))
})

test_that("best-effort cleanup never launches on a poisoned session", {
  .dsvert_dsi_clear_poisoned_sessions()
  on.exit(.dsvert_dsi_clear_poisoned_sessions(), add = TRUE)
  connection <- structure(
    list(sid = "poisoned-cleanup-session"),
    class = c("DSLiteConnection", "list"))
  conns <- list(site = connection)
  .dsvert_dsi_poison_sessions(.dsvert_dsi_job_session_key(connection))
  launches <- 0L

  cleaned <- testthat::with_mocked_bindings(
    .dsvert_cleanup_best_effort(
      conns, call("cleanupDS"), .aggregate = DSI::datashield.aggregate),
    datashield.sessions = function(...) invisible(NULL),
    dsIsAsync = function(...) list(aggregate = TRUE),
    dsAggregate = function(...) {
      launches <<- launches + 1L
      structure(list(), class = "test_job")
    },
    .package = "DSI")

  expect_false(cleaned)
  expect_identical(launches, 0L)
})

test_that("strict wrappers immediately propagate a poisoned DSI session", {
  conns <- list(site = list())
  poisoned <- .dsvert_dsi_poisoned_session_condition("site")
  aggregate <- function(...) stop(poisoned)

  strict <- tryCatch(
    .dsvert_aggregate_strict(
      conns, call("phaseDS"), .aggregate = aggregate),
    dsvert_dsi_poisoned_session = identity)
  expect_s3_class(strict, "dsvert_dsi_poisoned_session")
  expect_match(conditionMessage(strict), "fresh DSI login connections")

  cycle <- tryCatch(
    .dsvert_fanout_cycle(
      conns, list(site = call("exchangeDS")), .aggregate = aggregate),
    dsvert_dsi_poisoned_session = identity)
  expect_s3_class(cycle, "dsvert_dsi_poisoned_session")
  expect_match(conditionMessage(cycle), "fresh DSI login connections")
})

test_that("strict aggregate never exposes top-level connector details", {
  conns <- list(s1 = list())
  aggregate <- function(conns, expr, error, ...) {
    stop("SECRET connector failure")
  }

  condition <- tryCatch(
    dsVertClient:::.dsvert_aggregate_strict(
      conns, call("phaseDS"), operation = "safe phase",
      .aggregate = aggregate),
    error = identity)

  expect_false(grepl("SECRET", conditionMessage(condition), fixed = TRUE))
  expect_match(conditionMessage(condition), "partial or invalid site results")
})

test_that("unknown bracketed dsvert tokens are terminal and sanitized", {
  token <- "[dsvert_unknown_future:v9]"
  conns <- list(s1 = list())
  expr <- call(
    "psiPaddedRelayExchangeDS", request_json = "frozen",
    session_id = "12345678-1234-4234-9234-123456789abc",
    outbound_operation_id = "op_11111111111111111111111111111111",
    terminal_receipt_b64url = "")

  callback_attempts <- 0L
  callback <- function(conns, expr, error, ...) {
    callback_attempts <<- callback_attempts + 1L
    error("s1", paste(token, "SECRET callback"))
    list(s1 = NULL)
  }
  callback_condition <- tryCatch(
    .dsvert_aggregate_strict(conns, expr, .aggregate = callback),
    error = identity)
  expect_s3_class(callback_condition, "error")
  expect_match(conditionMessage(callback_condition), "transport failed")
  expect_false(grepl(token, conditionMessage(callback_condition), fixed = TRUE))
  expect_false(grepl(
    "SECRET", conditionMessage(callback_condition), fixed = TRUE))
  expect_identical(callback_attempts, 1L)

  top_attempts <- 0L
  remote_condition <- structure(
    list(message = paste(token, "SECRET top"), call = NULL),
    class = c("dsvert_unknown_future", "error", "condition"))
  top <- function(...) {
    top_attempts <<- top_attempts + 1L
    stop(remote_condition)
  }
  clock_values <- c(0, 0, 1)
  clock_calls <- 0L
  top_condition <- testthat::with_mocked_bindings(
    tryCatch(
      .dsvert_aggregate_strict(conns, expr, .aggregate = top),
      error = identity),
    .dsvert_dsi_retry_deadline_seconds = function() 0.5,
    .dsvert_monotonic_seconds = function() {
      clock_calls <<- clock_calls + 1L
      clock_values[[min(clock_calls, length(clock_values))]]
    },
    .dsvert_retry_sleep = function(...) NULL,
    .dsvert_retry_jitter = function() 1,
    .package = "dsVertClient")
  expect_s3_class(top_condition, "error")
  expect_false(inherits(top_condition, "dsvert_unknown_future"))
  expect_false(inherits(top_condition, "dsvert_retry_deadline_exceeded"))
  expect_match(conditionMessage(top_condition), "transport failed")
  expect_false(grepl(token, conditionMessage(top_condition), fixed = TRUE))
  expect_false(grepl("SECRET", conditionMessage(top_condition), fixed = TRUE))
  expect_identical(top_attempts, 1L)
})

test_that("phase-not-ready parsing is exact, typed, and sanitized", {
  tag <- "[dsvert_phase_not_ready:v1]"
  parsed <- .dsvert_client_parse_phase_not_ready(paste0(
    "[connector (400)] ", tag, " SECRET remote state"))

  expect_s3_class(parsed, "dsvert_phase_not_ready")
  expect_identical(parsed$code, "phase_not_ready")
  expect_identical(parsed$retryable, FALSE)
  expect_match(conditionMessage(parsed), tag, fixed = TRUE)
  expect_false(grepl("SECRET", conditionMessage(parsed), fixed = TRUE))
  expect_null(.dsvert_client_parse_phase_not_ready(
    "dsvert_phase_not_ready:v1"))
  expect_null(.dsvert_client_parse_phase_not_ready(
    "[dsvert_phase_not_ready:v10]"))
})

test_that("lifetime exhaustion parsing is exact, typed, and sanitized", {
  tag <- "[dsvert_dp_lifetime_budget_exhausted:v1]"
  parsed <- .dsvert_client_parse_dp_lifetime_budget_exhausted(paste0(
    "[connector (400)] ", tag, " SECRET remote state"))

  expect_s3_class(parsed, "dsvert_dp_lifetime_budget_exhausted")
  expect_identical(parsed$code, "dp_lifetime_budget_exhausted")
  expect_identical(parsed$retryable, FALSE)
  expect_identical(conditionMessage(parsed), tag)
  expect_false(grepl("SECRET", conditionMessage(parsed), fixed = TRUE))
  expect_null(.dsvert_client_parse_dp_lifetime_budget_exhausted(
    "dsvert_dp_lifetime_budget_exhausted:v1"))
  expect_null(.dsvert_client_parse_dp_lifetime_budget_exhausted(
    "[dsvert_dp_lifetime_budget_exhausted:v10]"))
})

test_that("strict aggregate rethrows only unmixed lifetime exhaustion once", {
  tag <- "[dsvert_dp_lifetime_budget_exhausted:v1]"
  conns <- list(s1 = list(), s2 = list())
  expr <- call(
    "psiPaddedRelayExchangeDS", request_json = "frozen",
    session_id = "12345678-1234-4234-9234-123456789abc",
    outbound_operation_id = "op_11111111111111111111111111111111",
    terminal_receipt_b64url = "")

  callback_attempts <- 0L
  callback <- function(conns, expr, error, ...) {
    callback_attempts <<- callback_attempts + 1L
    for (site in names(conns)) {
      error(site, paste0("[remote prefix] ", tag, " SECRET callback"))
    }
    stats::setNames(rep(list(NULL), length(conns)), names(conns))
  }
  callback_condition <- tryCatch(
    .dsvert_aggregate_strict(conns, expr, .aggregate = callback),
    error = identity)
  expect_s3_class(
    callback_condition, "dsvert_dp_lifetime_budget_exhausted")
  expect_identical(
    callback_condition$code, "dp_lifetime_budget_exhausted")
  expect_identical(callback_condition$retryable, FALSE)
  expect_identical(conditionMessage(callback_condition), tag)
  expect_identical(callback_attempts, 1L)

  top_attempts <- 0L
  top <- function(...) {
    top_attempts <<- top_attempts + 1L
    stop(paste0("[connector (400)] ", tag, " SECRET top"), call. = FALSE)
  }
  top_condition <- tryCatch(
    .dsvert_aggregate_strict(conns, expr, .aggregate = top),
    error = identity)
  expect_s3_class(top_condition, "dsvert_dp_lifetime_budget_exhausted")
  expect_identical(top_condition$code, "dp_lifetime_budget_exhausted")
  expect_identical(top_condition$retryable, FALSE)
  expect_identical(conditionMessage(top_condition), tag)
  expect_identical(top_attempts, 1L)

  peer_token <- paste0(
    "DSVERT_PEER_NOT_RECOGNIZED_V1|peer=s1|expected=unconfigured",
    "|observed=", strrep("a", 64L))
  other_tokens <- c(
    "[dsvert_phase_not_ready:v1]",
    "[dsvert_resource_backpressure:v1]",
    "[dsvert_resource_oversize:v1]",
    "[dsvert_retry_current_release_instance:new_release_instance]",
    "[dsvert_dp_public_failure:v1]",
    "[dsvert_unknown_future:v9]",
    peer_token)
  for (other_token in other_tokens) {
    mixed_attempts <- 0L
    mixed <- function(conns, expr, error, ...) {
      mixed_attempts <<- mixed_attempts + 1L
      for (site in names(conns)) error(site, paste(tag, other_token))
      stats::setNames(rep(list(NULL), length(conns)), names(conns))
    }
    mixed_condition <- tryCatch(
      .dsvert_aggregate_strict(conns, expr, .aggregate = mixed),
      error = identity)
    expect_s3_class(mixed_condition, "error")
    expect_false(inherits(
      mixed_condition, "dsvert_dp_lifetime_budget_exhausted"),
      info = other_token)
    expect_match(
      conditionMessage(mixed_condition), "transport failed",
      info = other_token)
    expect_identical(mixed_attempts, 1L, info = other_token)
  }

  cross_site_attempts <- 0L
  cross_site <- function(conns, expr, error, ...) {
    cross_site_attempts <<- cross_site_attempts + 1L
    error("s1", tag)
    error("s2", "[dsvert_phase_not_ready:v1]")
    list(s1 = NULL, s2 = NULL)
  }
  cross_site_condition <- tryCatch(
    .dsvert_aggregate_strict(conns, expr, .aggregate = cross_site),
    error = identity)
  expect_false(inherits(
    cross_site_condition, "dsvert_dp_lifetime_budget_exhausted"))
  expect_false(inherits(cross_site_condition, "dsvert_phase_not_ready"))
  expect_match(conditionMessage(cross_site_condition), "transport failed")
  expect_identical(cross_site_attempts, 1L)

  mixed_top_attempts <- 0L
  mixed_top <- function(...) {
    mixed_top_attempts <<- mixed_top_attempts + 1L
    stop(paste(tag, "[dsvert_phase_not_ready:v1]"), call. = FALSE)
  }
  mixed_top_condition <- tryCatch(
    .dsvert_aggregate_strict(conns, expr, .aggregate = mixed_top),
    error = identity)
  expect_false(inherits(
    mixed_top_condition, "dsvert_dp_lifetime_budget_exhausted"))
  expect_false(inherits(mixed_top_condition, "dsvert_phase_not_ready"))
  expect_match(conditionMessage(mixed_top_condition), "transport failed")
  expect_identical(mixed_top_attempts, 1L)
})

test_that("strict aggregate rethrows only an unmixed phase-not-ready tag", {
  tag <- "[dsvert_phase_not_ready:v1]"
  conns <- list(s1 = list(), s2 = list())
  expr <- call("phaseDS")

  callback_attempts <- 0L
  callback <- function(conns, expr, error, ...) {
    callback_attempts <<- callback_attempts + 1L
    for (site in names(conns)) {
      error(site, paste0("[remote prefix] ", tag, " SECRET callback"))
    }
    stats::setNames(rep(list(NULL), length(conns)), names(conns))
  }
  callback_condition <- tryCatch(
    .dsvert_aggregate_strict(conns, expr, .aggregate = callback),
    error = identity)
  expect_s3_class(callback_condition, "dsvert_phase_not_ready")
  expect_identical(callback_condition$code, "phase_not_ready")
  expect_identical(callback_condition$retryable, FALSE)
  expect_false(grepl(
    "SECRET", conditionMessage(callback_condition), fixed = TRUE))
  expect_identical(callback_attempts, 1L)

  top_attempts <- 0L
  top <- function(...) {
    top_attempts <<- top_attempts + 1L
    stop(paste0("[connector (400)] ", tag, " SECRET top"), call. = FALSE)
  }
  top_condition <- tryCatch(
    .dsvert_aggregate_strict(conns, expr, .aggregate = top),
    error = identity)
  expect_s3_class(top_condition, "dsvert_phase_not_ready")
  expect_identical(top_condition$code, "phase_not_ready")
  expect_false(grepl("SECRET", conditionMessage(top_condition), fixed = TRUE))
  expect_identical(top_attempts, 1L)

  mixed_attempts <- 0L
  mixed <- function(conns, expr, error, ...) {
    mixed_attempts <<- mixed_attempts + 1L
    error("s1", tag)
    error("s2", "SECRET unknown failure")
    list(s1 = NULL, s2 = NULL)
  }
  mixed_condition <- tryCatch(
    .dsvert_aggregate_strict(conns, expr, .aggregate = mixed),
    error = identity)
  expect_s3_class(mixed_condition, "error")
  expect_false(inherits(mixed_condition, "dsvert_phase_not_ready"))
  expect_match(conditionMessage(mixed_condition), "transport failed")
  expect_false(grepl("SECRET", conditionMessage(mixed_condition), fixed = TRUE))
  expect_identical(mixed_attempts, 1L)

  silent_null <- function(conns, expr, error, ...) {
    error("s1", tag)
    list(s1 = NULL, s2 = NULL)
  }
  silent_condition <- tryCatch(
    .dsvert_aggregate_strict(conns, expr, .aggregate = silent_null),
    error = identity)
  expect_false(inherits(silent_condition, "dsvert_phase_not_ready"))
  expect_match(conditionMessage(silent_condition), "transport failed")

  idempotent <- call(
    "psiPaddedRelayExchangeDS", request_json = "frozen",
    session_id = "12345678-1234-4234-9234-123456789abc",
    outbound_operation_id = "op_11111111111111111111111111111111",
    terminal_receipt_b64url = "")
  conflicting_attempts <- 0L
  conflicting <- function(conns, expr, error, ...) {
    conflicting_attempts <<- conflicting_attempts + 1L
    for (site in names(conns)) {
      error(site, paste(
        tag, "[dsvert_resource_backpressure:v1]"))
    }
    stats::setNames(rep(list(NULL), length(conns)), names(conns))
  }
  conflicting_condition <- tryCatch(
    .dsvert_aggregate_strict(conns, idempotent, .aggregate = conflicting),
    error = identity)
  expect_false(inherits(
    conflicting_condition, "dsvert_phase_not_ready"))
  expect_match(conditionMessage(conflicting_condition), "transport failed")
  expect_identical(conflicting_attempts, 1L)

  conflicting_top_attempts <- 0L
  conflicting_top <- function(...) {
    conflicting_top_attempts <<- conflicting_top_attempts + 1L
    stop(paste(tag, "[dsvert_resource_backpressure:v1]"), call. = FALSE)
  }
  conflicting_top_condition <- tryCatch(
    .dsvert_aggregate_strict(
      conns, idempotent, .aggregate = conflicting_top),
    error = identity)
  expect_false(inherits(
    conflicting_top_condition, "dsvert_phase_not_ready"))
  expect_match(
    conditionMessage(conflicting_top_condition), "transport failed")
  expect_identical(conflicting_top_attempts, 1L)

  retry_tag <- paste0(
    "[dsvert_retry_current_release_instance:new_release_instance]")
  conflicting_retry_attempts <- 0L
  conflicting_retry <- function(conns, expr, error, ...) {
    conflicting_retry_attempts <<- conflicting_retry_attempts + 1L
    for (site in names(conns)) error(site, paste(tag, retry_tag))
    stats::setNames(rep(list(NULL), length(conns)), names(conns))
  }
  conflicting_retry_condition <- tryCatch(
    .dsvert_aggregate_strict(
      conns, idempotent, .aggregate = conflicting_retry),
    error = identity)
  expect_false(inherits(
    conflicting_retry_condition, "dsvert_phase_not_ready"))
  expect_false(inherits(
    conflicting_retry_condition, "dsvert_release_instance_retry"))
  expect_match(
    conditionMessage(conflicting_retry_condition), "transport failed")
  expect_identical(conflicting_retry_attempts, 1L)

  conflicting_retry_top_attempts <- 0L
  conflicting_retry_top <- function(...) {
    conflicting_retry_top_attempts <<- conflicting_retry_top_attempts + 1L
    stop(paste(tag, retry_tag), call. = FALSE)
  }
  conflicting_retry_top_condition <- tryCatch(
    .dsvert_aggregate_strict(
      conns, idempotent, .aggregate = conflicting_retry_top),
    error = identity)
  expect_false(inherits(
    conflicting_retry_top_condition, "dsvert_phase_not_ready"))
  expect_false(inherits(
    conflicting_retry_top_condition, "dsvert_release_instance_retry"))
  expect_match(
    conditionMessage(conflicting_retry_top_condition), "transport failed")
  expect_identical(conflicting_retry_top_attempts, 1L)
})

test_that("strict aggregate cannot positionally misroute a per-site map", {
  conns <- list(s1 = list(), s2 = list())
  calls <- 0L
  aggregate <- function(conns, expr, ...) {
    calls <<- calls + 1L
    expect_identical(names(expr), names(conns))
    stats::setNames(rep(list(TRUE), length(conns)), names(conns))
  }
  expressions <- list(
    s2 = call("phaseDS", value = 2L),
    s1 = call("phaseDS", value = 1L))
  expect_identical(
    .dsvert_aggregate_strict(conns, expressions, .aggregate = aggregate),
    list(s1 = TRUE, s2 = TRUE))
  expect_identical(calls, 1L)

  expect_error(.dsvert_aggregate_strict(
    conns, unname(expressions), .aggregate = aggregate), "named exactly once")
  expect_error(.dsvert_aggregate_strict(
    conns, expressions["s1"], .aggregate = aggregate), "named exactly once")
  expect_identical(calls, 1L)
})

test_that("strict aggregate canonicalizes reorder and rejects mutated fan-in", {
  conns <- list(s1 = list(), s2 = list())
  reordered <- function(...) list(s2 = "two", s1 = "one")
  expect_identical(
    .dsvert_aggregate_strict(conns, call("phaseDS"),
                             .aggregate = reordered),
    list(s1 = "one", s2 = "two"))

  malformed <- list(
    missing = function(...) list(s1 = "one"),
    duplicated = function(...) structure(
      list("one", "substitution", "two"),
      names = c("s1", "s1", "s2")),
    omitted_with_attacker = function(...) list(s1 = "one", attacker = "two"))
  for (case in names(malformed)) {
    expect_error(
      .dsvert_aggregate_strict(
        conns, call("phaseDS"), .aggregate = malformed[[case]]),
      "malformed or misassociated|partial or invalid", info = case)
  }
})

test_that("NULL is accepted only under an explicit per-site contract", {
  conns <- list(s1 = list(), s2 = list())
  aggregate <- function(conns, expr, error, ...) {
    list(s1 = TRUE, s2 = NULL)
  }

  result <- dsVertClient:::.dsvert_aggregate_strict(
    conns, call("phaseDS"), operation = "nullable phase",
    allow_null = "s2", .aggregate = aggregate)
  expect_identical(result, list(s1 = TRUE, s2 = NULL))
})

test_that("logical TRUE contract rejects non-typed acknowledgements", {
  conns <- list(s1 = list(), s2 = list())
  aggregate <- function(conns, expr, error, ...) {
    list(s1 = TRUE, s2 = 1L)
  }

  expect_error(
    dsVertClient:::.dsvert_aggregate_strict(
      conns, call("phaseDS"), operation = "ack phase",
      result_contract = "logical_true", .aggregate = aggregate),
    "partial or invalid site results")
})

test_that("fan-out rejects missing, duplicate, or non-call expressions", {
  conns <- list(s1 = list(), s2 = list())
  aggregate <- function(...) stop("must not be called")

  expect_error(dsVertClient:::.dsvert_fanout_by_site(
    conns, list(s1 = call("phaseDS")), .aggregate = aggregate),
    "exactly once")
  expect_error(dsVertClient:::.dsvert_fanout_by_site(
    conns, structure(list(call("phaseDS"), call("phaseDS")),
                          names = c("s1", "s1")), .aggregate = aggregate),
    "exactly once")
  expect_error(dsVertClient:::.dsvert_fanout_by_site(
    conns, list(s1 = call("phaseDS"), s2 = 2L), .aggregate = aggregate),
    "must be an R call")
  expect_error(dsVertClient:::.dsvert_fanout_by_site(
    conns, list(s1 = call("phaseDS"), s2 = as.name("phaseDS")),
    .aggregate = aggregate),
    "must be an R call")
})

test_that("idempotent fan-out cycles classify only absent results as unavailable", {
  conns <- list(s1 = list(), s2 = list())
  expressions <- list(
    s1 = call("exchangeDS", payload = "one"),
    s2 = call("exchangeDS", payload = "two"))
  captured <- NULL
  aggregate <- function(conns, expr, async, error, errors.print) {
    captured <<- list(expr = expr, async = async, errors.print = errors.print)
    list(s2 = list(ack = 2), s1 = list(ack = 1))
  }

  cycle <- dsVertClient:::.dsvert_fanout_cycle(
    conns, expressions, operation = "opaque exchange", .aggregate = aggregate)

  expect_identical(cycle$state, "ok")
  expect_named(cycle$result, c("s1", "s2"))
  expect_identical(captured$expr, expressions)
  expect_true(captured$async)
  expect_false(captured$errors.print)

  absent <- function(conns, expr, async, error, errors.print) {
    list(s1 = list(ack = 1), s2 = NULL)
  }
  expect_identical(
    dsVertClient:::.dsvert_fanout_cycle(
      conns, expressions, .aggregate = absent)$state,
    "unavailable")

  transient <- function(...) stop("private connector text")
  condition <- dsVertClient:::.dsvert_fanout_cycle(
    conns, expressions, .aggregate = transient)
  expect_identical(condition, list(state = "unavailable", result = NULL))
})

test_that("idempotent fan-out cycles reject present outer protocol corruption", {
  conns <- list(s1 = list(), s2 = list())
  expressions <- list(s1 = call("exchangeDS"), s2 = call("exchangeDS"))
  wrong <- function(...) list(s1 = TRUE, attacker = TRUE)

  condition <- tryCatch(
    dsVertClient:::.dsvert_fanout_cycle(
      conns, expressions, operation = "opaque exchange", .aggregate = wrong),
    error = identity)
  expect_s3_class(condition, "error")
  expect_match(conditionMessage(condition), "malformed or misassociated")
  expect_false(grepl("attacker", conditionMessage(condition), fixed = TRUE))
})

test_that("oversize is terminal in one-cycle and assignment transports", {
  tag <- paste0(
    "[dsvert_resource_oversize:v1] resource_oversize: ",
    "fixed transport geometry cannot fit")
  conns <- list(s1 = list())
  expressions <- list(s1 = call("exchangeDS"))
  callback_failure <- function(conns, expr, async, error, errors.print) {
    error("s1", tag)
    list(s1 = NULL)
  }
  cycle <- tryCatch(.dsvert_fanout_cycle(
    conns, expressions, .aggregate = callback_failure), error = identity)
  expect_s3_class(cycle, "dsvert_resource_oversize")
  expect_false(cycle$retryable)

  top_failure <- function(...) stop(tag, call. = FALSE)
  cycle_top <- tryCatch(.dsvert_fanout_cycle(
    conns, expressions, .aggregate = top_failure), error = identity)
  expect_s3_class(cycle_top, "dsvert_resource_oversize")
  expect_false(cycle_top$retryable)

  assignment_failure <- function(
      conns, symbol, expr, async, success, error, errors.print) {
    error("s1", tag)
    invisible(NULL)
  }
  assigned <- tryCatch(.dsvert_assign_strict(
    conns, "value", expressions, .assign = assignment_failure),
    error = identity)
  expect_s3_class(assigned, "dsvert_resource_oversize")
  expect_false(assigned$retryable)
})

test_that("only tagged backpressure retries an explicit idempotent endpoint", {
  tag <- paste0(
    "[dsvert_resource_backpressure:v1] resource_backpressure: ",
    "relay session spool capacity is currently full")
  parsed <- .dsvert_client_parse_resource_backpressure(tag)
  expect_s3_class(parsed, "dsvert_resource_backpressure")
  expect_identical(parsed$code, "resource_backpressure")
  expect_true(parsed$retryable)
  expect_null(.dsvert_client_parse_resource_backpressure(
    "resource_backpressure without the fixed public tag"))

  conns <- list(s1 = list())
  expr <- call(
    "psiPaddedRelayExchangeDS", request_json = "frozen",
    session_id = "12345678-1234-4234-9234-123456789abc",
    outbound_operation_id = "op_11111111111111111111111111111111",
    terminal_receipt_b64url = "")
  calls <- 0L
  aggregate <- function(conns, expr, async, error, errors.print) {
    calls <<- calls + 1L
    if (calls == 1L) {
      error("s1", tag)
      return(list(s1 = NULL))
    }
    list(s1 = list(ok = TRUE))
  }
  result <- testthat::with_mocked_bindings(
    .dsvert_aggregate_strict(
      conns, expr, operation = "idempotent padded relay",
      .aggregate = aggregate),
    .dsvert_retry_sleep = function(...) NULL,
    .dsvert_retry_jitter = function() 1,
    .package = "dsVertClient")
  expect_identical(result, list(s1 = list(ok = TRUE)))
  expect_identical(calls, 2L)

  oversize_tag <- paste0(
    "[dsvert_resource_oversize:v1] resource_oversize: ",
    "relay payload cannot fit its fixed policy")
  parsed_oversize <- .dsvert_client_parse_resource_oversize(oversize_tag)
  expect_s3_class(parsed_oversize, "dsvert_resource_oversize")
  expect_identical(parsed_oversize$code, "resource_oversize")
  expect_false(parsed_oversize$retryable)
  expect_null(.dsvert_client_parse_resource_oversize(
    "resource_oversize without the fixed public tag"))

  for (failure_message in "signature or pinned contract mismatch") {
    calls <- 0L
    contract_failure <- function(conns, expr, async, error, errors.print) {
      calls <<- calls + 1L
      error("s1", failure_message)
      list(s1 = NULL)
    }
    expect_error(.dsvert_aggregate_strict(
      conns, expr, operation = "idempotent padded relay",
      .aggregate = contract_failure), "transport failed")
    expect_identical(calls, 1L, info = failure_message)
  }

  calls <- 0L
  terminal_oversize <- function(conns, expr, async, error, errors.print) {
    calls <<- calls + 1L
    error("s1", oversize_tag)
    list(s1 = NULL)
  }
  condition <- tryCatch(.dsvert_aggregate_strict(
    conns, expr, operation = "idempotent padded relay",
    .aggregate = terminal_oversize), error = identity)
  expect_s3_class(condition, "dsvert_resource_oversize")
  expect_identical(condition$code, "resource_oversize")
  expect_false(condition$retryable)
  expect_identical(calls, 1L)
})

test_that("oversized DSI expressions fail locally before connector side effects", {
  conns <- list(s1 = list())
  calls <- 0L
  aggregate <- function(...) {
    calls <<- calls + 1L
    list(s1 = TRUE)
  }
  withr::local_options(list(dsvert.dsi_max_expression_bytes = 32 * 1024L))

  condition <- tryCatch(
    dsVertClient:::.dsvert_aggregate_strict(
      conns, call("phaseDS", payload = strrep("A", 40 * 1024L)),
      .aggregate = aggregate),
    error = identity)
  expect_s3_class(condition, "dsvert_resource_oversize")
  expect_identical(condition$code, "resource_oversize")
  expect_false(condition$retryable)
  expect_identical(calls, 0L)

  expect_identical(
    dsVertClient:::.dsvert_aggregate_strict(
      conns, call("phaseDS", payload = strrep("A", 24 * 1024L)),
      result_contract = "logical_true", .aggregate = aggregate),
    list(s1 = TRUE))
  expect_identical(calls, 1L)
})

test_that("portable 480 KiB raw Base64 frames remain below DSLite boundary", {
  raw_bytes <- as.raw((seq_len(480 * 1024L) - 1L) %% 251L)
  encoded <- gsub("[\r\n]", "", jsonlite::base64_enc(raw_bytes))
  encoded <- chartr("+/", "-_", sub("=+$", "", encoded))
  expr <- call("exchangeDS", payload = encoded, offset = 0,
               total_bytes = as.numeric(length(raw_bytes)))
  size <- dsVertClient:::.dsvert_validate_dsi_expression_sizes(expr)

  expect_lt(unname(size), 768 * 1024L)
  expect_equal(nchar(encoded, type = "bytes"),
               ceiling(4 * length(raw_bytes) / 3))
})

test_that("portable 480 KiB raw frame crosses a real DSLite parser", {
  skip_if_not_installed("DSLite")
  server_name <- paste0("dsvert_dsi_boundary_", Sys.getpid())
  server <- DSLite::newDSLiteServer(tables = list(t = data.frame(x = 1)))
  server$aggregateMethod("dsvertBoundaryEchoDS", function(payload) {
    list(encoded_bytes = nchar(payload, type = "bytes"))
  })
  assign(server_name, server, envir = .GlobalEnv)
  withr::defer(rm(list = server_name, envir = .GlobalEnv))
  builder <- DSI::newDSLoginBuilder()
  builder$append(
    server = "site1", url = server_name, table = "t",
    driver = "DSLiteDriver")
  conns <- DSI::datashield.login(builder$build(), assign = FALSE)
  withr::defer(DSI::datashield.logout(conns))

  raw_bytes <- as.raw((seq_len(480 * 1024L) - 1L) %% 251L)
  encoded <- gsub("[\r\n]", "", jsonlite::base64_enc(raw_bytes))
  encoded <- chartr("+/", "-_", sub("=+$", "", encoded))
  result <- dsVertClient:::.dsvert_aggregate_strict(
    conns, call("dsvertBoundaryEchoDS", payload = encoded),
    operation = "DSLite request-boundary regression")

  expect_identical(
    as.numeric(result$site1$encoded_bytes), as.numeric(nchar(encoded)))
})

test_that("a 768 KiB source-window response crosses real DSLite separately", {
  skip_if_not_installed("DSLite")
  server_name <- paste0("dsvert_dsi_output_boundary_", Sys.getpid())
  server <- DSLite::newDSLiteServer(tables = list(t = data.frame(x = 1)))
  server$aggregateMethod("dsvertOutputBoundaryDS", function(bytes) {
    strrep("A", as.integer(bytes))
  })
  assign(server_name, server, envir = .GlobalEnv)
  withr::defer(rm(list = server_name, envir = .GlobalEnv))
  builder <- DSI::newDSLoginBuilder()
  builder$append(
    server = "site1", url = server_name, table = "t",
    driver = "DSLiteDriver")
  conns <- DSI::datashield.login(builder$build(), assign = FALSE)
  withr::defer(DSI::datashield.logout(conns))

  response_bytes <- 768L * 1024L
  expression <- call("dsvertOutputBoundaryDS", bytes = response_bytes)
  expect_lt(.dsvert_validate_dsi_expression_sizes(expression), 1024L)
  result <- .dsvert_aggregate_strict(
    conns, expression, operation = "DSLite output-boundary regression")
  expect_identical(nchar(result$site1, type = "bytes"), response_bytes)
})

test_that("a real old DSLite server declines sync negotiation without poison", {
  skip_if_not_installed("DSLite")
  .dsvert_dsi_clear_poisoned_sessions()
  withr::defer(.dsvert_dsi_clear_poisoned_sessions())
  server_name <- paste0("dsvert_old_sync_negotiation_", Sys.getpid())
  server <- DSLite::newDSLiteServer(tables = list(t = data.frame(x = 1)))
  server$aggregateMethod(
    "dsvertDPCapsuleSourceTicketDS",
    function(manifest_json) paste0("legacy:", manifest_json))
  assign(server_name, server, envir = .GlobalEnv)
  withr::defer(rm(list = server_name, envir = .GlobalEnv))
  builder <- DSI::newDSLoginBuilder()
  builder$append(
    server = "site1", url = server_name, table = "t",
    driver = "DSLiteDriver")
  conns <- DSI::datashield.login(builder$build(), assign = FALSE)
  withr::defer(DSI::datashield.logout(conns))

  legacy <- call("dsvertDPCapsuleSourceTicketDS", manifest_json = "m1")
  legacy_wire <- paste0("legacy:", .dsvert_dsi_text_encode("m1"))
  expect_identical(
    .dsvert_aggregate_strict(conns, legacy)$site1, legacy_wire)
  negotiation <- list(site1 = call(
    "dsvertDPCapsuleSourceTicketDS", manifest_json = "m1",
    transport_contract =
      .DSVERT_CLIENT_DP_CAPSULE_SOURCE_WINDOW_VERSION))
  expect_null(.dsvert_dp_capsule_source_optional_negotiation(
    conns, negotiation, DSI::datashield.aggregate))
  expect_identical(
    .dsvert_aggregate_strict(conns, legacy)$site1, legacy_wire)
})

test_that("ambiguous sync execute failures poison the affected DSI handle", {
  for (message in c(
      "connection reset after execute",
      "Internal server error: response failed after execute")) {
    .dsvert_dsi_clear_poisoned_sessions()
    withr::defer(.dsvert_dsi_clear_poisoned_sessions())
    key <- paste0("ambiguous-sync-handle-", digest::digest(
      message, algo = "sha256", serialize = FALSE))
    local_mocked_bindings(
      .dsvert_dsi_inspect_connection = function(connection) {
        list(kind = "armadillo",
             endpoint = list(scope = "https://example.test"))
      },
      .dsvert_dsi_job_session_key = function(connection) key,
      .package = "dsVertClient")
    condition <- testthat::with_mocked_bindings(
      tryCatch(
        .dsvert_dsi_direct_aggregate(
          list(site = list()), call("publicCapabilityDS"), async = FALSE,
          require_settled_sync_failure = TRUE),
        dsvert_dsi_poisoned_session = identity),
      datashield.sessions = function(...) invisible(NULL),
      dsIsAsync = function(...) list(aggregate = TRUE),
      dsAggregate = function(...) stop(simpleError(message)),
      dsFetch = function(...) stop("an ambiguous job must not be fetched"),
      .package = "DSI")

    expect_s3_class(condition, "dsvert_dsi_poisoned_session")
    expect_match(conditionMessage(condition), "fresh DSI login connections",
                 info = message)
    expect_true(.dsvert_dsi_session_is_poisoned(key), info = message)
  }
})

test_that("terminal Armadillo application rejection permits scalar fallback", {
  .dsvert_dsi_clear_poisoned_sessions()
  withr::defer(.dsvert_dsi_clear_poisoned_sessions())
  key <- "settled-sync-handle"
  failures <- character()
  local_mocked_bindings(
    .dsvert_dsi_inspect_connection = function(connection) {
      list(kind = "armadillo", endpoint = list(scope = "https://example.test"))
    },
    .dsvert_dsi_job_session_key = function(connection) key,
    .package = "dsVertClient")
  result <- testthat::with_mocked_bindings(
    .dsvert_dsi_direct_aggregate(
      list(site = list()), call("publicCapabilityDS"), async = FALSE,
      error = function(site, message) failures <<- c(failures, site),
      require_settled_sync_failure = TRUE),
    datashield.sessions = function(...) invisible(NULL),
    dsIsAsync = function(...) list(aggregate = TRUE),
    dsAggregate = function(...) stop(simpleError(paste0(
      "Command 'publicCapabilityDS' failed on armadillo: Error whilst ",
      "evaluating request -> unused argument ",
      "(transport_contract = window)"))),
    dsFetch = function(...) stop("a rejected command must not be fetched"),
    .package = "DSI")

  expect_identical(result, list(site = NULL))
  expect_identical(failures, "site")
  expect_false(.dsvert_dsi_session_is_poisoned(key))
})

test_that("strict assignments use one named asynchronous fan-out", {
  conns <- list(s1 = list(), s2 = list())
  values <- list(
    s1 = call("phaseAssignDS", value = 1L),
    s2 = call("phaseAssignDS", value = 2L))
  captured <- NULL
  assign_impl <- function(conns, symbol, expr, async, success, error,
                          errors.print) {
    captured <<- list(
      conns = conns, symbol = symbol, expr = expr, async = async,
      errors.print = errors.print)
    for (site in names(conns)) success(site)
    invisible(NULL)
  }

  expect_true(dsVertClient:::.dsvert_assign_strict(
    conns, "result", values, operation = "assignment phase",
    .assign = assign_impl))
  expect_identical(names(captured$conns), c("s1", "s2"))
  expect_identical(captured$symbol, "result")
  expect_identical(captured$expr, values)
  expect_true(captured$async)
  expect_false(captured$errors.print)
})

test_that("strict assignments mask partial connector errors", {
  conns <- list(s1 = list(), s2 = list())
  partial <- function(conns, symbol, expr, async, success, error,
                      errors.print) {
    success("s1")
    error("s2", "SECRET assignment detail")
    invisible(NULL)
  }
  condition <- tryCatch(
    dsVertClient:::.dsvert_assign_strict(
      conns, "result", call("phaseAssignDS"),
      operation = "assignment phase", .assign = partial),
    error = identity)

  expect_s3_class(condition, "error")
  expect_match(conditionMessage(condition), "partial assignment")
  expect_false(grepl("SECRET", conditionMessage(condition), fixed = TRUE))
})

test_that("strict assignments require one success callback per site", {
  conns <- list(s1 = list(), s2 = list())
  missing <- function(conns, symbol, expr, async, success, error,
                      errors.print) {
    success("s1")
    invisible(NULL)
  }
  duplicate <- function(conns, symbol, expr, async, success, error,
                        errors.print) {
    success("s1")
    success("s1")
    success("s2")
    invisible(NULL)
  }
  unknown <- function(conns, symbol, expr, async, success, error,
                      errors.print) {
    success("s1")
    success("unknown")
    success("s2")
    invisible(NULL)
  }
  for (implementation in list(missing, duplicate, unknown)) {
    expect_error(.dsvert_assign_strict(
      conns, "result", call("phaseAssignDS"),
      .assign = implementation), "partial assignment")
  }
})

test_that("strict assignments route distinct expressions through real DSLite", {
  skip_if_not_installed("DSLite")
  server_name <- paste0("dsvert_assign_fanout_", Sys.getpid())
  server <- DSLite::newDSLiteServer(tables = list(
    site1 = data.frame(x = 1), site2 = data.frame(x = 2)))
  server$assignMethod("dsvertAssignValueDS", function(value) {
    as.integer(value)
  })
  assign(server_name, server, envir = .GlobalEnv)
  withr::defer(rm(list = server_name, envir = .GlobalEnv))
  builder <- DSI::newDSLoginBuilder()
  builder$append(
    server = "site1", url = server_name, table = "site1",
    driver = "DSLiteDriver")
  builder$append(
    server = "site2", url = server_name, table = "site2",
    driver = "DSLiteDriver")
  conns <- DSI::datashield.login(builder$build(), assign = FALSE)
  withr::defer(DSI::datashield.logout(conns))

  expect_true(dsVertClient:::.dsvert_assign_strict(
    conns, "assigned_value", list(
      site1 = call("dsvertAssignValueDS", value = 11L),
      site2 = call("dsvertAssignValueDS", value = 22L)),
    operation = "real DSLite assignment fan-out"))
  observed <- vapply(conns, function(connection) {
    server$getSessionData(methods::slot(connection, "sid"), "assigned_value")
  }, integer(1L))
  expect_identical(unname(observed), c(11L, 22L))
})

test_that("strict assignments preserve the Opal S4 success path", {
  skip_if_not_installed("DSOpal")
  loadNamespace("DSOpal")
  opal <- new.env(parent = emptyenv())
  opal$url <- "https://opal.example.org"
  opal$config <- list(options = list())
  opal$rid <- "opal-assignment-session"
  class(opal) <- "opal"
  connection <- methods::new(
    "OpalConnection", name = "opal", opal = opal)
  conns <- list(opal = connection)
  launched <- 0L
  result <- testthat::with_mocked_bindings(
    .dsvert_assign_strict(
      conns, "assigned_value",
      list(opal = call("identity", 11L))),
    datashield.assign.expr = function(
        conns, symbol, expr, async, success, error, errors.print) {
      launched <<- launched + 1L
      success("opal")
      invisible(NULL)
    },
    .package = "DSI")

  expect_true(result)
  expect_identical(launched, 1L)
  expect_false(.dsvert_dsi_session_is_poisoned(
    .dsvert_dsi_job_session_key(connection)))
})

test_that("best-effort cleanup is async and never leaks or throws", {
  conns <- list(s1 = list(), s2 = list())
  captured <- NULL
  cleanup <- function(conns, expr, async, error, errors.print) {
    captured <<- list(expr = expr, async = async, errors.print = errors.print)
    error("s2", "SECRET cleanup detail")
    stop("SECRET top-level cleanup detail")
  }

  expect_false(dsVertClient:::.dsvert_cleanup_best_effort(
    conns, call("cleanupDS"), .aggregate = cleanup))
  expect_named(captured$expr, c("s1", "s2"))
  expect_true(captured$async)
  expect_false(captured$errors.print)
})
test_that("release-instance retry parser accepts only current safe actions", {
  fresh <- .dsvert_client_parse_release_instance_retry(paste(
    "server failed",
    "[dsvert_retry_current_release_instance:new_release_instance]"))
  unpublished <- .dsvert_client_parse_release_instance_retry(paste(
    "server failed",
    "[dsvert_retry_current_release_instance:retry_unpublished_instance]"))
  retired <- .dsvert_client_parse_release_instance_retry(
    "[dsvert_retry_current_release_instance:rematerialize_sticky]")
  expect_s3_class(fresh, "dsvert_release_instance_retry")
  expect_identical(fresh$retry_action, "new_release_instance")
  expect_s3_class(unpublished, "dsvert_release_instance_retry")
  expect_identical(
    unpublished$retry_action, "retry_unpublished_instance")
  expect_null(retired)
})
