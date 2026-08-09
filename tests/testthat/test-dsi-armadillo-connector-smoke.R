test_that("real Armadillo S4 methods cover async, sync, and poisoned sessions", {
  skip_if_not_installed("callr")
  skip_if_not_installed("DSMolgenisArmadillo", minimum_version = "3.0.2")
  skip_if_not_installed("httr")
  skip_if_not_installed("httpuv")

  expected <- list(
    version = "dsvert-armadillo-stub-v1",
    ok = TRUE,
    value = 42L)
  port <- httpuv::randomPort()
  server <- callr::r_bg(
    function(port, value) {
      execute_count <- 0L
      assignment_count <- 0L
      app <- list(call = function(request) {
        path <- request$PATH_INFO
        if (identical(path, "/ready") ||
            identical(path, "/actuator/info")) {
          return(list(
            status = 200L,
            headers = list("Content-Type" = "application/json"),
            body = "{\"status\":\"UP\"}"))
        }
        if (identical(path, "/counts")) {
          return(list(
            status = 200L,
            headers = list("Content-Type" = "application/json"),
            body = paste0(
              "{\"execute_count\":", execute_count,
              ",\"assignment_count\":", assignment_count, "}")))
        }
        if (identical(request$REQUEST_METHOD, "POST") &&
            identical(path, "/execute")) {
          execute_count <<- execute_count + 1L
          if (identical(execute_count, 4L)) {
            return(list(
              status = 500L,
              headers = list(
                "Content-Type" = "text/plain; charset=utf-8"),
              body = "accepted-but-response-ambiguous"))
          }
          return(list(status = 200L, headers = list(), body = raw()))
        }
        if (identical(request$REQUEST_METHOD, "POST") &&
            grepl("^/symbols/", path)) {
          assignment_count <<- assignment_count + 1L
          if (identical(assignment_count, 1L)) {
            return(list(
              status = 500L,
              headers = list(
                "Content-Type" = "text/plain; charset=utf-8"),
              body = "accepted-assignment-with-ambiguous-response"))
          }
          return(list(status = 200L, headers = list(), body = raw()))
        }
        if (identical(request$REQUEST_METHOD, "GET") &&
            identical(path, "/lastcommand")) {
          return(list(
            status = 200L,
            headers = list("Content-Type" = "application/json"),
            body = "{\"status\":\"COMPLETED\"}"))
        }
        if (identical(request$REQUEST_METHOD, "GET") &&
            identical(path, "/lastresult")) {
          return(list(
            status = 200L,
            headers = list("Content-Type" = "application/octet-stream"),
            body = serialize(value, NULL, version = 3L)))
        }
        list(status = 404L, headers = list(), body = raw())
      })
      listener <- httpuv::startServer("127.0.0.1", port, app)
      on.exit(httpuv::stopServer(listener), add = TRUE)
      repeat httpuv::service(100)
    },
    args = list(port = port, value = expected),
    supervise = TRUE,
    stdout = "|", stderr = "|")
  on.exit({
    if (server$is_alive()) server$kill()
    try(server$wait(timeout = 2000L), silent = TRUE)
  }, add = TRUE)

  endpoint <- paste0("http://127.0.0.1:", port)
  ready <- FALSE
  for (attempt in seq_len(100L)) {
    if (!server$is_alive()) {
      stop(paste(c(server$read_all_output(), server$read_all_error()),
                 collapse = "\n"), call. = FALSE)
    }
    response <- tryCatch(
      httr::GET(paste0(endpoint, "/ready"), httr::timeout(0.2)),
      error = function(e) NULL)
    if (!is.null(response) && identical(httr::status_code(response), 200L)) {
      ready <- TRUE
      break
    }
    Sys.sleep(0.02)
  }
  expect_true(ready)

  loadNamespace("DSMolgenisArmadillo")
  new_connection <- function() {
    methods::new(
      "ArmadilloConnection", name = "armadillo",
      handle = httr::handle(endpoint), user = "researcher",
      cookies = list(), token = "")
  }
  connection <- new_connection()
  conns <- list(armadillo = connection)
  testthat::local_mocked_bindings(
    .dsvert_dsi_test_loopback_allowed = function() TRUE,
    .package = "dsVertClient")
  .dsvert_dsi_clear_poisoned_sessions()
  withr::defer(.dsvert_dsi_clear_poisoned_sessions())

  cycle <- .dsvert_fanout_cycle(
    conns, list(armadillo = call("identity", 1L)),
    operation = "Armadillo connector smoke")
  expect_identical(cycle$state, "ok")
  expect_identical(cycle$result$armadillo, expected)

  strict <- .dsvert_aggregate_strict(
    conns, call("identity", 2L), operation = "Armadillo strict smoke")
  expect_identical(strict$armadillo, expected)

  synchronous <- .dsvert_dsi_direct_aggregate(
    conns, call("identity", 3L), async = FALSE,
    require_settled_sync_failure = TRUE)
  expect_identical(synchronous$armadillo, expected)

  ambiguous <- tryCatch(
    .dsvert_dsi_direct_aggregate(
      conns, call("identity", 4L), async = FALSE,
      require_settled_sync_failure = TRUE),
    dsvert_dsi_poisoned_session = identity)
  expect_s3_class(ambiguous, "dsvert_dsi_poisoned_session")
  session_key <- .dsvert_dsi_job_session_key(connection)
  expect_true(.dsvert_dsi_session_is_poisoned(session_key))

  blocked <- tryCatch(
    .dsvert_aggregate_strict(
      conns, call("identity", 5L), operation = "poisoned Armadillo smoke"),
    dsvert_dsi_poisoned_session = identity)
  expect_s3_class(blocked, "dsvert_dsi_poisoned_session")

  fresh_connection <- new_connection()
  fresh_key <- .dsvert_dsi_job_session_key(fresh_connection)
  expect_false(identical(fresh_key, session_key))
  expect_false(.dsvert_dsi_session_is_poisoned(fresh_key))
  fresh <- .dsvert_dsi_direct_aggregate(
    list(armadillo = fresh_connection), call("identity", 6L), async = FALSE,
    require_settled_sync_failure = TRUE)
  expect_identical(fresh$armadillo, expected)

  assignment_counts <- function() {
    response <- httr::GET(paste0(endpoint, "/counts"))
    httr::stop_for_status(response)
    httr::content(response, as = "parsed", type = "application/json")
  }
  .dsvert_dsi_poison_sessions(session_key)
  pre_poisoned <- tryCatch(
    .dsvert_assign_strict(
      conns, "aligned", list(armadillo = call("identity", 7L))),
    dsvert_dsi_poisoned_session = identity)
  expect_s3_class(pre_poisoned, "dsvert_dsi_poisoned_session")
  expect_identical(assignment_counts()$assignment_count, 0L)

  .dsvert_dsi_clear_poisoned_sessions()
  ambiguous_assignment <- tryCatch(
    .dsvert_assign_strict(
      conns, "aligned", list(armadillo = call("identity", 8L))),
    dsvert_dsi_poisoned_session = identity)
  expect_s3_class(
    ambiguous_assignment, "dsvert_dsi_poisoned_session")
  expect_true(.dsvert_dsi_session_is_poisoned(session_key))
  expect_identical(assignment_counts()$assignment_count, 1L)

  blocked_assignment <- tryCatch(
    .dsvert_assign_strict(
      conns, "aligned", list(armadillo = call("identity", 9L))),
    dsvert_dsi_poisoned_session = identity)
  expect_s3_class(blocked_assignment, "dsvert_dsi_poisoned_session")
  expect_identical(assignment_counts()$assignment_count, 1L)

  unstable <- testthat::with_mocked_bindings(
    tryCatch(.dsvert_assign_strict(
      list(armadillo = fresh_connection), "aligned",
      list(armadillo = call("identity", 10L))), error = identity),
    .dsvert_dsi_job_session_key = function(connection) NULL,
    .package = "dsVertClient")
  expect_s3_class(unstable, "error")
  expect_match(conditionMessage(unstable), "stable authenticated session")
  expect_identical(assignment_counts()$assignment_count, 1L)

  expect_true(.dsvert_assign_strict(
    list(armadillo = fresh_connection), "aligned",
    list(armadillo = call("identity", 11L))))
  expect_identical(assignment_counts()$assignment_count, 2L)
})
