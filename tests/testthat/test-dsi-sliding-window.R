test_that("portable DSI/DSLite cannot expose an aggregate sliding window", {
  skip_if_not_installed("DSLite")

  state <- new.env(parent = emptyenv())
  state$calls <- 0L
  state_name <- paste0("dsvert_window_probe_state_", Sys.getpid())
  assign(state_name, state, envir = .GlobalEnv)
  on.exit(rm(list = state_name, envir = .GlobalEnv), add = TRUE)
  server <- DSLite::newDSLiteServer(tables = list(t = data.frame(x = 1)))
  method_name <- paste0("dsvertWindowProbeDS", Sys.getpid())
  server$aggregateMethod(method_name, function(state_name, value) {
    state <- get(state_name, envir = .GlobalEnv, inherits = FALSE)
    state$calls <- state$calls + 1L
    value
  })
  object_name <- paste0("dsvert_window_probe_", Sys.getpid())
  assign(object_name, server, envir = .GlobalEnv)
  on.exit(rm(list = object_name, envir = .GlobalEnv), add = TRUE)

  builder <- DSI::newDSLoginBuilder()
  builder$append(
    server = "site", url = object_name, table = "t",
    driver = "DSLiteDriver")
  conns <- DSI::datashield.login(builder$build(), assign = FALSE)
  on.exit(DSI::datashield.logout(conns), add = TRUE)

  expect_false(DSI::dsIsAsync(conns[[1L]])$aggregate)
  handle <- DSI::dsAggregate(
    conns[[1L]], call(
      name = method_name, state_name = state_name, value = 1L), async = TRUE)
  # DSLite has executed the method before returning even when async=TRUE.
  expect_identical(state$calls, 1L)
  expect_true(DSI::dsIsCompleted(handle))
  expect_identical(DSI::dsFetch(handle), 1L)

  fetched <- DSI::datashield.aggregate(
    conns, call(
      name = method_name, state_name = state_name, value = 2L), async = TRUE,
    errors.print = FALSE)
  expect_identical(fetched, list(site = 2L))
  expect_false(inherits(fetched[[1L]], "DSResult"))
  expect_identical(state$calls, 2L)
})
