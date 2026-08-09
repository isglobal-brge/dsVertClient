opal_connection <- function(url = "https://opal.example.org",
                            options = list()) {
  structure(list(opal = list(url = url, config = list(options = options))),
            class = c("OpalConnection", "list"))
}

armadillo_connection <- function(url = "https://armadillo.example.org") {
  structure(list(handle = list(url = url)),
            class = c("ArmadilloConnection", "list"))
}

generic_connection <- function(url = "https://generic.example.org") {
  structure(list(url = url), class = c("UnknownDSIConnection", "list"))
}

test_that("outer DSI validation accepts in-process and verified Opal links", {
  dslite <- structure(list(), class = c("DSLiteConnection", "list"))
  expect_invisible(.dsvert_validate_dsi_transport_security(
    list(local = dslite)))
  expect_invisible(.dsvert_validate_dsi_transport_security(
    list(opal = opal_connection())))

  withr::local_options(list(
    httr_config = list(options = list(
      ssl_verifypeer = 0, ssl_verifyhost = 0))))
  expect_error(
    .dsvert_validate_dsi_transport_security(
      list(opal = opal_connection())),
    "certificate verification disabled")
  expect_invisible(.dsvert_validate_dsi_transport_security(list(
    opal = opal_connection(options = list(
      ssl_verifypeer = 1, ssl_verifyhost = 2)))))
  expect_error(
    .dsvert_validate_dsi_transport_security(
      list(armadillo = armadillo_connection())),
    "certificate verification disabled")
})

test_that("recognized Armadillo HTTPS is plug-and-play and generic stays denied", {
  expect_invisible(.dsvert_validate_dsi_transport_security(
    list(armadillo = armadillo_connection())))
  expect_error(
    .dsvert_validate_dsi_transport_security(
      list(generic = generic_connection())),
    "unsupported generic DSI connector")
})

test_that("the installed DataSHIELD surface has no unguarded base aliases", {
  skip_if_not_installed("DSLite")
  skip_if_not_installed("dsVert")

  configuration <- DSLite::defaultDSConfiguration(include = "dsVert")
  forbidden <- c("c", "list", "numeric", "character")
  expect_length(
    intersect(forbidden, configuration$AggregateMethods$name), 0L)

  server_name <- paste0("dsvert_alias_surface_", Sys.getpid())
  server <- DSLite::newDSLiteServer(tables = list(
    t = data.frame(canary = c(101, 202, 303))))
  server$config(configuration)
  assign(server_name, server, envir = .GlobalEnv)
  withr::defer(rm(list = server_name, envir = .GlobalEnv))

  builder <- DSI::newDSLoginBuilder()
  builder$append(
    server = "site", url = server_name, table = "t",
    driver = "DSLiteDriver")
  conns <- DSI::datashield.login(
    builder$build(), assign = TRUE, symbol = "D")
  withr::defer(DSI::datashield.logout(conns))

  expect_error(
    DSI::datashield.aggregate(
      conns, call(name = "c", quote(D$canary)),
      errors.print = FALSE),
    regexp = "DataSHIELD errors")
  expect_match(
    paste(unlist(DSI::datashield.errors(), use.names = FALSE),
          collapse = "\n"),
    "does not allow expression: c", fixed = TRUE)
})

test_that("Armadillo rejects global TLS downgrade and malformed httr config", {
  skip_if_not_installed("httr")
  connection <- armadillo_connection()

  withr::local_options(list(httr_config = NULL))
  expect_invisible(.dsvert_validate_dsi_transport_security(
    list(armadillo = connection)))

  withr::local_options(list(httr_config = httr::config()))
  expect_invisible(.dsvert_validate_dsi_transport_security(
    list(armadillo = connection)))

  for (config in list(
      httr::config(ssl_verifypeer = 0),
      httr::config(ssl_verifyhost = 0),
      httr::config(ssl_verifyhost = 1))) {
    withr::local_options(list(httr_config = config))
    expect_error(.dsvert_validate_dsi_transport_security(
      list(armadillo = connection)), "certificate verification disabled")
  }

  withr::local_options(list(httr_config = httr::config(
    ssl_verifypeer = 1, ssl_verifyhost = 2)))
  expect_invisible(.dsvert_validate_dsi_transport_security(
    list(armadillo = connection)))

  for (config in list(
      "not-an-httr-request",
      list(options = "not-a-list"),
      list(unrelated = TRUE))) {
    withr::local_options(list(httr_config = config))
    expect_error(.dsvert_validate_dsi_transport_security(
      list(armadillo = connection)), "cannot be inspected")
  }
})

test_that("remote HTTP is rejected without exposing endpoint credentials", {
  error <- expect_error(.dsvert_validate_dsi_transport_security(list(
    armadillo = armadillo_connection(
      "http://operator:secret@armadillo.example.org/path?token=hidden"))))
  expect_match(conditionMessage(error), "does not use verified HTTPS")
  expect_false(grepl("operator|secret|token|hidden", conditionMessage(error)))
})

test_that("endpoint identity retains port and base path without credentials", {
  first <- .dsvert_dsi_endpoint(
    "https://operator:secret@example.org:8443/opal-a?token=hidden")
  second <- .dsvert_dsi_endpoint("https://example.org:9443/opal-a")
  ipv6 <- .dsvert_dsi_endpoint("https://[::1]:8443/armadillo")
  expect_identical(first$host, "example.org")
  expect_identical(first$port, "8443")
  expect_identical(first$scope, "https://example.org:8443/opal-a")
  expect_identical(second$scope, "https://example.org:9443/opal-a")
  expect_identical(ipv6$scope, "https://[::1]:8443/armadillo")
  expect_false(grepl("operator|secret|token|hidden", first$scope))
  expect_null(.dsvert_dsi_endpoint("https://example.org:99999/path"))
})

test_that("DSMolgenisArmadillo S4 endpoint metadata is plug-and-play", {
  skip_if_not_installed("DSMolgenisArmadillo")
  loadNamespace("DSMolgenisArmadillo")
  skip_if_not_installed("httr")
  expect_true(
    utils::packageVersion("DSMolgenisArmadillo") >=
      numeric_version("3.0.2"))
  handle <- httr::handle("https://armadillo.example.org")
  connection <- methods::new(
    "ArmadilloConnection", name = "armadillo", handle = handle,
    user = "researcher", cookies = list(), token = "test-token")

  withr::local_options(list(dsvert.dsi_tls_attested = NULL))
  expect_null(getOption("dsvert.dsi_tls_attested"))
  expect_invisible(.dsvert_validate_dsi_transport_security(
    list(armadillo = connection)))
  lifecycle <- .dsvert_dsi_armadillo_handle_identity(connection)
  expect_match(lifecycle, "^curl-handle:h_[0-9a-f]{32}@")
  session_key <- .dsvert_dsi_job_session_key(connection)
  expect_match(session_key, "^[0-9a-f]{64}$")
  refreshed <- connection
  refreshed@token <- "fresh-test-token"
  expect_identical(
    .dsvert_dsi_armadillo_handle_identity(refreshed), lifecycle)
  expect_identical(
    .dsvert_dsi_job_session_key(refreshed), session_key)
  fresh_handle <- connection
  fresh_handle@handle <- httr::handle("https://armadillo.example.org")
  expect_false(identical(
    .dsvert_dsi_armadillo_handle_identity(fresh_handle), lifecycle))
  expect_false(identical(
    .dsvert_dsi_job_session_key(fresh_handle), session_key))
  expect_identical(DSI::dsIsAsync(connection)$aggregate, TRUE)
  for (generic in c(
      "dsAggregate", "dsIsCompleted", "dsFetch", "dsKeepAlive")) {
    signature <- if (generic %in% c("dsAggregate", "dsKeepAlive")) {
      "ArmadilloConnection"
    } else {
      "ArmadilloResult"
    }
    expect_s4_class(
      methods::selectMethod(generic, signature, optional = TRUE),
      "MethodDefinition")
  }
  connection@handle$url <- "http://armadillo.example.org"
  expect_error(.dsvert_validate_dsi_transport_security(
    list(armadillo = connection)), "does not use verified HTTPS")
})

test_that("plaintext loopback cannot be enabled by client configuration", {
  connection <- generic_connection("http://127.0.0.1:8080")
  withr::local_options(list(dsvert.dsi_allow_test_loopback = TRUE))
  withr::local_envvar(c(DSVERT_TEST_MODE = "true"))
  expect_error(
    .dsvert_validate_dsi_transport_security(list(local = connection)),
    "does not use verified HTTPS")
})

test_that("invalid transport policy and malformed connection sets fail closed", {
  expect_error(
    .dsvert_validate_dsi_transport_security(
      structure(list(opal_connection()), names = "")),
    "unique logical site names")
  expect_error(
    .dsvert_validate_dsi_transport_security(
      list(site = generic_connection("not-a-url"))),
    "no inspectable DSI endpoint")
})

test_that("strict aggregate validates real DSI calls but not injected tests", {
  insecure <- list(site = generic_connection("http://remote.example.org"))
  expect_error(
    .dsvert_validate_real_dsi_transport(
      insecure, DSI::datashield.aggregate),
    "does not use verified HTTPS")
  expect_invisible(.dsvert_validate_real_dsi_transport(
    insecure, function(...) list(site = TRUE)))
})

test_that("best-effort cleanup never crosses an insecure outer transport", {
  insecure <- list(site = generic_connection("http://remote.example.org"))

  expect_false(dsVertClient:::.dsvert_cleanup_best_effort(
    insecure, call("cleanupDS"),
    .aggregate = DSI::datashield.aggregate))
})
