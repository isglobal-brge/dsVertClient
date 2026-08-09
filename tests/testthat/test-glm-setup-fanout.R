test_that("GLM setup uses one DSI fan-out call per independent phase", {
  calls <- list()
  fake_aggregate <- function(conns, expr, ...) {
    calls[[length(calls) + 1L]] <<- list(
      servers = names(conns), expr = expr)

    if (is.call(expr)) {
      method <- as.character(expr[[1L]])
      if (identical(method, "glmRing63TransportInitDS")) {
        return(stats::setNames(lapply(names(conns), function(server) list(
          transport_pk = paste0("transport-", server),
          identity_pk = paste0("identity-", server),
          signature = paste0("signature-", server))), names(conns)))
      }
      if (identical(method, "mpcStoreTransportKeysDS")) {
        return(stats::setNames(as.list(rep(TRUE, length(conns))), names(conns)))
      }
      if (identical(method, "glmStandardizeDS")) {
        return(stats::setNames(lapply(names(conns), function(server) list(
          x_means = 0, x_sds = 1,
          y_mean = if (identical(server, "s1")) 0 else NULL,
          y_sd = if (identical(server, "s1")) 1 else NULL)), names(conns)))
      }
    }

    if (is.list(expr) && length(expr) == length(conns)) {
      return(stats::setNames(lapply(names(conns), function(server) {
        server_expr <- expr[[server]]
        expect_identical(as.character(server_expr[[1L]]), "glmStandardizeDS")
        list(
          x_means = stats::setNames(0, paste0("x_", server)),
          x_sds = stats::setNames(1, paste0("x_", server)),
          y_mean = if (identical(server, "s1")) 0 else NULL,
          y_sd = if (identical(server, "s1")) 1 else NULL)
      }), names(conns)))
    }
    stop("Unexpected aggregate expression in test")
  }

  datasources <- list(s1 = structure(list(), class = "fake_conn"),
                      s2 = structure(list(), class = "fake_conn"),
                      s3 = structure(list(), class = "fake_conn"))
  setup <- dsVertClient:::.glm_mpc_setup(
    datasources = datasources,
    server_names = names(datasources),
    server_list = names(datasources),
    non_label_servers = c("s2", "s3"),
    y_server = "s1", y_var = "y",
    x_vars = list(s1 = "x_s1", s2 = "x_s2", s3 = "x_s3"),
    data_name = "D", family = "gaussian",
    session_id = "01234567-89ab-4def-8123-456789abcdef",
    verbose = FALSE, missing = "fail",
    .aggregate = fake_aggregate)

  expect_length(calls, 3L)
  expect_true(all(vapply(calls, function(x)
    identical(x$servers, names(datasources)), logical(1L))))
  expect_true(is.call(calls[[1L]]$expr))
  expect_true(is.call(calls[[2L]]$expr))
  expect_true(is.list(calls[[3L]]$expr))
  expect_named(calls[[3L]]$expr, names(datasources))
  expect_identical(setup$missing_policy, "fail")
  expect_identical(setup$transport_pks,
                   list(s1 = "transport-s1", s2 = "transport-s2",
                        s3 = "transport-s3"))

  standardize_calls <- calls[[3L]]$expr
  expect_identical(standardize_calls$s1$y_var, "y")
  expect_null(standardize_calls$s2$y_var)
  expect_null(standardize_calls$s3$y_var)
  expect_true(all(vapply(standardize_calls, function(x)
    identical(x$missing, "fail"), logical(1L))))
})

test_that("shared transport setup takes two fan-out calls regardless of K", {
  calls <- list()
  fake_aggregate <- function(conns, expr, ...) {
    calls[[length(calls) + 1L]] <<- list(
      servers = names(conns), method = as.character(expr[[1L]]))
    if (identical(as.character(expr[[1L]]),
                  "glmRing63TransportInitDS")) {
      return(stats::setNames(lapply(names(conns), function(server) list(
        transport_pk = paste0("pk-", server),
        identity_pk = paste0("id-", server),
        signature = paste0("sig-", server))), names(conns)))
    }
    stats::setNames(as.list(rep(TRUE, length(conns))), names(conns))
  }

  datasources <- list(s1 = list(), s2 = list(), s3 = list(), s4 = list())
  pks <- dsVertClient:::.dsvert_setup_peer_transport(
    datasources, names(datasources), names(datasources),
    "01234567-89ab-4def-8123-456789abcdef",
    .aggregate = fake_aggregate)

  expect_length(calls, 2L)
  expect_identical(vapply(calls, `[[`, character(1L), "method"),
                   c("glmRing63TransportInitDS", "mpcStoreTransportKeysDS"))
  expect_true(all(vapply(calls, function(x)
    identical(x$servers, names(datasources)), logical(1L))))
  expect_identical(pks,
                   list(s1 = "pk-s1", s2 = "pk-s2", s3 = "pk-s3",
                        s4 = "pk-s4"))
})

test_that("GLM setup neither exposes connector errors nor retries a failed phase", {
  calls <- 0L
  failing_aggregate <- function(conns, expr, ...) {
    calls <<- calls + 1L
    stop("SECRET server-side diagnostic")
  }
  condition <- tryCatch(
    dsVertClient:::.glm_mpc_setup(
      datasources = list(s1 = list(), s2 = list()),
      server_names = c("s1", "s2"), server_list = c("s1", "s2"),
      non_label_servers = "s1", y_server = "s2", y_var = "y",
      x_vars = list(s1 = "x", s2 = character()), data_name = "D",
      family = "gaussian", session_id = "setup-no-retry",
      verbose = FALSE, .aggregate = failing_aggregate),
    error = identity)

  expect_s3_class(condition, "error")
  expect_match(conditionMessage(condition), "partial or invalid site results")
  expect_false(grepl("SECRET", conditionMessage(condition), fixed = TRUE))
  expect_identical(calls, 1L)
})
