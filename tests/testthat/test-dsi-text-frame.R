test_that("client DSI framing matches the server golden wire contract", {
  values <- c(
    "", "abc_DEF-09",
    paste0(strrep("A", 20), "\"", strrep("B", 20)),
    "{\"a\":\"b\\\\c\",\"n\":1}", "é")
  expected <- c(
    "DSV1R_", "DSV1R_abc_DEF-09",
    paste0("DSV1L_20-", strrep("A", 20), "2220-", strrep("B", 20)),
    "DSV1B_eyJhIjoiYlxcYyIsIm4iOjF9", "DSV1B_w6k")
  encoded <- vapply(
    values, .dsvert_dsi_text_encode, character(1L), USE.NAMES = FALSE)

  expect_identical(encoded, expected)
  expect_identical(vapply(
    encoded, .dsvert_dsi_text_decode, character(1L), USE.NAMES = FALSE),
    values)
})

test_that("DSV1 package generations fail closed in both directions and mixed", {
  raw <- "{\"version\":\"old-raw-client\"}"
  framed <- .dsvert_dsi_text_encode(raw)

  old_to_new <- tryCatch(
    .dsvert_dsi_text_decode(raw, "old client payload"),
    error = identity)
  new_to_old <- tryCatch(
    jsonlite::fromJSON(framed, simplifyVector = FALSE),
    error = identity)

  expect_s3_class(old_to_new, "error")
  expect_s3_class(new_to_old, "error")
  expect_identical(.dsvert_dsi_text_decode(framed), raw)

  first_data_free_phase <- c(
    current = TRUE,
    legacy = !inherits(new_to_old, "error"))
  protected_access_started <- FALSE
  if (all(first_data_free_phase)) protected_access_started <- TRUE
  expect_false(all(first_data_free_phase))
  expect_false(protected_access_started)
})

test_that("promoted DSI JSON map is exactly 49 methods and 97 formals", {
  expect_length(.DSVERT_DSI_TEXT_REMOTE_FORMALS, 49L)
  expect_identical(sum(lengths(.DSVERT_DSI_TEXT_REMOTE_FORMALS)), 97L)
  expect_false("exactGCExchangeDS" %in%
                 names(.DSVERT_DSI_TEXT_REMOTE_FORMALS))
  expect_false("psiPaddedFilterDS" %in%
                 names(.DSVERT_DSI_TEXT_REMOTE_FORMALS))

  skip_if_not_installed("dsVert")
  server <- asNamespace("dsVert")
  if (!all(vapply(names(.DSVERT_DSI_TEXT_REMOTE_FORMALS), exists,
                  logical(1L), envir = server, inherits = FALSE))) {
    skip("requires the paired current dsVert source namespace")
  }
  for (method in names(.DSVERT_DSI_TEXT_REMOTE_FORMALS)) {
    expect_true(exists(method, envir = server, inherits = FALSE),
                info = method)
    expect_setequal(
      intersect(names(formals(get(method, envir = server))),
                .DSVERT_DSI_TEXT_REMOTE_FORMALS[[method]]),
      .DSVERT_DSI_TEXT_REMOTE_FORMALS[[method]])
  }
})

test_that("rewriter frames only the explicit promoted arguments once", {
  for (method in names(.DSVERT_DSI_TEXT_REMOTE_FORMALS)) {
    payload <- as.character(jsonlite::toJSON(
      list(value = "quoted\\path"), auto_unbox = TRUE))
    arguments <- stats::setNames(
      rep(list(payload),
          length(.DSVERT_DSI_TEXT_REMOTE_FORMALS[[method]])),
      .DSVERT_DSI_TEXT_REMOTE_FORMALS[[method]])
    expression <- as.call(c(list(as.name(method)), arguments))
    framed <- .dsvert_dsi_text_frame_expressions(expression)
    for (argument in names(arguments)) {
      expect_true(.dsvert_dsi_text_is_framed(framed[[argument]]),
                  info = paste(method, argument))
      expect_identical(
        .dsvert_dsi_text_decode(framed[[argument]]), arguments[[argument]])
    }
    expect_error(
      .dsvert_dsi_text_frame_expressions(framed), "already framed")
  }

  exact <- call(
    "exactGCExchangeDS", session_id = "session", operation_id = "operation",
    delivery_payload = "safe_Base64url-")
  assign <- call("psiPaddedFilterDS", data = as.name("data"))
  expect_identical(.dsvert_dsi_text_frame_expressions(exact), exact)
  expect_identical(.dsvert_dsi_text_frame_expressions(assign), assign)
})

test_that("aggregate and cleanup boundaries frame before validation and send", {
  raw <- call(
    "psiPaddedRelayExchangeDS",
    request_json = "{\"routes\":[{\"peer\":\"p1\"}]}",
    session_id = "session")
  seen <- list()
  aggregate <- function(conns, expr, ...) {
    seen[[length(seen) + 1L]] <<- expr
    stats::setNames(list(TRUE), names(conns))
  }
  result <- .dsvert_aggregate_strict(
    list(site = list()), raw, result_contract = "logical_true",
    .aggregate = aggregate)
  expect_identical(result, list(site = TRUE))
  expect_true(.dsvert_dsi_text_is_framed(seen[[1L]]$request_json))
  expect_identical(
    .dsvert_dsi_text_decode(seen[[1L]]$request_json), raw$request_json)

  expect_error(
    .dsvert_transport_aggregate(
      aggregate, list(site = list()), raw, async = FALSE),
    "Unframed promoted DSI call")

  cleanup <- call(
    "exactGCCleanupDS", session_id = "session",
    cleanup_capability_json = as.character(jsonlite::toJSON(
      list(capability = "signed\\value"), auto_unbox = TRUE)))
  expect_true(.dsvert_cleanup_best_effort(
    list(site = list()), cleanup, .aggregate = aggregate))
  cleanup_seen <- seen[[2L]]$site
  expect_true(.dsvert_dsi_text_is_framed(
    cleanup_seen$cleanup_capability_json))
  expect_identical(
    .dsvert_dsi_text_decode(cleanup_seen$cleanup_capability_json),
    cleanup$cleanup_capability_json)
})

test_that("source W=8 measures framed expressions and shrinks adversarial B", {
  expression_bytes <- function(value) nchar(paste(
    deparse(value, width.cutoff = 500L), collapse = "\n"), type = "bytes")
  batch <- function(value) {
    json <- .dsvert_joint_dp_client_json(value)
    list(
      p1 = list(value = value, json = json),
      p2 = list(value = value, json = json))
  }

  typical <- rep(list(batch(list(
    ciphertext = strrep("A_", 5000L), signature = strrep("B_", 50L)))),
    8L)
  typical_prefix <- .dsvert_dp_capsule_source_accept_prefix(
    c("p1", "p2"), typical,
    capacity_bytes = c(p1 = 100000, p2 = 100000))
  typical_wire <- .dsvert_dsi_text_frame_expressions(
    typical_prefix$expressions)
  expect_identical(typical_prefix$count, 8L)
  expect_true(startsWith(typical_wire$p1$envelope_json, "DSV1L_"))
  expect_lte(expression_bytes(typical_wire$p1), 100000)

  adversarial <- rep(list(batch(list(x = as.list(seq_len(1000L))))), 8L)
  raw_w8 <- .dsvert_dp_capsule_source_accept_expressions(
    c("p1", "p2"), adversarial)
  framed_w8 <- .dsvert_dsi_text_frame_expressions(raw_w8)
  expect_true(startsWith(framed_w8$p1$envelope_json, "DSV1B_"))
  capacity <- expression_bytes(raw_w8$p1)
  expect_gt(expression_bytes(framed_w8$p1), capacity)
  adversarial_prefix <- .dsvert_dp_capsule_source_accept_prefix(
    c("p1", "p2"), adversarial,
    capacity_bytes = c(p1 = capacity, p2 = capacity))
  adversarial_wire <- .dsvert_dsi_text_frame_expressions(
    adversarial_prefix$expressions)
  expect_lt(adversarial_prefix$count, 8L)
  expect_lte(expression_bytes(adversarial_wire$p1), capacity)
})

test_that("real DSLite carries framed JSON for all seven promoted families", {
  skip_if_not_installed("DSLite")
  server <- DSLite::newDSLiteServer(tables = list(t = data.frame(x = 1)))
  register <- function(name, method) {
    environment(method) <- asNamespace("dsVertClient")
    server$aggregateMethod(name, method)
  }
  register(
    "dsvertDPCapsuleManifestSignDS",
    function(schema_json, workload_contract_json) list(
      schema = getFromNamespace(
        ".dsvert_dsi_text_decode", "dsVertClient")(schema_json),
      workload = getFromNamespace(
        ".dsvert_dsi_text_decode", "dsVertClient")(workload_contract_json)))
  register(
    "dsvertJointDPVectorAllocationCommitDS",
    function(first_prepare_json, second_prepare_json) list(
      first = getFromNamespace(
        ".dsvert_dsi_text_decode", "dsVertClient")(first_prepare_json),
      second = getFromNamespace(
        ".dsvert_dsi_text_decode", "dsVertClient")(second_prepare_json)))
  register(
    "dsvertJointDPVectorResultDS",
    function(manifest_json, first_prepare_json, second_prepare_json,
             session_id = NULL) list(
      manifest = getFromNamespace(
        ".dsvert_dsi_text_decode", "dsVertClient")(manifest_json),
      first = getFromNamespace(
        ".dsvert_dsi_text_decode", "dsVertClient")(first_prepare_json),
      second = getFromNamespace(
        ".dsvert_dsi_text_decode", "dsVertClient")(second_prepare_json)))
  register(
    "dsvertDPCapsuleSourceAcceptDS",
    function(envelope_json, transport_contract = "scalar") list(
      envelope = getFromNamespace(
        ".dsvert_dsi_text_decode", "dsVertClient")(envelope_json)))
  register(
    "dsvertDPCategoricalCrossFinalizeDS",
    function(manifest_json, analysis_id, session_id) list(
      manifest = getFromNamespace(
        ".dsvert_dsi_text_decode", "dsVertClient")(manifest_json)))
  register(
    "psiPaddedRelayExchangeDS",
    function(request_json, session_id, outbound_operation_id = "",
             terminal_receipt_b64url = "") list(
      request = getFromNamespace(
        ".dsvert_dsi_text_decode", "dsVertClient")(request_json)))
  register(
    "exactGCCleanupDS",
    function(session_id, cleanup_capability_json) list(
      capability = getFromNamespace(
        ".dsvert_dsi_text_decode", "dsVertClient")(
          cleanup_capability_json)))

  object <- paste0("dsvert_dsi_frame_", Sys.getpid())
  assign(object, server, envir = .GlobalEnv)
  on.exit(rm(list = object, envir = .GlobalEnv), add = TRUE)
  builder <- DSI::newDSLoginBuilder()
  builder$append(
    server = "site", url = object, table = "t", driver = "DSLiteDriver")
  conns <- DSI::datashield.login(builder$build(), assign = FALSE)
  on.exit(DSI::datashield.logout(conns), add = TRUE)

  one <- as.character(jsonlite::toJSON(
    list(a = "quoted path", n = 1), auto_unbox = TRUE))
  two <- "{\"values\":[1,2,3],\"ok\":true}"
  expressions <- list(
    manifest = call(
      "dsvertDPCapsuleManifestSignDS", schema_json = one,
      workload_contract_json = two),
    allocation = call(
      "dsvertJointDPVectorAllocationCommitDS",
      first_prepare_json = one, second_prepare_json = two),
    vector = call(
      "dsvertJointDPVectorResultDS", manifest_json = one,
      first_prepare_json = two, second_prepare_json = one,
      session_id = "session"),
    source = call(
      "dsvertDPCapsuleSourceAcceptDS", envelope_json = one,
      transport_contract = "scalar"),
    cross = call(
      "dsvertDPCategoricalCrossFinalizeDS", manifest_json = one,
      analysis_id = "analysis", session_id = "session"),
    psi = call(
      "psiPaddedRelayExchangeDS", request_json = one,
      session_id = "session"),
    cleanup = call(
      "exactGCCleanupDS", session_id = "session",
      cleanup_capability_json = one))

  expect_error(
    DSI::datashield.aggregate(conns, expressions$manifest),
    "DataSHIELD errors")
  for (family in names(expressions)) {
    result <- .dsvert_aggregate_strict(
      conns, expressions[[family]], operation = family,
      .aggregate = DSI::datashield.aggregate)[[1L]]
    expect_true(is.list(result), info = family)
    expect_true(all(unlist(result, use.names = FALSE) %in% c(one, two)),
                info = family)
  }

  long_json <- as.character(jsonlite::toJSON(list(
    payload = strrep("A_", 4096L)), auto_unbox = TRUE))
  long_expression <- call(
    "dsvertDPCapsuleManifestSignDS", schema_json = long_json,
    workload_contract_json = two)
  long_wire <- .dsvert_dsi_text_frame_expressions(long_expression)
  expect_true(startsWith(long_wire$schema_json, "DSV1L_"))
  long_result <- .dsvert_aggregate_strict(
    conns, long_expression, operation = "manifest L framing",
    .aggregate = DSI::datashield.aggregate)[[1L]]
  expect_identical(long_result$schema, long_json)
})
