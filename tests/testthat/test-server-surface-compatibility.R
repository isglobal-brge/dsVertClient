.collect_character_literals <- function(node) {
  if (missing(node)) return(character())
  values <- if (is.character(node)) node else character()
  if (is.call(node) || is.expression(node) || is.pairlist(node) ||
      is.list(node)) {
    values <- c(values, unlist(
      lapply(as.list(node), .collect_character_literals),
      use.names = FALSE
    ))
  }
  values
}

.collect_remote_constructions <- function(node, result = NULL) {
  if (is.null(result)) {
    result <- list(literal = character(), as_name = character(),
                   dynamic = character())
  }
  if (missing(node)) return(result)

  if (is.call(node)) {
    head <- as.character(node[[1L]])
    head <- head[[length(head)]]
    if (identical(head, "call")) {
      args <- as.list(node)[-1L]
      name_index <- which(names(args) == "name")
      if (length(name_index)) {
        value <- args[[name_index[[1L]]]]
        if (is.character(value) && length(value) == 1L) {
          result$literal <- c(result$literal, value)
        } else if (is.symbol(value)) {
          result$dynamic <- c(result$dynamic, as.character(value))
        } else {
          result$dynamic <- c(
            result$dynamic, paste(deparse(value), collapse = " "))
        }
      }
    }
    if (identical(head, "as.name") && length(node) == 2L &&
        is.character(node[[2L]]) && length(node[[2L]]) == 1L &&
        grepl("^[A-Za-z][A-Za-z0-9.]*DS$", node[[2L]])) {
      result$as_name <- c(result$as_name, node[[2L]])
    }
  }

  if (is.call(node) || is.expression(node) || is.pairlist(node) ||
      is.list(node)) {
    for (child in as.list(node)) {
      result <- .collect_remote_constructions(child, result)
    }
  }
  result
}

.companion_description <- function() {
  configured_root <- Sys.getenv("DSVERT_SERVER_SOURCE", unset = "")
  configured <- if (nzchar(configured_root)) {
    file.path(normalizePath(configured_root, mustWork = FALSE), "DESCRIPTION")
  } else {
    ""
  }
  if (file.exists(configured)) return(configured)

  source_tree <- testthat::test_path("..", "..", "..", "dsVert",
                                     "DESCRIPTION")
  if (file.exists(source_tree)) return(source_tree)

  installed <- system.file("DESCRIPTION", package = "dsVert")
  if (nzchar(installed) && file.exists(installed)) return(installed)
  testthat::skip("the companion dsVert DESCRIPTION is unavailable")
}

.companion_surface_inventory <- function() {
  configured_root <- Sys.getenv("DSVERT_SERVER_SOURCE", unset = "")
  configured <- if (nzchar(configured_root)) {
    file.path(normalizePath(configured_root, mustWork = FALSE), "inst", "docs",
              "remote_surface_classification.json")
  } else {
    ""
  }
  if (file.exists(configured)) return(configured)

  source_tree <- testthat::test_path(
    "..", "..", "..", "dsVert", "inst", "docs",
    "remote_surface_classification.json")
  if (file.exists(source_tree)) return(source_tree)

  installed <- system.file(
    "docs", "remote_surface_classification.json", package = "dsVert")
  if (nzchar(installed) && file.exists(installed)) return(installed)
  testthat::skip("the companion dsVert remote-surface inventory is unavailable")
}

test_that("an explicit companion server source takes precedence", {
  root <- withr::local_tempdir(pattern = "dsvert-companion-source-")
  dir.create(file.path(root, "inst", "docs"), recursive = TRUE)
  writeLines("Package: dsVert", file.path(root, "DESCRIPTION"))
  writeLines("{}", file.path(root, "inst", "docs",
                               "remote_surface_classification.json"))
  withr::local_envvar(c(DSVERT_SERVER_SOURCE = root))
  expected_root <- normalizePath(root, mustWork = TRUE)

  expect_identical(.companion_description(),
                   file.path(expected_root, "DESCRIPTION"))
  expect_identical(
    .companion_surface_inventory(),
    file.path(expected_root, "inst", "docs",
              "remote_surface_classification.json"))
})

.client_source_nodes <- function() {
  r_dir <- testthat::test_path("..", "..", "R")
  sources <- if (dir.exists(r_dir)) {
    list.files(r_dir, pattern = "[.]R$", full.names = TRUE)
  } else {
    character()
  }
  if (length(sources)) return(lapply(sources, parse))

  namespace <- asNamespace("dsVertClient")
  function_names <- lsf.str(envir = namespace, all = TRUE)
  lapply(function_names, function(name) {
    fun <- get(name, envir = namespace, inherits = FALSE)
    list(formals(fun), body(fun))
  })
}

.client_server_method_literals <- function() {
  unique(unlist(lapply(.client_source_nodes(), .collect_character_literals),
                use.names = FALSE))
}

.client_remote_constructions <- function() {
  parts <- lapply(.client_source_nodes(), .collect_remote_constructions)
  list(
    literal = sort(unique(unlist(lapply(parts, `[[`, "literal"),
                                 use.names = FALSE))),
    as_name = sort(unique(unlist(lapply(parts, `[[`, "as_name"),
                                 use.names = FALSE))),
    dynamic = sort(unique(unlist(lapply(parts, `[[`, "dynamic"),
                                 use.names = FALSE)))
  )
}

test_that("client DSI expressions are registered or locally quarantined", {
  fields <- read.dcf(.companion_description())
  method_fields <- intersect(c("AggregateMethods", "AssignMethods"),
                             colnames(fields))
  registration <- unique(unlist(lapply(method_fields, function(field) {
    trimws(strsplit(fields[1L, field], ",", fixed = TRUE)[[1L]])
  }), use.names = FALSE))
  registered <- sub("^[^=]+=", "", registration)
  registered_ds <- sort(registered[
    grepl("^[A-Za-z][A-Za-z0-9.]*DS$", registered)])

  inventory <- jsonlite::read_json(
    .companion_surface_inventory(), simplifyVector = FALSE)
  expect_identical(
    inventory$schema, "dsvert-remote-surface-classification-v3")
  dynamic_allowed <- unlist(
    inventory$client_ast_resolution$bounded_dynamic_call$allowed_endpoints,
    use.names = FALSE)
  constructions <- .client_remote_constructions()

  expect_identical(
    constructions$dynamic,
    inventory$client_ast_resolution$bounded_dynamic_call$symbol)
  expect_setequal(
    constructions$as_name,
    unlist(inventory$client_ast_resolution$as_call_as_name,
           use.names = FALSE))
  resolved <- sort(unique(c(
    constructions$literal, constructions$as_name, dynamic_allowed)))
  retired <- sort(unique(unlist(
    inventory$retired_registered_surface, use.names = FALSE)))
  locally_quarantined <- sort(unique(unlist(
    inventory$client_ast_resolution$locally_quarantined_unregistered_endpoints,
    use.names = FALSE)))

  expect_false("dsvertHistogramDS" %in% constructions$literal)
  expect_false("dsvertHistogramDS" %in% registered_ds)
  expect_false("dsvertContingencyDS" %in% constructions$literal)
  expect_false("dsvertContingencyDS" %in% registered_ds)

  unregistered <- setdiff(resolved, registered_ds)
  expect_true(all(unregistered %in%
                    sort(unique(c(retired, locally_quarantined)))))
  expect_setequal(setdiff(retired, unregistered), "dsvertImputeColumnDS")
  expect_setequal(setdiff(registered_ds, resolved), character())
  expect_identical(
    length(constructions$literal),
    as.integer(inventory$counts$direct_literal_call_endpoints))
  expect_identical(
    length(constructions$as_name),
    as.integer(inventory$counts$as_call_as_name_endpoints))
  expect_identical(
    length(dynamic_allowed),
    as.integer(inventory$counts$bounded_dynamic_call_endpoints))

  classified <- unlist(
    inventory$registered_endpoint_classes$production_safe_purpose_bound,
    use.names = FALSE)
  expect_setequal(classified, registered_ds)
  expect_identical(anyDuplicated(classified), 0L)
  expect_length(retired, 105L)
  expect_setequal(
    locally_quarantined,
    c("dsvertFormalFinalizerHandoffSourceDS",
      "dsvertFormalGLMControlSourceDS", "dsvertJointDPCapsuleStatusDS"))
  expect_identical(anyDuplicated(retired), 0L)
  expect_false(any(retired %in% registered_ds))
  expect_length(
    inventory$client_ast_resolution$registered_without_product_call_builder,
    0L)
  expect_length(
    inventory$client_ast_resolution$reachable_unregistered_endpoints, 0L)
})

test_that("the retired consumer audit names only preserved public frontdoors", {
  inventory <- jsonlite::read_json(
    .companion_surface_inventory(), simplifyVector = FALSE)
  routes <- inventory$retired_surface_consumer_audit$
    quarantined_frontdoor_routes
  documented <- unlist(routes, use.names = FALSE)
  expect_identical(anyDuplicated(documented), 0L)
  expect_true(all(documented %in% getNamespaceExports("dsVertClient")))
  expect_setequal(
    names(routes),
    c("legacy_glm", "cox", "negative_binomial", "multinomial", "ordinal",
      "lmm", "gee", "ipw", "mi"))
})

test_that("the only dynamic call name is closed before DSI expression construction", {
  inventory <- jsonlite::read_json(
    .companion_surface_inventory(), simplifyVector = FALSE)
  allowed <- unlist(
    inventory$client_ast_resolution$bounded_dynamic_call$allowed_endpoints,
    use.names = FALSE)
  store_blob <- get(
    ".dsvert_store_blob", envir = asNamespace("dsVertClient"),
    inherits = FALSE)

  expect_identical(eval(formals(store_blob)$store_function), allowed[[1L]])
  body_ds_literals <- unique(.collect_character_literals(body(store_blob)))
  body_ds_literals <- body_ds_literals[
    grepl("^[A-Za-z][A-Za-z0-9.]*DS$", body_ds_literals)]
  expect_setequal(body_ds_literals, allowed)

  args <- list(
    blob = "opaque", key = "slot", conn = list(),
    session_id = "session", .aggregate = function(...) NULL)
  expect_error(
    do.call(store_blob, c(args, list(store_function = "arbitraryDS"))),
    "must name a supported idempotent blob endpoint", fixed = TRUE)
  for (endpoint in allowed) {
    expect_error(
      do.call(store_blob, c(args, list(store_function = endpoint))),
      "conn must contain exactly one named DataSHIELD target", fixed = TRUE)
  }
})

test_that("retired redundant triple-consume endpoints have no client caller", {
  literals <- .client_server_method_literals()

  expect_false(any(c(
    "k2BeaverVecmulConsumeTripleDS",
    "k2StoreGradTripleDS"
  ) %in% literals))
})
