# Parse the deliberately narrow formula contract supported by the current
# vertical MPC design. Transformations and interactions require an explicit,
# disclosure-reviewed design-matrix protocol; treating their term labels as
# column names silently changes the requested model.
.dsvert_formula_identifier <- function(value) {
  is.character(value) && length(value) == 1L && !is.na(value) &&
    grepl("^[A-Za-z][A-Za-z0-9_.]{0,127}$", value)
}

.dsvert_formula_reference <- function(expression, qualified, what) {
  if (is.symbol(expression)) {
    column <- as.character(expression)
    if (.dsvert_formula_identifier(column)) {
      return(list(
        reference = column, server = NULL, column = column,
        qualified = FALSE))
    }
  }
  if (isTRUE(qualified) && is.call(expression) &&
      length(expression) == 3L && identical(expression[[1L]], quote(`$`)) &&
      is.symbol(expression[[2L]]) && is.symbol(expression[[3L]])) {
    server <- as.character(expression[[2L]])
    column <- as.character(expression[[3L]])
    if (.dsvert_formula_identifier(server) &&
        .dsvert_formula_identifier(column)) {
      return(list(
        reference = paste0(server, "$", column), server = server,
        column = column, qualified = TRUE))
    }
  }
  stop(
    what, " must be a plain column", if (isTRUE(qualified)) {
      " or an exact server$column reference"
    } else {
      ""
    },
    call. = FALSE)
}

.dsvert_plain_formula <- function(formula, qualified = FALSE) {
  if (!is.logical(qualified) || length(qualified) != 1L || is.na(qualified)) {
    stop("qualified must be one non-missing logical value", call. = FALSE)
  }
  if (is.character(formula)) {
    if (length(formula) != 1L || is.na(formula) || !nzchar(formula) ||
        nchar(formula, type = "bytes") > 10000L) {
      stop("formula must be one non-empty formula string", call. = FALSE)
    }
    formula <- tryCatch(
      stats::as.formula(formula, env = baseenv()),
      error = function(e) stop("Invalid formula syntax", call. = FALSE)
    )
  }
  if (!inherits(formula, "formula") || length(formula) != 3L) {
    stop("formula must be a two-sided formula", call. = FALSE)
  }
  response <- tryCatch(
    .dsvert_formula_reference(
      formula[[2L]], qualified = qualified, what = "The response"),
    error = function(error) stop(
      "The response must be one pre-existing numeric column; transformed ",
      "or multivariate responses are not supported", call. = FALSE))
  if ("." %in% all.names(formula[[3L]], functions = FALSE)) {
    stop(
      "Only additive pre-existing numeric columns are supported in formula; ",
      "the data-dependent '.' expansion is not available across servers.",
      call. = FALSE
    )
  }
  valid_qualified_ast <- function(expression) {
    if (!is.call(expression)) return(TRUE)
    if (identical(expression[[1L]], quote(`$`))) {
      return(!is.null(tryCatch(
        .dsvert_formula_reference(
          expression, qualified = TRUE, what = "A formula term"),
        error = function(error) NULL)))
    }
    arguments <- as.list(expression)[-1L]
    all(vapply(arguments, valid_qualified_ast, logical(1L)))
  }
  if (isTRUE(qualified) && !valid_qualified_ast(formula[[3L]])) {
    stop(
      "Only additive pre-existing numeric columns are supported in formula. ",
      "A qualified column must be written exactly as server$column.",
      call. = FALSE)
  }

  parsed_terms <- tryCatch(
    stats::terms(formula),
    error = function(e) stop("Invalid formula syntax", call. = FALSE)
  )
  labels <- attr(parsed_terms, "term.labels")
  orders <- attr(parsed_terms, "order")

  parse_reference <- function(label) {
    expr <- tryCatch(parse(text = label, keep.source = FALSE),
                     error = function(e) NULL)
    if (length(expr) != 1L) return(NULL)
    tryCatch(
      .dsvert_formula_reference(
        expr[[1L]], qualified = qualified, what = "A predictor"),
      error = function(error) NULL)
  }
  predictors <- lapply(labels, parse_reference)
  if (any(vapply(predictors, is.null, logical(1L))) ||
      any(orders != 1L)) {
    stop(
      "Only additive pre-existing numeric columns are supported in formula. ",
      "Pre-encode factors/transformations locally and precompute within-site ",
      "terms; cross-site interactions require a dedicated MPC protocol.",
      call. = FALSE
    )
  }

  list(
    formula = formula,
    response = response$reference,
    predictors = vapply(predictors, `[[`, character(1L), "reference"),
    intercept = identical(as.integer(attr(parsed_terms, "intercept")), 1L),
    references = c(list(response), predictors)
  )
}

# Resolve formula columns using public schema metadata only. This helper does
# not discover metadata or contact DataSHIELD servers.
.dsvert_resolve_formula_schema <- function(
    formula, schema, response_kinds = "numeric",
    predictor_kinds = "numeric") {
  expected <- c("server", "column", "kind", "role")
  if (!is.data.frame(schema) || !setequal(names(schema), expected) ||
      nrow(schema) < 1L || any(!vapply(schema[expected], is.character,
                                       logical(1L))) ||
      anyNA(schema[expected]) ||
      any(!vapply(schema[c("server", "column")], function(values) {
        all(vapply(values, .dsvert_formula_identifier, logical(1L)))
      }, logical(1L))) ||
      any(!schema$kind %in% c("identifier", "numeric", "categorical")) ||
      any(!schema$role %in% c("data", "id")) ||
      any((schema$kind == "identifier") != (schema$role == "id"))) {
    stop(paste(
      "schema must be non-empty public metadata with character columns",
      "server, column, kind and role"), call. = FALSE)
  }
  normalize_kinds <- function(value, argument) {
    if (!is.character(value) || !length(value) || anyNA(value) ||
        any(!value %in% c("numeric", "categorical")) ||
        anyDuplicated(value)) {
      stop(argument,
           " must contain unique public kinds: numeric or categorical",
           call. = FALSE)
    }
    unname(value)
  }
  response_kinds <- normalize_kinds(response_kinds, "response_kinds")
  predictor_kinds <- normalize_kinds(predictor_kinds, "predictor_kinds")
  schema <- schema[expected]
  physical_key <- paste(schema$server, schema$column, sep = "\r")
  if (anyDuplicated(physical_key)) {
    stop("Public schema has duplicated server/column metadata",
         call. = FALSE)
  }

  parsed <- .dsvert_plain_formula(formula, qualified = TRUE)
  references <- parsed$references
  uses <- c("response", rep("predictor", length(references) - 1L))
  resolved <- lapply(seq_along(references), function(index) {
    reference <- references[[index]]
    matches <- schema$column == reference$column
    if (isTRUE(reference$qualified)) {
      if (!reference$server %in% schema$server) {
        stop("Unknown formula server: ", reference$server, call. = FALSE)
      }
      matches <- matches & schema$server == reference$server
      if (!any(matches)) {
        if (reference$column %in% schema$column) {
          stop(
            "Column ", reference$column, " is not owned by server ",
            reference$server, call. = FALSE)
        }
        stop("Formula variable is missing from public schema: ",
             reference$reference, call. = FALSE)
      }
    } else if (!any(matches)) {
      stop("Formula variable is missing from public schema: ",
           reference$column, call. = FALSE)
    } else if (sum(matches) != 1L) {
      stop(
        "Formula variable has multiple owners; qualify it as server$column: ",
        reference$column, call. = FALSE)
    }
    row <- schema[which(matches), , drop = FALSE]
    if (identical(row$role[[1L]], "id")) {
      stop(
        "Identifier columns cannot be used as a ", uses[[index]], ": ",
        reference$reference, call. = FALSE)
    }
    permitted <- if (identical(uses[[index]], "response")) {
      response_kinds
    } else {
      predictor_kinds
    }
    if (!row$kind[[1L]] %in% permitted) {
      stop(
        "Formula ", uses[[index]], " has incompatible public kind '",
        row$kind[[1L]], "': ", reference$reference,
        "; expected ", paste(permitted, collapse = " or "),
        call. = FALSE)
    }
    data.frame(
      use = uses[[index]], reference = reference$reference,
      server = row$server, column = row$column, kind = row$kind,
      role = row$role, stringsAsFactors = FALSE)
  })
  bindings <- do.call(rbind, resolved)
  predictor_rows <- bindings$use == "predictor"
  canonical_predictors <- bindings$reference[predictor_rows]
  canonical <- paste(
    bindings$reference[[1L]], "~",
    paste(c(if (isTRUE(parsed$intercept)) "1" else "0",
            canonical_predictors), collapse = " + "))
  list(
    formula = parsed$formula,
    response = parsed$response,
    predictors = parsed$predictors,
    intercept = parsed$intercept,
    canonical = canonical,
    bindings = bindings,
    owners = stats::setNames(bindings$server, bindings$reference)
  )
}
