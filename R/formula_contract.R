# Parse the deliberately narrow formula contract supported by the current
# vertical MPC design. Transformations and interactions require an explicit,
# disclosure-reviewed design-matrix protocol; treating their term labels as
# column names silently changes the requested model.
.dsvert_plain_formula <- function(formula) {
  if (is.character(formula)) {
    if (length(formula) != 1L || is.na(formula) || !nzchar(formula) ||
        nchar(formula, type = "bytes") > 10000L) {
      stop("formula must be one non-empty formula string", call. = FALSE)
    }
    formula <- tryCatch(
      stats::as.formula(formula),
      error = function(e) stop("Invalid formula syntax", call. = FALSE)
    )
  }
  if (!inherits(formula, "formula") || length(formula) != 3L) {
    stop("formula must be a two-sided formula", call. = FALSE)
  }
  if (!is.symbol(formula[[2L]])) {
    stop("The response must be one pre-existing numeric column; transformed ",
         "or multivariate responses are not supported", call. = FALSE)
  }
  if ("." %in% all.names(formula[[3L]], functions = FALSE)) {
    stop(
      "Only additive pre-existing numeric columns are supported in formula; ",
      "the data-dependent '.' expansion is not available across servers.",
      call. = FALSE
    )
  }

  parsed_terms <- tryCatch(
    stats::terms(formula),
    error = function(e) stop("Invalid formula syntax", call. = FALSE)
  )
  labels <- attr(parsed_terms, "term.labels")
  orders <- attr(parsed_terms, "order")

  is_plain_column <- function(label) {
    expr <- tryCatch(parse(text = label, keep.source = FALSE),
                     error = function(e) NULL)
    length(expr) == 1L && is.symbol(expr[[1L]]) &&
      identical(as.character(expr[[1L]]), label)
  }
  plain <- vapply(labels, is_plain_column, logical(1L))
  if (any(!plain) || any(orders != 1L)) {
    stop(
      "Only additive pre-existing numeric columns are supported in formula. ",
      "Pre-encode factors/transformations locally and precompute within-site ",
      "terms; cross-site interactions require a dedicated MPC protocol.",
      call. = FALSE
    )
  }

  list(
    formula = formula,
    response = as.character(formula[[2L]]),
    predictors = labels,
    intercept = identical(as.integer(attr(parsed_terms, "intercept")), 1L)
  )
}
