.nb_collect_calls <- function(node, targets) {
  found <- list()
  walk <- function(value) {
    if (!is.call(value)) return(invisible(NULL))
    head <- value[[1L]]
    if (is.symbol(head) && as.character(head) %in% targets) {
      found[[length(found) + 1L]] <<- value
    }
    lapply(as.list(value)[-1L], walk)
    invisible(NULL)
  }
  walk(node)
  found
}

.nb_missing_formal <- function(value) {
  identical(value, quote(expr = ))
}

test_that("NB affine and scale helpers receive exactly their declared contract", {
  helpers <- c(
    ".nb_fullreg_nd_affine_key", ".nb_fullreg_nd_scale_key")
  callers <- list(
    .nb_fullreg_nd_gradient_from_residual,
    .nb_fullreg_nd_beta_score_fisher)
  calls <- unlist(lapply(callers, function(fun) {
    .nb_collect_calls(body(fun), helpers)
  }), recursive = FALSE)
  expect_gte(length(calls), 6L)

  for (invocation in calls) {
    helper_name <- as.character(invocation[[1L]])
    helper <- get(helper_name, envir = asNamespace("dsVertClient"))
    declared <- names(formals(helper))
    supplied <- names(as.list(invocation)[-1L])
    expect_true(all(nzchar(supplied)),
      info = paste(helper_name, "must use explicit named forwarding"))
    expect_true(all(supplied %in% declared),
      info = paste(helper_name, "received undeclared arguments:",
                   paste(setdiff(supplied, declared), collapse = ", ")))
    required <- names(Filter(.nb_missing_formal, formals(helper)))
    expect_true(all(required %in% supplied),
      info = paste(helper_name, "omitted required arguments:",
                   paste(setdiff(required, supplied), collapse = ", ")))
  }
})
