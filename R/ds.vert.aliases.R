#' User-facing ds.vert.* API aliases
#'
#' These wrappers provide a compact, formula-style public surface while keeping
#' the historical CamelCase functions as compatibility backends. Wrappers that
#' have distinct K=2 and K>=3 implementations dispatch from the number of active
#' DataSHIELD connections.
#'
#' @name ds.vert.aliases
NULL

.dsvert_datasources <- function(datasources = NULL) {
  if (is.null(datasources)) datasources <- DSI::datashield.connections_find()
  if (length(datasources) == 0L) {
    stop("No DataSHIELD connections found", call. = FALSE)
  }
  datasources
}

.dsvert_k <- function(datasources = NULL) {
  length(.dsvert_datasources(datasources))
}

.dsvert_set_frontdoor <- function(x, frontdoor, route = NULL, K = NULL) {
  if (is.list(x)) {
    x$frontdoor <- frontdoor
    if (!is.null(route)) x$route <- route
    if (!is.null(K)) x$k_mode <- if (K == 2L) "K=2" else "K>=3"
  }
  x
}

.dsvert_route_result <- function(x, frontdoor, route, datasources = NULL) {
  K <- tryCatch(length(.dsvert_datasources(datasources)),
                error = function(e) NULL)
  .dsvert_set_frontdoor(x, frontdoor, route = route, K = K)
}

#' @rdname ds.vert.aliases
#' @export
ds.vert.align <- function(data_name, id_col, newobj = "D_aligned",
                          datasources = NULL, ...) {
  datasources <- .dsvert_datasources(datasources)
  out <- ds.psiAlign(data_name = data_name, id_col = id_col, newobj = newobj,
                     datasources = datasources, ...)
  .dsvert_set_frontdoor(out, "ds.vert.align", "ds.psiAlign",
                        length(datasources))
}

#' @rdname ds.vert.aliases
#' @export
ds.vert.is_aligned <- function(newobj = "DA", datasources = NULL, ...) {
  datasources <- .dsvert_datasources(datasources)
  out <- ds.isPsiAligned(newobj = newobj, datasources = datasources, ...)
  .dsvert_set_frontdoor(out, "ds.vert.is_aligned", "ds.isPsiAligned",
                        length(datasources))
}

#' @rdname ds.vert.aliases
#' @export
ds.vert.desc <- function(data_name, datasources = NULL, ...) {
  datasources <- .dsvert_datasources(datasources)
  out <- ds.vertDesc(data_name = data_name, datasources = datasources, ...)
  .dsvert_set_frontdoor(out, "ds.vert.desc", "ds.vertDesc",
                        length(datasources))
}

#' @rdname ds.vert.aliases
#' @export
ds.vert.cor <- function(data_name, variables = NULL,
                        datasources = NULL, ...) {
  datasources <- .dsvert_datasources(datasources)
  out <- ds.vertCor(data_name = data_name, variables = variables,
                    datasources = datasources, ...)
  .dsvert_set_frontdoor(out, "ds.vert.cor", "ds.vertCor",
                        length(datasources))
}

#' @rdname ds.vert.aliases
#' @export
ds.vert.pca <- function(data_name = NULL, variables = NULL,
                        datasources = NULL, ...) {
  datasources <- .dsvert_datasources(datasources)
  out <- ds.vertPCA(data_name = data_name, variables = variables,
                    datasources = datasources, ...)
  .dsvert_set_frontdoor(out, "ds.vert.pca", "ds.vertPCA",
                        length(datasources))
}

#' @rdname ds.vert.aliases
#' @export
ds.vert.chisq <- function(data_name, var1, var2,
                          datasources = NULL, ...) {
  datasources <- .dsvert_datasources(datasources)
  out <- ds.vertChisq(data_name = data_name, var1 = var1, var2 = var2,
                      datasources = datasources, ...)
  .dsvert_set_frontdoor(out, "ds.vert.chisq", "ds.vertChisq",
                        length(datasources))
}

#' @rdname ds.vert.aliases
#' @export
ds.vert.fisher <- function(data_name, var1, var2,
                           datasources = NULL, ...) {
  datasources <- .dsvert_datasources(datasources)
  out <- ds.vertFisher(data_name = data_name, var1 = var1, var2 = var2,
                       datasources = datasources, ...)
  .dsvert_set_frontdoor(out, "ds.vert.fisher", "ds.vertFisher",
                        length(datasources))
}

#' @rdname ds.vert.aliases
#' @export
ds.vert.chisq_cross <- function(data, var1, var2,
                                datasources = NULL, ...) {
  datasources <- .dsvert_datasources(datasources)
  out <- ds.vertChisqCross(data = data, var1 = var1, var2 = var2,
                           datasources = datasources, ...)
  .dsvert_set_frontdoor(out, "ds.vert.chisq_cross", "ds.vertChisqCross",
                        length(datasources))
}

#' @rdname ds.vert.aliases
#' @export
ds.vert.glm <- function(formula, data = NULL, datasources = NULL, ...) {
  datasources <- .dsvert_datasources(datasources)
  out <- ds.vertGLM(formula = formula, data = data,
                    datasources = datasources, ...)
  .dsvert_set_frontdoor(out, "ds.vert.glm", "ds.vertGLM",
                        length(datasources))
}

#' @rdname ds.vert.aliases
#' @export
ds.vert.cox <- function(formula, data = NULL, datasources = NULL, ...) {
  datasources <- .dsvert_datasources(datasources)
  out <- ds.vertCox(formula = formula, data = data,
                    datasources = datasources, ...)
  .dsvert_set_frontdoor(out, "ds.vert.cox", "ds.vertCox",
                        length(datasources))
}

#' @rdname ds.vert.aliases
#' @export
ds.vert.coxph <- function(formula, data = NULL, ...) {
  out <- ds.vert.cox(formula = formula, data = data, ...)
  if (is.list(out)) out$frontdoor <- "ds.vert.coxph"
  out
}

#' @rdname ds.vert.aliases
#' @export
ds.vert.nb <- function(formula, data = NULL,
                       method = c("auto", "fast", "mom", "profile",
                                  "accurate"),
                       datasources = NULL, ...) {
  method <- match.arg(method)
  datasources <- .dsvert_datasources(datasources)
  K <- length(datasources)
  route <- switch(method,
    auto = "ds.vertNBMoMTheta",
    fast = "ds.vertNBMoMTheta",
    mom = "ds.vertNBMoMTheta",
    profile = "ds.vertNB",
    accurate = "ds.vertNBFullRegTheta")
  fit <- switch(route,
    ds.vertNBMoMTheta = ds.vertNBMoMTheta(
      formula = formula, data = data, datasources = datasources, ...),
    ds.vertNB = ds.vertNB(
      formula = formula, data = data, datasources = datasources, ...),
    ds.vertNBFullRegTheta = ds.vertNBFullRegTheta(
      formula = formula, data = data, variant = "full_reg_nd",
      datasources = datasources, ...))
  .dsvert_set_frontdoor(fit, "ds.vert.nb", route, K)
}

#' @rdname ds.vert.aliases
#' @export
ds.vert.multinom <- function(formula, data = NULL,
                             datasources = NULL, ...) {
  datasources <- .dsvert_datasources(datasources)
  out <- ds.vertMultinom(formula = formula, data = data,
                         datasources = datasources, ...)
  .dsvert_set_frontdoor(out, "ds.vert.multinom",
                        "ds.vertMultinomJointNewton", length(datasources))
}

#' @rdname ds.vert.aliases
#' @export
ds.vert.ordinal <- function(formula, data = NULL,
                            datasources = NULL, ...) {
  datasources <- .dsvert_datasources(datasources)
  out <- ds.vertOrdinal(formula = formula, data = data,
                        datasources = datasources, ...)
  .dsvert_set_frontdoor(out, "ds.vert.ordinal",
                        "ds.vertOrdinalJointNewton", length(datasources))
}

#' @rdname ds.vert.aliases
#' @export
ds.vert.lmm <- function(formula, data = NULL, cluster_col,
                        max_iter = 30L, inner_iter = 50L,
                        max_outer = 30L, tol = NULL, ring = NULL,
                        verbose = TRUE, datasources = NULL, ...) {
  extra <- list(...)
  if (length(extra)) {
    arg_names <- names(extra)
    arg_names[!nzchar(arg_names)] <- "<unnamed>"
    stop("unused argument(s): ", paste(arg_names, collapse = ", "),
         call. = FALSE)
  }
  datasources <- .dsvert_datasources(datasources)
  K <- length(datasources)
  if (K >= 3L) {
    fit <- ds.vertLMM.k3(formula = formula, data = data,
                         cluster_col = cluster_col,
                         max_outer = max_outer,
                         tol = tol %||% 1e-4,
                         ring = ring %||% "ring127",
                         verbose = verbose,
                         datasources = datasources)
    .dsvert_set_frontdoor(fit, "ds.vert.lmm", "ds.vertLMM.k3", K)
  } else {
    fit <- ds.vertLMM(formula = formula, data = data,
                      cluster_col = cluster_col,
                      max_iter = max_iter,
                      inner_iter = inner_iter,
                      tol = tol %||% 1e-4,
                      ring = ring %||% "ring63",
                      verbose = verbose,
                      datasources = datasources)
    .dsvert_set_frontdoor(fit, "ds.vert.lmm", "ds.vertLMM", K)
  }
}

#' @rdname ds.vert.aliases
#' @export
ds.vert.gee <- function(formula, data = NULL, datasources = NULL, ...) {
  datasources <- .dsvert_datasources(datasources)
  out <- ds.vertGEE(formula = formula, data = data,
                    datasources = datasources, ...)
  .dsvert_set_frontdoor(out, "ds.vert.gee", "ds.vertGEE",
                        length(datasources))
}

#' @rdname ds.vert.aliases
#' @export
ds.vert.glmm <- function(formula, data = NULL, cluster_col,
                         method = c("laplace", "pql"),
                         datasources = NULL, ...) {
  method <- match.arg(method)
  datasources <- .dsvert_datasources(datasources)
  if (identical(method, "laplace")) {
    out <- ds.vertGLMMLaplace(formula = formula, data = data,
                              cluster_col = cluster_col,
                              datasources = datasources, ...)
    .dsvert_set_frontdoor(out, "ds.vert.glmm", "ds.vertGLMMLaplace",
                          length(datasources))
  } else {
    out <- ds.vertGLMM(formula = formula, data = data,
                       cluster_col = cluster_col,
                       datasources = datasources, ...)
    .dsvert_set_frontdoor(out, "ds.vert.glmm", "ds.vertGLMM",
                          length(datasources))
  }
}

#' @rdname ds.vert.aliases
#' @export
ds.vert.glmer <- function(formula, data = NULL, cluster_col,
                          datasources = NULL, ...) {
  datasources <- .dsvert_datasources(datasources)
  out <- ds.vertGLMMLaplace(formula = formula, data = data,
                            cluster_col = cluster_col,
                            datasources = datasources, ...)
  .dsvert_set_frontdoor(out, "ds.vert.glmer", "ds.vertGLMMLaplace",
                        length(datasources))
}

#' @rdname ds.vert.aliases
#' @export
ds.vert.ipw <- function(outcome_formula, propensity_formula, data = NULL,
                        datasources = NULL, ...) {
  datasources <- .dsvert_datasources(datasources)
  out <- ds.vertIPW(outcome_formula = outcome_formula,
                    propensity_formula = propensity_formula,
                    data = data, datasources = datasources, ...)
  .dsvert_set_frontdoor(out, "ds.vert.ipw", "ds.vertIPW",
                        length(datasources))
}

#' @rdname ds.vert.aliases
#' @export
ds.vert.mi <- function(formula, data = NULL, impute_columns = NULL,
                       datasources = NULL, ...) {
  datasources <- .dsvert_datasources(datasources)
  out <- ds.vertMI(formula = formula, data = data,
                   impute_columns = impute_columns,
                   datasources = datasources, ...)
  .dsvert_set_frontdoor(out, "ds.vert.mi", "ds.vertMI",
                        length(datasources))
}

#' @rdname ds.vert.aliases
#' @export
ds.vert.lasso <- function(fit, lambda_1, ...) {
  .dsvert_route_result(ds.vertLASSO(fit = fit, lambda_1 = lambda_1, ...),
                       "ds.vert.lasso", "ds.vertLASSO")
}

#' @rdname ds.vert.aliases
#' @export
ds.vert.lasso_iter <- function(formula, data = NULL, ...) {
  .dsvert_route_result(ds.vertLASSOIter(formula = formula, data = data, ...),
                       "ds.vert.lasso_iter", "ds.vertLASSOIter")
}

#' @rdname ds.vert.aliases
#' @export
ds.vert.lasso_proximal <- function(fit, lambda, ...) {
  .dsvert_route_result(ds.vertLASSOProximal(fit = fit, lambda = lambda, ...),
                       "ds.vert.lasso_proximal", "ds.vertLASSOProximal")
}

#' @rdname ds.vert.aliases
#' @export
ds.vert.lasso_1step <- function(fit, lambda, ...) {
  .dsvert_route_result(ds.vertLASSO1Step(fit = fit, lambda = lambda, ...),
                       "ds.vert.lasso_1step", "ds.vertLASSO1Step")
}

#' @rdname ds.vert.aliases
#' @export
ds.vert.lasso_cv <- function(fit, lambda_grid = NULL, ...) {
  .dsvert_route_result(ds.vertLASSOCV(fit = fit, lambda_grid = lambda_grid,
                                      ...),
                       "ds.vert.lasso_cv", "ds.vertLASSOCV")
}

#' @rdname ds.vert.aliases
#' @export
ds.vert.lr <- function(reduced, full) {
  .dsvert_set_frontdoor(ds.vertLR(reduced = reduced, full = full),
                        "ds.vert.lr", "ds.vertLR")
}

#' @rdname ds.vert.aliases
#' @export
ds.vert.confint <- function(fit, parm = NULL, level = 0.95) {
  .dsvert_set_frontdoor(ds.vertConfint(fit = fit, parm = parm, level = level),
                        "ds.vert.confint", "ds.vertConfint")
}

#' @rdname ds.vert.aliases
#' @export
ds.vert.wald <- function(fit, parm, null = 0) {
  .dsvert_set_frontdoor(ds.vertWald(fit = fit, parm = parm, null = null),
                        "ds.vert.wald", "ds.vertWald")
}

#' @rdname ds.vert.aliases
#' @export
ds.vert.contrast <- function(fit, K, m = NULL) {
  .dsvert_set_frontdoor(ds.vertContrast(fit = fit, K = K, m = m),
                        "ds.vert.contrast", "ds.vertContrast")
}
