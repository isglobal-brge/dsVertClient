#' Public ds.vert.* API aliases
#'
#' These wrappers preserve the compact formula-style API names. Availability
#' is defined by \code{ds.vertMethodStatus()} and the selected signed capsule,
#' not by the presence of an alias. Aliases for quarantined families raise a
#' typed \code{dsvert_route_unavailable} condition before any DSI call.
#' \code{ds.vert.glm()} can delegate an explicit Gaussian
#' \code{dp_analysis_id} to \code{ds.vertDPGaussian()}; its no-id legacy route
#' is unavailable. A binomial/Poisson \code{analysis_id} opens one signed
#' finite likelihood-grid Synopsis; \code{formal_analysis_id} reads a
#' completed custodian-signed public certificate; a
#' \code{fresh_formal_analysis_id} runs only the corresponding
#' custodian-configured durable formal analysis. Neither runs the retired
#' iterative route.
#' \code{ds.vert.multinom(..., frequency = x)} is an intercept-only log-odds
#' postprocess from one validated sticky Frequency release. With a signed
#' \code{analysis_id}, it selects additive softmax coefficients from one finite
#' sticky-DP likelihood grid; covariance and sampling inference remain unavailable.
#' \code{ds.vert.ordinal(..., frequency = x)} likewise returns only
#' intercept-only cumulative-logit thresholds from the same kind of release;
#' its complete ordered category domain must be explicit. It can likewise read
#' the artifact by \code{server}. With a signed \code{analysis_id}, it selects
#' additive cumulative-logit coefficients and thresholds from one finite
#' sticky-DP likelihood grid; covariance and sampling inference remain unavailable.
#' \code{ds.vert.nb(..., frequency = x)} returns an intercept-only NB2
#' method-of-moments fit for a bounded non-negative integer frequency domain.
#' With a signed \code{analysis_id}, it also selects additive coefficients and
#' dispersion from one finite sticky-DP likelihood grid; it has no covariance
#' or sampling inference and cannot create another release.
#' \code{ds.vert.lmm(..., analysis_id = x)} reads the matching signed
#' random-intercept, finite ML/REML, or finite one-or-more-random-slope grid.
#' It has no arbitrary random-effect formula, continuous optimisation, standard
#' errors or classical inference.
#' \code{ds.vert.glmm(..., analysis_id = x)} reads either that binary
#' \code{outcome ~ 1} population-average moment artifact or a custodian-signed
#' finite random-intercept likelihood grid for additive bare covariates; it is
#' not the retired PQL route and has no unconstrained optimisation or inference.
#' No alias re-enables a retired remote endpoint or weakens the signed-artifact
#' and custodian-owned policy gates of an available backend.
#' For \code{ds.vert.pca()}, an authenticated \code{cor_result} can be supplied
#' through \code{...}; it is local post-processing of that signed release and
#' does not require a DataSHIELD connection or create a new DP release.
#' \code{ds.vert.align()} returns a credential-free
#' \code{ds.vertFederation}. Pass that object as \code{data_name} (or
#' \code{data} for Gaussian GLM) to reuse the aligned symbol. Each consumer
#' re-attests the same sites and PSI contract; it does not rerun PSI. The
#' federation also caches the custodian-published, policy-only column kind and
#' role catalogue. Formula frontdoors use it to resolve unique names and to
#' require explicit \code{server$column} qualification for collisions.
#'
#' @param data_name,data,formula,datasources Aligned data-frame symbol, reusable
#'   \code{ds.vertFederation}, model \code{formula}, and DataSHIELD connections.
#'   Model formulas may qualify ambiguous columns as
#'   \code{site_name$column}; the expression is parsed, never evaluated.
#' @param server Required source-owner when an intercept-only NB2, multinomial,
#'   or ordinal alias resolves its canonical signed Frequency artifact itself.
#' @param id_col,newobj Record-identifier column and output symbol for alignment.
#'   \code{data_name} and \code{id_col} may each be one string broadcast to all
#'   sites or a complete named per-site character/list map.
#' @param variables,var1,var2 Column selections for descriptive / bivariate routes.
#' @param cluster_col Grouping column for the mixed-model routes.
#' @param analysis_id,dp_analysis_id Custodian-configured signed random-intercept,
#'   binomial/Poisson, Gaussian, NB2, multinomial, ordinal, or durable Cox
#'   analysis id. \code{analysis_id} is required by \code{ds.vert.lmm()} and
#'   \code{ds.vert.glmm()}, enables finite-grid NB2, multinomial and ordinal
#'   aliases, and selects the profile Cox analysis; \code{dp_analysis_id}
#'   enables only independent Gaussian GEE post-processing.
#' @param precision,method,ring,verbose Binomial-sigmoid precision preset,
#'   estimator/route selector, fixed-point ring, and progress flag. For
#'   code{ds.vert.glmm()}, code{method = "auto"} selects the signed
#'   intercept-only moment or additive finite-grid artifact from the formula;
#'   the explicit values must agree with that formula.
#' @param formal_analysis_id,fresh_formal_analysis_id Custodian-owned selector
#'   for an already completed formal Cox profile, binomial/Poisson GLM
#'   certificate, or discrete-time public certificate. The fresh GLM selector
#'   is accepted by \code{ds.vert.glm()} and \code{ds.vert.gee()} only to run
#'   its configured durable analysis. For Cox, \code{analysis_id} is the
#'   standard profile selector and \code{fresh_formal_analysis_id} remains an
#'   alias. The GEE adapter remains limited to the GLM independence-working
#'   point-estimate scope, and the discrete selector remains a distinct
#'   fixed-grid pooled-logistic estimand.
#' @param max_iter,inner_iter,max_outer,tol Iteration caps and convergence
#'   tolerance for the iterative fits.
#' @param outcome_formula,propensity_formula Outcome and propensity models (IPW).
#' @param impute_columns,m Columns to impute and number of imputations (MI).
#' @param lambda,lambda_1,lambda_grid Penalty or penalty grid for the LASSO routes.
#' @param fit,reduced,full,parm,level,type,null,K Inference-helper inputs
#'   (fitted object, nested models, parameter, confidence level, interval type,
#'   null value, class count). For code{ds.vert.confint()},
#'   code{type = "mechanism"} requests the signed Gaussian DP-mechanism
#'   region; it is not a sampling confidence interval.
#' @param ... Further arguments forwarded to the backend.
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
    x$k_mode <- if (is.null(K) || length(K) != 1L || is.na(K) ||
                      !is.numeric(K) || !is.finite(K) || K != floor(K) ||
                      K < 1L) {
      "unknown"
    } else if (K == 1L) {
      "K=1"
    } else if (K == 2L) {
      "K=2"
    } else {
      "K>=3"
    }
  }
  x
}

.dsvert_route_result <- function(x, frontdoor, route, datasources = NULL) {
  K <- tryCatch(length(.dsvert_datasources(datasources)),
                error = function(e) NULL)
  .dsvert_set_frontdoor(x, frontdoor, route = route, K = K)
}

.dsvert_arg_family <- function(args, default = "gaussian") {
  family <- args$family %||% default
  if (is.character(family) && length(family) >= 1L) return(family[[1L]])
  if (inherits(family, "family") && is.character(family$family)) {
    return(family$family[[1L]])
  }
  default
}

.dsvert_precision_intervals <- function() {
  as.integer(getOption("dsvert.frontdoor_binomial_sigmoid_intervals", 150L))
}

.dsvert_apply_binomial_precision <- function(args,
                                             precision = c("auto", "high",
                                                           "fast"),
                                             family = NULL,
                                             force_binomial = FALSE) {
  precision <- match.arg(precision)
  if (identical(precision, "fast")) return(args)
  if (is.null(family)) family <- .dsvert_arg_family(args)
  if ((isTRUE(force_binomial) || identical(family, "binomial")) &&
      is.null(args$binomial_sigmoid_intervals)) {
    args$binomial_sigmoid_intervals <- .dsvert_precision_intervals()
  }
  args
}

.dsvert_add_policy <- function(x, method = NULL, precision = NULL) {
  if (is.list(x)) {
    if (!is.null(method)) x$method_frontdoor <- method
    if (!is.null(precision)) x$precision_frontdoor <- precision
  }
  x
}

#' @rdname ds.vert.aliases
#' @export
ds.vert.align <- function(data_name, id_col, newobj = "D_aligned",
                          datasources = NULL, ...) {
  datasources <- .dsvert_datasources(datasources)
  data_name <- .dsvert_site_character(data_name, datasources, "data_name")
  id_col <- .dsvert_site_character(id_col, datasources, "id_col")
  attestation <- ds.psiAlign(
    data_name = data_name, id_col = id_col, newobj = newobj,
    datasources = datasources, ...)
  attestation <- .dsvert_validate_psi_padded_attestation(attestation)
  public_schema <- .dsvert_federation_public_schema(
    symbol = newobj,
    datasources = datasources,
    id_columns = id_col,
    attestation = attestation)
  .dsvert_new_federation(
    symbol = newobj,
    sites = names(datasources),
    source_symbols = data_name,
    id_columns = id_col,
    attestation = attestation,
    public_schema = public_schema)
}

#' @rdname ds.vert.aliases
#' @export
ds.vert.is_aligned <- function(newobj = "DA", datasources = NULL, ...) {
  datasources <- .dsvert_datasources(datasources)
  if (inherits(newobj, "ds.vertFederation")) {
    status <- tryCatch(
      .dsvert_federation_status(newobj, datasources),
      error = function(e) NULL)
    out <- list(aligned = !is.null(status), n_common = NA_integer_)
    return(.dsvert_set_frontdoor(
      out, "ds.vert.is_aligned", "ds.isPsiAligned", length(datasources)))
  }
  out <- ds.isPsiAligned(newobj = newobj, datasources = datasources, ...)
  .dsvert_set_frontdoor(out, "ds.vert.is_aligned", "ds.isPsiAligned",
                        length(datasources))
}

#' @rdname ds.vert.aliases
#' @export
ds.vert.desc <- function(data_name, datasources = NULL, ...) {
  datasources <- .dsvert_datasources(datasources)
  resolved <- .dsvert_resolve_federation(data_name, datasources)
  out <- ds.vertDesc(
    data_name = resolved$value, datasources = resolved$datasources, ...)
  if (inherits(out, "data.frame")) return(out)
  .dsvert_set_frontdoor(out, "ds.vert.desc", "ds.vertDesc",
                        length(datasources))
}

#' @rdname ds.vert.aliases
#' @export
ds.vert.cor <- function(data_name, variables = NULL,
                        datasources = NULL, ...) {
  datasources <- .dsvert_datasources(datasources)
  resolved <- .dsvert_resolve_federation(data_name, datasources)
  out <- ds.vertCor(
    data_name = resolved$value, variables = variables,
    datasources = resolved$datasources, ...)
  .dsvert_set_frontdoor(out, "ds.vert.cor", "ds.vertCor",
                        length(datasources))
}

#' @rdname ds.vert.aliases
#' @export
ds.vert.pca <- function(data_name = NULL, variables = NULL,
                        datasources = NULL, ...) {
  arguments <- list(...)
  if (!is.null(arguments$cor_result)) {
    out <- ds.vertPCA(
      data_name = data_name, variables = variables,
      datasources = datasources, ...)
    return(.dsvert_set_frontdoor(out, "ds.vert.pca", "ds.vertPCA"))
  }
  datasources <- .dsvert_datasources(datasources)
  resolved <- .dsvert_resolve_federation(data_name, datasources)
  out <- ds.vertPCA(
    data_name = resolved$value, variables = variables,
    datasources = resolved$datasources, ...)
  .dsvert_set_frontdoor(out, "ds.vert.pca", "ds.vertPCA",
                        length(datasources))
}

#' @rdname ds.vert.aliases
#' @export
ds.vert.chisq <- function(data_name, var1 = NULL, var2 = NULL,
                          datasources = NULL, ...) {
  if (inherits(data_name, "ds.vertFederation")) {
    datasources <- .dsvert_datasources(datasources)
    resolved <- .dsvert_resolve_federation(data_name, datasources)
    data_name <- resolved$value
    datasources <- resolved$datasources
  }
  existing_release <- inherits(data_name, "ds.vertDPContingency")
  if (!existing_release) datasources <- .dsvert_datasources(datasources)
  out <- ds.vertChisq(data_name = data_name, var1 = var1, var2 = var2,
                      datasources = datasources, ...)
  .dsvert_set_frontdoor(out, "ds.vert.chisq", "ds.vertChisq",
                        if (existing_release) NULL else length(datasources))
}

#' @rdname ds.vert.aliases
#' @export
ds.vert.fisher <- function(data_name, var1 = NULL, var2 = NULL,
                           datasources = NULL, ...) {
  if (inherits(data_name, "ds.vertFederation")) {
    datasources <- .dsvert_datasources(datasources)
    resolved <- .dsvert_resolve_federation(data_name, datasources)
    data_name <- resolved$value
    datasources <- resolved$datasources
  }
  existing_release <- inherits(data_name, "ds.vertDPContingency")
  if (!existing_release) datasources <- .dsvert_datasources(datasources)
  out <- ds.vertFisher(data_name = data_name, var1 = var1, var2 = var2,
                       datasources = datasources, ...)
  .dsvert_set_frontdoor(out, "ds.vert.fisher", "ds.vertFisher",
                        if (existing_release) NULL else length(datasources))
}

#' @rdname ds.vert.aliases
#' @export
ds.vert.chisq_cross <- function(data, var1 = NULL, var2 = NULL,
                                datasources = NULL, ...) {
  if (inherits(data, "ds.vertFederation")) {
    datasources <- .dsvert_datasources(datasources)
    resolved <- .dsvert_resolve_federation(data, datasources)
    data <- resolved$value
    datasources <- resolved$datasources
  }
  existing_release <- inherits(data, "ds.vertDPContingency")
  if (!existing_release) datasources <- .dsvert_datasources(datasources)
  out <- ds.vertChisqCross(data = data, var1 = var1, var2 = var2,
                           datasources = datasources, ...)
  .dsvert_set_frontdoor(out, "ds.vert.chisq_cross", "ds.vertChisqCross",
                        if (existing_release) NULL else length(datasources))
}

#' @rdname ds.vert.aliases
#' @export
ds.vert.glm <- function(formula, data = NULL,
                        precision = c("auto", "high", "fast"),
                        datasources = NULL, ...) {
  dots_call <- match.call(expand.dots = FALSE)$...
  dots_names <- if (is.null(dots_call)) character() else names(dots_call)
  grid <- "analysis_id" %in% dots_names
  formal <- "formal_analysis_id" %in% dots_names
  fresh <- "fresh_formal_analysis_id" %in% dots_names
  gaussian_dp <- "dp_analysis_id" %in% dots_names
  if (!grid && !formal && !fresh && !gaussian_dp) {
    .dsvert_block_retired_remote_route("legacy_glm")
  }
  precision <- match.arg(precision)
  if (grid || formal || fresh) {
    if (!missing(precision) && !identical(precision, "auto")) {
      stop(paste(
        "precision is server-owned for formal GLM selectors;",
        "omit it or use precision='auto'"), call. = FALSE)
    }
    # Direct promise forwarding is intentional: the sealed formal adapter
    # rejects unsupported argument names and blocks before forcing a
    # datasource expression or discovering any DSI connection.
    out <- ds.vertGLM(
      formula = formula, data = data, datasources = datasources, ...)
    out <- .dsvert_set_frontdoor(
      out, "ds.vert.glm",
      if (grid) "ds.vertGLM.finite_grid" else if (fresh) {
        "ds.vertGLM.fresh_formal"
      } else {
        "ds.vertGLM.formal"
      },
      if (is.null(datasources)) NULL else length(datasources))
    return(.dsvert_add_policy(out, precision = "server-owned"))
  }
  dots <- list(...)
  datasources <- .dsvert_datasources(datasources)
  args <- .dsvert_apply_binomial_precision(
    c(list(formula = formula, data = data, datasources = datasources),
      dots),
    precision = precision)
  out <- do.call(ds.vertGLM, args)
  out <- .dsvert_set_frontdoor(out, "ds.vert.glm", "ds.vertGLM",
                               length(datasources))
  .dsvert_add_policy(out, precision = precision)
}

#' @rdname ds.vert.aliases
#' @export
ds.vert.cox <- function(formula, data = NULL,
                        method = c("profile", "discrete"),
                        datasources = NULL, analysis_id = NULL,
                        formal_analysis_id = NULL,
                        fresh_formal_analysis_id = NULL, ...) {
  if (sum(!vapply(list(analysis_id, formal_analysis_id,
                       fresh_formal_analysis_id), is.null, logical(1L))) > 1L) {
    stop(paste(
      "analysis_id, formal_analysis_id and fresh_formal_analysis_id are",
      "mutually exclusive"),
         call. = FALSE)
  }
  if (!is.null(analysis_id)) {
    selected <- if (missing(method)) "profile" else match.arg(method)
    if (!identical(selected, "profile")) {
      stop("analysis_id does not accept method='discrete'", call. = FALSE)
    }
    out <- ds.vertCox(
      formula = formula, data = data, datasources = datasources,
      analysis_id = analysis_id, ...)
    return(.dsvert_set_frontdoor(
      out, "ds.vert.cox", "ds.vertCox.profile.analysis",
      if (is.null(datasources)) NULL else length(datasources)))
  }
  if (!is.null(fresh_formal_analysis_id)) {
    selected <- if (missing(method)) "profile" else match.arg(method)
    if (!identical(selected, "profile")) {
      stop("fresh_formal_analysis_id does not accept method='discrete'",
           call. = FALSE)
    }
    out <- ds.vertCox(
      formula = formula, data = data, datasources = datasources,
      fresh_formal_analysis_id = fresh_formal_analysis_id, ...)
    return(.dsvert_set_frontdoor(
      out, "ds.vert.cox", "ds.vertCox.profile.fresh",
      if (is.null(datasources)) NULL else length(datasources)))
  }
  if (!is.null(formal_analysis_id)) {
    selected <- if (missing(method)) "profile" else match.arg(method)
    out <- if (identical(selected, "discrete")) {
      ds.vertCoxDiscreteNonDisclosive(
        formula = formula, data = data, datasources = datasources,
        formal_analysis_id = formal_analysis_id, ...)
    } else {
      ds.vertCox(
        formula = formula, data = data, datasources = datasources,
        formal_analysis_id = formal_analysis_id, ...)
    }
    return(.dsvert_set_frontdoor(
      out, "ds.vert.cox", paste0("ds.vertCox.", selected, ".formal"),
      if (is.null(datasources)) NULL else length(datasources)))
  }
  .dsvert_block_retired_remote_route("cox")
  method <- match.arg(method)
  datasources <- .dsvert_datasources(datasources)
  if (identical(method, "profile")) {
    out <- ds.vertCox(formula = formula, data = data,
                      datasources = datasources, ...)
    route <- "ds.vertCoxProfileNonDisclosive"
  } else {
    out <- ds.vertCoxDiscreteNonDisclosive(
      formula = formula, data = data, target = "discrete_logit",
      datasources = datasources, ...)
    route <- "ds.vertCoxDiscreteNonDisclosive"
  }
  out <- .dsvert_set_frontdoor(out, "ds.vert.cox", route,
                               length(datasources))
  .dsvert_add_policy(out, method = method)
}

#' @rdname ds.vert.aliases
#' @export
ds.vert.coxph <- function(formula, data = NULL, method = "profile",
                          analysis_id = NULL,
                          formal_analysis_id = NULL,
                          fresh_formal_analysis_id = NULL, ...) {
  if (sum(!vapply(list(analysis_id, formal_analysis_id,
                       fresh_formal_analysis_id), is.null, logical(1L))) > 1L) {
    stop(paste(
      "analysis_id, formal_analysis_id and fresh_formal_analysis_id are",
      "mutually exclusive"),
         call. = FALSE)
  }
  if (!is.null(analysis_id)) {
    if (!identical(method, "profile")) {
      stop("analysis_id does not accept method", call. = FALSE)
    }
    out <- ds.vert.cox(
      formula = formula, data = data, analysis_id = analysis_id, ...)
    if (is.list(out)) out$frontdoor <- "ds.vert.coxph"
    return(out)
  }
  if (!is.null(fresh_formal_analysis_id)) {
    if (!identical(method, "profile")) {
      stop("fresh_formal_analysis_id does not accept method", call. = FALSE)
    }
    out <- ds.vert.cox(
      formula = formula, data = data,
      fresh_formal_analysis_id = fresh_formal_analysis_id, ...)
    if (is.list(out)) out$frontdoor <- "ds.vert.coxph"
    return(out)
  }
  if (!is.null(formal_analysis_id)) {
    if (!identical(method, "profile")) {
      stop("formal_analysis_id does not accept method", call. = FALSE)
    }
    out <- ds.vert.cox(
      formula = formula, data = data,
      formal_analysis_id = formal_analysis_id, ...)
    if (is.list(out)) out$frontdoor <- "ds.vert.coxph"
    return(out)
  }
  .dsvert_block_retired_remote_route("cox")
  method <- match.arg(method, "profile")
  out <- ds.vert.cox(formula = formula, data = data, method = method, ...)
  if (is.list(out)) out$frontdoor <- "ds.vert.coxph"
  out
}

#' @rdname ds.vert.aliases
#' @export
ds.vert.nb <- function(formula, data = NULL, server = NULL,
                       method = "accurate",
                       datasources = NULL, ...) {
  extras <- list(...)
  if (!is.null(extras$analysis_id)) {
    if (!is.null(server) || !identical(method, "accurate")) {
      stop("The signed NB2 grid does not accept server or a non-default method",
           call. = FALSE)
    }
    out <- do.call(ds.vertNBFullRegTheta, c(
      list(formula = formula, data = data, datasources = datasources), extras))
    return(.dsvert_set_frontdoor(out, "ds.vert.nb",
                                 "ds.vertNBFullRegTheta", NULL))
  }
  if (!is.null(extras$frequency)) {
    if ("datasources" %in% names(match.call())[-1L] || !is.null(server)) {
      stop("Frequency-backed NB2 does not accept datasources or server", call. = FALSE)
    }
    if (!identical(method, "accurate")) {
      stop("Frequency-backed NB2 does not accept method", call. = FALSE)
    }
    out <- do.call(ds.vertNBFullRegTheta, c(
      list(formula = formula, data = data), extras))
    return(.dsvert_set_frontdoor(out, "ds.vert.nb",
                                 "ds.vertNBFullRegTheta", NULL))
  }
  if (!is.null(server)) {
    if (!identical(method, "accurate")) {
      stop("Frequency-backed NB2 does not accept method", call. = FALSE)
    }
    out <- do.call(ds.vertNBFullRegTheta, c(
      list(formula = formula, data = data, server = server,
           datasources = datasources), extras))
    return(.dsvert_set_frontdoor(out, "ds.vert.nb",
                                 "ds.vertNBFullRegTheta", NULL))
  }
  .dsvert_block_retired_remote_route("negative_binomial")
  method <- match.arg(method, "accurate")
  datasources <- .dsvert_datasources(datasources)
  K <- length(datasources)
  route <- "ds.vertNBFullRegTheta"
  fit <- ds.vertNBFullRegTheta(
    formula = formula, data = data, variant = "full_reg_nd",
    datasources = datasources, ...)
  fit <- .dsvert_set_frontdoor(fit, "ds.vert.nb", route, K)
  .dsvert_add_policy(fit, method = method)
}

#' @rdname ds.vert.aliases
#' @export
ds.vert.multinom <- function(formula, data = NULL, server = NULL,
                             datasources = NULL, ...) {
  extras <- list(...)
  if (!is.null(extras$analysis_id)) {
    if (!is.null(server)) {
      stop("The signed multinomial grid does not accept server", call. = FALSE)
    }
    out <- do.call(ds.vertMultinom, c(
      list(formula = formula, data = data, datasources = datasources), extras))
    return(.dsvert_set_frontdoor(out, "ds.vert.multinom",
                                 "ds.vertMultinom", NULL))
  }
  if (!is.null(extras$frequency)) {
    if ("datasources" %in% names(match.call())[-1L] || !is.null(server)) {
      stop("Frequency-backed multinomial does not accept datasources or server", call. = FALSE)
    }
    out <- do.call(ds.vertMultinom, c(
      list(formula = formula, data = data), extras))
    return(.dsvert_set_frontdoor(out, "ds.vert.multinom",
                                 "ds.vertMultinom", NULL))
  }
  if (!is.null(server)) {
    out <- do.call(ds.vertMultinom, c(
      list(formula = formula, data = data, server = server,
           datasources = datasources), extras))
    return(.dsvert_set_frontdoor(out, "ds.vert.multinom",
                                 "ds.vertMultinom", NULL))
  }
  .dsvert_block_retired_remote_route("multinomial")
  datasources <- .dsvert_datasources(datasources)
  out <- ds.vertMultinom(formula = formula, data = data,
                         datasources = datasources, ...)
  .dsvert_set_frontdoor(out, "ds.vert.multinom",
                        "ds.vertMultinomJointNewton", length(datasources))
}

#' @rdname ds.vert.aliases
#' @export
ds.vert.ordinal <- function(formula, data = NULL, server = NULL,
                            datasources = NULL, analysis_id = NULL, ...) {
  extras <- list(...)
  if (!is.null(analysis_id)) {
    if (!is.null(server) || !is.null(extras$frequency)) {
      stop("The signed ordinal grid does not accept server or frequency",
           call. = FALSE)
    }
    out <- do.call(ds.vertOrdinal, c(
      list(formula = formula, data = data, datasources = datasources,
           analysis_id = analysis_id), extras))
    return(.dsvert_set_frontdoor(out, "ds.vert.ordinal",
                                 "ds.vertOrdinal", NULL))
  }
  if (!is.null(extras$frequency)) {
    if ("datasources" %in% names(match.call())[-1L] || !is.null(server)) {
      stop("Frequency-backed ordinal does not accept datasources or server", call. = FALSE)
    }
    out <- do.call(ds.vertOrdinal, c(
      list(formula = formula, data = data), extras))
    return(.dsvert_set_frontdoor(out, "ds.vert.ordinal",
                                 "ds.vertOrdinal", NULL))
  }
  if (!is.null(server)) {
    out <- do.call(ds.vertOrdinal, c(
      list(formula = formula, data = data, server = server,
           datasources = datasources), extras))
    return(.dsvert_set_frontdoor(out, "ds.vert.ordinal",
                                 "ds.vertOrdinal", NULL))
  }
  .dsvert_block_retired_remote_route("ordinal")
  datasources <- .dsvert_datasources(datasources)
  out <- ds.vertOrdinal(formula = formula, data = data,
                        datasources = datasources, ...)
  .dsvert_set_frontdoor(out, "ds.vert.ordinal",
                        "ds.vertOrdinalJointNewton", length(datasources))
}

#' @rdname ds.vert.aliases
#' @export
ds.vert.lmm <- function(formula, data = NULL, cluster_col,
                        analysis_id = NULL,
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
  fit <- ds.vertLMM(formula = formula, data = data,
                    cluster_col = cluster_col,
                    analysis_id = analysis_id,
                    reml = FALSE,
                    max_iter = max_iter,
                    inner_iter = inner_iter,
                    tol = tol %||% 1e-4,
                    ring = ring %||% "ring63",
                    verbose = verbose,
                    datasources = datasources)
  .dsvert_set_frontdoor(fit, "ds.vert.lmm", "ds.vertLMM", NULL)
}

#' @rdname ds.vert.aliases
#' @export
ds.vert.gee <- function(formula, data = NULL,
                        precision = c("auto", "high", "fast"),
                        datasources = NULL, formal_analysis_id = NULL,
                        fresh_formal_analysis_id = NULL,
                        dp_analysis_id = NULL, analysis_id = NULL, ...) {
  selected_analysis_ids <- sum(!vapply(
    list(formal_analysis_id, fresh_formal_analysis_id, dp_analysis_id,
         analysis_id),
    is.null, logical(1L)))
  if (selected_analysis_ids > 1L) {
    stop(paste(
      "formal_analysis_id, fresh_formal_analysis_id, dp_analysis_id and",
      "analysis_id are",
      "mutually exclusive"),
         call. = FALSE)
  }
  if (selected_analysis_ids == 1L) {
    explicit_arguments <- names(match.call())[-1L]
    if ("precision" %in% explicit_arguments) {
      stop(paste(
        "signed-artifact GEE does not accept legacy controls:",
        "precision"), call. = FALSE)
    }
    arguments <- list(formula = formula, data = data, datasources = datasources)
    if (!is.null(formal_analysis_id)) {
      arguments$formal_analysis_id <- formal_analysis_id
    } else if (!is.null(fresh_formal_analysis_id)) {
      arguments$fresh_formal_analysis_id <- fresh_formal_analysis_id
    } else if (!is.null(dp_analysis_id)) {
      arguments$dp_analysis_id <- dp_analysis_id
    } else {
      arguments$analysis_id <- analysis_id
    }
    out <- do.call(ds.vertGEE, c(arguments, list(...)))
    return(.dsvert_set_frontdoor(out, "ds.vert.gee", "ds.vertGEE",
                                 length(datasources)))
  }
  .dsvert_block_retired_remote_route("gee")
  precision <- match.arg(precision)
  datasources <- .dsvert_datasources(datasources)
  args <- .dsvert_apply_binomial_precision(
    c(list(formula = formula, data = data, datasources = datasources),
      list(...)),
    precision = precision)
  out <- do.call(ds.vertGEE, args)
  out <- .dsvert_set_frontdoor(out, "ds.vert.gee", "ds.vertGEE",
                               length(datasources))
  .dsvert_add_policy(out, precision = precision)
}

#' @rdname ds.vert.aliases
#' @export
ds.vert.glmm <- function(formula, data = NULL, cluster_col,
                         analysis_id = NULL,
                         method = c("auto", "moment", "finite_grid"),
                         datasources = NULL, ...) {
  method <- match.arg(method)
  terms <- if (inherits(formula, "formula")) {
    attr(stats::terms(formula), "term.labels")
  } else character()
  resolved_method <- if (length(terms)) "finite_grid" else "moment"
  if (!identical(method, "auto") && !identical(method, resolved_method)) {
    stop("method does not match the GLMM formula scope", call. = FALSE)
  }
  datasources <- .dsvert_datasources(datasources)
  out <- ds.vertGLMM(formula = formula, data = data,
                     cluster_col = cluster_col, analysis_id = analysis_id,
                     datasources = datasources, ...)
  out <- .dsvert_set_frontdoor(out, "ds.vert.glmm", "ds.vertGLMM",
                               length(datasources))
  .dsvert_add_policy(out, method = resolved_method)
}

#' @rdname ds.vert.aliases
#' @export
ds.vert.ipw <- function(outcome_formula, propensity_formula, data = NULL,
                        precision = c("auto", "high", "fast"),
                        datasources = NULL, ...) {
  precision <- match.arg(precision)
  if (!identical(precision, "auto")) {
    stop("intercept-only IPW has no precision control", call. = FALSE)
  }
  datasources <- .dsvert_datasources(datasources)
  args <- c(list(outcome_formula = outcome_formula,
                 propensity_formula = propensity_formula,
                 data = data, datasources = datasources),
            list(...))
  out <- do.call(ds.vertIPW, args)
  out <- .dsvert_set_frontdoor(out, "ds.vert.ipw", "ds.vertIPW",
                               length(datasources))
  .dsvert_add_policy(out, precision = precision)
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
ds.vert.lasso_iter <- function(formula, data = NULL,
                               method = c("auto", "accurate", "fast"),
                               ...) {
  method <- match.arg(method)
  args <- c(list(formula = formula, data = data), list(...))
  if (is.null(args$exact_non_gaussian)) {
    args$exact_non_gaussian <- !identical(method, "fast")
  }
  route <- if (identical(args$family %||% "gaussian", "gaussian")) {
    "ds.vertLASSOIter(signed-gaussian-synopsis)"
  } else if (isTRUE(args$exact_non_gaussian)) {
    "ds.vertLASSOIter(aggregate-score-unavailable)"
  } else {
    "ds.vertLASSOIter(one-step-surrogate-unavailable)"
  }
  out <- .dsvert_route_result(do.call(ds.vertLASSOIter, args),
                              "ds.vert.lasso_iter", route)
  .dsvert_add_policy(out, method = method)
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
ds.vert.lr <- function(reduced, full,
                       type = c("sampling", "mechanism"), level = 0.95) {
  .dsvert_set_frontdoor(ds.vertLR(
    reduced = reduced, full = full, type = type, level = level),
                        "ds.vert.lr", "ds.vertLR")
}

#' @rdname ds.vert.aliases
#' @export
ds.vert.confint <- function(fit, parm = NULL, level = 0.95,
                            type = c("sampling", "mechanism")) {
  .dsvert_set_frontdoor(ds.vertConfint(
    fit = fit, parm = parm, level = level, type = type),
                        "ds.vert.confint", "ds.vertConfint")
}

#' @rdname ds.vert.aliases
#' @export
ds.vert.wald <- function(fit, parm, null = 0,
                         type = c("sampling", "mechanism"), level = 0.95) {
  .dsvert_set_frontdoor(ds.vertWald(
    fit = fit, parm = parm, null = null, type = type, level = level),
                        "ds.vert.wald", "ds.vertWald")
}

#' @rdname ds.vert.aliases
#' @export
ds.vert.contrast <- function(fit, K, m = NULL,
                             type = c("sampling", "mechanism"),
                             level = 0.95) {
  .dsvert_set_frontdoor(ds.vertContrast(
    fit = fit, K = K, m = m, type = type, level = level),
                        "ds.vert.contrast", "ds.vertContrast")
}
