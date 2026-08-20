.dsvert_effective_glm_ring <- function(family, requested_ring) {
  if (!identical(family, "gaussian")) return(127L)
  as.integer(requested_ring)
}

.dsvert_glm_fit_statistics <- function(family, deviance, n_obs,
                                        n_parameters, y_sd = NULL,
                                        weights_active = FALSE,
                                        offset_active = FALSE,
                                        intercept_included = TRUE) {
  deviance <- as.numeric(deviance)[1L]
  n_obs <- as.integer(n_obs)[1L]
  n_parameters <- as.integer(n_parameters)[1L]

  deviance_type <- if (!is.finite(deviance)) {
    "unavailable"
  } else if (identical(family, "gaussian")) {
    if (isTRUE(weights_active)) "weighted_gaussian_rss" else "gaussian_rss"
  } else {
    paste0("canonical_", family)
  }

  valid_gaussian_null <- identical(family, "gaussian") &&
    !isTRUE(weights_active) && !isTRUE(offset_active) &&
    isTRUE(intercept_included) && length(y_sd) == 1L &&
    is.numeric(y_sd) && is.finite(y_sd) && y_sd > 0 &&
    is.finite(n_obs) && n_obs > 1L
  null_deviance <- if (valid_gaussian_null) {
    (n_obs - 1L) * y_sd^2
  } else {
    NA_real_
  }
  pseudo_r2 <- if (is.finite(deviance) && is.finite(null_deviance) &&
                   null_deviance > 0) {
    1 - deviance / null_deviance
  } else {
    NA_real_
  }

  aic <- NA_real_
  aic_type <- "unavailable"
  if (is.finite(deviance) && is.finite(n_obs) && n_obs > 0L &&
      is.finite(n_parameters) && n_parameters >= 0L) {
    if (identical(family, "gaussian") && !isTRUE(weights_active) &&
        deviance > 0) {
      # Gaussian ML estimates one additional scale/dispersion parameter.
      aic <- n_obs * (log(2 * pi) + 1 + log(deviance / n_obs)) +
        2 * (n_parameters + 1L)
      aic_type <- "gaussian_ml_including_dispersion"
    } else if (identical(family, "binomial") &&
               !isTRUE(weights_active)) {
      aic <- deviance + 2 * n_parameters
      aic_type <- "binomial_canonical"
    }
  }

  list(
    deviance_type = deviance_type,
    null_deviance = null_deviance,
    pseudo_r2 = pseudo_r2,
    aic = aic,
    aic_type = aic_type
  )
}

.dsvert_glm_alignment_metadata <- function(status, n_obs) {
  unavailable <- function() list(
    alignment_attested = FALSE,
    alignment_manifest_hash = NULL,
    cohort_id = NULL
  )
  if (!is.list(status) || !isTRUE(status$aligned) ||
      length(status$n_common) != 1L || !is.numeric(status$n_common) ||
      !is.finite(status$n_common) || !isTRUE(status$n_common == n_obs) ||
      !is.list(status$manifests) || !length(status$manifests)) {
    return(unavailable())
  }
  hashes <- vapply(status$manifests, function(manifest) {
    if (!is.list(manifest) || !is.character(manifest$hash) ||
        length(manifest$hash) != 1L || is.na(manifest$hash)) return(NA_character_)
    manifest$hash
  }, character(1L))
  if (anyNA(hashes) || length(unique(hashes)) != 1L ||
      !grepl("^[0-9a-f]{64}$", hashes[[1L]])) {
    return(unavailable())
  }
  list(
    alignment_attested = TRUE,
    alignment_manifest_hash = hashes[[1L]],
    cohort_id = hashes[[1L]]
  )
}

#' @title Sticky-Synopsis Gaussian GLM frontdoor
#' @description This public frontdoor has one available analysis route:
#'   an explicit \code{dp_analysis_id} with \code{family = "gaussian"}
#'   delegates to \code{ds.vertDPGaussian()} and returns that signed,
#'   contribution-bounded sticky joint-DP Synopsis estimand. A call without
#'   \code{dp_analysis_id} raises a typed \code{dsvert_route_unavailable}
#'   condition before any DSI call. An explicit \code{formal_analysis_id} for
#'   binomial or Poisson also fails before DSI until its durable worker and
#'   single common joint-DP opening are promoted. The retained iterative
#'   Ring/Beaver code below the local gate is unreachable through this
#'   frontdoor and carries no public disclosure, accuracy, or availability
#'   claim.
#'
#' @details
#' \strong{Available route.} Supply an additive formula, an aligned data name,
#' \code{family = "gaussian"} and a custodian-configured
#' \code{dp_analysis_id}. The signed artifact owns clipping bounds, the
#' complete-case cohort, contribution caps, privacy parameters and variable
#' ownership. The adapter accepts only arguments that describe that bounded
#' Gaussian estimand; it never falls back to the retired iterative GLM.
#'
#' \strong{Unavailable routes.} Default/no-id calls and all legacy iterative
#' routes stop locally with zero DSI calls. The \code{formal_analysis_id}
#' selector is reserved for binomial/logit and Poisson/log models but is also
#' sealed locally in this release. No binomial or Poisson fit is therefore
#' returned by this function.
#'
#' @param formula,data,x_vars,y_server Additive model specification, aligned
#'   data name and optional signed-artifact ownership checks for the Gaussian
#'   Synopsis route.
#' @param family Must be \code{"gaussian"} with \code{dp_analysis_id}.
#'   Binomial and Poisson are not available through the public frontdoor.
#' @param lambda,no_intercept,data_name,y_var,missing Gaussian Synopsis
#'   estimand selectors. \code{lambda} is the explicit non-negative ridge
#'   penalty; \code{missing}, when supplied, must be
#'   \code{"complete_case_capsule"}.
#' @param verbose,datasources Progress flag and DataSHIELD connections used
#'   only after the Gaussian signed-artifact request has passed local checks.
#' @param dp_analysis_id Custodian-configured signed bounded Gaussian artifact
#'   id. This is required for the available route.
#' @param formal_analysis_id Reserved custodian-configured binomial/Poisson
#'   selector. Supplying it returns a typed
#'   \code{dsvert_formal_glm_frontdoor_unavailable} condition before DSI.
#' @param max_iter,tol,log_n,offset,weights,ring,binomial_sigmoid_intervals,eta_privacy,keep_session,std_mode,start,compute_se,compute_deviance,gradient_only,numeric_backend
#'   Retained legacy arguments. They are rejected when explicitly supplied to
#'   the Gaussian Synopsis adapter, and the no-id legacy route is unavailable.
#' @return With a valid Gaussian \code{dp_analysis_id}, a
#'   \code{ds.vertDPGaussian} object containing bounded noisy sufficient-
#'   statistic regression output and no classical standard errors, p-values,
#'   individual fitted values, residuals or scores. All other routes raise a
#'   typed condition before DSI and return no fitted object.
#' @seealso \code{\link{ds.vertDPGaussian}},
#'   \code{\link{ds.vertMethodStatus}}
#' @examples
#' \dontrun{
#' fit <- ds.vertGLM(
#'   y ~ x1 + x2, data = "D", family = "gaussian",
#'   dp_analysis_id = "custodian-gaussian-analysis")
#' }
#' @importFrom DSI datashield.aggregate datashield.connections_find
#' @export
ds.vertGLM <- function(formula, data = NULL, x_vars = NULL, y_server = NULL,
                       family = "gaussian", max_iter = 100, tol = 1e-4,
                       lambda = 0, log_n = 12,
                       offset = NULL, weights = NULL,
                       # Ring63 (frac_bits=20, default, back-compat) or
                       # Ring127 (frac_bits=50, STRICT-capable per
                       # Catrina-Saxena 2^-fracbits scaling). Ring127
                       # routes through the Uint128 Go primitives that
                       # already exist for task #116 Cox/LMM. Used by
                       # IPW/#98 for STRICT closure; other families may
                       # opt in as the Ring127 regression suite expands.
                       ring = 63L,
                       binomial_sigmoid_intervals = NULL,
                       verbose = TRUE, datasources = NULL,
                       eta_privacy = "auto",
                       # Keep the MPC session alive on the servers after
                       # the fit returns. Exposes fit$session_id and
                       # fit$transport_pks + fit$server_list so follow-on
                       # helpers (e.g. ds.vertLMM's cluster-residual pass,
                       # ds.vertGEE's sandwich meat) can reuse the
                       # already-aligned shares without re-running PSI +
                       # transport-keys + standardisation. Caller is
                       # responsible for eventually calling
                       # mpcCleanupDS(session_id) on every server.
                       keep_session = FALSE,
                       # Suppress the auto-added intercept. Useful when
                       # the caller is supplying a pre-transformed
                       # design matrix in which one of the predictor
                       # columns already encodes the intercept (e.g.
                       # ds.vertLMM's cluster-mean-centred GLS fit where
                       # "1 - lambda_i" replaces the constant).
                       no_intercept = FALSE,
                       # "full" (center+scale), "scale_only" (sd only,
                       # preserves column means), or "none" (raw).
                       # ds.vertLMM's closed-form GLS path uses
                       # "scale_only" + no_intercept=TRUE.
                       std_mode = "full",
                       # Internal fixed-evaluation path for follow-on
                       # methods such as LMM variance components.
                       start = NULL,
                       compute_se = TRUE,
                       compute_deviance = TRUE,
                       gradient_only = FALSE,
                       # Legacy positional args for backward compatibility
                       data_name = NULL, y_var = NULL,
                       missing = "fail",
                       numeric_backend = "auto",
                       dp_analysis_id = NULL,
                       formal_analysis_id = NULL) {
  call_matched <- match.call()

  if (!is.null(dp_analysis_id) && !is.null(formal_analysis_id)) {
    stop("dp_analysis_id and formal_analysis_id are mutually exclusive",
         call. = FALSE)
  }

  if (!is.null(formal_analysis_id)) {
    return(.dsvert_formal_glm_frontdoor_adapter(
      explicit_arguments = names(call_matched)[-1L],
      formula = if (missing(formula)) NULL else formula,
      data = data, family = family, verbose = verbose,
      datasources = datasources, analysis_id = formal_analysis_id))
  }

  if (!is.null(dp_analysis_id)) {
    return(.dsvert_dp_gaussian_glm_adapter(
      explicit_arguments = names(call_matched)[-1L],
      formula = if (missing(formula)) NULL else formula,
      data = data, x_vars = x_vars, y_server = y_server,
      family = family, lambda = lambda, no_intercept = no_intercept,
      data_name = data_name, y_var = y_var, missing = missing,
      verbose = verbose, datasources = datasources,
      analysis_id = dp_analysis_id))
  }

  .dsvert_block_retired_remote_route("legacy_glm")

  if (!is.character(missing) || length(missing) != 1L ||
      !missing %in% c("fail", "mean_impute")) {
    stop("missing must be one of 'fail' or 'mean_impute'", call. = FALSE)
  }
  if (!is.character(family) || length(family) != 1L || is.na(family) ||
      !family %in% c("gaussian", "binomial", "poisson")) {
    stop("family must be 'gaussian', 'binomial', or 'poisson'",
         call. = FALSE)
  }
  if (!is.numeric(max_iter) || length(max_iter) != 1L ||
      !is.finite(max_iter) || max_iter < 0 || max_iter != floor(max_iter)) {
    stop("max_iter must be one non-negative integer", call. = FALSE)
  }
  max_iter <- as.integer(max_iter)
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0) {
    stop("tol must be one finite positive number", call. = FALSE)
  }
  if (!is.numeric(lambda) || length(lambda) != 1L ||
      !is.finite(lambda) || lambda < 0) {
    stop("lambda must be one finite non-negative number", call. = FALSE)
  }
  if (!is.numeric(ring) || length(ring) != 1L || !is.finite(ring) ||
      ring != floor(ring) || !as.integer(ring) %in% c(63L, 127L)) {
    stop("ring must be 63 or 127", call. = FALSE)
  }
  requested_ring <- as.integer(ring)
  ring <- requested_ring
  if (!is.character(numeric_backend) || length(numeric_backend) != 1L ||
      is.na(numeric_backend) ||
      !tolower(numeric_backend) %in%
        c("auto", "ring63", "ring127", "exact_gc", "multiprecision")) {
    stop("numeric_backend must be auto, ring63, ring127, exact_gc, or ",
         "multiprecision", call. = FALSE)
  }
  numeric_backend <- tolower(numeric_backend)
  if (!is.character(eta_privacy) || length(eta_privacy) != 1L ||
      is.na(eta_privacy) ||
      !eta_privacy %in% c("auto", "k2_beaver", "secure_agg")) {
    stop("eta_privacy must be 'auto', 'k2_beaver', or 'secure_agg'",
         call. = FALSE)
  }
  if (!is.character(std_mode) || length(std_mode) != 1L || is.na(std_mode) ||
      !std_mode %in% c("full", "scale_only", "none")) {
    stop("std_mode must be 'full', 'scale_only', or 'none'", call. = FALSE)
  }
  logical_args <- list(
    verbose = verbose, keep_session = keep_session,
    no_intercept = no_intercept, compute_se = compute_se,
    compute_deviance = compute_deviance, gradient_only = gradient_only
  )
  invalid_logical <- names(logical_args)[!vapply(
    logical_args,
    function(x) is.logical(x) && length(x) == 1L && !is.na(x),
    logical(1L)
  )]
  if (length(invalid_logical)) {
    stop(paste(invalid_logical, collapse = ", "),
         " must be non-missing logical scalars", call. = FALSE)
  }
  if (!is.null(start) &&
      (!is.numeric(start) || !length(start) || any(!is.finite(start)))) {
    stop("start must be NULL or a non-empty finite numeric vector",
         call. = FALSE)
  }
  for (argument in c("offset", "weights")) {
    value <- get(argument, inherits = FALSE)
    if (!is.null(value) &&
        (!is.character(value) || length(value) != 1L || is.na(value) ||
         !nzchar(value))) {
      stop(argument, " must be NULL or one non-empty column name",
           call. = FALSE)
    }
  }
  if (identical(family, "gaussian") && !is.null(offset)) {
    stop(
      "Gaussian offsets are not implemented on the standardized MPC scale; ",
      "use an explicit offset-adjusted response or a supported Poisson offset",
      call. = FALSE
    )
  }

  # ===========================================================================
  # Parse formula or legacy arguments
  # ===========================================================================
  if (!missing(formula)) {
    if (inherits(formula, "formula")) {
      formula_spec <- .dsvert_plain_formula(formula)
      y_var <- formula_spec$response
      if (is.null(x_vars)) x_vars <- formula_spec$predictors
      if (!formula_spec$intercept) no_intercept <- TRUE
      data_name <- data
    } else if (is.character(formula) && grepl("~", formula)) {
      formula_spec <- .dsvert_plain_formula(formula)
      y_var <- formula_spec$response
      if (is.null(x_vars)) x_vars <- formula_spec$predictors
      if (!formula_spec$intercept) no_intercept <- TRUE
      data_name <- data
    } else if (is.character(formula) && !grepl("~", formula)) {
      # Legacy: first arg is data_name (backward compat)
      data_name <- formula
      y_var <- data
      # x_vars stays as passed
    }
  }
  if (is.null(data_name) && !is.null(data)) data_name <- data

  # ===========================================================================
  # Input Validation + Smart Auto-Detection
  # ===========================================================================

  # data can be: "DA" (same name on all servers) or list(s1="tableA", s2="tableB")
  # For now we require all servers use the same data frame name (standard in DataSHIELD).
  if (is.list(data_name) && !is.null(names(data_name))) {
    # Named list: validate and extract. For future use.
    # Currently DataSHIELD requires the same assign name on all servers.
    stop("Named list for data is not yet supported. Use ds.psiAlign() to align data first, ",
         "then pass the aligned name (e.g. 'DA').", call. = FALSE)
  }
  if (!is.character(data_name) || length(data_name) != 1)
    stop("data must be a single character string (the name of the aligned data frame on all servers).",
         call. = FALSE)
  if (!is.character(y_var) || length(y_var) != 1)
    stop("y_var must be a single character string", call. = FALSE)
  if (!is.null(x_vars) && !is.character(x_vars) && !is.list(x_vars)) {
    stop("x_vars must be NULL, a character vector, or a named list",
         call. = FALSE)
  }
  if (is.list(x_vars) &&
      (is.null(names(x_vars)) || any(!nzchar(names(x_vars))) ||
       anyDuplicated(names(x_vars)) ||
       any(!vapply(x_vars, is.character, logical(1L))))) {
    stop("x_vars must be a named list mapping server names to variable vectors",
         call. = FALSE)
  }
  if (!is.null(binomial_sigmoid_intervals)) {
    binomial_sigmoid_intervals <- as.integer(binomial_sigmoid_intervals)
    if (length(binomial_sigmoid_intervals) != 1L ||
        !is.finite(binomial_sigmoid_intervals) ||
        binomial_sigmoid_intervals < 10L) {
      stop("binomial_sigmoid_intervals must be NULL or an integer >= 10",
           call. = FALSE)
    }
    if (identical(family, "binomial")) {
      old_binomial_intervals <- getOption(
        "dsvert.glm_num_intervals_binomial", NULL)
      on.exit(options(
        dsvert.glm_num_intervals_binomial = old_binomial_intervals),
        add = TRUE)
      options(dsvert.glm_num_intervals_binomial =
                binomial_sigmoid_intervals)
    }
  }
  effective_binomial_sigmoid_intervals <- if (identical(family, "binomial")) {
    suppressWarnings(as.integer(getOption(
      "dsvert.glm_num_intervals_binomial",
      getOption("dsvert.glm_num_intervals", 100L))[[1L]]))
  } else {
    NA_integer_
  }
  if (identical(family, "binomial") &&
      (!is.finite(effective_binomial_sigmoid_intervals) ||
       effective_binomial_sigmoid_intervals < 10L)) {
    effective_binomial_sigmoid_intervals <- 100L
  }

  if (is.null(datasources))
    datasources <- DSI::datashield.connections_find()
  if (length(datasources) == 0)
    stop("No DataSHIELD connections found", call. = FALSE)

  # Authenticate and validate every participating custodian's public numeric
  # policy before querying schema, row counts, alignment state, or creating MPC
  # state. There is no unattested compatibility route.
  numeric_policies <- .dsvert_require_numeric_policies(datasources)

  # Auto-detect: query servers for their columns and map variables automatically
  col_results <- NULL
  if (is.null(x_vars) || is.character(x_vars)) {
    user_x_vars <- x_vars  # NULL = use all available, character = specific vars
    if (verbose) message("[Auto-detect] Querying server columns...")
    col_results <- .dsvert_aggregate_strict(
      conns = datasources,
      expr = call(name = "dsvertColNamesDS", data_name = data_name),
      operation = "GLM column discovery")
    server_names <- names(datasources)

    # Build column map: which server has which variable (exclude IDs)
    col_map <- list()
    all_available <- character(0)
    for (srv in server_names) {
      cols <- setdiff(col_results[[srv]]$columns, c("id", "patient_id"))
      col_map[[srv]] <- cols
      all_available <- c(all_available, cols)
    }
    all_available <- unique(all_available)

    # Validate: all requested variables must exist somewhere
    if (!is.null(user_x_vars)) {
      missing_variables <- setdiff(user_x_vars, all_available)
      if (length(missing_variables) > 0)
        stop("Variables not found on any server: ",
             paste(missing_variables, collapse = ", "),
             "\n  Available: ", paste(all_available, collapse = ", "), call. = FALSE)
    }
    if (!y_var %in% all_available)
      stop("Response variable '", y_var, "' not found on any server.\n",
           "  Available columns per server:\n",
           paste(sapply(server_names, function(s)
             paste0("    ", s, ": ", paste(col_map[[s]], collapse = ", "))),
             collapse = "\n"), call. = FALSE)

    # Find y_server automatically
    if (is.null(y_server)) {
      y_servers <- server_names[sapply(server_names, function(s) y_var %in% col_map[[s]])]
      y_server <- y_servers[1]
      if (length(y_servers) > 1 && verbose)
        message("  Note: '", y_var, "' found on multiple servers (",
                paste(y_servers, collapse = ", "), "). Using: ", y_server)
      if (verbose) message("  y_var '", y_var, "' found on: ", y_server)
    }

    # Build x_vars automatically
    x_vars <- list()
    for (srv in server_names) {
      feats <- col_map[[srv]]
      # Always exclude y_var from features (even on non-label servers)
      feats <- setdiff(feats, y_var)
      # If user specified specific vars, filter to those
      if (!is.null(user_x_vars))
        feats <- intersect(feats, user_x_vars)
      x_vars[[srv]] <- feats
    }

    # Remove response-only vars from feature lists
    # Keep server in x_vars even if it has 0 features (label server)
    if (verbose) {
      for (srv in server_names) {
        role <- if (srv == y_server) " (label)" else ""
        if (length(x_vars[[srv]]) > 0)
          message("  ", srv, role, ": ", paste(x_vars[[srv]], collapse = ", "))
        else
          message("  ", srv, role, ": (response only)")
      }
    }
  }

  if (!is.list(x_vars) || is.null(names(x_vars)))
    stop("x_vars must be a named list mapping server names to variable vectors",
         call. = FALSE)
  if (is.null(y_server))
    stop("y_server must be specified: the server holding '", y_var, "'",
         call. = FALSE)
  if (!y_server %in% names(x_vars))
    stop("y_server '", y_server, "' must be in x_vars", call. = FALSE)
  if (anyDuplicated(names(x_vars)) || any(!nzchar(names(x_vars))) ||
      any(!vapply(x_vars, is.character, logical(1L)))) {
    stop("x_vars must map unique non-empty server names to character vectors",
         call. = FALSE)
  }
  predictor_names <- unlist(x_vars, use.names = FALSE)
  if (anyNA(predictor_names) || any(!nzchar(predictor_names)) ||
      anyDuplicated(predictor_names)) {
    stop("Every predictor must be a unique non-empty column across servers",
         call. = FALSE)
  }

  n_partitions_check <- length(x_vars)
  non_label_count <- n_partitions_check - 1L
  eta_privacy <- .dsvert_select_eta_privacy(
    eta_privacy, n_partitions = n_partitions_check
  )

  use_secure_agg <- (eta_privacy == "secure_agg")
  use_k2_beaver <- (eta_privacy == "k2_beaver")

  if (isTRUE(gradient_only)) {
    if (is.null(start)) {
      stop("gradient_only requires a standardized `start` vector",
           call. = FALSE)
    }
    compute_deviance <- FALSE
    max_iter <- max(1L, as.integer(max_iter))
  }

  # Adaptive log_n for large n: max_slots = 2^(log_n-1)
  # n <= 4096 -> log_n=13, n <= 8192 -> log_n=14, n <= 16384 -> log_n=15
  # Note: this is set here before n_obs is known; adjusted later after getObsCountDS

  # ===========================================================================
  # Setup (datasources already resolved above for auto-detect)
  # ===========================================================================

  server_names <- names(datasources)
  if (!all(names(x_vars) %in% server_names)) {
    missing_servers <- setdiff(names(x_vars), server_names)
    stop("Unknown server(s) in x_vars: ",
         paste(missing_servers, collapse = ", "),
         call. = FALSE)
  }

  server_list <- names(x_vars)
  non_label_servers <- setdiff(server_list, y_server)
  n_partitions <- length(x_vars)

  # Offset/weight discovery needs the same public column catalogue. Query all
  # participating sites once, never sequentially until a match is found.
  if ((!is.null(offset) || !is.null(weights)) && is.null(col_results)) {
    col_results <- .dsvert_aggregate_strict(
      conns = datasources[server_list],
      expr = call(name = "dsvertColNamesDS", data_name = data_name),
      operation = "GLM column discovery")
  }

  # Get observation count (lightweight sync call before helpers are available)
  first_conn <- which(server_names == server_list[1])
  count_result <- .dsvert_aggregate_strict(
    conns = datasources[first_conn],
    expr = call(name = "getObsCountDS", data_name),
    operation = "GLM observation-count preflight"
  )
  if (is.list(count_result) && length(count_result) == 1)
    count_result <- count_result[[1]]
  n_obs <- count_result$n_obs

  # Reuse an existing authenticated PSI row-order manifest when available.
  # This does not execute PSI; older/non-PSI deployments simply produce no
  # attestation, in which case cohort-sensitive LR inference fails closed.
  alignment_status <- .psi_alignment_status(
    data_name, datasources[server_list]
  )
  alignment_metadata <- .dsvert_glm_alignment_metadata(
    alignment_status, n_obs
  )

  # Adaptive log_n: ensure max_slots >= n_obs
  max_slots <- 2^(log_n - 1)
  if (n_obs > max_slots) {
    new_log_n <- ceiling(log2(n_obs)) + 1
    if (new_log_n > 15) {
      warning(sprintf("n_obs=%d exceeds max for log_n=15 (%d slots). Using log_n=15 with observation batching.",
                       n_obs, 2^14))
      new_log_n <- 15L
    }
    if (verbose) message(sprintf("  [Adaptive] log_n bumped %d->%d for n=%d observations", log_n, new_log_n, n_obs))
    log_n <- as.integer(new_log_n)
  }

  n_vars_total <- sum(sapply(x_vars, length))

  # Numeric preflight runs before a session or any MPC state is created.  All
  # magnitude/error bounds come from the custodians; no analyst argument can
  # relax them.  The current R orchestrator has adapters only for ring63/127,
  # so an attested future backend is still refused until its adapter exists.
  numeric_workload <- .dsvert_numeric_glm_workload(
    n_obs = n_obs,
    n_predictors = n_vars_total,
    family = family,
    max_iter = max_iter,
    compute_se = isTRUE(compute_se),
    compute_deviance = isTRUE(compute_deviance),
    weights_active = !is.null(weights),
    offset_active = !is.null(offset))
  numeric_certificate <- .dsvert_numeric_preflight_from_policies(
    numeric_policies[server_list], numeric_workload,
    requested_backend = numeric_backend,
    requested_ring = requested_ring)
  if (!numeric_certificate$effective_backend %in% c("ring63", "ring127")) {
    unavailable_certificate <- numeric_certificate
    unavailable_certificate$status <- "numeric_backend_unavailable"
    unavailable_certificate$reason <- "client_backend_adapter_unavailable"
    .dsvert_stop_numeric(
      "numeric_backend_unavailable",
      paste0("The selected ", numeric_certificate$effective_backend,
             " backend has no integrated ds.vertGLM client adapter"),
      unavailable_certificate, "client_backend_adapter_unavailable")
  }
  effective_ring <- as.integer(numeric_certificate$ring_bits)
  ring <- effective_ring

  if (verbose) {
    message(sprintf("=== Encrypted-Label BCD-IRLS for %s GLM ===", family))
    message(sprintf("Observations: %d, Variables: %d, Partitions: %d",
                    n_obs, n_vars_total, n_partitions))
    message(sprintf("Label server: %s (holds '%s')", y_server, y_var))
    if (length(non_label_servers) > 0)
      message(sprintf("Non-label servers: %s",
                      paste(non_label_servers, collapse = ", ")))
  }

  # Generate session_id for all protocol phases (crypto state isolation)
  session_id <- .dsvert_uuid4()

  # Guaranteed cleanup on exit (even if error occurs mid-protocol).
  # Cleanup remains best-effort because it runs only after the protocol can no
  # longer advance. It uses one async fan-out per cleanup phase and suppresses
  # connector diagnostics.
  # Skipped when keep_session = TRUE so follow-on helpers can reuse the
  # MPC session (caller assumes cleanup responsibility).
  on.exit({
    if (!isTRUE(keep_session)) {
      .dsvert_cleanup_best_effort(
        datasources[server_list],
        call(name = "mpcCleanupDS", session_id = session_id))
      .dsvert_cleanup_best_effort(
        datasources[server_list], call(name = "mpcGcDS"))
    }
    .dsvert_reset_chunk_size()
  }, add = TRUE)

  # ===========================================================================
  # Phase 0-1: Transport key setup + standardize features
  #   (delegated to .glm_mpc_setup in ds.vertGLM.setup.R)
  # ===========================================================================
  setup <- .glm_mpc_setup(
    datasources      = datasources,
    server_names     = server_names,
    server_list      = server_list,
    non_label_servers = non_label_servers,
    y_server         = y_server,
    y_var            = y_var,
    x_vars           = x_vars,
    data_name        = data_name,
    family           = family,
    session_id       = session_id,
    verbose          = verbose,
    # no_intercept=TRUE: skip y-standardisation (no mean-shift).
    standardize_y_override = if (isTRUE(no_intercept)) FALSE else NULL,
    std_mode = std_mode,
    missing = missing,
    numeric_ring = ring
  )

  # Unpack setup results
  transport_pks <- setup$transport_pks
  x_means       <- setup$x_means
  x_sds         <- setup$x_sds
  y_mean        <- setup$y_mean
  y_sd          <- setup$y_sd
  std_data      <- setup$std_data
  standardize_y <- setup$standardize_y
  missing_policy <- setup$missing_policy
  numeric_attestations <- setup$numeric_attestations
  .dsAgg        <- setup$.dsAgg
  .sendBlob     <- setup$.sendBlob

  if (length(numeric_certificate$policy_ids)) {
    input_bindings <- vapply(server_list, function(server) {
      .dsvert_numeric_validate_attestation(
        attestation = numeric_attestations[[server]],
        certificate = numeric_certificate,
        server = server,
        kind = "glm_standardized_input",
        session_id = session_id,
        data_name = data_name,
        variables = unique(c(
          x_vars[[server]], if (server == y_server) y_var else NULL)),
        family = family,
        ring = ring,
        n = n_obs)
    }, character(1L))
    numeric_certificate <- .dsvert_numeric_attach_attestations(
      numeric_certificate, input_bindings, all_inputs = TRUE)
  }

  if (max_iter < 1L && !isTRUE(compute_se) &&
      !isTRUE(compute_deviance) && !isTRUE(gradient_only) &&
      is.null(start) && is.null(offset) && is.null(weights)) {
    all_names <- unlist(x_vars[server_list], use.names = FALSE)
    all_x_means <- unlist(x_means[server_list], use.names = FALSE)
    all_x_sds <- unlist(x_sds[server_list], use.names = FALSE)
    all_coefs <- c("(Intercept)" = 0, setNames(rep(0, length(all_names)),
                                               all_names))
    n_parameters <- length(all_coefs) - as.integer(isTRUE(no_intercept))
    fit_statistics <- .dsvert_glm_fit_statistics(
      family = family, deviance = NA_real_, n_obs = n_obs,
      n_parameters = n_parameters, y_sd = y_sd,
      weights_active = FALSE, offset_active = FALSE,
      intercept_included = !isTRUE(no_intercept)
    )
    numeric_certificate <- .dsvert_numeric_finalize_certificate(
      numeric_certificate, all_coefs, converged = TRUE)
    result <- list(
      coefficients = all_coefs,
      std_errors = rep(NA_real_, length(all_coefs)),
      statistic_values = rep(NA_real_, length(all_coefs)),
      t_values = if (identical(family, "gaussian")) {
        rep(NA_real_, length(all_coefs))
      } else {
        NULL
      },
      z_values = rep(NA_real_, length(all_coefs)),
      p_values = rep(NA_real_, length(all_coefs)),
      coefficient_reference = if (identical(family, "gaussian")) "t" else "normal",
      covariance = NULL,
      covariance_information = NULL,
      covariance_unscaled = NULL,
      iterations = 0L,
      converged = TRUE,
      family = family,
      n_obs = n_obs,
      n_vars = length(all_coefs),
      n_parameters = n_parameters,
      df_residual = n_obs - n_parameters,
      lambda = lambda,
      deviance = NA_real_,
      deviance_type = fit_statistics$deviance_type,
      null_deviance = fit_statistics$null_deviance,
      pseudo_r2 = fit_statistics$pseudo_r2,
      aic = fit_statistics$aic,
      aic_type = fit_statistics$aic_type,
      y_server = y_server,
      eta_privacy = eta_privacy,
      missing_policy = missing_policy,
      data_name = data_name,
      y_var = y_var,
      x_vars = x_vars,
      offset = offset,
      weights = weights,
      ring = effective_ring,
      requested_ring = requested_ring,
      effective_ring = effective_ring,
      numeric_certificate = numeric_certificate,
      alignment_attested = alignment_metadata$alignment_attested,
      alignment_manifest_hash = alignment_metadata$alignment_manifest_hash,
      cohort_id = alignment_metadata$cohort_id,
      x_means = setNames(all_x_means, all_names),
      x_sds = setNames(all_x_sds, all_names),
      y_sd = if (exists("y_sd", inherits = FALSE)) y_sd else NULL,
      y_mean = if (exists("y_mean", inherits = FALSE)) y_mean else NULL,
      hessian_std = NULL,
      gradient_std = NULL,
      gradient = NULL,
      call = call_matched)
    class(result) <- c("ds.glm", "list")
    return(result)
  }

  # ===========================================================================
  # Offset registration (Poisson / NB rate regression).
  # ===========================================================================
  # When the caller passes offset = "colname", we auto-detect which
  # server holds that column via the col_results from setup (or by
  # querying once here) and register the offset on that server only.
  # The server-side k2SetOffsetDS caches the FP-encoded offset; the
  # modified k2ComputeEtaShareDS picks it up each iteration. Offsets
  # never leave their home server.
  if (!is.null(offset)) {
    if (!is.character(offset) || length(offset) != 1L) {
      stop("offset must be a single character string (column name)",
           call. = FALSE)
    }
    offset_srv <- NULL
    for (.srv in server_list) {
      site_catalog <- col_results[[.srv]]
      cols <- if (is.list(site_catalog) &&
                  is.character(site_catalog$columns)) {
        site_catalog$columns
      } else {
        character()
      }
      if (!is.null(cols) && offset %in% cols) {
        offset_srv <- .srv
        break
      }
    }
    if (is.null(offset_srv)) {
      stop("Offset column '", offset, "' not found on any server",
           call. = FALSE)
    }
    if (verbose) message(sprintf("Registering offset '%s' on server %s",
                                  offset, offset_srv))
    .ci <- which(server_names == offset_srv)
    offset_result <- .dsAgg(datasources[.ci], call(name = "k2SetOffsetDS",
      data_name = data_name,
      offset_column = offset,
      numeric_family = family,
      session_id = session_id))
    if (is.list(offset_result) && length(offset_result) == 1L) {
      offset_result <- offset_result[[1L]]
    }
    if (length(numeric_certificate$policy_ids)) {
      offset_binding <- .dsvert_numeric_validate_attestation(
        offset_result$numeric_attestation, numeric_certificate,
        server = offset_srv, kind = "glm_offset",
        session_id = session_id, data_name = data_name,
        variables = offset, family = family, ring = ring, n = n_obs)
      numeric_certificate <- .dsvert_numeric_attach_attestations(
        numeric_certificate, offset_binding)
    }
  }

  # ===========================================================================
  # Per-patient weights registration (for IPW / survey-weighted regression).
  # ===========================================================================
  # Weights live plaintext only on the server that holds the weights column.
  # The holder splits w and sqrt(w) into additive shares for the two DCF
  # parties; weighted gradients are then computed by share-domain Beaver
  # multiplication of w_share * residual_share. This avoids transferring
  # patient-level weights to a peer that could combine them with local
  # covariates to reconstruct hidden treatment/outcome information.
  weights_active <- FALSE
  if (!is.null(weights)) {
    if (!is.character(weights) || length(weights) != 1L) {
      stop("weights must be a single character string (column name)",
           call. = FALSE)
    }
    weights_srv <- NULL
    for (.srv in server_list) {
      site_catalog <- col_results[[.srv]]
      cols <- if (is.list(site_catalog) &&
                  is.character(site_catalog$columns)) {
        site_catalog$columns
      } else {
        character()
      }
      if (!is.null(cols) && weights %in% cols) {
        weights_srv <- .srv
        break
      }
    }
    if (is.null(weights_srv)) {
      stop("Weights column '", weights, "' not found on any server",
           call. = FALSE)
    }
    if (verbose) message(sprintf("Registering weights '%s' on server %s",
                                  weights, weights_srv))
    weights_ci <- which(server_names == weights_srv)
    weights_ring <- ring
    if (use_secure_agg) {
      # K>=3 uses two DCF parties: fusion server plus coordinator.
      fusion_srv <- .k3_select_fusion_server(server_list, y_server, x_vars)
      dcf_weight_parties <- c(fusion_srv, y_server)
      dcf_role <- if (weights_srv == dcf_weight_parties[1L]) {
        "dcf0"
      } else if (weights_srv == dcf_weight_parties[2L]) {
        "dcf1"
      } else {
        "dealer"
      }

      setres <- .dsAgg(datasources[weights_ci], call(name = "k2ShareWeightsDS",
        data_name = data_name,
        weights_column = weights,
        dcf0_pk = transport_pks[[dcf_weight_parties[1L]]],
        dcf1_pk = transport_pks[[dcf_weight_parties[2L]]],
        dcf_role = dcf_role,
        ring = weights_ring,
        numeric_family = family,
        session_id = session_id))
      if (is.list(setres) && length(setres) == 1L) setres <- setres[[1]]

      if (dcf_role == "dealer") {
        .sendBlob(setres$dcf0_blob, setres$dcf0_transfer,
                  which(server_names == dcf_weight_parties[1L]))
        .sendBlob(setres$dcf1_blob, setres$dcf1_transfer,
                  which(server_names == dcf_weight_parties[2L]))
        .sendBlob(setres$dcf0_sqrt_blob, setres$dcf0_sqrt_transfer,
                  which(server_names == dcf_weight_parties[1L]))
        .sendBlob(setres$dcf1_sqrt_blob, setres$dcf1_sqrt_transfer,
                  which(server_names == dcf_weight_parties[2L]))
        for (srv in dcf_weight_parties) {
          .dsAgg(datasources[which(server_names == srv)],
            call(name = "k2ReceiveWeightSharesDS", numeric_family = family,
                 peer_name = weights_srv,
                 session_id = session_id))
        }
      } else {
        peer_srv <- if (dcf_role == "dcf0") dcf_weight_parties[2L] else dcf_weight_parties[1L]
        peer_ci <- which(server_names == peer_srv)
        .sendBlob(setres$peer_blob, setres$peer_transfer, peer_ci)
        .sendBlob(setres$peer_sqrt_blob, setres$peer_sqrt_transfer, peer_ci)
        .dsAgg(datasources[peer_ci], call(name = "k2ReceiveWeightSharesDS",
          numeric_family = family, peer_name = weights_srv,
          session_id = session_id))
      }
    } else {
      # K=2 DCF party is the other label/non-label server.
      peer_srv <- if (weights_srv == y_server) non_label_servers[1] else y_server
      peer_ci <- which(server_names == peer_srv)
      dcf0_srv <- y_server
      dcf1_srv <- non_label_servers[1]
      dcf_role <- if (weights_srv == dcf0_srv) "dcf0" else "dcf1"
      setres <- .dsAgg(datasources[weights_ci], call(name = "k2ShareWeightsDS",
        data_name = data_name,
        weights_column = weights,
        dcf0_pk = transport_pks[[dcf0_srv]],
        dcf1_pk = transport_pks[[dcf1_srv]],
        dcf_role = dcf_role,
        ring = weights_ring,
        numeric_family = family,
        session_id = session_id))
      if (is.list(setres) && length(setres) == 1L) setres <- setres[[1]]
      .sendBlob(setres$peer_blob, setres$peer_transfer, peer_ci)
      .sendBlob(setres$peer_sqrt_blob, setres$peer_sqrt_transfer, peer_ci)
      .dsAgg(datasources[peer_ci], call(name = "k2ReceiveWeightSharesDS",
        numeric_family = family, peer_name = weights_srv,
        session_id = session_id))
    }
    if (length(numeric_certificate$policy_ids)) {
      weight_binding <- .dsvert_numeric_validate_attestation(
        setres$numeric_attestation, numeric_certificate,
        server = weights_srv, kind = "glm_weights",
        session_id = session_id, data_name = data_name,
        variables = weights, family = family, ring = ring, n = n_obs)
      numeric_certificate <- .dsvert_numeric_attach_attestations(
        numeric_certificate, weight_binding)
    }
    weights_active <- TRUE
  }

  # ===========================================================================
  # Phase 3: Iterative Ring63 Beaver (on standardized scale)
  # ===========================================================================
  label_intercept <- !standardize_y

  coordinator <- y_server
  coordinator_conn <- which(server_names == coordinator)

  betas <- list()
  for (server in server_list) {
    p <- length(x_vars[[server]])
    if (server == coordinator && label_intercept && !use_k2_beaver) p <- p + 1
    betas[[server]] <- rep(0, p)
  }

  converged <- FALSE
  final_iter <- 0

  # Encrypted eta blobs (opaque to client): server_name -> base64url string
  encrypted_etas <- list()

  if (use_k2_beaver) {
    # =========================================================================
    # Phase 3 (K=2 Strict v2): Chebyshev Beaver MPC
    # =========================================================================
    # Secure polynomial evaluation of sigmoid/exp on secret shares.
    # Neither party sees eta, mu, residuals, or weights in plaintext.
    # Only p_k gradient scalars + 2 intercept scalars revealed per iteration.

    if (verbose) message("\n[Phase 3] BCD iterations (K=2 strict Chebyshev Beaver)...")

    nl <- non_label_servers[1]
    nl_conn <- which(server_names == nl)

    loop_result <- .k2_strict_loop(
      datasources = datasources,
      server_names = server_names,
      server_list = server_list,
      coordinator = coordinator,
      coordinator_conn = coordinator_conn,
      non_label_servers = non_label_servers,
      nl = nl,
      nl_conn = nl_conn,
      x_vars = x_vars,
      y_var = y_var,
      std_data = std_data,
      transport_pks = transport_pks,
      session_id = session_id,
      family = family,
      lambda = lambda,
      max_iter = max_iter,
      tol = tol,
      n_obs = n_obs,
      verbose = verbose,
      .dsAgg = .dsAgg,
      .sendBlob = .sendBlob,
      weights_active = isTRUE(weights_active),
      offset_active = !is.null(offset),
      no_intercept  = isTRUE(no_intercept),
      ring = ring,
      compute_se = isTRUE(compute_se),
      compute_deviance = isTRUE(compute_deviance),
      gradient_only = isTRUE(gradient_only),
      start = start,
      numeric_certificate = numeric_certificate
    )

    betas <- loop_result$betas
    converged <- loop_result$converged
    final_iter <- loop_result$iterations
    # Store the loop's intercept for destandardization
    k2_loop_intercept <- loop_result$intercept

  } else if (use_secure_agg) {
    # ALL families: Ring63 Beaver gradient (Gaussian=identity link, others=DCF wide spline)
    if (verbose) message("\n[Phase 3] Ring", ring, " Beaver Gradient (K=",
                         length(server_list), " servers, family=", family, ")...")
    k3_result <- .k3_ring63_gradient_loop(
      datasources = datasources, server_list = server_list,
      server_names = server_names, x_vars = x_vars,
      coordinator = coordinator, coordinator_conn = coordinator_conn,
      non_label_servers = non_label_servers, transport_pks = transport_pks,
      std_data = std_data, y_var = y_var, family = family,
      betas = betas, n_obs = n_obs, lambda = lambda,
      session_id = session_id,
      max_iter = max_iter, tol = tol, verbose = verbose,
      label_intercept = label_intercept,
      .dsAgg = .dsAgg, .sendBlob = .sendBlob,
      weights_active = isTRUE(weights_active),
      no_intercept = isTRUE(no_intercept),
      start = start,
      compute_se = isTRUE(compute_se),
      compute_deviance = isTRUE(compute_deviance),
      gradient_only = isTRUE(gradient_only),
      ring = ring,
      numeric_certificate = numeric_certificate)
    betas <- k3_result$betas
    converged <- k3_result$converged
    final_iter <- k3_result$final_iter

  } else {
    stop("Unsupported eta_privacy mode: ", eta_privacy, call. = FALSE)
  }
  # ===========================================================================
  # Phase 4: Unstandardize coefficients
  # ===========================================================================
  all_coefs_std <- numeric()
  all_x_means <- numeric()
  all_x_sds <- numeric()
  all_names <- character()
  beta_0_from_label <- 0
  # For K=2 beaver (wide spline), use the intercept from the loop
  if (use_k2_beaver && exists("k2_loop_intercept")) {
    beta_0_from_label <- k2_loop_intercept
  }
  if (use_secure_agg && exists("k3_result") &&
      !is.null(k3_result$intercept)) {
    beta_0_from_label <- k3_result$intercept
  }

  for (server in server_list) {
    server_beta <- betas[[server]]
    if (server == coordinator && label_intercept && !use_k2_beaver) {
      beta_0_from_label <- server_beta[1]
      server_beta <- server_beta[-1]
    }
    all_coefs_std <- c(all_coefs_std, server_beta)
    all_x_means <- c(all_x_means, x_means[[server]])
    all_x_sds <- c(all_x_sds, x_sds[[server]])
    all_names <- c(all_names, x_vars[[server]])
  }

  if (standardize_y && !is.null(y_sd)) {
    all_coefs_orig <- all_coefs_std * y_sd / all_x_sds
    # IPW fix (2026-04-21 PM): under weighted fit, the loop's alpha_std
    # absorbs the (ybar_W - ybar)/sigma_y shift because y is centered by the
    # UNWEIGHTED mean but the weighted-score optimum has non-zero
    # mean residual in standardized space. Add the beta_0_from_label *
    # sigma_y term to unstandardize correctly. Under unweighted fit alpha_std
    # converges to approx 0 (weighted-by-unity mean of centered y is 0), so
    # the term is approx 0 and back-compat is preserved -- verified by the
    # Ring63 w_unit probe staying at max|Deltabeta| = 1.12e-4 STRICT.
    intercept <- beta_0_from_label * y_sd + y_mean -
                 sum(all_coefs_orig * all_x_means)
  } else {
    all_coefs_orig <- all_coefs_std / all_x_sds
    intercept <- beta_0_from_label - sum(all_coefs_orig * all_x_means)
  }
  if (isTRUE(no_intercept)) {
    # The caller is supplying a design matrix that already contains a
    # column encoding the intercept (e.g. "1 - lambda_i" for LMM GLS).
    # Report intercept = 0 to avoid double-counting.
    intercept <- 0
  }

  all_coefs <- c(intercept, all_coefs_orig)
  names(all_coefs) <- c("(Intercept)", all_names)
  n_vars_total <- length(all_coefs)
  n_parameters <- n_vars_total - as.integer(isTRUE(no_intercept))
  df_residual <- n_obs - n_parameters

  # ===========================================================================
  # Phase 5: Deviance (server-side, on ORIGINAL-scale data)
  # ===========================================================================
  # Deviance is computed on the coordinator (label server) using stored
  # eta values. The client relays the final encrypted etas from non-label
  # servers so the coordinator can reconstruct eta_total server-side.
  # The client NEVER sees the n-length eta_total vector.

  if (use_secure_agg) {
    # Secure deviance already computed in the Ring63 loop via Beaver dot-product.
    # No individual eta values are revealed. Only the scalar Sum(mu-y)^2 is returned.
    if (!is.null(k3_result$deviance)) {
      if (verbose) message(sprintf("\n[Phase 5] Secure deviance (Beaver): %.4f", k3_result$deviance))
    }
  } else if (use_k2_beaver) {
    # K=2: secure deviance already computed in loop
    if (!is.null(loop_result$deviance)) {
      if (verbose) message(sprintf("\n[Phase 5] Secure deviance (K=2 Beaver): %.4f", loop_result$deviance))
    }
  }

  # Deviance is computed in share-domain: Gaussian weighted/unweighted RSS,
  # binomial canonical deviance, or Poisson canonical deviance. Weighted
  # non-Gaussian fits deliberately return NA until their canonical expression
  # is implemented.
  deviance <- NA; null_deviance <- NA
  if (use_secure_agg && exists("k3_result") && !is.null(k3_result$deviance)) {
    deviance <- k3_result$deviance
  } else if (use_k2_beaver && exists("loop_result") && !is.null(loop_result$deviance)) {
    deviance <- loop_result$deviance
  }
  # Gaussian deviance is computed in standardized space -- destandardize
  if (family == "gaussian" && !is.null(y_sd) && !is.na(deviance)) {
    deviance <- deviance * y_sd^2
  }
  fit_statistics <- .dsvert_glm_fit_statistics(
    family = family, deviance = deviance, n_obs = n_obs,
    n_parameters = n_parameters, y_sd = y_sd,
    weights_active = isTRUE(weights_active),
    offset_active = !is.null(offset),
    intercept_included = !isTRUE(no_intercept)
  )
  null_deviance <- fit_statistics$null_deviance
  pseudo_r2 <- fit_statistics$pseudo_r2
  deviance_type <- fit_statistics$deviance_type
  aic <- fit_statistics$aic
  aic_type <- fit_statistics$aic_type

  if (verbose)
    message(sprintf("\nDeviance: %.4f, Null deviance: %.4f, Pseudo R2: %.4f",
                    deviance, null_deviance, pseudo_r2))

  # ===========================================================================
  # Standard Errors + P-values (exact, via finite-difference Hessian)
  # Computed in STANDARDIZED space, then destandardized via Jacobian.
  # ===========================================================================
  inv_H <- NULL
  if (use_secure_agg && exists("k3_result")) inv_H <- k3_result$inv_hessian
  if (use_k2_beaver && exists("loop_result")) inv_H <- loop_result$inv_hessian

  empty_inference <- stats::setNames(
    rep(NA_real_, length(all_coefs)), names(all_coefs))
  std_errors <- empty_inference
  statistic_values <- empty_inference
  t_values <- NULL
  z_values <- empty_inference
  p_values <- empty_inference
  covariance <- NULL
  covariance_information <- NULL

  if (isTRUE(compute_se) &&
      !is.null(inv_H) && !is.null(attr(inv_H, "raw_hessian"))) {
    H_raw <- attr(inv_H, "raw_hessian")
    # K>=3 finite-difference Hessian is generated in the protocol's beta
    # order: coordinator features first, then the remaining servers. The
    # public coefficients are reported in `server_list` order below. Align
    # the Hessian before covariance/SE and downstream Gram reconstruction.
    if (is.null(rownames(H_raw)) && isTRUE(use_secure_agg)) {
      hess_feature_order <- unlist(
        x_vars[c(coordinator, setdiff(server_list, coordinator))],
        use.names = FALSE)
      hess_order <- c("(Intercept)", hess_feature_order)
      if (length(hess_order) == nrow(H_raw)) {
        dimnames(H_raw) <- list(hess_order, hess_order)
      }
    }
    if (!is.null(rownames(H_raw))) {
      target_order <- names(all_coefs)
      perm <- match(target_order, rownames(H_raw))
      if (all(!is.na(perm)) && length(perm) == nrow(H_raw)) {
        H_raw <- H_raw[perm, perm, drop = FALSE]
      }
    }
    # Fisher = n x (Hessian - lambdaI) where Hessian = X_std^T W X_std / n + lambdaI
    H_adj <- H_raw - lambda * diag(nrow(H_raw))
    fisher_std <- n_obs * H_adj
    if (any(!is.finite(fisher_std))) {
      unavailable_certificate <- numeric_certificate
      unavailable_certificate$status <- "numeric_backend_unavailable"
      unavailable_certificate$reason <- "non_finite_information_matrix"
      .dsvert_stop_numeric(
        "numeric_backend_unavailable",
        "The protected information matrix contains a non-finite numeric result",
        unavailable_certificate, "non_finite_information_matrix")
    }
    identifiability_reason <- if (identical(family, "binomial") &&
                                    identical(as.numeric(lambda), 0)) {
      "separation_or_singular_information"
    } else {
      "singular_information_matrix"
    }
    cov_std <- .dsvert_solve_identifiable(
      fisher_std,
      context = paste0("The protected ", family, " model"),
      reason = identifiability_reason,
      symmetric = TRUE,
      certificate = numeric_certificate)
    numeric_certificate <- .dsvert_numeric_mark_estimable(
      numeric_certificate,
      reason = "positive_definite_full_rank_fisher_information")

    if (!is.null(cov_std)) {
      # Destandardize: construct Jacobian J where theta_orig = J x theta_std + const
      p_feat <- length(all_x_sds)
      J <- diag(p_feat + 1)  # (intercept + features)

      if (standardize_y && !is.null(y_sd)) {
        # Gaussian: beta_orig_j = beta_std_j x y_sd / x_sd_j
        for (jj in seq_len(p_feat)) {
          J[jj + 1, jj + 1] <- y_sd / all_x_sds[jj]
          J[1, jj + 1] <- -y_sd * all_x_means[jj] / all_x_sds[jj]
        }
        J[1, 1] <- y_sd
      } else {
        # Binomial/Poisson: beta_orig_j = beta_std_j / x_sd_j
        for (jj in seq_len(p_feat)) {
          J[jj + 1, jj + 1] <- 1.0 / all_x_sds[jj]
          J[1, jj + 1] <- -all_x_means[jj] / all_x_sds[jj]
        }
      }

      # Transform covariance to original space
      cov_orig <- J %*% cov_std %*% t(J)
      dimnames(cov_orig) <- list(names(all_coefs), names(all_coefs))

      # `cov_std` is the inverse Fisher matrix. For Gaussian models with
      # standardized y, transforming by J adds the y_sd^2 scale, but ordinary
      # glm() inference uses sigma_hat^2 * (X'WX)^-1. Keep the unscaled bread
      # available for sandwich users and scale the public covariance by the
      # residual variance estimate when deviance is available.
      cov_info_orig <- cov_orig
      if (family == "gaussian") {
        y_scale <- if (standardize_y && !is.null(y_sd) &&
                        is.finite(y_sd) && y_sd > 0) {
          y_sd^2
        } else {
          1
        }
        cov_info_orig <- cov_orig / y_scale
        if (!is.na(deviance) && is.finite(deviance)) {
          covariance_df_residual <- max(df_residual, 1L)
          sigma2_hat <- deviance / covariance_df_residual
          cov_orig <- cov_info_orig * sigma2_hat
          dimnames(cov_orig) <- list(names(all_coefs), names(all_coefs))
        }
      }
      dimnames(cov_info_orig) <- list(names(all_coefs), names(all_coefs))

      se_orig <- sqrt(pmax(diag(cov_orig), 0))
      std_errors <- se_orig
      names(std_errors) <- names(all_coefs)

      statistic_values <- all_coefs / std_errors
      statistic_values[!is.finite(statistic_values)] <- NA
      z_values <- statistic_values  # historical field; see coefficient_reference
      if (identical(family, "gaussian")) {
        t_values <- statistic_values
        p_values <- if (is.finite(df_residual) && df_residual > 0L) {
          2 * stats::pt(-abs(statistic_values), df = df_residual)
        } else {
          rep(NA_real_, length(statistic_values))
        }
      } else {
        p_values <- 2 * stats::pnorm(-abs(statistic_values))
      }
      names(statistic_values) <- names(z_values) <- names(p_values) <-
        names(all_coefs)
      if (!is.null(t_values)) names(t_values) <- names(all_coefs)
    }
    # Expose the full covariance matrix (original scale) so downstream
    # client-side inference (multi-coef Wald, GEE sandwich, CI on linear
    # contrasts) can reuse it without another MPC round.
    if (exists("cov_orig", inherits = FALSE)) covariance <- cov_orig
    if (exists("cov_info_orig", inherits = FALSE)) {
      covariance_information <- cov_info_orig
    }

    if (verbose && any(!is.na(std_errors))) {
      statistic_label <- if (identical(family, "gaussian")) "t value" else "z value"
      p_label <- if (identical(family, "gaussian")) "Pr(>|t|)" else "Pr(>|z|)"
      message("\nCoefficients:")
      message(sprintf("  %-15s %10s %10s %10s %10s", "", "Estimate",
                      "Std.Error", statistic_label, p_label))
      for (nm in names(all_coefs)) {
        sig <- if (!is.na(p_values[nm]) && p_values[nm] < 0.001) "***"
               else if (!is.na(p_values[nm]) && p_values[nm] < 0.01) "**"
               else if (!is.na(p_values[nm]) && p_values[nm] < 0.05) "*"
               else ""
        message(sprintf("  %-15s %10.4f %10.4f %10.3f %10.6f %s",
          nm, all_coefs[nm], std_errors[nm], statistic_values[nm],
          p_values[nm], sig))
      }
    }
  }

  gradient_std <- NULL
  if (use_secure_agg && exists("k3_result") &&
      !is.null(k3_result$gradient_std)) {
    gradient_std <- k3_result$gradient_std
  }
  if (use_k2_beaver && exists("loop_result") &&
      !is.null(loop_result$gradient_std)) {
    gradient_std <- loop_result$gradient_std
  }
  gradient_original <- NULL
  if (!is.null(gradient_std)) {
    target_order <- names(all_coefs)
    if (!is.null(names(gradient_std))) {
      perm <- match(target_order, names(gradient_std))
      if (all(!is.na(perm))) gradient_std <- gradient_std[perm]
    }
    names(gradient_std) <- target_order
    if (!standardize_y) {
      gradient_original <- gradient_std
      slope_names <- setdiff(target_order, "(Intercept)")
      if ("(Intercept)" %in% target_order && length(slope_names) > 0L) {
        g0 <- gradient_std["(Intercept)"]
        for (nm in slope_names) {
          gradient_original[nm] <- all_x_means[match(nm, all_names)] * g0 +
            all_x_sds[match(nm, all_names)] * gradient_std[nm]
        }
      }
      names(gradient_original) <- target_order
    }
  }

  # ===========================================================================
  # Assemble Result
  # ===========================================================================
  numeric_certificate <- .dsvert_numeric_finalize_certificate(
    numeric_certificate,
    coefficients = all_coefs,
    converged = converged,
    returned_numeric = list(
      deviance = deviance,
      covariance = covariance,
      covariance_information = covariance_information,
      std_errors = std_errors,
      gradient = gradient_original))
  result <- list(
    coefficients = all_coefs,
    std_errors = std_errors,
    statistic_values = statistic_values,
    t_values = t_values,
    z_values = z_values,
    p_values = p_values,
    coefficient_reference = if (identical(family, "gaussian")) "t" else "normal",
    covariance = covariance,
    covariance_information = covariance_information,
    covariance_unscaled = covariance_information,
    iterations = final_iter,
    converged = converged,
    family = family,
    binomial_sigmoid_intervals = effective_binomial_sigmoid_intervals,
    n_obs = n_obs,
    n_vars = n_vars_total,
    n_parameters = n_parameters,
    df_residual = df_residual,
    lambda = lambda,
    estimand = if (isTRUE(as.numeric(lambda) > 0)) {
      paste0("explicit_ridge_penalized_", family, "_glm")
    } else {
      paste0("unpenalized_", family, "_glm")
    },
    deviance = deviance,
    deviance_type = deviance_type,
    null_deviance = null_deviance,
    pseudo_r2 = pseudo_r2,
    aic = aic,
    aic_type = aic_type,
    y_server = y_server,
    eta_privacy = eta_privacy,
    missing_policy = missing_policy,
    data_name = data_name,
    y_var = y_var,
    x_vars = x_vars,
    offset = offset,
    weights = weights,
    ring = effective_ring,
    requested_ring = requested_ring,
    effective_ring = effective_ring,
    numeric_certificate = numeric_certificate,
    alignment_attested = alignment_metadata$alignment_attested,
    alignment_manifest_hash = alignment_metadata$alignment_manifest_hash,
    cohort_id = alignment_metadata$cohort_id,
    x_means = setNames(all_x_means, all_names),
    x_sds   = setNames(all_x_sds,   all_names),
    y_sd    = if (exists("y_sd", inherits = FALSE)) y_sd else NULL,
    y_mean  = if (exists("y_mean", inherits = FALSE)) y_mean else NULL,
    hessian_std = if (exists("H_raw", inherits = FALSE)) H_raw else NULL,
    gradient_std = gradient_std,
    gradient = gradient_original,
    call = call_matched
  )

  if (isTRUE(keep_session)) {
    result$session_id <- session_id
    result$transport_pks <- transport_pks
    result$server_list <- server_list
    result$x_vars <- x_vars
    result$y_var <- y_var
    result$data_name <- data_name
    result$std_data <- std_data
    result$standardize_y <- standardize_y
    result$ring <- ring
  }

  # Cleanup handled by on.exit()

  class(result) <- c("ds.glm", "list")
  return(result)
}

#' @title Print Method for ds.glm Objects
#' @description Prints a summary of GLM results.
#' @param x A ds.glm object
#' @param ... Additional arguments (ignored)
#' @export
print.ds.glm <- function(x, ...) {
  cat("\nVertically Partitioned GLM (Block Coordinate Descent)\n")
  cat("=======================================================\n\n")

  cat("Call:\n")
  print(x$call)
  cat("\n")

  cat("Family:", x$family, "\n")
  cat("Observations:", x$n_obs, "\n")
  cat("Predictors:", x$n_vars, "\n")
  cat("Regularization (lambda):", x$lambda, "\n")
  if (!is.null(x$estimand)) cat("Estimand:", x$estimand, "\n")
  if (!is.null(x$effective_ring)) {
    cat("MPC ring: requested", x$requested_ring %||% x$effective_ring,
        ", effective", x$effective_ring, "\n")
  }
  if (!is.null(x$numeric_certificate)) {
    cat("Numeric preflight:", x$numeric_certificate$status,
        "[", x$numeric_certificate$effective_backend, "]",
        "| execution certified:",
        if (isTRUE(x$numeric_certificate$numerically_certified)) "yes" else "no",
        "\n")
  }
  if (!is.null(x$y_server))
    cat("Label server:", x$y_server, "\n")
  cat("Iterations:", x$iterations, "\n")
  cat("Converged:", x$converged, "\n\n")

  cat("Coefficients:\n")
  print(round(x$coefficients, 6))

  invisible(x)
}

#' @title Summary Method for ds.glm Objects
#' @description Prints detailed summary including deviance and fit statistics.
#' @param object A ds.glm object
#' @param ... Additional arguments (ignored)
#' @export
summary.ds.glm <- function(object, ...) {
  cat("\nVertically Partitioned GLM - Summary\n")
  cat("====================================\n\n")

  cat("Call:\n")
  print(object$call)
  cat("\n")

  cat("Family:", object$family, "\n")
  cat("Observations:", object$n_obs, "\n")
  cat("Predictors:", object$n_vars, "\n")
  cat("Regularization (lambda):", object$lambda, "\n")
  if (!is.null(object$effective_ring)) {
    cat("MPC ring: requested", object$requested_ring %||% object$effective_ring,
        ", effective", object$effective_ring, "\n")
  }
  if (!is.null(object$numeric_certificate)) {
    cat("Numeric preflight:", object$numeric_certificate$status,
        "[", object$numeric_certificate$effective_backend, "]",
        "| execution certified:",
        if (isTRUE(object$numeric_certificate$numerically_certified)) "yes" else "no",
        "\n")
    if (is.finite(object$numeric_certificate$total_numeric_error_max)) {
      cat("Certified per-path arithmetic error bound:",
          format(object$numeric_certificate$total_numeric_error_max,
                 scientific = TRUE), "\n")
    }
  }
  if (!is.null(object$y_server))
    cat("Label server:", object$y_server, "\n")
  cat("\n")

  cat("Convergence:\n")
  cat("  Iterations:", object$iterations, "\n")
  cat("  Converged:", object$converged, "\n\n")

  cat("Deviance:\n")
  if (length(object$null_deviance) == 1L &&
      is.finite(object$null_deviance)) {
    cat("  Null deviance:    ", sprintf("%.4f", object$null_deviance),
        " on", object$n_obs - 1, "degrees of freedom\n")
  } else {
    cat("  Null deviance:     unavailable for this fit contract\n")
  }
  residual_df <- object$df_residual %||% (object$n_obs - object$n_vars)
  cat("  Residual deviance:", sprintf("%.4f", object$deviance),
      " on", residual_df, "degrees of freedom\n\n")

  cat("Model Fit:\n")
  r2_label <- if (identical(object$family, "gaussian")) {
    "R-squared"
  } else {
    "Pseudo R-squared"
  }
  cat(" ", r2_label, ":", sprintf("%.4f", object$pseudo_r2), "\n")
  cat("  AIC:", sprintf("%.4f", object$aic),
      "[", object$aic_type %||% "unspecified", "]\n\n")

  cat("Coefficients:\n")
  coef_df <- data.frame(Estimate = object$coefficients)
  print(round(coef_df, 6))

  invisible(object)
}

#' @title Coefficients Method for ds.glm Objects
#' @description Extract coefficients from a ds.glm object.
#' @param object A ds.glm object
#' @param ... Additional arguments (ignored)
#' @return Named numeric vector of coefficients
#' @export
coef.ds.glm <- function(object, ...) {
  object$coefficients
}
