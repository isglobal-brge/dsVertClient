#' @title Generalized Linear Model for Vertically Partitioned Data
#' @description Fits a GLM across vertically partitioned data using Ring63
#'   Beaver MPC with DCF wide spline for the link function. The system
#'   auto-detects which server holds each variable. Only p-dimensional
#'   aggregate gradients are revealed to the client per iteration.
#'   No observation-level data is ever disclosed.
#'
#' @param formula A formula (e.g. \code{npreg ~ age + bmi + glu}) or
#'   character string (\code{"npreg ~ age + bmi + glu"}). Can also be
#'   a data_name string for backward compatibility.
#' @param data Character string. Name of the (aligned) data frame on
#'   each server.
#' @param x_vars Optional. Character vector of predictor names, or a
#'   named list mapping server names to variable vectors. If NULL and
#'   formula is used, extracted from the formula. If NULL and no formula,
#'   all available columns (minus y and IDs) are used.
#' @param y_server Character string. Name of the server holding the response.
#' @param family Character string. GLM family: "gaussian", "binomial",
#'   or "poisson". Default is "gaussian".
#' @param max_iter Integer. Maximum L-BFGS iterations. Default is 100.
#' @param tol Numeric. Convergence tolerance on coefficient change.
#'   Default is 1e-4.
#' @param lambda Numeric. L2 regularization parameter. Default is 1e-4.
#' @param log_n Integer. Legacy parameter (ignored). Kept for backward
#'   compatibility.
#' @param verbose Logical. Print progress messages. Default is TRUE.
#' @param datasources DataSHIELD connection object or list of connections.
#' @param eta_privacy Character. \code{"auto"} (default) selects
#'   \code{"k2_beaver"} for K=2 or \code{"secure_agg"} for K>=3.
#'
#' @return A list with class "ds.glm" containing:
#'   \itemize{
#'     \item \code{coefficients}: Named coefficient vector (original scale)
#'     \item \code{std_errors}: Standard errors (finite-difference Hessian)
#'     \item \code{covariance}: Model covariance matrix on the original scale
#'     \item \code{covariance_information}: Inverse Fisher/bread matrix
#'       without Gaussian residual-variance scaling, for sandwich methods
#'     \item \code{z_values}: z-statistics (coef / SE)
#'     \item \code{p_values}: Two-sided p-values
#'     \item \code{iterations}: Number of iterations
#'     \item \code{converged}: Logical
#'     \item \code{family}: Family used
#'     \item \code{n_obs}: Number of observations
#'     \item \code{deviance}: Residual sum of squares
#'     \item \code{pseudo_r2}: 1 - deviance/null_deviance
#'   }
#'
#' @details
#' \subsection{Protocol}{
#' All computation uses Ring63 fixed-point arithmetic with Beaver MPC:
#' \enumerate{
#'   \item \strong{Transport keys}: X25519 keypairs on all servers (~0.5s)
#'   \item \strong{Standardize}: Each server standardizes its features
#'   \item \strong{Input sharing}: Features split into additive Ring63 shares
#'     between 2 DCF parties. Non-DCF servers contribute shares.
#'   \item \strong{L-BFGS loop}: Per iteration:
#'     \itemize{
#'       \item Compute eta shares (Ring63 matrix-vector)
#'       \item DCF wide spline for sigmoid/exp (binomial/Poisson) or
#'         identity link (Gaussian)
#'       \item Beaver matvec for gradient (server-generated triples)
#'       \item Client aggregates Ring63 shares -> p gradient scalars
#'       \item L-BFGS quasi-Newton update
#'     }
#'   \item \strong{SE}: p+1 gradient evaluations (finite-difference Hessian)
#'   \item \strong{Deviance}: Beaver dot-product for residual sum of squares
#'   \item \strong{Unstandardize}: Coefficients + SE via Jacobian transform
#' }
#' }
#'
#' \subsection{Security}{
#' No observation-level data is disclosed. The client sees only p-dimensional
#' aggregate gradients per iteration. Beaver triples are generated server-side
#' (never seen by client). Dealer rotation for K>=4 ensures the analyst must
#' compromise (K-1)/K servers to extract data.
#' }
#'
#' @references
#' van Kesteren, E.J. et al. (2019). Privacy-preserving generalized linear
#' models using distributed block coordinate descent. arXiv:1911.03183.
#'
#' Mouchet, C. et al. (2021). "Multiparty Homomorphic Encryption from
#' Ring-Learning-With-Errors". \emph{Proceedings on Privacy Enhancing Technologies}.
#'
#' @seealso \code{\link{ds.vertCor}} for correlation analysis,
#'   \code{\link{ds.vertPCA}} for PCA analysis
#'
#' @examples
#' \dontrun{
#' # Simplest: formula interface (auto-detects everything)
#' model <- ds.vertGLM(npreg ~ age + bmi + glu + bp + skin,
#'                      data = "DA", family = "gaussian")
#'
#' # String formula
#' model <- ds.vertGLM("diabetes ~ age + bmi + glu",
#'                      data = "DA", family = "binomial")
#'
#' # Auto-detect all features
#' model <- ds.vertGLM("DA", "npreg", family = "poisson")
#'
#' # Manual server mapping (legacy)
#' model <- ds.vertGLM("DA", "npreg",
#'   list(s1 = c("age", "bmi"), s2 = c("glu", "bp")),
#'   y_server = "s2", family = "gaussian")
#' }
#'
#' @param data_name Internal alias for \code{data}; if both are supplied,
#'   \code{data_name} takes precedence (back-compat path).
#' @param y_var Internal alias for the LHS of \code{formula}; if both
#'   are supplied, \code{y_var} takes precedence (back-compat path).
#' @param offset Optional numeric vector or column name on the outcome
#'   server added to the linear predictor (e.g. \code{log(person_years)}
#'   for a Poisson rate model).
#' @param weights Optional numeric vector or column name of per-row
#'   weights (e.g. inverse-probability weights for IPW).
#' @param ring Integer (63 or 127). Selects the MPC ring / fracBits
#'   pipeline; Ring127 (fracBits=50) is STRICT-capable per
#'   Catrina-Saxena. Default 63L for back-compat.
#' @param keep_session Logical. If TRUE, leave the MPC session alive
#'   on the servers and expose \code{session_id}, \code{transport_pks},
#'   and \code{server_list} on the returned fit so follow-on helpers
#'   (LMM cluster residuals, GEE sandwich meat, etc.) can reuse the
#'   already-aligned shares. Caller must invoke \code{mpcCleanupDS}
#'   eventually.
#' @param no_intercept Logical. Suppress the auto-added intercept.
#'   Useful when the design matrix already encodes one (e.g. cluster-
#'   mean-centred GLS fit).
#' @param std_mode Character. Standardisation mode: \code{"full"}
#'   (default) standardises both X and y; alternative modes (e.g.
#'   \code{"x_only"}) skip y standardisation for offset/weights paths.
#' @param start Optional named coefficient vector for internal K>=3 fixed
#'   evaluation paths. When supplied with \code{max_iter = 0}, the secure
#'   loop evaluates residual/deviance shares at the supplied coefficients
#'   without optimisation.
#' @param compute_se Logical. Compute finite-difference Hessian/standard
#'   errors. Internal fixed-evaluation callers can set FALSE to avoid
#'   extra MPC rounds when only residual shares are needed.
#' @param compute_deviance Logical. Compute the secure aggregate deviance.
#'   Internal score-only callers can set FALSE to avoid the extra MPC pass.
#' @param gradient_only Logical. Internal diagnostic/optimizer path. If TRUE,
#'   evaluate and return the aggregate score at \code{start} without taking
#'   an optimizer step. If \code{compute_se = TRUE}, the aggregate finite-
#'   difference Hessian is also returned. The returned score/Hessian are
#'   low-dimensional aggregates only; no observation-level eta/probability/
#'   residual vector is opened.
#' @importFrom DSI datashield.aggregate datashield.connections_find
#' @export
ds.vertGLM <- function(formula, data = NULL, x_vars = NULL, y_server = NULL,
                       family = "gaussian", max_iter = 100, tol = 1e-4,
                       lambda = 1e-4, log_n = 12,
                       offset = NULL, weights = NULL,
                       # Ring63 (frac_bits=20, default, back-compat) or
                       # Ring127 (frac_bits=50, STRICT-capable per
                       # Catrina-Saxena 2^-fracbits scaling). Ring127
                       # routes through the Uint128 Go primitives that
                       # already exist for task #116 Cox/LMM. Used by
                       # IPW/#98 for STRICT closure; other families may
                       # opt in as the Ring127 regression suite expands.
                       ring = 63L,
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
                       data_name = NULL, y_var = NULL) {
  call_matched <- match.call()

  # ===========================================================================
  # Parse formula or legacy arguments
  # ===========================================================================
  if (!missing(formula)) {
    if (inherits(formula, "formula")) {
      # R formula object: npreg ~ age + bmi + ped
      f_terms <- terms(formula)
      y_var <- as.character(attr(f_terms, "variables")[[2]])
      x_from_formula <- attr(f_terms, "term.labels")
      if (is.null(x_vars)) x_vars <- x_from_formula
      data_name <- data
    } else if (is.character(formula) && grepl("~", formula)) {
      # String formula: "npreg ~ age + bmi + ped"
      f <- as.formula(formula)
      f_terms <- terms(f)
      y_var <- as.character(attr(f_terms, "variables")[[2]])
      x_from_formula <- attr(f_terms, "term.labels")
      if (is.null(x_vars)) x_vars <- x_from_formula
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
  if (!family %in% c("gaussian", "binomial", "poisson"))
    stop("family must be 'gaussian', 'binomial', or 'poisson'",
         call. = FALSE)

  if (is.null(datasources))
    datasources <- DSI::datashield.connections_find()
  if (length(datasources) == 0)
    stop("No DataSHIELD connections found", call. = FALSE)

  # Auto-detect: query servers for their columns and map variables automatically
  if (is.null(x_vars) || is.character(x_vars)) {
    user_x_vars <- x_vars  # NULL = use all available, character = specific vars
    if (verbose) message("[Auto-detect] Querying server columns...")
    col_results <- DSI::datashield.aggregate(datasources,
      call(name = "dsvertColNamesDS", data_name = data_name))
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
      missing <- setdiff(user_x_vars, all_available)
      if (length(missing) > 0)
        stop("Variables not found on any server: ", paste(missing, collapse = ", "),
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

  n_partitions_check <- length(x_vars)
  non_label_count <- n_partitions_check - 1

  # Route to the appropriate protocol:
  #   K=2: Beaver MPC with wide spline DCF (non-disclosive, ~1e-3 precision)
  #   K>=3: Secure aggregation with Enc(r) + L-BFGS (~8e-3 precision)
  if (eta_privacy == "auto") {
    if (non_label_count >= 2) {
      eta_privacy <- "secure_agg"
    } else if (non_label_count >= 1) {
      eta_privacy <- "k2_beaver"
    } else {
      stop("Need at least 2 servers (1 label + 1 non-label)", call. = FALSE)
    }
  }

  use_secure_agg <- (eta_privacy == "secure_agg")
  use_k2_beaver <- (eta_privacy == "k2_beaver")

  if (use_secure_agg && non_label_count < 2)
    stop("secure_agg requires >= 3 servers", call. = FALSE)
  if (use_k2_beaver && non_label_count != 1)
    stop("K=2 mode requires exactly 2 servers", call. = FALSE)
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
    missing <- setdiff(names(x_vars), server_names)
    stop("Unknown server(s) in x_vars: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }

  server_list <- names(x_vars)
  non_label_servers <- setdiff(server_list, y_server)
  n_partitions <- length(x_vars)

  # Get observation count (lightweight sync call before helpers are available)
  first_conn <- which(server_names == server_list[1])
  count_result <- DSI::datashield.aggregate(
    conns = datasources[first_conn],
    expr = call(name = "getObsCountDS", data_name)
  )
  if (is.list(count_result) && length(count_result) == 1)
    count_result <- count_result[[1]]
  n_obs <- count_result$n_obs

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
  session_id <- local({
    hex <- sample(c(0:9, letters[1:6]), 32, replace = TRUE)
    hex[13] <- "4"  # UUID v4
    hex[17] <- sample(c("8","9","a","b"), 1)  # variant 1
    paste0(
      paste(hex[1:8], collapse = ""), "-",
      paste(hex[9:12], collapse = ""), "-",
      paste(hex[13:16], collapse = ""), "-",
      paste(hex[17:20], collapse = ""), "-",
      paste(hex[21:32], collapse = "")
    )
  })

  # Guaranteed cleanup on exit (even if error occurs mid-protocol).
  # Uses DSI::datashield.aggregate directly (not .dsAgg) so cleanup works
  # even if .glm_mpc_setup() fails before returning the closure.
  # Skipped when keep_session = TRUE so follow-on helpers can reuse the
  # MPC session (caller assumes cleanup responsibility).
  on.exit({
    if (!isTRUE(keep_session)) {
      for (.srv in server_list) {
        .ci <- which(server_names == .srv)
        tryCatch(
          DSI::datashield.aggregate(conns = datasources[.ci],
            expr = call(name = "mpcCleanupDS", session_id = session_id)),
          error = function(e) NULL)
        tryCatch(
          DSI::datashield.aggregate(conns = datasources[.ci],
            expr = call(name = "mpcGcDS")),
          error = function(e) NULL)
      }
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
    std_mode = std_mode
  )

  # Unpack setup results
  transport_pks <- setup$transport_pks
  x_means       <- setup$x_means
  x_sds         <- setup$x_sds
  y_mean        <- setup$y_mean
  y_sd          <- setup$y_sd
  std_data      <- setup$std_data
  standardize_y <- setup$standardize_y
  .dsAgg        <- setup$.dsAgg
  .sendBlob     <- setup$.sendBlob

  if (max_iter < 1L && !isTRUE(compute_se) &&
      !isTRUE(compute_deviance) && !isTRUE(gradient_only) &&
      is.null(start) && is.null(offset) && is.null(weights)) {
    all_names <- unlist(x_vars[server_list], use.names = FALSE)
    all_x_means <- unlist(x_means[server_list], use.names = FALSE)
    all_x_sds <- unlist(x_sds[server_list], use.names = FALSE)
    all_coefs <- c("(Intercept)" = 0, setNames(rep(0, length(all_names)),
                                               all_names))
    result <- list(
      coefficients = all_coefs,
      std_errors = rep(NA_real_, length(all_coefs)),
      z_values = rep(NA_real_, length(all_coefs)),
      p_values = rep(NA_real_, length(all_coefs)),
      covariance = NULL,
      covariance_information = NULL,
      covariance_unscaled = NULL,
      iterations = 0L,
      converged = TRUE,
      family = family,
      n_obs = n_obs,
      n_vars = length(all_coefs),
      lambda = lambda,
      deviance = NA_real_,
      null_deviance = if (family == "gaussian") (n_obs - 1) * (y_sd %||% 1)^2 else NA,
      pseudo_r2 = NA_real_,
      aic = NA_real_,
      y_server = y_server,
      eta_privacy = eta_privacy,
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
      .ci <- which(server_names == .srv)
      cols <- tryCatch(
        DSI::datashield.aggregate(datasources[.ci],
          call(name = "dsvertColNamesDS", data_name = data_name))[[1]]$columns,
        error = function(e) NULL)
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
    .dsAgg(datasources[.ci], call(name = "k2SetOffsetDS",
      data_name = data_name,
      offset_column = offset,
      session_id = session_id))
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
      .ci <- which(server_names == .srv)
      cols <- tryCatch(
        DSI::datashield.aggregate(datasources[.ci],
          call(name = "dsvertColNamesDS", data_name = data_name))[[1]]$columns,
        error = function(e) NULL)
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
        session_id = session_id))
      if (is.list(setres) && length(setres) == 1L) setres <- setres[[1]]

      if (dcf_role == "dealer") {
        .sendBlob(setres$dcf0_blob, "k2_peer_weight_share",
                  which(server_names == dcf_weight_parties[1L]))
        .sendBlob(setres$dcf1_blob, "k2_peer_weight_share",
                  which(server_names == dcf_weight_parties[2L]))
        .sendBlob(setres$dcf0_sqrt_blob, "k2_peer_sqrt_weight_share",
                  which(server_names == dcf_weight_parties[1L]))
        .sendBlob(setres$dcf1_sqrt_blob, "k2_peer_sqrt_weight_share",
                  which(server_names == dcf_weight_parties[2L]))
        for (srv in dcf_weight_parties) {
          .dsAgg(datasources[which(server_names == srv)],
            call(name = "k2ReceiveWeightSharesDS", session_id = session_id))
        }
      } else {
        peer_srv <- if (dcf_role == "dcf0") dcf_weight_parties[2L] else dcf_weight_parties[1L]
        peer_ci <- which(server_names == peer_srv)
        .sendBlob(setres$peer_blob, "k2_peer_weight_share", peer_ci)
        .sendBlob(setres$peer_sqrt_blob, "k2_peer_sqrt_weight_share", peer_ci)
        .dsAgg(datasources[peer_ci], call(name = "k2ReceiveWeightSharesDS",
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
        session_id = session_id))
      if (is.list(setres) && length(setres) == 1L) setres <- setres[[1]]
      .sendBlob(setres$peer_blob, "k2_peer_weight_share", peer_ci)
      .sendBlob(setres$peer_sqrt_blob, "k2_peer_sqrt_weight_share", peer_ci)
      .dsAgg(datasources[peer_ci], call(name = "k2ReceiveWeightSharesDS",
        session_id = session_id))
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
      no_intercept  = isTRUE(no_intercept),
      ring = ring,
      compute_se = isTRUE(compute_se),
      compute_deviance = isTRUE(compute_deviance),
      gradient_only = isTRUE(gradient_only),
      start = start
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
      ring = ring)
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

  # Deviance: from secure Beaver Sumr^2 (available in both K=2 and K>=3)
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
  # Null deviance: for Gaussian = Sum(y-ybar)^2 = (n-1)*var(y). In std space, var(y_std)=1.
  null_deviance <- if (family == "gaussian") (n_obs - 1) * (y_sd %||% 1)^2 else NA
  pseudo_r2 <- if (!is.na(null_deviance) && null_deviance > 0) 1 - (deviance / null_deviance) else NA
  aic <- if (!is.na(deviance)) deviance + 2 * n_vars_total else NA

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

  std_errors <- rep(NA, n_vars_total)
  z_values <- rep(NA, n_vars_total)
  p_values <- rep(NA, n_vars_total)
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
    cov_std <- tryCatch(solve(fisher_std), error = function(e) NULL)

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
          df_resid <- max(n_obs - n_vars_total, 1L)
          sigma2_hat <- deviance / df_resid
          cov_orig <- cov_info_orig * sigma2_hat
          dimnames(cov_orig) <- list(names(all_coefs), names(all_coefs))
        }
      }
      dimnames(cov_info_orig) <- list(names(all_coefs), names(all_coefs))

      se_orig <- sqrt(pmax(diag(cov_orig), 0))
      std_errors <- se_orig
      names(std_errors) <- names(all_coefs)

      z_values <- all_coefs / std_errors
      z_values[!is.finite(z_values)] <- NA
      p_values <- 2 * stats::pnorm(-abs(z_values))
      names(z_values) <- names(p_values) <- names(all_coefs)
    }
    # Expose the full covariance matrix (original scale) so downstream
    # client-side inference (multi-coef Wald, GEE sandwich, CI on linear
    # contrasts) can reuse it without another MPC round.
    if (exists("cov_orig", inherits = FALSE)) covariance <- cov_orig
    if (exists("cov_info_orig", inherits = FALSE)) {
      covariance_information <- cov_info_orig
    }

    if (verbose && any(!is.na(std_errors))) {
      message("\nCoefficients:")
      message(sprintf("  %-15s %10s %10s %10s %10s", "", "Estimate", "Std.Error", "z value", "Pr(>|z|)"))
      for (nm in names(all_coefs)) {
        sig <- if (!is.na(p_values[nm]) && p_values[nm] < 0.001) "***"
               else if (!is.na(p_values[nm]) && p_values[nm] < 0.01) "**"
               else if (!is.na(p_values[nm]) && p_values[nm] < 0.05) "*"
               else ""
        message(sprintf("  %-15s %10.4f %10.4f %10.3f %10.6f %s",
          nm, all_coefs[nm], std_errors[nm], z_values[nm], p_values[nm], sig))
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
  result <- list(
    coefficients = all_coefs,
    std_errors = std_errors,
    z_values = z_values,
    p_values = p_values,
    covariance = covariance,
    covariance_information = covariance_information,
    covariance_unscaled = covariance_information,
    iterations = final_iter,
    converged = converged,
    family = family,
    n_obs = n_obs,
    n_vars = n_vars_total,
    lambda = lambda,
    deviance = deviance,
    null_deviance = null_deviance,
    pseudo_r2 = pseudo_r2,
    aic = aic,
    y_server = y_server,
    eta_privacy = eta_privacy,
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
  if (!is.null(object$y_server))
    cat("Label server:", object$y_server, "\n")
  cat("\n")

  cat("Convergence:\n")
  cat("  Iterations:", object$iterations, "\n")
  cat("  Converged:", object$converged, "\n\n")

  cat("Deviance:\n")
  cat("  Null deviance:    ", sprintf("%.4f", object$null_deviance),
      " on", object$n_obs - 1, "degrees of freedom\n")
  cat("  Residual deviance:", sprintf("%.4f", object$deviance),
      " on", object$n_obs - object$n_vars, "degrees of freedom\n\n")

  cat("Model Fit:\n")
  cat("  Pseudo R-squared (McFadden):", sprintf("%.4f", object$pseudo_r2), "\n")
  cat("  AIC:", sprintf("%.4f", object$aic), "\n\n")

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
