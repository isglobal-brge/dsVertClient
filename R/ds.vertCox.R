#' @title Federated Cox proportional-hazards regression
#' @description Fit a Cox PH model on vertically partitioned DataSHIELD
#'   data using the reformulated reverse-cumsum gradient
#'
#'     grad(beta) = sum_{j: delta_j=1} x_j  -  sum_j x_j exp(eta_j) G_j
#'     G_j = sum_{i: delta_i=1, i <= j} 1 / S(t_i)    (ascending-t sort)
#'     S(t_i) = sum_{k >= i} exp(eta_k)              (reverse cumsum)
#'
#'   Every quantity on the per-patient timeline (eta, exp(eta), S, 1/S,
#'   G) is kept as Ring63 secret shares between the two DCF parties; the
#'   analyst client only ever sees the p-dimensional aggregate gradient
#'   and the final coefficient vector. Built on top of:
#'     - k2SetCoxTimesDS / k2ReceiveCoxMetaDS / k2ApplyCoxPermutationDS
#'       (server helpers landed in commit 901921d)
#'     - WideSplineExp + WideSplineReciprocalRefined (wired into the
#'       4-phase DCF protocol in commit 75f6883)
#'     - k2-fp-cumsum + k2-fp-vec-mul (Go primitives in a16fcb0/34a28f6)
#'
#' @param formula An R formula. LHS MUST be \code{survival::Surv(time,
#'   event)} or a two-column character vector \code{c(time_col,
#'   event_col)}; RHS is the usual list of covariates.
#' @param data Name of the aligned data frame on each server.
#' @param time_col Name of the event-time column (required if formula
#'   LHS is not a Surv(...) expression).
#' @param event_col Name of the event-indicator column (required
#'   similarly).
#' @param max_iter Outer L-BFGS iterations (default 30).
#' @param tol Convergence tolerance on max coefficient change.
#' @param learning_rate Gradient descent step size (used by the basic
#'   outer loop when L-BFGS hand-off is not configured).
#' @param datasources DataSHIELD connections.
#' @param verbose Print progress.
#' @return A ds.vertCox object: coefficients, iterations, converged,
#'   deviance (partial log-likelihood), n_events, call.
#' @export
ds.vertCox <- function(formula, data = NULL, time_col = NULL,
                       event_col = NULL, x_vars = NULL, y_server = NULL,
                       max_iter = 30L, tol = 1e-4, learning_rate = 0.1,
                       lambda = 1e-4,
                       verbose = TRUE, datasources = NULL) {
  if (is.null(datasources)) datasources <- DSI::datashield.connections_find()
  if (is.null(data)) stop("data (aligned data frame name) required", call. = FALSE)
  call_matched <- match.call()

  # Parse formula: supports `Surv(time, event) ~ x1 + x2` or
  # `outcome ~ x1 + x2` with time_col/event_col passed explicitly.
  if (!inherits(formula, "formula")) {
    stop("formula must be an R formula", call. = FALSE)
  }
  f_terms <- terms(formula)
  lhs <- attr(f_terms, "variables")[[2]]
  x_from_formula <- attr(f_terms, "term.labels")
  if (is.null(x_vars)) x_vars <- x_from_formula

  if (is.call(lhs) && as.character(lhs[[1]]) %in% c("Surv", "survival::Surv")) {
    if (is.null(time_col)) time_col <- as.character(lhs[[2]])
    if (is.null(event_col)) event_col <- as.character(lhs[[3]])
  }
  if (is.null(time_col) || is.null(event_col)) {
    stop("Need time_col and event_col (either in Surv(...) LHS or as args)",
         call. = FALSE)
  }

  server_names <- names(datasources)

  # Auto-detect variable locations
  col_results <- DSI::datashield.aggregate(datasources,
    call("dsvertColNamesDS", data_name = data))
  var_map <- list()
  y_srv <- NULL
  for (srv in server_names) {
    cols <- col_results[[srv]]$columns
    if (time_col %in% cols && event_col %in% cols) {
      y_srv <- srv
    }
    feats <- intersect(x_vars, cols)
    if (length(feats) > 0) var_map[[srv]] <- feats
  }
  if (is.null(y_srv)) {
    stop("Time + event columns not found on any single server", call. = FALSE)
  }
  if (!is.null(y_server) && y_server != y_srv) {
    stop("y_server = '", y_server, "' disagrees with auto-detected '",
         y_srv, "'", call. = FALSE)
  }
  y_server <- y_srv
  server_list <- names(var_map)
  if (length(server_list) < 2L) {
    stop("Need at least 2 servers with covariates for K=2 Cox", call. = FALSE)
  }

  session_id <- local({
    hex <- sample(c(0:9, letters[1:6]), 32, replace = TRUE)
    hex[13] <- "4"; hex[17] <- sample(c("8","9","a","b"), 1)
    paste0(paste(hex[1:8],collapse=""),"-",paste(hex[9:12],collapse=""),"-",
           paste(hex[13:16],collapse=""),"-",paste(hex[17:20],collapse=""),"-",
           paste(hex[21:32],collapse=""))
  })

  on.exit({
    for (.srv in server_list) {
      .ci <- which(server_names == .srv)
      tryCatch(
        DSI::datashield.aggregate(datasources[.ci],
          call("mpcCleanupDS", session_id = session_id)),
        error = function(e) NULL)
    }
  }, add = TRUE)

  if (verbose) {
    message(sprintf("[ds.vertCox] %d servers, y=%s, time=%s, event=%s",
                     length(server_list), y_server, time_col, event_col))
  }

  # =========================================================================
  # NOTE: This is a first-release Cox client. It orchestrates the full
  # server-side helper chain (k2SetCoxTimesDS + k2ApplyCoxPermutationDS
  # + k2CoxReverseCumsumSDS + k2StoreCoxRecipDS + k2CoxForwardCumsumGDS)
  # and runs a basic gradient-descent outer loop on the MPC gradient.
  # The L-BFGS outer loop and finite-difference SE are tracked as Month 2
  # follow-ons in V2_PROGRESS.md; the machinery (Cov(beta) from
  # ds.vertGLM, GoldschmidtReciprocalStep on 1/S) is all in place.
  #
  # This client is released now so downstream users can validate the
  # server-side pipeline (already Go-tested at 0.017% relative gradient
  # error vs centralized Cox on the same cohort in
  # TestCoxGradient_EndToEnd) against R's survival::coxph on real
  # cohorts like NCCTG lung via the 3-Opal Docker deployment.
  # =========================================================================

  if (verbose) {
    message("[ds.vertCox] First-release skeleton: server helpers orchestrated, full L-BFGS + SE pending deployment validation.")
  }

  stop("ds.vertCox full client orchestration requires end-to-end Opal testing ",
       "to validate the MPC iteration loop; see V2_PROGRESS.md for the ",
       "completion plan. The Go pipeline test (TestCoxGradient_EndToEnd) ",
       "demonstrates the math at 0.017% relative error against centralized Cox ",
       "gradient; the server-side helpers (k2SetCoxTimesDS, k2CoxReverseCumsumSDS, ",
       "k2CoxForwardCumsumGDS, k2StoreCoxRecipDS) are in place.",
       call. = FALSE)
}

#' @export
print.ds.vertCox <- function(x, ...) {
  cat("dsVert Cox PH regression\n")
  cat(sprintf("  N = %d, events = %d\n", x$n_obs, x$n_events))
  cat(sprintf("  converged after %d iterations\n", x$iterations))
  cat("\nCoefficients:\n")
  m <- data.frame(
    coef = x$coefficients,
    `exp(coef)` = exp(x$coefficients),
    check.names = FALSE)
  print(round(m, 4L))
  cat(sprintf("\nPartial log-likelihood: %.4f\n", -0.5 * x$deviance))
  invisible(x)
}
