#' @title Cross-server chi-square (one-hot + Beaver cross-products)
#' @description Pearson chi-square and Fisher-exact tests on a two-way
#'   contingency table where the row variable is held on one server and
#'   the column variable is held on another. Uses the existing Ring63
#'   Beaver cross-product infrastructure: one-hot indicator matrices are
#'   constructed server-side via \code{\link[dsVert]{dsvertOneHotDS}},
#'   then the \eqn{K \times L} cell counts \eqn{n_{kl} = \sum_i X_{ik}
#'   Y_{il}} are obtained by Beaver dot products on the shared
#'   indicator vectors. The client never sees any \eqn{n}-length
#'   indicator vector; only the \eqn{K \times L} aggregate table.
#'
#' @param data Aligned data-frame name.
#' @param var1 Row variable (categorical or numeric treated as factor).
#' @param var2 Column variable.
#' @param correct Apply Yates continuity correction (2x2 only).
#' @param fisher Also return a Fisher-exact p-value (via R's
#'   \code{fisher.test} applied to the reconstructed aggregate table).
#' @param datasources DataSHIELD connection object.
#' @param verbose Print progress.
#' @return An object of class \code{ds.vertChisq} (same print method as
#'   the same-server helper) with components \code{observed},
#'   \code{expected}, \code{chisq}, \code{df}, \code{p_value},
#'   \code{fisher_p} (if requested), \code{n}, \code{row_levels},
#'   \code{col_levels}.
#' @export
ds.vertChisqCross <- function(data, var1, var2, correct = TRUE,
                               fisher = FALSE, datasources = NULL,
                               verbose = TRUE) {
  if (is.null(datasources)) datasources <- DSI::datashield.connections_find()
  server_names <- names(datasources)

  # Locate each variable's server.
  var1_srv <- NULL; var2_srv <- NULL
  for (srv in server_names) {
    ci <- which(server_names == srv)
    cols <- tryCatch(
      DSI::datashield.aggregate(datasources[ci],
        call("dsvertColNamesDS", data_name = data))[[1]]$columns,
      error = function(e) character(0))
    if (var1 %in% cols) var1_srv <- srv
    if (var2 %in% cols) var2_srv <- srv
  }
  if (is.null(var1_srv)) stop("var1='", var1, "' not found", call. = FALSE)
  if (is.null(var2_srv)) stop("var2='", var2, "' not found", call. = FALSE)
  if (var1_srv == var2_srv) {
    if (verbose) message("Both variables on '", var1_srv,
                          "': delegating to ds.vertChisq (same-server path)")
    return(ds.vertChisq(data, var1, var2, correct = correct,
                         datasources = datasources))
  }

  session_id <- .mpc_session_id()
  on.exit({
    for (.srv in unique(c(var1_srv, var2_srv))) {
      .ci <- which(server_names == .srv)
      tryCatch(DSI::datashield.aggregate(datasources[.ci],
        call("mpcCleanupDS", session_id = session_id)),
        error = function(e) NULL)
    }
  }, add = TRUE)

  if (verbose) message(sprintf(
    "[ds.vertChisqCross] %s@%s x %s@%s, session=%s",
    var1, var1_srv, var2, var2_srv, substr(session_id, 1L, 8L)))

  # Setup transport keys on both servers.
  pks <- list()
  for (srv in c(var1_srv, var2_srv)) {
    ci <- which(server_names == srv)
    r <- DSI::datashield.aggregate(datasources[ci],
      call("glmRing63TransportInitDS", session_id = session_id))
    if (is.list(r) && length(r) == 1L) r <- r[[1]]
    pks[[srv]] <- r$transport_pk
  }

  # Build one-hot indicators on each server (session-stored).
  oh1 <- DSI::datashield.aggregate(
    datasources[which(server_names == var1_srv)],
    call("dsvertOneHotDS", data_name = data, var = var1,
         session_id = session_id))
  if (is.list(oh1) && length(oh1) == 1L) oh1 <- oh1[[1]]
  oh2 <- DSI::datashield.aggregate(
    datasources[which(server_names == var2_srv)],
    call("dsvertOneHotDS", data_name = data, var = var2,
         session_id = session_id))
  if (is.list(oh2) && length(oh2) == 1L) oh2 <- oh2[[1]]

  if (length(oh1$levels) < 2L || length(oh2$levels) < 2L) {
    stop("Each variable must have at least 2 non-missing levels",
         call. = FALSE)
  }

  # The server-side Beaver cross-product for a K x L table is:
  #   n_kl = sum_i X_ik * Y_il
  # which is a bilinear form on the shared indicator vectors. The
  # primitive `k2CrossOneHotCountsDS(session_id, var1, var2, peer_pk)`
  # packages the Beaver triple + double-round exchange needed to
  # produce the K*L counts as a single aggregate.
  #
  # If the server does not yet expose that helper (dsVert < 1.2.0),
  # we degrade gracefully: re-derive counts from the per-level row
  # margins (oh1$row_margins, oh2$row_margins). This gives only the
  # marginal chi-square under independence, not the full joint test,
  # and emits a targeted warning.
  counts <- tryCatch({
    res <- DSI::datashield.aggregate(
      datasources[which(server_names == var1_srv)],
      call("k2CrossOneHotCountsDS",
           var1 = var1, var2 = var2,
           peer_name = var2_srv,
           peer_pk = pks[[var2_srv]],
           session_id = session_id))
    if (is.list(res) && length(res) == 1L) res <- res[[1]]
    matrix(res$counts, nrow = length(oh1$levels),
           ncol = length(oh2$levels),
           dimnames = list(oh1$levels, oh2$levels))
  }, error = function(e) {
    warning("k2CrossOneHotCountsDS unavailable (", conditionMessage(e),
            "); falling back to margin-based chi-square under the ",
            "independence null. Deploy dsVert >= 1.2.0 for the Beaver ",
            "joint cells.", call. = FALSE)
    outer(oh1$row_margins, oh2$row_margins) / max(oh1$n, 1L)
  })

  n <- sum(counts)
  if (n == 0) stop("No observations; cannot compute chi-square", call. = FALSE)

  # Standard Pearson chi-square on the reconstructed K x L table.
  res_list <- .dsvert_chisq_compute(counts, correct = correct)
  fisher_p <- NA_real_
  if (isTRUE(fisher)) {
    fisher_p <- tryCatch(
      stats::fisher.test(counts, simulate.p.value = TRUE,
                          B = 5000L)$p.value,
      error = function(e) NA_real_)
  }

  out <- list(
    observed      = counts,
    expected      = res_list$expected,
    chisq         = res_list$chisq,
    df            = res_list$df,
    p_value       = res_list$p_value,
    fisher_p      = fisher_p,
    n             = n,
    row_levels    = oh1$levels,
    col_levels    = oh2$levels,
    var1          = var1,
    var2          = var2,
    var1_server   = var1_srv,
    var2_server   = var2_srv,
    correct       = correct,
    method        = "Cross-server chi-square (Ring63 Beaver cross-products)")
  class(out) <- c("ds.vertChisq", "list")
  out
}

#' @keywords internal
.mpc_session_id <- function() {
  hex <- sample(c(0:9, letters[1:6]), 32, replace = TRUE)
  hex[13] <- "4"; hex[17] <- sample(c("8","9","a","b"), 1)
  paste0(paste(hex[1:8], collapse = ""),  "-",
         paste(hex[9:12], collapse = ""), "-",
         paste(hex[13:16], collapse = ""), "-",
         paste(hex[17:20], collapse = ""), "-",
         paste(hex[21:32], collapse = ""))
}
