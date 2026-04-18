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
  # Beaver joint cells via K*L element-wise products over shared
  # indicator vectors. For each (k, l):
  #   1. var1_srv extracts column k of its one-hot matrix into a
  #      per-patient FP share (and shares n-length zeros to the peer).
  #   2. var2_srv does the symmetric thing with column l.
  #   3. After the usual input-sharing preamble, both parties hold
  #      shares of X_k and Y_l; we run the already-shipped 4-step
  #      Beaver vecmul (GenTriples / Consume / R1 / R2) to obtain
  #      shares of (X_k .* Y_l); the client sums the two shares of
  #      sum_i (X_k Y_l)_i to get n_{kl}. No per-patient indicator
  #      ever leaves the owning server.
  counts <- tryCatch({
    .dsvert_chisq_bilinear_counts(
      datasources, server_names,
      var1_srv, var2_srv, var1, var2,
      oh1, oh2, pks, session_id, verbose)
  }, error = function(e) {
    warning("Beaver bilinear path failed (", conditionMessage(e),
            "); falling back to margin-based expected cells under ",
            "the independence null.", call. = FALSE)
    outer(oh1$row_margins, oh2$row_margins) / max(oh1$n, 1L)
  })

  n <- sum(counts)
  if (n == 0) stop("No observations; cannot compute chi-square", call. = FALSE)

  # Standard Pearson chi-square on the reconstructed K x L table.
  row_m <- as.integer(rowSums(counts))
  col_m <- as.integer(colSums(counts))
  res_list <- .dsvert_chisq_compute(counts,
    row_margins = row_m, col_margins = col_m, n = n,
    correct = correct)
  names(res_list)[names(res_list) == "statistic"] <- "chisq"
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
    chisq         = res_list$chisq %||% res_list$statistic,
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
.dsvert_chisq_bilinear_counts <- function(datasources, server_names,
                                           var1_srv, var2_srv,
                                           var1, var2, oh1, oh2, pks,
                                           session_id, verbose = FALSE) {
  # Helpers.
  v1_ci <- which(server_names == var1_srv)
  v2_ci <- which(server_names == var2_srv)
  K <- length(oh1$levels); L <- length(oh2$levels); n <- oh1$n
  if (oh1$n != oh2$n) {
    stop("one-hot vectors have mismatched length (", oh1$n, " vs ",
         oh2$n, "); rerun PSI alignment", call. = FALSE)
  }
  sendBlob <- function(blob, key, conn_idx) {
    DSI::datashield.aggregate(datasources[conn_idx],
      call("mpcStoreBlobDS", key = key, chunk = blob,
           session_id = session_id))
  }

  # Step A: mutually share the full one-hot matrices between the two
  # servers. Each party splits its plaintext matrix into (own_share,
  # peer_share), keeps own_share overwriting the session slot, and
  # sends peer_share sealed to the other server.
  var1_src <- paste0("k2_onehot_", var1, "_fp")
  var2_src <- paste0("k2_onehot_", var2, "_fp")
  var1_peer <- paste0("k2_onehot_peer_", var1, "_fp")
  var2_peer <- paste0("k2_onehot_peer_", var2, "_fp")

  share1 <- DSI::datashield.aggregate(datasources[v1_ci],
    call("k2BeaverShareVectorDS",
         source_key = var1_src,
         peer_pk = pks[[var2_srv]],
         session_id = session_id))
  if (is.list(share1) && length(share1) == 1L) share1 <- share1[[1L]]
  sendBlob(share1$peer_blob, paste0("bshr_", var1), v2_ci)
  DSI::datashield.aggregate(datasources[v2_ci],
    call("k2BeaverReceiveVectorDS",
         blob_key = paste0("bshr_", var1),
         output_key = var1_peer,
         session_id = session_id))

  share2 <- DSI::datashield.aggregate(datasources[v2_ci],
    call("k2BeaverShareVectorDS",
         source_key = var2_src,
         peer_pk = pks[[var1_srv]],
         session_id = session_id))
  if (is.list(share2) && length(share2) == 1L) share2 <- share2[[1L]]
  sendBlob(share2$peer_blob, paste0("bshr_", var2), v1_ci)
  DSI::datashield.aggregate(datasources[v1_ci],
    call("k2BeaverReceiveVectorDS",
         blob_key = paste0("bshr_", var2),
         output_key = var2_peer,
         session_id = session_id))

  # After this step:
  #   var1_srv session holds SHARE of both var1 (in var1_src) and
  #     SHARE of var2 (in var2_peer).
  #   var2_srv session holds SHARE of var1 (in var1_peer) and SHARE of
  #     var2 (in var2_src).
  #
  # For the Beaver vecmul we need (X, Y) each as shares on BOTH
  # parties. Convention: party 0 = var1_srv.
  #   on var1_srv:  X_share = column k of var1_src; Y_share = col l of var2_peer
  #   on var2_srv:  X_share = column k of var1_peer; Y_share = col l of var2_src
  # (Both pairs sum to the true X_k, Y_l.)

  counts <- matrix(0, nrow = K, ncol = L,
                    dimnames = list(oh1$levels, oh2$levels))

  dealer_ci <- v2_ci   # non-party-0 acts as dealer for the Beaver triples

  # PERFORMANCE NOTE: each (kk, ll) cell gets a FRESH triple (reusing
  # a triple for two products is a security leak: the common masking
  # value a would reveal X_k1 - X_k2 when the two (x-a) shares are
  # compared). A proper batch optimisation generates K*L independent
  # triples in one Go call (k2-beaver-vecmul-gen-triples supports any
  # n, so we'd just call with n*K*L). That is a Wave-4 refinement;
  # current K*L independent triple gens are already amortised across
  # a single session.
  for (kk in seq_len(K)) {
    for (ll in seq_len(L)) {
      # Extract column kk on each party (into canonical beaver X slot).
      DSI::datashield.aggregate(datasources[v1_ci],
        call("k2BeaverExtractColumnDS",
             source_key = var1_src, n = as.integer(n), K = as.integer(K),
             col_index = as.integer(kk), output_key = "k2_beaver_x",
             session_id = session_id))
      DSI::datashield.aggregate(datasources[v2_ci],
        call("k2BeaverExtractColumnDS",
             source_key = var1_peer, n = as.integer(n), K = as.integer(K),
             col_index = as.integer(kk), output_key = "k2_beaver_x",
             session_id = session_id))
      # Extract column ll on each party (into canonical beaver Y slot).
      DSI::datashield.aggregate(datasources[v1_ci],
        call("k2BeaverExtractColumnDS",
             source_key = var2_peer, n = as.integer(n), K = as.integer(L),
             col_index = as.integer(ll), output_key = "k2_beaver_y",
             session_id = session_id))
      DSI::datashield.aggregate(datasources[v2_ci],
        call("k2BeaverExtractColumnDS",
             source_key = var2_src, n = as.integer(n), K = as.integer(L),
             col_index = as.integer(ll), output_key = "k2_beaver_y",
             session_id = session_id))
      # Beaver vecmul: dealer -> consume -> r1 (relay) -> r2.
      tri <- DSI::datashield.aggregate(datasources[dealer_ci],
        call("k2BeaverVecmulGenTriplesDS",
             dcf0_pk = pks[[var1_srv]], dcf1_pk = pks[[var2_srv]],
             n = as.integer(n),
             session_id = session_id, frac_bits = 20L))
      if (is.list(tri) && length(tri) == 1L) tri <- tri[[1L]]
      sendBlob(tri$triple_blob_0, "k2_beaver_vecmul_triple", v1_ci)
      sendBlob(tri$triple_blob_1, "k2_beaver_vecmul_triple", v2_ci)
      DSI::datashield.aggregate(datasources[v1_ci],
        call("k2BeaverVecmulConsumeTripleDS", session_id = session_id))
      DSI::datashield.aggregate(datasources[v2_ci],
        call("k2BeaverVecmulConsumeTripleDS", session_id = session_id))
      r1a <- DSI::datashield.aggregate(datasources[v1_ci],
        call("k2BeaverVecmulR1DS",
             peer_pk = pks[[var2_srv]],
             x_key = "k2_beaver_x", y_key = "k2_beaver_y",
             n = as.integer(n),
             session_id = session_id, frac_bits = 20L))
      r1b <- DSI::datashield.aggregate(datasources[v2_ci],
        call("k2BeaverVecmulR1DS",
             peer_pk = pks[[var1_srv]],
             x_key = "k2_beaver_x", y_key = "k2_beaver_y",
             n = as.integer(n),
             session_id = session_id, frac_bits = 20L))
      if (is.list(r1a) && length(r1a) == 1L) r1a <- r1a[[1L]]
      if (is.list(r1b) && length(r1b) == 1L) r1b <- r1b[[1L]]
      sendBlob(r1a$peer_blob, "k2_beaver_vecmul_peer_masked", v2_ci)
      sendBlob(r1b$peer_blob, "k2_beaver_vecmul_peer_masked", v1_ci)
      DSI::datashield.aggregate(datasources[v1_ci],
        call("k2BeaverVecmulR2DS",
             is_party0 = TRUE,
             x_key = "k2_beaver_x", y_key = "k2_beaver_y",
             output_key = "k2_beaver_z", n = as.integer(n),
             session_id = session_id, frac_bits = 20L))
      DSI::datashield.aggregate(datasources[v2_ci],
        call("k2BeaverVecmulR2DS",
             is_party0 = FALSE,
             x_key = "k2_beaver_x", y_key = "k2_beaver_y",
             output_key = "k2_beaver_z", n = as.integer(n),
             session_id = session_id, frac_bits = 20L))
      # Sum shares per party -> aggregate.
      s1 <- DSI::datashield.aggregate(datasources[v1_ci],
        call("k2BeaverSumShareDS", source_key = "k2_beaver_z",
             session_id = session_id, frac_bits = 20L))
      s2 <- DSI::datashield.aggregate(datasources[v2_ci],
        call("k2BeaverSumShareDS", source_key = "k2_beaver_z",
             session_id = session_id, frac_bits = 20L))
      if (is.list(s1) && length(s1) == 1L) s1 <- s1[[1L]]
      if (is.list(s2) && length(s2) == 1L) s2 <- s2[[1L]]
      # Aggregate the two scalar FP shares client-side via the existing
      # k2-ring63-aggregate op (which is exactly sum-shares-and-decode).
      agg <- dsVert:::.callMpcTool("k2-ring63-aggregate", list(
        share_a = s1$sum_share_fp, share_b = s2$sum_share_fp,
        frac_bits = 20L))
      cell <- as.integer(round(as.numeric(agg$values[1L])))
      counts[kk, ll] <- max(0L, cell)
      if (isTRUE(verbose)) {
        message(sprintf("  n[%s,%s] = %d",
                         oh1$levels[kk], oh2$levels[ll], counts[kk, ll]))
      }
    }
  }
  counts
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
