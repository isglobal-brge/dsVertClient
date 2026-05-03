#' Share-domain weighted residual helper for GLM/IPW.
#'
#' The DCF parties hold additive shares of both the residual vector and the
#' weight vector. This helper runs one Beaver vecmul round to install a share
#' of w * residual (or sqrt(w) * residual for Gaussian deviance) as the active
#' residual consumed by the existing GLM gradient/deviance machinery.
#'
#' @keywords internal
#' @noRd
.glm_apply_shared_weight_residual <- function(
    datasources, dcf_parties, dcf_conns, dealer_conn, transport_pks,
    session_id, n_obs, .dsAgg, .sendBlob,
    weight_key = "k2_weights_share_fp",
    output_key = "k2_weighted_residual_share_fp",
    ring = 63L) {
  ring <- as.integer(ring)
  if (!ring %in% c(63L, 127L)) stop("ring must be 63 or 127", call. = FALSE)
  frac_bits <- if (ring == 127L) 50L else 20L
  n_int <- as.integer(n_obs)

  for (ci in dcf_conns) {
    .dsAgg(datasources[ci], call(name = "k2PrepareWeightedResidualShareDS",
                                 session_id = session_id))
  }

  tri <- .dsAgg(datasources[dealer_conn],
    call(name = "k2BeaverVecmulGenTriplesDS",
         dcf0_pk = transport_pks[[dcf_parties[1L]]],
         dcf1_pk = transport_pks[[dcf_parties[2L]]],
         n = as.numeric(n_int), session_id = session_id,
         frac_bits = frac_bits, ring = ring))
  if (is.list(tri) && length(tri) == 1L) tri <- tri[[1L]]
  .sendBlob(tri$triple_blob_0, "k2_beaver_vecmul_triple", dcf_conns[1L])
  .sendBlob(tri$triple_blob_1, "k2_beaver_vecmul_triple", dcf_conns[2L])

  for (ci in dcf_conns) {
    .dsAgg(datasources[ci], call(name = "k2BeaverVecmulConsumeTripleDS",
                                 session_id = session_id))
  }

  r1 <- vector("list", 2L)
  for (i in seq_along(dcf_parties)) {
    ci <- dcf_conns[i]
    peer <- dcf_parties[3L - i]
    r <- .dsAgg(datasources[ci], call(name = "k2BeaverVecmulR1DS",
      peer_pk = transport_pks[[peer]],
      x_key = weight_key,
      y_key = "k2_weight_residual_share_fp",
      n = as.numeric(n_int), session_id = session_id,
      frac_bits = frac_bits, ring = ring))
    if (is.list(r) && length(r) == 1L) r <- r[[1L]]
    r1[[i]] <- r
  }
  .sendBlob(r1[[1L]]$peer_blob, "k2_beaver_vecmul_peer_masked", dcf_conns[2L])
  .sendBlob(r1[[2L]]$peer_blob, "k2_beaver_vecmul_peer_masked", dcf_conns[1L])

  for (i in seq_along(dcf_parties)) {
    .dsAgg(datasources[dcf_conns[i]], call(name = "k2BeaverVecmulR2DS",
      is_party0 = (i == 1L),
      x_key = weight_key,
      y_key = "k2_weight_residual_share_fp",
      output_key = output_key,
      n = as.numeric(n_int), session_id = session_id,
      frac_bits = frac_bits, ring = ring))
  }

  for (ci in dcf_conns) {
    .dsAgg(datasources[ci], call(name = "k2FinalizeWeightedResidualShareDS",
                                 input_key = output_key,
                                 session_id = session_id))
  }
  invisible(NULL)
}
