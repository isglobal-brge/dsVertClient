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
    datasources, dcf_parties, dcf_conns, dealer_conn = NULL, transport_pks,
    session_id, n_obs, .dsAgg, .sendBlob,
    weight_key = "k2_weights_share_fp",
    output_key = "k2_weighted_residual_share_fp",
    ring = 63L) {
  ring <- as.integer(ring)
  if (!ring %in% c(63L, 127L)) stop("ring must be 63 or 127", call. = FALSE)
  frac_bits <- if (ring == 127L) 50L else 20L
  n_int <- as.integer(n_obs)

  if (identical(ring, 127L)) {
    exact_product <- if (
      identical(weight_key, "k2_weights_share_fp") &&
      identical(output_key, "k2_weighted_residual_share_fp")) {
      "weight"
    } else if (
      identical(weight_key, "k2_sqrt_weights_share_fp") &&
      identical(output_key, "k2_sqrt_weighted_residual_share_fp")) {
      "sqrt_weight"
    } else {
      stop("Ring127 weighted residual requested an unapproved exact purpose.",
           call. = FALSE)
    }
    all_names <- names(datasources) %||% dcf_parties
    selected_names <- all_names[dcf_conns]
    conns <- datasources[dcf_conns]
    names(conns) <- selected_names
    .dsvert_setup_exact_gc_transport(
      datasources = datasources, server_names = all_names,
      servers = dcf_conns, session_id = session_id,
      .aggregate = .dsAgg)
    prepared <- .dsvert_aggregate_strict(
      conns,
      call(name = "k2PrepareWeightedResidualShareDS",
           exact_product = exact_product, session_id = session_id),
      operation = "GLM producer-bound weighted residual preparation",
      .aggregate = .dsAgg)
    manifests <- stats::setNames(lapply(selected_names, function(server) {
      value <- prepared[[server]]
      if (!is.list(value) || !isTRUE(value$stored) ||
          !is.list(value$exact_vecmul_manifest)) {
        stop("A GLM peer did not mint its exact multiplication manifest.",
             call. = FALSE)
      }
      value$exact_vecmul_manifest
    }), selected_names)
    .dsvert_exact_gc_vecmul_run(
      datasources = datasources,
      server_names = all_names,
      servers = dcf_conns, session_id = session_id, total_n = n_int,
      x_key = weight_key, y_key = "k2_weight_residual_share_fp",
      output_key = output_key, input_manifests = manifests,
      transport_ready = TRUE,
      .aggregate = .dsAgg)
    for (ci in dcf_conns) {
      .dsAgg(datasources[ci], call(
        name = "k2FinalizeWeightedResidualShareDS",
        input_key = output_key, session_id = session_id))
    }
    return(invisible(NULL))
  }

  for (ci in dcf_conns) {
    .dsAgg(datasources[ci], call(name = "k2PrepareWeightedResidualShareDS",
                                 session_id = session_id))
  }

  .ot_beaver_prepare_vecmul(
    datasources = datasources,
    party_conns = dcf_conns,
    party_names = dcf_parties,
    transport_pks = transport_pks,
    session_id = session_id,
    n = n_int,
    ring = ring,
    .dsAgg = .dsAgg,
    .sendBlob = .sendBlob)

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
  .sendBlob(r1[[1L]]$peer_blob, r1[[1L]]$peer_transfer, dcf_conns[2L])
  .sendBlob(r1[[2L]]$peer_blob, r1[[2L]]$peer_transfer, dcf_conns[1L])

  for (i in seq_along(dcf_parties)) {
    .dsAgg(datasources[dcf_conns[i]], call(name = "k2BeaverVecmulR2DS",
      is_party0 = (i == 1L),
      peer_name = dcf_parties[3L - i],
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
