#' @title GLM Setup: Transport Keys + Standardization
#' @description Initializes transport keys on all servers and standardizes
#'   features. Pure Ring63 MPC.
#' @return List with transport_pks, x_means, x_sds, y_mean, y_sd,
#'   std_data, standardize_y, .dsAgg, .sendBlob
#' @keywords internal

.glm_mpc_setup <- function(datasources, server_names, server_list,
                           non_label_servers, y_server, y_var, x_vars,
                           data_name, family, session_id, verbose,
                           standardize_y_override = NULL,
                           std_mode = "full", missing = "fail",
                           numeric_ring = 63L,
                           .aggregate = DSI::datashield.aggregate) {

  # =========================================================================
  # Helpers (closures capturing datasources, session_id)
  # =========================================================================
  .dsAgg <- function(conns, expr, ...) {
    if (length(list(...))) {
      stop("Additional aggregate controls are unavailable in a protected MPC phase",
           call. = FALSE)
    }
    .dsvert_aggregate_strict(
      conns = conns, expr = expr, operation = "GLM MPC protocol phase",
      .aggregate = .aggregate)
  }

  .sendBlob <- function(blob, contract, conn_idx) {
    .dsvert_store_transfer_or_legacy(
      blob, contract, datasources[conn_idx], session_id,
      producer_conns = datasources,
      .aggregate = .aggregate)
  }

  # =========================================================================
  # Phase 0: Transport Key Setup (Ring63) + Identity Verification
  # =========================================================================
  transport_pks <- list()
  identity_info <- list()

  if (length(non_label_servers) > 0) {
    if (verbose) message("\n[Phase 0] Transport key setup (", length(server_list), " servers)...")
    t0_key <- proc.time()[[3]]

    tk_results <- .dsAgg(
      conns = datasources[server_list],
      expr = call(name = "glmRing63TransportInitDS", session_id = session_id))
    for (server in server_list) {
      tk_result <- tk_results[[server]]
      transport_pks[[server]] <- tk_result$transport_pk
      if (!is.null(tk_result$identity_pk))
        identity_info[[server]] <- list(
          identity_pk = tk_result$identity_pk, signature = tk_result$signature)
    }

    pk_sorted <- transport_pks[sort(names(transport_pks))]
    id_sorted <- if (length(identity_info) > 0) identity_info[sort(names(identity_info))] else NULL
    .to_b64url <- function(x) gsub("\\+","-",gsub("/","_",gsub("=+$","",x,perl=TRUE),fixed=TRUE),fixed=TRUE)
    .json_to_b64url <- function(x) .to_b64url(gsub("\n","",jsonlite::base64_enc(charToRaw(jsonlite::toJSON(x, auto_unbox = TRUE))),fixed=TRUE))
    pk_b64 <- .json_to_b64url(pk_sorted)
    id_b64 <- if (!is.null(id_sorted)) .json_to_b64url(id_sorted) else ""
    .dsvert_aggregate_strict(
      conns = datasources[server_list],
      expr = call(name = "mpcStoreTransportKeysDS",
                  transport_keys_b64 = pk_b64, identity_info_b64 = id_b64,
                  session_id = session_id),
      operation = "GLM pinned-key binding",
      result_contract = "logical_true", .aggregate = .aggregate)

    if (verbose) message(sprintf("  [Key Setup] Transport keys exchanged (%d servers, %.1fs)",
                                   length(server_list), proc.time()[[3]] - t0_key))
  }

  # =========================================================================
  # Phase 1: Standardize features
  # =========================================================================
  if (verbose) message("\n[Phase 1] Standardizing features across ", length(server_list), " servers...")
  t0_std <- proc.time()[[3]]
  std_data <- paste0(data_name, "_std")
  # Override can disable y-standardisation (ds.vertLMM's GLS fit: the
  # design matrix already encodes the intercept via (1-lambda_i), so
  # mean-shifting y breaks the no-auto-intercept OLS on standardised
  # features). When the override is FALSE we ALSO skip x-standardisation
  # so the K=2 loop sees the raw design matrix.
  standardize_y <- if (!is.null(standardize_y_override))
    isTRUE(standardize_y_override) else (family == "gaussian")
  skip_std <- !is.null(standardize_y_override) && !isTRUE(standardize_y_override)
  # std_mode: "full" / "scale_only" / "none". "scale_only" keeps the
  # design well-conditioned for L-BFGS without subtracting the column
  # means (required for ds.vertLMM's no-const closed-form GLS fit).
  .std_mode <- if (skip_std) "none" else std_mode

  x_means <- list()
  x_sds <- list()
  y_mean <- NULL
  y_sd <- NULL
  numeric_attestations <- list()

  std_exprs <- stats::setNames(lapply(server_list, function(server) {
    y_arg <- if (server == y_server && standardize_y) y_var else NULL

    srv_x <- x_vars[[server]]
    if (length(srv_x) == 0) srv_x <- NULL
    call(name = "glmStandardizeDS",
         data_name = data_name, output_name = std_data,
         x_vars = srv_x, y_var = y_arg,
         session_id = session_id,
         skip_standardize = isTRUE(skip_std),
         mode = .std_mode,
         missing = missing,
         numeric_family = family,
         numeric_ring = as.integer(numeric_ring),
         numeric_y_var = if (server == y_server) y_var else NULL)
  }), server_list)
  std_results <- .dsAgg(
    conns = datasources[server_list], expr = std_exprs)

  for (server in server_list) {
    std_result <- std_results[[server]]
    x_means[[server]] <- std_result$x_means
    x_sds[[server]] <- std_result$x_sds
    if (!is.null(std_result$y_mean)) {
      y_mean <- std_result$y_mean
      y_sd <- std_result$y_sd
    }
    numeric_attestations[[server]] <- std_result$numeric_attestation
  }

  if (verbose) {
    total_feats <- sum(sapply(x_vars, length))
    message(sprintf("  [Standardize] %d total features standardized (y %s, %.1fs)",
                    total_feats, if (standardize_y) "standardized" else "raw",
                    proc.time()[[3]] - t0_std))
  }

  # =========================================================================
  # Return
  # =========================================================================
  list(
    transport_pks  = transport_pks,
    x_means        = x_means,
    x_sds          = x_sds,
    y_mean         = y_mean,
    y_sd           = y_sd,
    std_data       = std_data,
    standardize_y  = standardize_y,
    missing_policy = missing,
    numeric_attestations = numeric_attestations,
    .dsAgg         = .dsAgg,
    .sendBlob      = .sendBlob
  )
}
