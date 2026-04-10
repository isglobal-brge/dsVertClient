#' @title GLM Setup: Transport Keys + Standardization
#' @description Initializes transport keys on all servers and standardizes
#'   features. Pure Ring63 MPC.
#' @return List with transport_pks, x_means, x_sds, y_mean, y_sd,
#'   std_data, standardize_y, .dsAgg, .sendBlob
#' @keywords internal

.glm_mpc_setup <- function(datasources, server_names, server_list,
                           non_label_servers, y_server, y_var, x_vars,
                           data_name, family, n_obs, log_n, log_scale,
                           generate_rlk, use_secure_agg, use_k2_beaver,
                           reuse_session, session_id, verbose) {

  # =========================================================================
  # Helpers (closures capturing datasources, session_id)
  # =========================================================================
  .dsAgg <- function(conns, expr, ...) {
    tryCatch(
      DSI::datashield.aggregate(conns = conns, expr = expr, ...),
      error = function(e) {
        msg <- conditionMessage(e)
        fn_name <- if (is.call(expr)) as.character(expr[[1]]) else "?"
        srv_name <- tryCatch(names(conns)[1], error = function(x) "?")
        message("  [dsAgg ERROR] ", fn_name, " on ", srv_name, ": ", msg)
        ds_errs <- tryCatch(DSI::datashield.errors(), error = function(x) NULL)
        if (!is.null(ds_errs) && length(ds_errs) > 0)
          for (nm in names(ds_errs)) message("    ", nm, ": ", ds_errs[[nm]])
        if (grepl("500|NullPointer|Internal Server Error", msg)) {
          Sys.sleep(2)
          DSI::datashield.aggregate(conns = conns, expr = expr, ...)
        } else stop(e)
      })
  }

  .sendBlob <- function(blob, key, conn_idx) {
    .dsvert_adaptive_send(blob, function(chunk_str, chunk_idx, n_chunks) {
      if (n_chunks == 1L) {
        .dsAgg(conns = datasources[conn_idx],
          expr = call("mpcStoreBlobDS", key = key, chunk = chunk_str,
                      session_id = session_id))
      } else {
        .dsAgg(conns = datasources[conn_idx],
          expr = call("mpcStoreBlobDS", key = key, chunk = chunk_str,
                      chunk_index = chunk_idx, n_chunks = n_chunks,
                      session_id = session_id))
      }
    })
  }

  # =========================================================================
  # Phase 0: Transport Key Setup (Ring63)
  # =========================================================================
  transport_pks <- list()

  if (length(non_label_servers) > 0) {
    if (verbose) message("\n[Phase 0] Transport key setup (", length(server_list), " servers)...")
    t0_key <- proc.time()[[3]]

    for (server in server_list) {
      conn_idx <- which(server_names == server)
      tk_result <- .dsAgg(conns = datasources[conn_idx],
        expr = call("glmRing63TransportInitDS", session_id = session_id))
      if (is.list(tk_result)) tk_result <- tk_result[[1]]
      transport_pks[[server]] <- tk_result$transport_pk
    }

    pk_sorted <- transport_pks[sort(names(transport_pks))]
    for (server in server_list) {
      conn_idx <- which(server_names == server)
      .dsAgg(conns = datasources[conn_idx],
        expr = call("mpcStoreTransportKeysDS",
                    transport_keys = pk_sorted, session_id = session_id))
    }

    if (verbose) message(sprintf("  [Key Setup] Transport keys exchanged (%d servers, %.1fs)",
                                   length(server_list), proc.time()[[3]] - t0_key))
  }

  # =========================================================================
  # Phase 1: Standardize features
  # =========================================================================
  if (verbose) message("\n[Phase 1] Standardizing features across ", length(server_list), " servers...")
  t0_std <- proc.time()[[3]]
  std_data <- paste0(data_name, "_std")
  standardize_y <- (family == "gaussian")

  x_means <- list()
  x_sds <- list()
  y_mean <- NULL
  y_sd <- NULL

  for (server in server_list) {
    conn_idx <- which(server_names == server)
    y_arg <- if (server == y_server && standardize_y) y_var else NULL

    std_result <- .dsAgg(conns = datasources[conn_idx],
      expr = call("glmStandardizeDS",
                  data_name = data_name, output_name = std_data,
                  x_vars = x_vars[[server]], y_var = y_arg,
                  session_id = session_id))
    if (is.list(std_result) && length(std_result) == 1)
      std_result <- std_result[[1]]

    x_means[[server]] <- std_result$x_means
    x_sds[[server]] <- std_result$x_sds
    if (!is.null(std_result$y_mean)) {
      y_mean <- std_result$y_mean
      y_sd <- std_result$y_sd
    }
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
    .dsAgg         = .dsAgg,
    .sendBlob      = .sendBlob
  )
}
