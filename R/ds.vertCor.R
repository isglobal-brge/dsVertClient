#' @title Ring63 Privacy-Preserving Correlation for Vertically Partitioned Data
#' @description Computes Pearson correlation matrix using Ring63 Beaver MPC.
#'   No CKKS, no smudging noise. Only aggregate scalars disclosed.
#' @param data_name Character. Aligned data frame name.
#' @param variables Named list: server -> variable names.
#' @param datasources DataSHIELD connections.
#' @return List with correlation matrix, variable names, n_obs.
#' @export
ds.vertCor <- function(data_name, variables, log_n = 12, log_scale = 40,
                       datasources = NULL, reuse_mhe = FALSE) {
  if (!is.list(variables) || is.null(names(variables)))
    stop("variables must be a named list", call. = FALSE)
  if (length(variables) < 2)
    stop("At least 2 servers required", call. = FALSE)

  session_id <- local({
    hex <- sample(c(0:9, letters[1:6]), 32, replace = TRUE)
    hex[13] <- "4"; hex[17] <- sample(c("8","9","a","b"), 1)
    paste0(paste(hex[1:8],collapse=""),"-",paste(hex[9:12],collapse=""),"-",
           paste(hex[13:16],collapse=""),"-",paste(hex[17:20],collapse=""),"-",
           paste(hex[21:32],collapse=""))
  })

  if (is.null(datasources)) datasources <- DSI::datashield.connections_find()
  server_names <- names(datasources)
  server_list <- names(variables)

  on.exit({
    for (.s in server_list) {
      .ci <- which(server_names == .s)
      tryCatch(DSI::datashield.aggregate(datasources[.ci],
        call("mheCleanupDS", session_id = session_id)), error = function(e) NULL)
    }
  })

  .dsAgg <- function(conns, expr) {
    tryCatch(DSI::datashield.aggregate(conns = conns, expr = expr),
      error = function(e) { message("  [ERROR] ", conditionMessage(e)); stop(e) })
  }
  .sendBlob <- function(blob, key, conn_idx) {
    .dsvert_adaptive_send(blob, function(chunk_str, chunk_idx, n_chunks) {
      if (n_chunks == 1L) {
        .dsAgg(datasources[conn_idx], call("mheStoreBlobDS", key = key,
          chunk = chunk_str, session_id = session_id))
      } else {
        .dsAgg(datasources[conn_idx], call("mheStoreBlobDS", key = key,
          chunk = chunk_str, chunk_index = chunk_idx, n_chunks = n_chunks,
          session_id = session_id))
      }
    })
  }
  .to_b64url <- function(x) gsub("+","-",gsub("/","_",gsub("=+$","",x,perl=TRUE),fixed=TRUE),fixed=TRUE)
  .b64url_to_b64 <- function(x) {
    x <- gsub("-","+",gsub("_","/",x,fixed=TRUE),fixed=TRUE)
    pad <- nchar(x)%%4; if(pad==2) x<-paste0(x,"=="); if(pad==3) x<-paste0(x,"="); x
  }

  p_total <- sum(sapply(server_list, function(s) length(variables[[s]])))
  all_vars <- unlist(variables[server_list])
  message(sprintf("=== Ring63 Correlation (p=%d, %d servers) ===", p_total, length(server_list)))

  # Phase 0: Transport key setup
  t0 <- proc.time()[[3]]
  message("[Phase 0] Transport keys...")
  for (server in server_list) {
    ci <- which(server_names == server)
    r <- .dsAgg(datasources[ci], call("glmRing63TransportInitDS", session_id = session_id))
    if (is.list(r)) r <- r[[1]]
  }
  transport_pks <- list()
  for (server in server_list) {
    ci <- which(server_names == server)
    r <- .dsAgg(datasources[ci], call("glmRing63TransportInitDS", session_id = session_id))
    if (is.list(r)) r <- r[[1]]
    transport_pks[[server]] <- r$transport_pk
  }
  pk_sorted <- transport_pks[sort(names(transport_pks))]
  for (server in server_list) {
    ci <- which(server_names == server)
    .dsAgg(datasources[ci], call("mheStoreTransportKeysDS",
      transport_keys = pk_sorted, session_id = session_id))
  }

  # Phase 1: Standardize + get observation count
  message("[Phase 1] Standardizing...")
  # Get n_obs
  r_count <- .dsAgg(datasources[which(server_names == server_list[1])],
    call("getObsCountDS", data_name = data_name))
  if (is.list(r_count)) r_count <- r_count[[1]]
  n_obs <- as.integer(r_count$n_obs)

  for (server in server_list) {
    ci <- which(server_names == server)
    .dsAgg(datasources[ci], call("glmStandardizeDS",
      data_name = data_name, output_name = "dsvert_cor_std",
      x_vars = variables[[server]], y_var = NULL, session_id = session_id))
  }

  # Phase 2: Local correlations (within-server, plaintext)
  message("[Phase 2] Local correlations...")
  local_cors <- list()
  for (server in server_list) {
    ci <- which(server_names == server)
    r <- .dsAgg(datasources[ci], call("localCorDS",
      data_name = data_name, variables = variables[[server]],
      session_id = session_id))
    if (is.list(r) && length(r) == 1) r <- r[[1]]
    local_cors[[server]] <- r$correlation
  }

  # Phase 3: Input sharing (features to 2 computation parties)
  message("[Phase 3] Input sharing...")
  fusion <- server_list[1]; coord <- server_list[length(server_list)]
  dcf_parties <- c(fusion, coord)
  dcf_conns <- sapply(dcf_parties, function(s) which(server_names == s))
  non_dcf <- setdiff(server_list, dcf_parties)

  # DCF parties share with each other
  for (server in dcf_parties) {
    ci <- which(server_names == server)
    peer <- dcf_parties[dcf_parties != server]
    peer_pk <- .to_b64url(transport_pks[[peer]])
    r <- .dsAgg(datasources[ci], call("k2ShareInputDS",
      data_name = "dsvert_cor_std", x_vars = variables[[server]],
      y_var = NULL, peer_pk = peer_pk, session_id = session_id))
    if (is.list(r) && length(r) == 1) r <- r[[1]]
    peer_ci <- which(server_names == peer)
    .sendBlob(r$encrypted_x_share, "k2_peer_x_share", peer_ci)
  }
  # Receive peer shares
  for (di in seq_along(dcf_parties)) {
    ci <- dcf_conns[di]
    peer_p <- length(variables[[dcf_parties[3 - di]]])
    .dsAgg(datasources[ci], call("k2ReceiveShareDS",
      peer_p = as.integer(peer_p), session_id = session_id))
  }
  # Non-DCF servers share with both DCF parties
  for (server in non_dcf) {
    ci <- which(server_names == server)
    fusion_pk <- .to_b64url(transport_pks[[fusion]])
    r <- .dsAgg(datasources[ci], call("k2ShareInputDS",
      data_name = "dsvert_cor_std", x_vars = variables[[server]],
      y_var = NULL, peer_pk = fusion_pk, session_id = session_id))
    if (is.list(r) && length(r) == 1) r <- r[[1]]
    .sendBlob(r$encrypted_x_share, paste0("k2_extra_x_share_", server), dcf_conns[1])
    coord_pk <- .to_b64url(transport_pks[[coord]])
    r2 <- .dsAgg(datasources[ci], call("glmRing63ExportOwnShareDS",
      peer_pk = coord_pk, session_id = session_id))
    if (is.list(r2) && length(r2) == 1) r2 <- r2[[1]]
    .sendBlob(r2$encrypted_own_share, paste0("k2_extra_x_share_", server), dcf_conns[2])
  }
  for (server in non_dcf) {
    for (di in seq_along(dcf_parties)) {
      .dsAgg(datasources[dcf_conns[di]], call("glmRing63ReceiveExtraShareDS",
        extra_key = paste0("k2_extra_x_share_", server),
        extra_p = as.integer(length(variables[[server]])),
        session_id = session_id))
    }
  }

  # Initialize X_full via zero-beta eta computation
  p_coord <- length(variables[[coord]])
  p_fusion <- length(variables[[fusion]])
  p_extras <- p_total - p_coord - p_fusion

  for (di in seq_along(dcf_parties)) {
    ci <- dcf_conns[di]
    is_coord <- (dcf_parties[di] == coord)
    b_coord <- rep(0, p_coord)
    if (is_coord) {
      b_nl <- c(rep(0, p_fusion), rep(0, p_extras))
    } else {
      b_nl <- c(rep(0, p_extras), rep(0, p_fusion))
    }
    .dsAgg(datasources[ci], call("k2ComputeEtaShareDS",
      beta_coord = b_coord, beta_nl = b_nl,
      intercept = 0, is_coordinator = is_coord, session_id = session_id))
    # Reorder fusion's X_full to canonical order
    if (!is_coord && p_extras > 0) {
      .dsAgg(datasources[ci], call("glmRing63ReorderXFullDS",
        p_coord = as.integer(p_coord), p_fusion = as.integer(p_fusion),
        p_extras = as.integer(p_extras), session_id = session_id))
    }
  }

  # Set y_share = zeros (no response for correlation)
  for (di in seq_along(dcf_parties)) {
    ci <- dcf_conns[di]
    .dsAgg(datasources[ci], call("glmRing63CorSetZeroYDS", session_id = session_id))
  }

  # Phase 4: Cross-server correlations via Beaver
  message(sprintf("[Phase 4] Beaver cross-products (%d columns)...", p_total))
  dealer <- if (length(non_dcf) > 0) non_dcf[1] else fusion
  dealer_conn <- which(server_names == dealer)

  xtx <- matrix(0, p_total, p_total)
  for (j in seq_len(p_total)) {
    # Set column j as "mu" (residual = mu - 0 = col_j)
    # Use blob to pass column index (avoids DataSHIELD integer parsing issue)
    for (di in seq_along(dcf_parties)) {
      ci <- dcf_conns[di]
      .sendBlob(paste0(j - 1L, ",", p_total), "cor_col_params", ci)
      .dsAgg(datasources[ci], call("glmRing63CorSetColDS",
        from_storage = TRUE, session_id = session_id))
    }
    # Beaver matvec: X^T × col_j
    grad_t <- .dsAgg(datasources[dealer_conn],
      call("glmRing63GenGradTriplesDS",
           dcf0_pk = transport_pks[[dcf_parties[1]]],
           dcf1_pk = transport_pks[[dcf_parties[2]]],
           n = as.integer(n_obs), p = as.integer(p_total),
           session_id = session_id))
    if (is.list(grad_t)) grad_t <- grad_t[[1]]
    .sendBlob(grad_t$grad_blob_0, "k2_grad_triple_fp", dcf_conns[1])
    .sendBlob(grad_t$grad_blob_1, "k2_grad_triple_fp", dcf_conns[2])

    r1 <- list()
    for (di in seq_along(dcf_parties)) {
      ci <- dcf_conns[di]
      peer <- dcf_parties[3 - di]
      .dsAgg(datasources[ci], call("k2StoreGradTripleDS", session_id = session_id))
      r <- .dsAgg(datasources[ci], call("k2GradientR1DS",
        peer_pk = transport_pks[[peer]], session_id = session_id))
      if (is.list(r) && length(r) == 1) r <- r[[1]]
      r1[[di]] <- r
    }
    .sendBlob(r1[[1]]$encrypted_r1, "k2_grad_peer_r1", dcf_conns[2])
    .sendBlob(r1[[2]]$encrypted_r1, "k2_grad_peer_r1", dcf_conns[1])

    r2 <- list()
    for (di in seq_along(dcf_parties)) {
      ci <- dcf_conns[di]
      r <- .dsAgg(datasources[ci], call("k2GradientR2DS",
        party_id = as.integer(di - 1), session_id = session_id))
      if (is.list(r) && length(r) == 1) r <- r[[1]]
      r2[[di]] <- r
    }

    agg <- dsVert:::.callMheTool("k2-ring63-aggregate", list(
      share_a = r2[[1]]$gradient_fp, share_b = r2[[2]]$gradient_fp, frac_bits = 20L))
    xtx[, j] <- agg$values
    message(sprintf("  Column %d/%d done", j, p_total))
  }

  # Phase 5: Assembly
  corr <- xtx / (n_obs - 1)
  rownames(corr) <- colnames(corr) <- all_vars

  # Overlay local correlations (more precise for within-server)
  for (server in server_list) {
    idx <- which(all_vars %in% variables[[server]])
    if (length(idx) > 1 && !is.null(local_cors[[server]])) {
      lc <- local_cors[[server]]
      if (is.matrix(lc)) corr[idx, idx] <- lc
    }
  }

  elapsed <- proc.time()[[3]] - t0
  message(sprintf("Correlation complete (%.0fs)", elapsed))

  structure(list(
    correlation = corr,
    var_names = all_vars,
    n_obs = n_obs,
    method = "Ring63-Beaver",
    servers = server_list,
    local_correlations = local_cors
  ), class = "ds.cor")
}
