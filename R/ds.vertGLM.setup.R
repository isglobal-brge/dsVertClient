# ===========================================================================
# .glm_mhe_setup() — Phase 0 (MHE key setup), Phase 1 (standardize),
#                     Phase 2 (encrypt y)
#
# Extracted from ds.vertGLM() for readability.  All protocol logic is
# identical to the original monolithic version.
#
# Returns a list with everything the BCD loop and post-processing need:
#   transport_pks, cpk (K>=3 non-Gaussian only), x_means, x_sds,
#   y_mean, y_sd, std_data, standardize_y, .dsAgg, .sendBlob
# ===========================================================================

.glm_mhe_setup <- function(datasources, server_names, server_list,
                           non_label_servers, y_server, y_var, x_vars,
                           data_name, family, n_obs, log_n, log_scale,
                           generate_rlk, use_secure_agg, use_k2_beaver,
                           reuse_mhe, session_id, verbose) {

  # =========================================================================
  # Helpers (closures capturing datasources, session_id)
  # =========================================================================
  .need_async <- FALSE

  .heavy_fns <- c("mheInitDS", "mheCombineDS", "mheEncryptRawDS",
    "glmHEEncryptEtaDS", "glmHELinkStepDS", "glmHEGradientEncDS",
    "mhePartialDecryptBatchWrappedDS", "mheFuseBatchDS",
    "mheRLKAggregateR1DS", "mheRLKRound2DS",
    "mhePartialDecryptWrappedDS", "mheFuseServerDS")

  .dsAgg <- function(conns, expr, ...) {
    if (.need_async && length(conns) == 1 && is.call(expr) &&
        as.character(expr[[1]]) %in% .heavy_fns) {
      conn <- conns[[1]]
      if (inherits(conn, "OpalConnection")) {
        opal <- conn@opal
        script <- asNamespace("DSOpal")$.deparse(expr)
        cmd_id <- opalr::opal.post(opal, "datashield", "session",
          opal$rid, "aggregate", query = list(async = "true"),
          body = script, contentType = "application/x-rscript",
          acceptType = "application/octet-stream")
        repeat {
          Sys.sleep(2)
          cmd <- opalr::opal.command(opal, cmd_id, wait = FALSE)
          if (cmd$status %in% c("COMPLETED", "FAILED", "CANCELED"))
            break
        }
        if (cmd$status != "COMPLETED") {
          msg <- cmd$error %||% cmd$status
          opalr::opal.command_rm(opal, cmd_id)
          stop("Server-side command failed on '", names(conns)[1],
               "': ", msg, call. = FALSE)
        }
        result <- opalr::opal.command_result(opal, cmd_id)
        opalr::opal.command_rm(opal, cmd_id)
        return(setNames(list(result), names(conns)[1]))
      }
    }
    # Retry once on transient 500 errors (Opal OrientDB/session bugs)
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
        } else {
          stop(e)
        }
      }
    )
  }

  # Send CT chunks to a server (reusable helper, adaptive chunk size)
  .sendCTChunks <- function(ct_str, conn_idx) {
    .dsvert_adaptive_send(ct_str, function(chunk_str, chunk_idx, n_chunks) {
      .dsAgg(
        conns = datasources[conn_idx],
        expr = call("mheStoreCTChunkDS",
                    chunk_index = chunk_idx,
                    chunk = chunk_str,
                    session_id = session_id)
      )
    })
  }

  # Relay a wrapped share to the fusion server in chunks (adaptive chunk size)
  .sendWrappedShare <- function(wrapped_share, party_id, fusion_conn_idx) {
    chunk_size <- .dsvert_get_chunk_size()
    n_chars <- nchar(wrapped_share)
    n_chunks <- ceiling(n_chars / chunk_size)
    for (ch in seq_len(n_chunks)) {
      start <- (ch - 1) * chunk_size + 1
      end <- min(ch * chunk_size, n_chars)
      .dsAgg(
        conns = datasources[fusion_conn_idx],
        expr = call("mheStoreWrappedShareDS",
                    party_id = as.integer(party_id),
                    share_data = substr(wrapped_share, start, end),
                    session_id = session_id)
      )
    }
  }

  # Send a blob to a server with adaptive chunking and fallback
  .sendBlob <- function(blob, key, conn_idx) {
    .dsvert_adaptive_send(blob, function(chunk_str, chunk_idx, n_chunks) {
      if (n_chunks == 1L) {
        .dsAgg(
          conns = datasources[conn_idx],
          expr = call("mheStoreBlobDS", key = key, chunk = chunk_str,
                      session_id = session_id)
        )
      } else {
        .dsAgg(
          conns = datasources[conn_idx],
          expr = call("mheStoreBlobDS",
                      key = key,
                      chunk = chunk_str,
                      chunk_index = chunk_idx,
                      n_chunks = n_chunks,
                      session_id = session_id)
        )
      }
    })
  }

  # =========================================================================
  # Phase 0: MHE Key Setup + Transport Keys
  # =========================================================================
  transport_pks <- list()
  cpk <- NULL

  if (use_k2_beaver && length(non_label_servers) > 0) {
    # K=2 wide spline: generate transport keys on SERVERS via lightweight mheInitDS.
    # This creates the X25519 keypair server-side (needed for k2ReceiveShareDS decrypt).
    if (verbose) message("\n[Phase 0] Transport key setup (K=2 wide spline)...")
    t0_key <- proc.time()[[3]]
    crp_k2 <- NULL
    gkg_k2 <- NULL
    if (verbose) message("  [Key Setup] Generating X25519 keypairs on ", length(server_list), " servers...")
    for (server in server_list) {
      conn_idx <- which(server_names == server)
      party_id <- as.integer(which(server_list == server) - 1L)
      tk_result <- .dsAgg(
        conns = datasources[conn_idx],
        expr = call("mheInitDS",
                    party_id = party_id,
                    crp = crp_k2, gkg_seed = gkg_k2,
                    num_obs = as.integer(n_obs),
                    log_n = 14L, log_scale = 40L,
                    generate_rlk = FALSE,
                    session_id = session_id)
      )
      if (is.list(tk_result)) tk_result <- tk_result[[1]]
      transport_pks[[server]] <- tk_result$transport_pk
      if (is.null(crp_k2)) {
        crp_k2 <- tk_result$crp
        gkg_k2 <- tk_result$gkg_seed
      }
    }
    # Store transport keys on each server (for peer PK lookup in k2ShareInputDS)
    if (verbose) message("  [Key Setup] Distributing transport public keys...")
    for (server in server_list) {
      conn_idx <- which(server_names == server)
      pk_sorted <- transport_pks[sort(names(transport_pks))]
      .dsAgg(
        conns = datasources[conn_idx],
        expr = call("mheStoreTransportKeysDS",
                    transport_keys = pk_sorted,
                    session_id = session_id)
      )
    }
    if (verbose) message(sprintf("  [Key Setup] Transport keys exchanged for %d servers (%.1fs)",
                                   length(server_list), proc.time()[[3]] - t0_key))
  }

  if (length(non_label_servers) > 0 && !use_k2_beaver) {
    if (verbose) message("\n[Phase 0] MHE key setup (", length(server_list), " servers, logN=", log_n, ", logScale=", log_scale, ")...")
    t0_mhe <- proc.time()[[3]]

    mhe_reused <- FALSE

    # --- MHE Context Reuse Check ---
    # NOTE: Use '_' as separator -- Opal's R expression parser rejects many
    # special characters ('|', '/', etc.) inside string arguments passed via
    # datashield.aggregate(). Only alphanumeric, '_', ',', ':', '-' are safe.
    context_id <- paste0(
      "peers_", paste(sort(server_list), collapse = ","),
      "_logn_", log_n,
      "_logscale_", log_scale,
      "_numobs_", n_obs,
      "_rlk_", tolower(as.character(generate_rlk))
    )

    if (isTRUE(reuse_mhe)) {
      reuse_results <- list()
      all_reusable <- TRUE
      for (server in server_list) {
        conn_idx <- which(server_names == server)
        res <- .dsAgg(
          conns = datasources[conn_idx],
          expr = call("mheReuseContextDS",
                      context_id = context_id,
                      session_id = session_id)
        )
        if (is.list(res) && !is.null(res[[1]])) res <- res[[1]]
        reuse_results[[server]] <- res
        if (!isTRUE(res$reusable)) {
          all_reusable <- FALSE
          break
        }
      }

      if (all_reusable) {
        mhe_reused <- TRUE
        for (server in server_list) {
          transport_pks[[server]] <- reuse_results[[server]]$transport_pk
        }
        if (verbose) message("  MHE keys reused from cached context (skipped)")
      }
    }

    if (!mhe_reused) {
      # ---- Coin-Tossing CRP (distributed randomness) ----
      if (verbose) message("  [Key Setup] Coin-toss CRP: generating commitments (", length(server_list), " servers)...")

      commit_results <- DSI::datashield.aggregate(
        conns = datasources[server_list],
        expr = call("mheCoinTossCommitDS", session_id = session_id)
      )
      commitments <- character(length(server_list))
      for (i in seq_along(server_list)) {
        r <- commit_results[[server_list[i]]]
        if (is.list(r) && length(r) == 1 && is.list(r[[1]])) r <- r[[1]]
        commitments[i] <- r$commitment
      }

      reveal_results <- DSI::datashield.aggregate(
        conns = datasources[server_list],
        expr = call("mheCoinTossRevealDS", session_id = session_id)
      )
      contributions <- character(length(server_list))
      for (i in seq_along(server_list)) {
        r <- reveal_results[[server_list[i]]]
        if (is.list(r) && length(r) == 1 && is.list(r[[1]])) r <- r[[1]]
        contributions[i] <- r$contribution
      }

      for (i in seq_along(server_list)) {
        ci <- which(server_names == server_list[i])
        DSI::datashield.aggregate(
          conns = datasources[ci],
          expr = call("mheCoinTossDeriveCRPDS",
                      contributions = unname(contributions),
                      commitments = unname(commitments),
                      log_n = as.integer(log_n),
                      log_scale = as.integer(log_scale),
                      party_id = as.integer(i - 1L),
                      session_id = session_id)
        )
      }

      if (verbose) message("  [Key Setup] CRP derived from ", length(server_list), "-party coin-toss (commit-reveal-derive)")

      # ---- Key Setup (all from coin-toss CRP) ----
      if (verbose) message("  [Key Setup] Generating per-party secret keys + PK shares...")
      conn_idx <- which(server_names == server_list[1])
      result0 <- .dsAgg(
        conns = datasources[conn_idx],
        expr = call("mheInitDS",
                    party_id = 0L, from_storage = TRUE,
                    num_obs = as.integer(n_obs),
                    log_n = as.integer(log_n),
                    log_scale = as.integer(log_scale),
                    generate_rlk = generate_rlk,
                    session_id = session_id)
      )
      if (is.list(result0)) result0 <- result0[[1]]
      crp <- result0$crp
      gkg_seed <- result0$gkg_seed

      pk_shares <- list()
      gkg_shares <- list()
      rlk_r1_shares <- list()
      pk_shares[[server_list[1]]] <- result0$public_key_share
      gkg_shares[[server_list[1]]] <- result0$galois_key_shares
      transport_pks[[server_list[1]]] <- result0$transport_pk
      if (generate_rlk && !is.null(result0$rlk_round1_share)) {
        rlk_r1_shares[[server_list[1]]] <- result0$rlk_round1_share
      }

      # Send CRP, GKG seed, and party_id to all non-party0 servers (sequential
      # blob transfer, then parallel mheInitDS calls).
      other_servers <- server_list[-1]
      other_conn_idxs <- integer(length(other_servers))
      for (i in seq_along(other_servers)) {
        server <- other_servers[i]
        ci <- which(server_names == server)
        other_conn_idxs[i] <- ci
        .sendBlob(crp, "crp", ci)
        .sendBlob(gkg_seed, "gkg_seed", ci)
        .sendBlob(as.character(i), "party_id", ci)
      }

      # Launch mheInitDS on all non-party0 servers in parallel
      init_expr <- call("mheInitDS",
                        party_id = 0L,
                        from_storage = TRUE,
                        num_obs = as.integer(n_obs),
                        log_n = as.integer(log_n),
                        log_scale = as.integer(log_scale),
                        generate_rlk = generate_rlk,
                        session_id = session_id)

      init_results <- .dsAgg(
        conns = datasources[other_servers],
        expr = init_expr
      )

      for (i in seq_along(other_servers)) {
        server <- other_servers[i]
        result <- init_results[[server]]
        pk_shares[[server]] <- result$public_key_share
        gkg_shares[[server]] <- result$galois_key_shares
        transport_pks[[server]] <- result$transport_pk
        if (generate_rlk && !is.null(result$rlk_round1_share)) {
          rlk_r1_shares[[server]] <- result$rlk_round1_share
        }
      }

      if (verbose) message("  [Key Setup] PK shares collected from ", length(server_list), " servers")

      # Two-round RLK generation protocol (binomial/Poisson HE-Link only)
      if (generate_rlk && length(rlk_r1_shares) > 0) {
        if (verbose) message("  [Key Setup] Generating collective relinearization key (2-round protocol)...")

        # Round 1 aggregation: send all R1 shares to coordinator (party 0)
        combine_conn <- which(server_names == server_list[1])
        for (i in seq_along(server_list)) {
          .sendBlob(rlk_r1_shares[[server_list[i]]],
                    paste0("rlk_r1_", i - 1), combine_conn)
        }

        agg_r1_result <- .dsAgg(
          conns = datasources[combine_conn],
          expr = call("mheRLKAggregateR1DS",
                      from_storage = TRUE,
                      n_parties = as.integer(length(server_list)),
                      session_id = session_id)
        )
        if (is.list(agg_r1_result)) agg_r1_result <- agg_r1_result[[1]]
        agg_r1 <- agg_r1_result$aggregated_round1

        # Distribute aggregated R1 to ALL servers for parallel round 2
        for (srv in server_list) {
          srv_conn <- which(server_names == srv)
          .sendBlob(agg_r1, "rlk_agg_r1", srv_conn)
        }

        # Round 2: all servers generate R2 shares in PARALLEL
        # Party 0 has agg_r1 in session; others got it via blob
        all_r2 <- DSI::datashield.aggregate(
          conns = datasources[server_list],
          expr = call("mheRLKRound2DS",
                      from_storage = TRUE,
                      session_id = session_id)
        )
        rlk_r2_shares <- list()
        for (srv in server_list) {
          r2 <- all_r2[[srv]]
          if (is.list(r2) && !is.null(r2$rlk_round2_share))
            rlk_r2_shares[[srv]] <- r2$rlk_round2_share
          else if (!is.null(r2[[1]]$rlk_round2_share))
            rlk_r2_shares[[srv]] <- r2[[1]]$rlk_round2_share
        }

        # Store aggregated R1 and R2 shares on coordinator for mheCombineDS
        .sendBlob(agg_r1, "rlk_agg_r1", combine_conn)
        for (i in seq_along(server_list)) {
          .sendBlob(rlk_r2_shares[[server_list[i]]],
                    paste0("rlk_r2_", i - 1), combine_conn)
        }
      }

      # Store PK shares, CRP, GKG seed, and GKG shares on the combining server
      # via chunked blob storage. Cryptographic objects can be several MB each,
      # exceeding R's expression parser stack limit if passed as call arguments.
      conn_idx <- which(server_names == server_list[1])

      for (i in seq_along(server_list)) {
        .sendBlob(pk_shares[[server_list[i]]], paste0("pk_", i - 1), conn_idx)
      }
      .sendBlob(crp, "crp", conn_idx)
      .sendBlob(gkg_seed, "gkg_seed", conn_idx)

      n_gkg_shares <- length(gkg_shares[[server_list[1]]])
      for (i in seq_along(server_list)) {
        shares <- gkg_shares[[server_list[i]]]
        for (j in seq_along(shares)) {
          .sendBlob(shares[j], paste0("gkg_", i - 1, "_", j - 1), conn_idx)
        }
      }

      if (verbose) message("  [Key Setup] Combining PK shares + Galois keys on coordinator...")
      combined <- .dsAgg(
        conns = datasources[conn_idx],
        expr = call("mheCombineDS",
                    from_storage = TRUE,
                    n_parties = as.integer(length(server_list)),
                    n_gkg_shares = as.integer(n_gkg_shares),
                    num_obs = as.integer(n_obs),
                    log_n = as.integer(log_n),
                    log_scale = as.integer(log_scale),
                    session_id = session_id)
      )
      if (is.list(combined)) combined <- combined[[1]]
      cpk <- combined$collective_public_key
      galois_keys <- combined$galois_keys
      relin_key <- combined$relin_key  # Non-NULL when RLK was generated

      # Distribute CPK, Galois keys, and RLK to other servers (parallel blob send)
      other_srv <- server_list[-1]
      if (verbose) message("  [Key Setup] Distributing CPK + Galois keys to ", length(other_srv), " servers (parallel)...")
      # Send CPK to all non-party0 servers simultaneously (per-chunk parallel)
      .dsvert_adaptive_send(cpk, function(chunk_str, chunk_idx, n_chunks) {
        if (n_chunks == 1L) {
          DSI::datashield.aggregate(datasources[other_srv],
            call("mheStoreBlobDS", key="cpk", chunk=chunk_str, session_id=session_id))
        } else {
          DSI::datashield.aggregate(datasources[other_srv],
            call("mheStoreBlobDS", key="cpk", chunk=chunk_str,
                 chunk_index=chunk_idx, n_chunks=n_chunks, session_id=session_id))
        }
      })
      # Galois keys: parallel per-key across servers
      if (!is.null(galois_keys)) {
        for (gk_i in seq_along(galois_keys)) {
          .dsvert_adaptive_send(galois_keys[gk_i], function(chunk_str, chunk_idx, n_chunks) {
            if (n_chunks == 1L) {
              DSI::datashield.aggregate(datasources[other_srv],
                call("mheStoreBlobDS", key=paste0("gk_",gk_i-1), chunk=chunk_str, session_id=session_id))
            } else {
              DSI::datashield.aggregate(datasources[other_srv],
                call("mheStoreBlobDS", key=paste0("gk_",gk_i-1), chunk=chunk_str,
                     chunk_index=chunk_idx, n_chunks=n_chunks, session_id=session_id))
            }
          })
        }
      }
      # RLK: parallel
      if (!is.null(relin_key) && nzchar(relin_key)) {
        .dsvert_adaptive_send(relin_key, function(chunk_str, chunk_idx, n_chunks) {
          if (n_chunks == 1L) {
            DSI::datashield.aggregate(datasources[other_srv],
              call("mheStoreBlobDS", key="rk", chunk=chunk_str, session_id=session_id))
          } else {
            DSI::datashield.aggregate(datasources[other_srv],
              call("mheStoreBlobDS", key="rk", chunk=chunk_str,
                   chunk_index=chunk_idx, n_chunks=n_chunks, session_id=session_id))
          }
        })
      }
      # Store CPK on all non-party0 servers in parallel
      .dsAgg(
        conns = datasources[server_list[-1]],
        expr = call("mheStoreCPKDS", from_storage = TRUE,
                    session_id = session_id)
      )

      # Register context for future reuse (parallel across all servers)
      .dsAgg(
        conns = datasources[server_list],
        expr = call("mheRegisterContextDS",
                    context_id = context_id,
                    session_id = session_id)
      )
    }  # end if (!mhe_reused)

    # Distribute transport keys for share-wrapping and secure routing
    # (always needed - fresh transport keys even when MHE keys are reused)
    fusion_server <- server_list[1]  # Party 0
    fusion_conn_idx <- which(server_names == fusion_server)
    non_fusion_servers <- setdiff(server_list, fusion_server)

    tk_map <- list(fusion = transport_pks[[fusion_server]])
    for (s in server_list) tk_map[[s]] <- transport_pks[[s]]
    # Distribute transport keys to all servers in parallel (same key map)
    .dsAgg(
      conns = datasources[server_list],
      expr = call("mheStoreTransportKeysDS", transport_keys = tk_map,
                  session_id = session_id)
    )

    if (verbose) message(sprintf("  [Key Setup] Complete (%.1fs)", proc.time()[[3]] - t0_mhe))

  }

  # =========================================================================
  # Phase 1: Standardize features (and y for Gaussian) for fast BCD
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

    std_result <- .dsAgg(
      conns = datasources[conn_idx],
      expr = call("glmStandardizeDS",
                  data_name = data_name,
                  output_name = std_data,
                  x_vars = x_vars[[server]],
                  y_var = y_arg,
                  session_id = session_id)
    )
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
  # Phase 2: Encrypt y and distribute (only if non-label servers exist)
  # =========================================================================
  if (length(non_label_servers) > 0 && !use_k2_beaver) {
    if (verbose) message("\n[Phase 2] Encrypting response variable (n=", n_obs, ")...")
    t0_enc <- proc.time()[[3]]

    # Encrypt STANDARDIZED y (Gaussian) or raw y (non-Gaussian)
    enc_data <- if (standardize_y) std_data else data_name

    conn_idx <- which(server_names == y_server)
    enc_result <- .dsAgg(
      conns = datasources[conn_idx],
      expr = call("mheEncryptRawDS",
                  data_name = enc_data, y_var = y_var,
                  store_local = TRUE,
                  session_id = session_id)
    )
    if (is.list(enc_result)) enc_result <- enc_result[[1]]
    ct_y <- enc_result$encrypted_y
    if (verbose) message(sprintf("  [Encrypt] y encrypted on %s (%d chars)",
                                   y_server, nchar(ct_y)))
    for (server in non_label_servers) {
      conn_idx <- which(server_names == server)
      n_chunks_sent <- .dsvert_adaptive_send(ct_y, function(chunk_str, chunk_idx, n_chunks) {
        .dsAgg(
          conns = datasources[conn_idx],
          expr = call("mheStoreEncChunkDS",
                      col_index = 1L,
                      chunk_index = chunk_idx,
                      chunk = chunk_str,
                      session_id = session_id)
        )
      })
      .dsAgg(
        conns = datasources[conn_idx],
        expr = call("mheAssembleEncColumnDS",
                    col_index = 1L, n_chunks = as.integer(n_chunks_sent),
                    session_id = session_id)
      )
      if (verbose)
        message(sprintf("  [Encrypt] ct(y) -> %s (%d chunks, %.1fs)",
                        server, n_chunks_sent, proc.time()[[3]] - t0_enc))
    }
  }

  # =========================================================================
  # Return all state needed by the BCD loop and post-processing
  # =========================================================================
  list(
    transport_pks  = transport_pks,
    cpk            = cpk,
    x_means        = x_means,
    x_sds          = x_sds,
    y_mean         = y_mean,
    y_sd           = y_sd,
    std_data       = std_data,
    standardize_y  = standardize_y,
    .dsAgg         = .dsAgg,
    .sendBlob      = .sendBlob,
    .sendCTChunks  = .sendCTChunks,
    .sendWrappedShare = .sendWrappedShare
  )
}
