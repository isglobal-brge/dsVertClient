#' @title K=2 MPC Backend Client Orchestration
#' @description Client-side orchestration for the dedicated two-party MPC
#'   backend (binomial and Poisson GLMs). The client acts as control plane
#'   only: prepares jobs, exchanges endpoint metadata, triggers execution,
#'   and collects safe outputs. The data plane is a direct authenticated
#'   sidecar channel between the two servers.
#'
#' @name k2-mpc-client
NULL

#' Internal: K=2 MPC GLM orchestration
#'
#' Called from ds.vertGLM when K=2 and family is binomial or poisson.
#'
#' @keywords internal
.ds.vertGLM_k2_mpc <- function(data_name, y_var, x_vars, y_server,
                                 family, lambda, max_iter, tol,
                                 log_n, log_scale, verbose,
                                 datasources, session_id, std_data,
                                 server_list, server_names,
                                 coordinator, coordinator_conn,
                                 non_label_servers, n_obs,
                                 x_means, x_sds, y_mean, y_sd,
                                 standardize_y) {

  nl <- non_label_servers[1]
  nl_conn <- which(server_names == nl)

  # =========================================================================
  # Step 1: Build canonical manifest for peer validation
  # =========================================================================
  manifest_parts <- c(
    paste("servers", paste(sort(server_list), collapse = ",")),
    paste("family", family),
    paste("n_obs", n_obs),
    paste("y_server", y_server),
    paste("lambda", lambda)
  )
  manifest <- digest::digest(paste(manifest_parts, collapse = "|"), algo = "sha256")
  if (verbose) message(sprintf("  Manifest hash: %s", substr(manifest, 1, 16)))

  # =========================================================================
  # Step 2: Prepare job on both servers
  # =========================================================================
  if (verbose) message("\n[Phase 3a] Preparing K=2 MPC job on both servers...")

  job_id <- session_id  # reuse session_id as job_id

  .dsAgg <- function(conns, expr, ...) {
    DSI::datashield.aggregate(conns, expr, ...)
  }

  prep_results <- list()
  for (server in server_list) {
    conn_idx <- which(server_names == server)
    role <- if (server == coordinator) "label" else "nonlabel"

    r <- .dsAgg(
      conns = datasources[conn_idx],
      expr = call("k2MpcPrepareDS",
                  data_name = std_data,
                  y_var = if (role == "label") y_var else NULL,
                  x_vars = x_vars[[server]],
                  family = family,
                  role = role,
                  job_id = job_id,
                  manifest = manifest,
                  lambda = lambda,
                  max_iter = as.integer(max_iter),
                  tol = tol,
                  frac_bits = 20L,
                  session_id = session_id)
    )
    if (is.list(r) && length(r) == 1) r <- r[[1]]
    prep_results[[server]] <- r

    if (verbose)
      message(sprintf("  %s (%s): port %d, n=%d, p=%d",
                      server, role, r$port, r$n_obs, r$p_local))
  }

  # =========================================================================
  # Step 3: Exchange endpoint metadata between servers
  # =========================================================================
  if (verbose) message("\n[Phase 3b] Exchanging sidecar endpoints...")

  for (server in server_list) {
    conn_idx <- which(server_names == server)
    peer <- setdiff(server_list, server)

    .dsAgg(
      conns = datasources[conn_idx],
      expr = call("k2MpcStorePeerDS",
                  job_id = job_id,
                  peer_info = prep_results[[peer]],
                  session_id = session_id)
    )
  }
  if (verbose) message("  Peer endpoints exchanged and validated")

  # =========================================================================
  # Step 4: Launch MPC training on both servers
  # =========================================================================
  if (verbose)
    message(sprintf("\n[Phase 3c] Running MPC %s training (max %d iterations)...",
                    family, max_iter))

  # Both servers need to run simultaneously (the MPC protocol is interactive).
  # In DataSHIELD, we can't run two blocking calls in parallel on different
  # servers. The workaround: launch both using async dispatch.
  #
  # For non-Opal backends (DSLite), this won't work. For real Opals,
  # we use opalr's async command submission.

  run_results <- list()

  # Check if connections support async
  is_opal <- all(sapply(datasources, function(c) inherits(c, "OpalConnection")))

  if (is_opal) {
    # Async dispatch: submit both, then poll for completion
    cmd_ids <- list()
    for (server in server_list) {
      conn_idx <- which(server_names == server)
      conn <- datasources[[conn_idx]]
      opal <- conn@opal

      script <- paste0('dsVert::k2MpcRunDS(job_id = "', job_id,
                       '", session_id = "', session_id, '")')

      cmd_id <- opalr::opal.post(opal, "datashield", "session",
                                  opal$rid, "aggregate",
                                  query = list(async = "true"),
                                  body = script,
                                  contentType = "application/x-rscript")
      cmd_ids[[server]] <- list(opal = opal, cmd_id = cmd_id$id)
      if (verbose) message(sprintf("  %s: sidecar launched (async cmd %s)", server, cmd_id$id))
    }

    # Poll for completion
    for (server in server_list) {
      opal <- cmd_ids[[server]]$opal
      cid <- cmd_ids[[server]]$cmd_id

      if (verbose) message(sprintf("  Waiting for %s...", server))

      # Poll with timeout
      timeout <- max_iter * 60  # seconds
      t0 <- proc.time()[["elapsed"]]
      repeat {
        status <- opalr::opal.get(opal, "datashield", "session",
                                   opal$rid, "command", cid)
        if (status$status == "COMPLETED") break
        if (status$status == "FAILED") {
          stop(sprintf("MPC training failed on %s: %s", server,
                       status$error %||% "unknown error"), call. = FALSE)
        }
        if (proc.time()[["elapsed"]] - t0 > timeout) {
          stop(sprintf("MPC training timed out on %s after %ds", server, timeout),
               call. = FALSE)
        }
        Sys.sleep(2)
      }

      # Get result
      result <- opalr::opal.get(opal, "datashield", "session",
                                 opal$rid, "command", cid, "result")
      run_results[[server]] <- result
      if (verbose)
        message(sprintf("  %s: completed (%d iters, conv=%s, %.0fs)",
                        server, result$iterations, result$converged,
                        result$runtime_seconds))
    }
  } else {
    # Synchronous fallback (for DSLite testing — won't work for real MPC)
    warning("K=2 MPC requires Opal connections for async sidecar dispatch. ",
            "DSLite cannot run interactive 2-party protocols.", call. = FALSE)

    for (server in server_list) {
      conn_idx <- which(server_names == server)
      r <- .dsAgg(
        conns = datasources[conn_idx],
        expr = call("k2MpcRunDS", job_id = job_id, session_id = session_id)
      )
      if (is.list(r) && length(r) == 1) r <- r[[1]]
      run_results[[server]] <- r
    }
  }

  # =========================================================================
  # Step 5: Assemble final model
  # =========================================================================
  if (verbose) message("\n[Phase 4] Assembling model...")

  label_result <- run_results[[coordinator]]
  nl_result <- run_results[[nl]]

  # Collect standardized coefficients
  all_coefs_std <- numeric()
  all_names <- character()
  beta_0 <- 0

  for (server in server_list) {
    r <- run_results[[server]]
    betas_server <- r$beta

    if (server == coordinator && !is.null(r$intercept)) {
      beta_0 <- r$intercept
    }

    all_coefs_std <- c(all_coefs_std, betas_server)
    all_names <- c(all_names, x_vars[[server]])
  }

  # Un-standardize
  all_x_means <- numeric()
  all_x_sds <- numeric()
  for (server in server_list) {
    all_x_means <- c(all_x_means, x_means[[server]])
    all_x_sds <- c(all_x_sds, x_sds[[server]])
  }

  if (standardize_y && !is.null(y_sd)) {
    all_coefs_orig <- all_coefs_std * y_sd / all_x_sds
    intercept <- y_mean - sum(all_coefs_orig * all_x_means)
  } else {
    all_coefs_orig <- all_coefs_std / all_x_sds
    intercept <- beta_0 - sum(all_coefs_orig * all_x_means)
  }

  all_coefs <- c(intercept, all_coefs_orig)
  names(all_coefs) <- c("(Intercept)", all_names)

  # =========================================================================
  # Step 6: Cleanup
  # =========================================================================
  for (server in server_list) {
    conn_idx <- which(server_names == server)
    tryCatch(
      .dsAgg(datasources[conn_idx],
             expr = call("k2MpcCleanupDS", job_id = job_id,
                         session_id = session_id)),
      error = function(e) NULL
    )
  }

  # =========================================================================
  # Build result object
  # =========================================================================
  converged <- label_result$converged && nl_result$converged
  iterations <- max(label_result$iterations, nl_result$iterations)

  deviance <- if (!is.null(label_result$deviance)) label_result$deviance else NA
  null_deviance <- if (!is.null(label_result$null_deviance)) label_result$null_deviance else NA

  structure(
    list(
      coefficients = all_coefs,
      converged = converged,
      iterations = iterations,
      deviance = deviance,
      null_deviance = null_deviance,
      pseudo_r2 = if (!is.na(deviance) && !is.na(null_deviance))
        1 - deviance / null_deviance else NA,
      n_obs = n_obs,
      family = family,
      eta_privacy = "k2_mpc",
      aic = if (!is.na(deviance)) deviance + 2 * length(all_coefs) else NA,
      call = match.call()
    ),
    class = "ds.glm"
  )
}
