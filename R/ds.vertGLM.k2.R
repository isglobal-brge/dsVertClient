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
  # Set sidecar host options if configured (for Docker/multi-host deployments)
  k2_hosts <- getOption("dsvert.k2_mpc_hosts", NULL)
  if (!is.null(k2_hosts)) {
    for (server in server_list) {
      if (server %in% names(k2_hosts)) {
        conn_idx <- which(server_names == server)
        .dsAgg(datasources[conn_idx],
          call("options", dsvert.k2_mpc_host = k2_hosts[[server]]))
      }
    }
  }

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

    peer_r <- prep_results[[peer]]
    .dsAgg(
      conns = datasources[conn_idx],
      expr = call("k2MpcStorePeerDS",
                  job_id = job_id,
                  peer_host = peer_r$host,
                  peer_port = peer_r$port,
                  peer_role = peer_r$role,
                  peer_manifest = peer_r$manifest_hash,
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

  # Sequential launch: label (background) → nonlabel (foreground)
  # Label starts sidecar in background (listening), returns immediately.
  # Nonlabel starts sidecar foreground (connects to label), blocks until done.
  # Then we poll label for its results.

  # Step 1: Label launches sidecar in background
  label_run <- .dsAgg(
    conns = datasources[coordinator_conn],
    expr = call("k2MpcRunDS", job_id = job_id, session_id = session_id)
  )
  if (is.list(label_run) && length(label_run) == 1) label_run <- label_run[[1]]
  if (verbose) message("  Label sidecar launched in background")

  # Step 2: Nonlabel runs sidecar (connects to label, blocks until complete)
  nl_run <- .dsAgg(
    conns = datasources[nl_conn],
    expr = call("k2MpcRunDS", job_id = job_id, session_id = session_id)
  )
  if (is.list(nl_run) && length(nl_run) == 1) nl_run <- nl_run[[1]]
  if (verbose)
    message(sprintf("  Nonlabel completed: %d iters, conv=%s",
                    nl_run$iterations, nl_run$converged))
  run_results[[nl]] <- nl_run

  # Step 3: Poll label for results (sidecar should be done by now)
  label_result <- .dsAgg(
    conns = datasources[coordinator_conn],
    expr = call("k2MpcResultDS", job_id = job_id, session_id = session_id)
  )
  if (is.list(label_result) && length(label_result) == 1) label_result <- label_result[[1]]
  if (verbose)
    message(sprintf("  Label completed: %d iters, conv=%s",
                    label_result$iterations, label_result$converged))
  run_results[[coordinator]] <- label_result

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
