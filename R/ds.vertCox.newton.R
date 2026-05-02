#' @title One-step Newton Cox at beta=0 (client orchestrator)
#' @description
#'   Called by ds.vertCox with one_step_newton = TRUE.
#'   Pre-conditions (caller has already:)
#'     * standardised X via .glm_mpc_setup()
#'     * shared X (k2ShareInputDS / k2ReceiveShareDS)
#'     * registered Cox times (k2SetCoxTimesDS / k2ReceiveCoxMetaDS)
#'     * applied permutation (k2ApplyCoxPermutationDS)
#'   so ss$k2_x_share_fp, ss$k2_peer_x_share_fp, ss$k2_cox_delta_fp are
#'   all sorted by time on both parties.
#'
#'   Does:
#'     1. dsvertCoxNewtonPrepDS on both parties (local extract + cumsum).
#'     2. dsvertCoxNewtonGradDS on both parties -> aggregate grad vector.
#'     3. For each Fisher pair (j <= k), run two Beaver vecmul protocols
#'        (one for the X*X term, one for the S*S term), then scalar
#'        aggregate.
#'     4. Solve beta = solve(Fisher, grad). Returns std-scale beta.
#'
#' @keywords internal
.ds_vertCox_newton_one_step <- function(datasources, server_names,
                                         server_list, y_server, nl,
                                         session_id, n_obs, transport_pks,
                                         p_coord, p_nl,
                                         .dsAgg = NULL, .sendBlob = NULL,
                                         verbose = TRUE, ring = 63L) {
  if (is.null(.dsAgg))    .dsAgg    <- DSI::datashield.aggregate
  if (is.null(.sendBlob)) {
    stop(".sendBlob must be provided by caller (uses adaptive chunked ",
         "relay via mpcStoreBlobDS).", call. = FALSE)
  }
  ring <- as.integer(ring)
  if (!ring %in% c(63L, 127L)) stop("ring must be 63 or 127", call. = FALSE)
  ring_tag <- if (ring == 127L) "ring127" else "ring63"
  frac_bits <- if (ring == 127L) 50L else 20L
  dealer_ci <- which(server_names == nl)
  single <- function(r) if (is.list(r) && length(r) == 1L) r[[1L]] else r

  # ---- 1. Prep on each party (independent, no blob relay). ----
  # Canonical column order: [p_coord (y_server's covariates) | p_nl (peer's)].
  # Each party maps its local "own"/"peer" share matrices to this global
  # indexing so grad_j / Fisher_jk refer to the SAME column on both sides.
  p_total <- as.integer(p_coord) + as.integer(p_nl)
  for (s in server_list) {
    ci <- which(server_names == s)
    is_coord <- (s == y_server)
    single(.dsAgg(datasources[ci],
      call(name = "dsvertCoxNewtonPrepDS", session_id = session_id,
           is_coordinator = is_coord,
           p_coord = as.integer(p_coord),
           p_nl = as.integer(p_nl))))
  }
  if (verbose) {
    message(sprintf("[ds.vertCox] Newton prep done. n=%d, p=%d",
                    n_obs, p_total))
  }

  # ---- 2. grad(0) -- per-party scalar share p-vector, then aggregate. ----
  grad_per_srv <- list()
  for (s in server_list) {
    ci <- which(server_names == s)
    r <- single(.dsAgg(datasources[ci],
      call(name = "dsvertCoxNewtonGradDS", session_id = session_id)))
    grad_per_srv[[s]] <- r
  }
  grad_vec <- numeric(p_total)
  for (j in seq_len(p_total)) {
    key <- sprintf("grad_%d", j)
    agg <- dsVert:::.callMpcTool("k2-ring63-aggregate", list(
      share_a = grad_per_srv[[y_server]][[key]],
      share_b = grad_per_srv[[nl]][[key]],
      frac_bits = frac_bits, ring = ring_tag))
    grad_vec[j] <- as.numeric(agg$values)[1L]
  }
  if (verbose) {
    message(sprintf("[ds.vertCox] Newton grad(0) ||.||=%.4g",
                    sqrt(sum(grad_vec^2))))
    message("[ds.vertCox] Newton grad(0) values: ",
            paste(sprintf("%.3g", grad_vec), collapse = ", "))
  }

  # ---- 3. Fisher(0) via Beaver vecmul on (X_j, X_k) and (S_j, S_k). ----
  Fisher <- matrix(0, p_total, p_total)

  .one_beaver_product <- function(which_vec, weight_key, j, k) {
    # Returns the scalar of Sum (weight) * (A_j * A_k) where A = X or S.
    # Steps: load pair -> gen triple -> consume -> R1 -> R2 -> scalar.
    for (s in server_list) {
      ci <- which(server_names == s)
      .dsAgg(datasources[ci], call(name = "dsvertCoxNewtonLoadPairDS",
        j = as.integer(j), k = as.integer(k), which = which_vec,
        session_id = session_id))
    }
    tri <- single(.dsAgg(datasources[dealer_ci],
      call(name = "k2BeaverVecmulGenTriplesDS",
           dcf0_pk = transport_pks[[y_server]],
           dcf1_pk = transport_pks[[nl]],
           n = as.integer(n_obs),
           session_id = session_id, frac_bits = frac_bits,
           ring = ring)))
    .sendBlob(tri$triple_blob_0, "k2_beaver_vecmul_triple",
              which(server_names == y_server))
    .sendBlob(tri$triple_blob_1, "k2_beaver_vecmul_triple",
              which(server_names == nl))
    for (s in server_list) {
      ci <- which(server_names == s)
      .dsAgg(datasources[ci],
        call(name = "k2BeaverVecmulConsumeTripleDS", session_id = session_id))
    }
    r1 <- list()
    for (s in server_list) {
      ci <- which(server_names == s); peer <- setdiff(server_list, s)
      r <- single(.dsAgg(datasources[ci], call(name = "k2BeaverVecmulR1DS",
        peer_pk = transport_pks[[peer]],
        x_key = "cox_n_Beaver_A_fp",
        y_key = "cox_n_Beaver_B_fp",
        n = as.integer(n_obs),
        session_id = session_id, frac_bits = frac_bits,
        ring = ring)))
      r1[[s]] <- r
    }
    .sendBlob(r1[[y_server]]$peer_blob, "k2_beaver_vecmul_peer_masked",
              which(server_names == nl))
    .sendBlob(r1[[nl]]$peer_blob, "k2_beaver_vecmul_peer_masked",
              which(server_names == y_server))
    for (s in server_list) {
      ci <- which(server_names == s)
      .dsAgg(datasources[ci], call(name = "k2BeaverVecmulR2DS",
        is_party0 = (s == y_server),
        x_key = "cox_n_Beaver_A_fp",
        y_key = "cox_n_Beaver_B_fp",
        output_key = "cox_n_Beaver_Z_fp",
        n = as.integer(n_obs),
        session_id = session_id, frac_bits = frac_bits,
        ring = ring))
    }
    scalar_per_srv <- list()
    for (s in server_list) {
      ci <- which(server_names == s)
      r <- single(.dsAgg(datasources[ci],
        call(name = "dsvertCoxNewtonFisherScalarDS",
             weight_key = weight_key, session_id = session_id)))
      scalar_per_srv[[s]] <- r
    }
    agg <- dsVert:::.callMpcTool("k2-ring63-aggregate", list(
      share_a = scalar_per_srv[[y_server]]$scalar_share_fp,
      share_b = scalar_per_srv[[nl]]$scalar_share_fp,
      frac_bits = frac_bits, ring = ring_tag))
    as.numeric(agg$values)[1L]
  }

  n_pairs <- p_total * (p_total + 1L) %/% 2L
  if (verbose) {
    message(sprintf(
      "[ds.vertCox] Newton Fisher: %d unique pairs x 2 terms = %d Beaver rounds",
      n_pairs, 2L * n_pairs))
  }
  pair_idx <- 0L
  for (j in seq_len(p_total)) {
    for (k in j:p_total) {
      pair_idx <- pair_idx + 1L
      t0 <- proc.time()[[3L]]
      term1 <- .one_beaver_product("X", "W1cum", j, k)
      term2 <- .one_beaver_product("S", "W2", j, k)
      Fisher[j, k] <- term1 - term2
      if (k != j) Fisher[k, j] <- Fisher[j, k]
      if (verbose) {
        message(sprintf(
          "  pair %d/%d (%d,%d)  F=%.4g  (T1=%.4g, T2=%.4g, %.1fs)",
          pair_idx, n_pairs, j, k, Fisher[j, k], term1, term2,
          proc.time()[[3L]] - t0))
      }
    }
  }

  # Symmetrise defensively (FP quantisation asymmetry at ~1e-6 level).
  Fisher <- 0.5 * (Fisher + t(Fisher))

  # ---- 4. Solve beta = Fisher^{-1} grad. ----
  # Tiny ridge for stability; dominated by Fisher diagonal magnitude.
  ridge <- 1e-8 * max(abs(diag(Fisher)))
  Fisher_reg <- Fisher + diag(ridge, p_total)
  beta_std <- tryCatch(solve(Fisher_reg, grad_vec),
                        error = function(e) {
                          message("[ds.vertCox] Fisher near-singular, ",
                                  "falling back to ridge-stabilised solve")
                          solve(Fisher_reg + diag(1e-4 * max(abs(diag(Fisher))),
                                                   p_total), grad_vec)
                        })

  list(beta_std = beta_std, grad = grad_vec, fisher = Fisher,
       p_total = p_total, n_obs = n_obs)
}
