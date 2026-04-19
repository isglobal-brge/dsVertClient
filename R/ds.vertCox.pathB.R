#' @title Cox Path B client orchestrator
#' @description
#'   Iterative Newton with Fisher(beta_k) via Beaver on session-live
#'   mu, G, S shares. See docs/error_bounds/cox_newton_onestep.md and
#'   docs/acceptance/path_b_targets.md for the acceptance gates that
#'   govern this code.
#'
#'   Pre-conditions (same as Newton one-step):
#'     - standardised X shared (k2ShareInputDS)
#'     - Cox times + delta broadcast (k2SetCoxTimesDS / k2ReceiveCoxMetaDS)
#'     - permutation applied (k2ApplyCoxPermutationDS)
#'     - one-step Newton prep has populated cox_n_Xc_<j>_fp per canonical j
#'     - one-step Newton fit has produced an initial beta_1 at caller
#'
#'   Per iter: 1 call to `.cox_score_round` (populates mu, G, S shares) +
#'     (p) Beaver vecmuls for X_j * mu
#'     +1 Beaver vecmul for 1/S^2
#'     + 2 * p(p+1)/2 Beaver vecmuls for Term1 (X*X then *muG)
#'     + 2 * p(p+1)/2 Beaver vecmuls for Term2 (T*T then *1/S^2)
#'     = p + 1 + 2*p(p+1) Beaver vecmuls.
#'   For p=5: 5 + 1 + 60 = 66 Beaver vecmuls per iter. Max 5 iters per
#'   P3 budget.
#'
#' @keywords internal
.ds_vertCox_path_b_fisher <- function(beta_std, datasources, server_names,
                                       server_list, y_server, nl, session_id,
                                       n_obs, p_total, transport_pks,
                                       .cox_score_round,
                                       .dsAgg, .sendBlob, verbose = FALSE) {
  dealer_ci <- which(server_names == nl)
  single <- function(r) if (is.list(r) && length(r) == 1L) r[[1L]] else r

  # ---- Run DCF pipeline at beta_k. Side-effect populates mu, G, 1/S,
  #      mu*G shares on both parties under canonical key names.
  neg_grad_k <- .cox_score_round(beta_std)

  # Helper: one Beaver vecmul of session slot x_key × y_key → output slot.
  # Returns nothing; output lives in `output_key` on both parties.
  .beaver_mul <- function(x_key, y_key, output_key) {
    tri <- single(.dsAgg(datasources[dealer_ci],
      call("k2BeaverVecmulGenTriplesDS",
           dcf0_pk = transport_pks[[y_server]],
           dcf1_pk = transport_pks[[nl]],
           n = as.integer(n_obs),
           session_id = session_id, frac_bits = 20L)))
    .sendBlob(tri$triple_blob_0, "k2_beaver_vecmul_triple",
              which(server_names == y_server))
    .sendBlob(tri$triple_blob_1, "k2_beaver_vecmul_triple",
              which(server_names == nl))
    for (s in server_list) {
      ci <- which(server_names == s)
      .dsAgg(datasources[ci],
        call("k2BeaverVecmulConsumeTripleDS", session_id = session_id))
    }
    r1 <- list()
    for (s in server_list) {
      ci <- which(server_names == s); peer <- setdiff(server_list, s)
      r <- single(.dsAgg(datasources[ci], call("k2BeaverVecmulR1DS",
        peer_pk = transport_pks[[peer]],
        x_key = x_key, y_key = y_key,
        n = as.integer(n_obs),
        session_id = session_id, frac_bits = 20L)))
      r1[[s]] <- r
    }
    .sendBlob(r1[[y_server]]$peer_blob, "k2_beaver_vecmul_peer_masked",
              which(server_names == nl))
    .sendBlob(r1[[nl]]$peer_blob, "k2_beaver_vecmul_peer_masked",
              which(server_names == y_server))
    for (s in server_list) {
      ci <- which(server_names == s)
      .dsAgg(datasources[ci], call("k2BeaverVecmulR2DS",
        is_party0 = (s == y_server),
        x_key = x_key, y_key = y_key,
        output_key = output_key,
        n = as.integer(n_obs),
        session_id = session_id, frac_bits = 20L))
    }
    invisible(NULL)
  }

  # Helper: aggregate scalar-share-FP from both parties into a plaintext scalar.
  .scalar_reveal <- function(input_key, weight_key = NULL) {
    per_srv <- list()
    for (s in server_list) {
      ci <- which(server_names == s)
      r <- single(.dsAgg(datasources[ci],
        call("dsvertCoxPathBScalarDS",
             input_key = input_key,
             weight_key = weight_key,
             session_id = session_id)))
      per_srv[[s]] <- r
    }
    agg <- dsVert:::.callMpcTool("k2-ring63-aggregate", list(
      share_a = per_srv[[y_server]]$scalar_share_fp,
      share_b = per_srv[[nl]]$scalar_share_fp,
      frac_bits = 20L))
    as.numeric(agg$values)[1L]
  }

  # ---- Step 1: T_j(i) = rev_cumsum(X_j * mu) for each canonical j.
  # We compute share of X_j * mu via Beaver, then apply strata-aware
  # reverse cumsum locally (linear op on shares).
  for (j in seq_len(p_total)) {
    Xkey <- sprintf("cox_n_Xc_%d_fp", j)    # from Newton prep
    tmp_key <- sprintf("cox_pb_Xmu_%d_fp", j)
    Tkey <- sprintf("cox_pb_T_%d_fp", j)
    .beaver_mul(Xkey, "k2_cox_mu_share_fp", tmp_key)
    for (s in server_list) {
      ci <- which(server_names == s)
      .dsAgg(datasources[ci], call("dsvertCoxPathBCumsumDS",
        input_key = tmp_key, output_key = Tkey,
        reverse = TRUE, session_id = session_id))
    }
  }

  # ---- Step 2: share of 1/S^2 via Beaver((1/S) * (1/S)).
  .beaver_mul("k2_cox_recip_S_share_fp", "k2_cox_recip_S_share_fp",
              "cox_pb_recipS2_fp")

  # ---- Step 3: per pair (j <= k), compute Fisher entries.
  Fisher <- matrix(0, p_total, p_total)
  for (j in seq_len(p_total)) {
    for (k in j:p_total) {
      Xj <- sprintf("cox_n_Xc_%d_fp", j)
      Xk <- sprintf("cox_n_Xc_%d_fp", k)
      Tj <- sprintf("cox_pb_T_%d_fp", j)
      Tk <- sprintf("cox_pb_T_%d_fp", k)

      # Term1 = sum_m X_mj X_mk (muG)_m
      .beaver_mul(Xj, Xk, "cox_pb_XX_fp")
      .beaver_mul("cox_pb_XX_fp", "k2_cox_mu_g_share_fp",
                  "cox_pb_XXmuG_fp")
      T1 <- .scalar_reveal("cox_pb_XXmuG_fp", weight_key = NULL)

      # Term2 = sum_i delta_i T_j(i) T_k(i) / S(t_i)^2
      .beaver_mul(Tj, Tk, "cox_pb_TT_fp")
      .beaver_mul("cox_pb_TT_fp", "cox_pb_recipS2_fp",
                  "cox_pb_TTS2_fp")
      T2 <- .scalar_reveal("cox_pb_TTS2_fp",
                            weight_key = "k2_cox_delta_fp")

      Fisher[j, k] <- T1 - T2
      if (k != j) Fisher[k, j] <- Fisher[j, k]
    }
  }
  Fisher <- 0.5 * (Fisher + t(Fisher))  # symmetrise (FP noise)

  # grad_total(beta_k) = -n * (neg_grad - lambda * beta) converted
  # back from the per-obs normalisation that .cox_score_round uses.
  grad_total <- as.numeric(-n_obs * (neg_grad_k - 0 * beta_std))
  # (we pass lambda=0 here because the Newton path handles ridge
  # separately at the outer level; the neg_grad already includes
  # whatever lambda was in effect.)

  list(fisher = Fisher, grad = grad_total, neg_grad = neg_grad_k)
}
