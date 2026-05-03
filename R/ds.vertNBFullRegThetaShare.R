#' @title Non-disclosive K=2 share-domain NB full-reg theta MLE orchestrator
#' @description Internal orchestration helper invoked by
#'   \code{ds.vertNBFullRegTheta(variant = "full_reg_nd")}. Per Newton-theta
#'   iter, drives the share-domain pipeline that closes D-INV-4 (no
#'   per-patient eta^nl reveal at label) by keeping eta, mu, log(mu+theta),
#'   1/(theta+mu) and (y+theta)*1/(theta+mu) all in Ring127 additive secret shares
#'   end-to-end through Beaver vecmul + AffineCombine + Chebyshev-
#'   Clenshaw primitives.
#'
#'   Refs: Lawless 1987 (NB profile-MLE theta score); Venables-Ripley 2002
#'   Sec.7.4 (\code{glm.nb} Newton); Catrina-Saxena 2010 Sec.3.3 (multiplicative
#'   depth ULP); Beaver 1991 (precomputed multiplication triples);
#'   Demmler-Schneider-Zohner ABY 2015 Sec.III.B (K=2 OT-Beaver
#'   dishonest-majority threat model); Trefethen ATAP Sec.8 (Bernstein-
#'   ellipse Chebyshev rel error).
#' @keywords internal
#' @noRd

.nb_fullreg_nd_session_setup <- function(formula, data, base_fit,
                                          datasources, server_names,
                                          y_srv, nl_srv,
                                          y_ci, nl_ci,
                                          x_label, x_nl,
                                          beta_label, beta_nl, int_val,
                                          y_var_char, session_id,
                                          .dsAgg, .sendBlob, verbose = FALSE) {
  # Initialise Ed25519 transport on label so it can receive NL's eta^nl share
  # blob (and later any blob NL relays during Newton iters).
  init_y <- .dsAgg(datasources[y_ci],
    call(name = "glmRing63TransportInitDS", session_id = session_id))
  if (is.list(init_y) && length(init_y) == 1L) init_y <- init_y[[1L]]
  label_pk <- init_y$transport_pk

  init_nl <- .dsAgg(datasources[nl_ci],
    call(name = "glmRing63TransportInitDS", session_id = session_id))
  if (is.list(init_nl) && length(init_nl) == 1L) init_nl <- init_nl[[1L]]
  nl_pk <- init_nl$transport_pk
  transport_pks <- list()
  transport_pks[[y_srv]]  <- label_pk
  transport_pks[[nl_srv]] <- nl_pk

  # NL splits eta^nl into Ring127 additive shares; transports peer-share
  # blob to label.
  share_r <- .dsAgg(datasources[nl_ci],
    call(name = "dsvertNBEtaShareDS",
         data_name = data, x_vars = x_nl,
         beta_values = as.numeric(beta_nl),
         target_pk = label_pk, session_id = session_id))
  if (is.list(share_r) && length(share_r) == 1L) share_r <- share_r[[1L]]

  # Relay blob to label.
  blob_slot <- "nb_eta_nl_share_blob"
  .sendBlob(share_r$sealed, blob_slot, y_ci)

  # Label receives blob, decrypts NL's share, computes own eta_label + beta_0
  # plaintext, FP-encodes, adds via k2-fp-add -> label's Ring127 share of
  # eta_total. Y is cached at label for later psi(y+theta) computation.
  recv_r <- .dsAgg(datasources[y_ci],
    call(name = "dsvertNBEtaTotalReceiveDS",
         data_name = data, y_var = y_var_char,
         x_vars_label = x_label,
         beta_values_label = as.numeric(beta_label),
         beta_intercept = as.numeric(int_val),
         peer_eta_share_blob_key = blob_slot,
         session_id = session_id))
  if (is.list(recv_r) && length(recv_r) == 1L) recv_r <- recv_r[[1L]]

  # Sanity: NL's k2_nb_eta_share_fp slot is populated by dsvertNBEtaShareDS;
  # confirm via the helper.
  conf_r <- .dsAgg(datasources[nl_ci],
    call(name = "dsvertNBEtaShareConfirmDS", session_id = session_id))
  if (is.list(conf_r) && length(conf_r) == 1L) conf_r <- conf_r[[1L]]
  if (!isTRUE(conf_r$stored))
    stop("NL eta_total share confirmation failed", call. = FALSE)
  conf_l <- .dsAgg(datasources[y_ci],
    call(name = "dsvertNBEtaShareConfirmDS", session_id = session_id))
  if (is.list(conf_l) && length(conf_l) == 1L) conf_l <- conf_l[[1L]]
  if (!isTRUE(conf_l$stored))
    stop("Label eta_total share confirmation failed", call. = FALSE)

  if (isTRUE(verbose))
    message(sprintf("[NBFullRegND] session setup OK n=%d (D-INV-4 closed: eta^nl in shares)",
                     as.integer(recv_r$n)))

  list(transport_pks = transport_pks, n = as.integer(recv_r$n),
       label_pk = label_pk, nl_pk = nl_pk)
}

.nb_fullreg_nd_transport_setup <- function(datasources, server_names,
                                           server_list, session_id,
                                           .dsAgg) {
  transport_pks <- list()
  identity_info <- list()
  for (server in server_list) {
    ci <- which(server_names == server)
    r <- .dsAgg(datasources[ci],
      call(name = "glmRing63TransportInitDS", session_id = session_id))
    if (is.list(r) && length(r) == 1L) r <- r[[1L]]
    transport_pks[[server]] <- r$transport_pk
    if (!is.null(r$identity_pk)) {
      identity_info[[server]] <- list(
        identity_pk = r$identity_pk,
        signature = r$signature)
    }
  }

  pk_sorted <- transport_pks[sort(names(transport_pks))]
  id_sorted <- if (length(identity_info) > 0L)
    identity_info[sort(names(identity_info))] else NULL
  json_to_b64url <- function(x) {
    .to_b64url(gsub("\n", "", jsonlite::base64_enc(charToRaw(
      jsonlite::toJSON(x, auto_unbox = TRUE))), fixed = TRUE))
  }
  pk_b64 <- json_to_b64url(pk_sorted)
  id_b64 <- if (!is.null(id_sorted)) json_to_b64url(id_sorted) else ""

  for (server in server_list) {
    ci <- which(server_names == server)
    .dsAgg(datasources[ci],
      call(name = "mpcStoreTransportKeysDS",
           transport_keys_b64 = pk_b64,
           identity_info_b64 = id_b64,
           session_id = session_id))
  }
  transport_pks
}

.nb_fullreg_nd_session_setup_k3 <- function(data, x_vars, beta_all,
                                            y_var_char, y_srv,
                                            datasources, server_names,
                                            server_list, session_id,
                                            .dsAgg, .sendBlob,
                                            verbose = FALSE) {
  if (!exists(".k3_select_fusion_server", mode = "function")) {
    stop("K>=3 NB full_reg_nd requires the K>=3 DCF input-sharing helpers",
         call. = FALSE)
  }
  fusion_server <- .k3_select_fusion_server(server_list, y_srv, x_vars)
  dcf_parties <- c(y_srv, fusion_server)
  dcf_conns <- vapply(dcf_parties, function(s) which(server_names == s),
                      integer(1L))
  non_dcf_servers <- setdiff(server_list, dcf_parties)

  transport_pks <- .nb_fullreg_nd_transport_setup(
    datasources = datasources, server_names = server_names,
    server_list = server_list, session_id = session_id,
    .dsAgg = .dsAgg)

  # Share original, unstandardised X and y into Ring127 additive shares held
  # only by the two DCF parties. All non-DCF servers split their columns once:
  # one share goes to fusion, the complement to the outcome/coordinator party.
  for (server in dcf_parties) {
    ci <- which(server_names == server)
    peer_dcf <- setdiff(dcf_parties, server)
    srv_x <- x_vars[[server]]
    if (length(srv_x) == 0L) srv_x <- NULL
    r <- .dsAgg(datasources[ci], call(name = "k2ShareInputDS",
      data_name = data,
      x_vars = srv_x,
      y_var = if (identical(server, y_srv)) y_var_char else NULL,
      peer_pk = .to_b64url(transport_pks[[peer_dcf]]),
      ring = 127L,
      session_id = session_id))
    if (is.list(r) && length(r) == 1L) r <- r[[1L]]
    peer_ci <- which(server_names == peer_dcf)
    .sendBlob(r$encrypted_x_share, "k2_peer_x_share", peer_ci)
    if (!is.null(r$encrypted_y_share)) {
      .sendBlob(r$encrypted_y_share, "k2_peer_y_share", peer_ci)
    }
  }

  for (server in non_dcf_servers) {
    srv_x <- x_vars[[server]]
    if (length(srv_x) == 0L) next
    ci <- which(server_names == server)
    r <- .dsAgg(datasources[ci], call(name = "k2ShareInputDS",
      data_name = data,
      x_vars = srv_x,
      y_var = NULL,
      peer_pk = .to_b64url(transport_pks[[fusion_server]]),
      ring = 127L,
      session_id = session_id))
    if (is.list(r) && length(r) == 1L) r <- r[[1L]]
    .sendBlob(r$encrypted_x_share, paste0("k2_extra_x_share_", server),
              which(server_names == fusion_server))

    r2 <- .dsAgg(datasources[ci], call(name = "glmRing63ExportOwnShareDS",
      peer_pk = .to_b64url(transport_pks[[y_srv]]),
      session_id = session_id))
    if (is.list(r2) && length(r2) == 1L) r2 <- r2[[1L]]
    .sendBlob(r2$encrypted_own_share, paste0("k2_extra_x_share_", server),
              which(server_names == y_srv))
  }

  for (i in seq_along(dcf_parties)) {
    dcf_ci <- dcf_conns[[i]]
    peer_srv <- dcf_parties[[3L - i]]
    .dsAgg(datasources[dcf_ci], call(name = "k2ReceiveShareDS",
      peer_p = as.integer(length(x_vars[[peer_srv]])),
      session_id = session_id))
  }

  for (server in non_dcf_servers) {
    extra_p <- length(x_vars[[server]])
    if (extra_p == 0L) next
    for (dcf_ci in dcf_conns) {
      .dsAgg(datasources[dcf_ci], call(name = "glmRing63ReceiveExtraShareDS",
        extra_key = paste0("k2_extra_x_share_", server),
        extra_p = as.integer(extra_p),
        session_id = session_id))
    }
  }

  beta_map <- list()
  idx <- 1L
  if (length(x_vars[[y_srv]]) > 0L) {
    beta_map[[y_srv]] <- seq.int(idx, idx + length(x_vars[[y_srv]]) - 1L)
    idx <- idx + length(x_vars[[y_srv]])
  }
  for (srv in server_list) {
    if (identical(srv, y_srv)) next
    p_s <- length(x_vars[[srv]])
    if (p_s == 0L) next
    beta_map[[srv]] <- seq.int(idx, idx + p_s - 1L)
    idx <- idx + p_s
  }

  feature_order <- c(x_vars[[y_srv]],
                     unlist(x_vars[setdiff(server_list, y_srv)],
                            use.names = FALSE))
  beta_vec <- unname(as.numeric(beta_all[feature_order]))
  beta_slice <- function(srv) {
    idx_srv <- beta_map[[srv]]
    if (is.null(idx_srv) || length(idx_srv) == 0L) NULL else beta_vec[idx_srv]
  }
  int_val <- as.numeric(beta_all[["(Intercept)"]])

  for (server in dcf_parties) {
    ci <- which(server_names == server)
    is_coord <- identical(server, y_srv)
    if (is_coord) {
      b_coord <- beta_slice(y_srv)
      b_nl <- c(beta_slice(fusion_server))
      for (ns in non_dcf_servers) b_nl <- c(b_nl, beta_slice(ns))
    } else {
      b_coord <- beta_slice(y_srv)
      b_nl <- c()
      for (ns in non_dcf_servers) b_nl <- c(b_nl, beta_slice(ns))
      b_nl <- c(b_nl, beta_slice(server))
    }
    .dsAgg(datasources[ci], call(name = "k2ComputeEtaShareDS",
      beta_coord = b_coord,
      beta_nl = b_nl,
      intercept = if (is_coord) int_val else 0.0,
      is_coordinator = is_coord,
      output_key = "k2_nb_eta_share_fp",
      session_id = session_id))
  }

  n_r <- .dsAgg(datasources[which(server_names == y_srv)],
    call(name = "getObsCountDS", data))
  if (is.list(n_r) && length(n_r) == 1L) n_r <- n_r[[1L]]
  n_obs <- as.integer(n_r$n_obs)

  if (isTRUE(verbose)) {
    message(sprintf("[NBFullRegND] K>=3 Ring127 setup OK n=%d DCF=(%s,%s)",
                    n_obs, y_srv, fusion_server))
  }
  list(transport_pks = transport_pks,
       n = n_obs,
       dcf_parties = dcf_parties,
       fusion_server = fusion_server,
       dealer_servers = if (length(non_dcf_servers) > 0L)
         non_dcf_servers else fusion_server)
}


# Reveal a Ring127 scalar share-sum from both servers + decode to float.
# Each server returns its local k2-fp-sum of the share at `key`; client
# combines via k2-ring63-aggregate (which routes to ring127 modular add).
.nb_fullreg_nd_reveal_sum <- function(key, datasources, server_list,
                                      server_names, y_server, nl,
                                      session_id, .dsAgg) {
  per_srv <- list()
  for (server in server_list) {
    ci <- which(server_names == server)
    r <- .dsAgg(datasources[ci],
      call(name = "k2BeaverSumShareDS", source_key = key,
           frac_bits = 50L, ring = "ring127", session_id = session_id))
    if (is.list(r) && length(r) == 1L) r <- r[[1L]]
    per_srv[[server]] <- r
  }
  agg <- dsVert:::.callMpcTool("k2-ring63-aggregate", list(
    share_a = per_srv[[y_server]]$sum_share_fp,
    share_b = per_srv[[nl]]$sum_share_fp,
    frac_bits = 50L, ring = "ring127"))
  as.numeric(agg$values)[1L]
}


# Note: prior `.nb_fullreg_nd_log_scale` plaintext-rescale heuristic was
# removed when the NR-LOG path landed (commits 2026-04-29). The single-
# scale [1, 10] core extended-log primitive bottlenecked at MARGINAL
# 11.6x sigma-probe ratio because most NB operating-range elements fell
# outside [1, 10]. NR-LOG on wide-Chebyshev [0.1, 1000] seed + 5 NR
# iters drives rel to ULP regardless of operating range.

# Per-Newton-iter share-domain score evaluation. Returns scalar score +
# deriv built from share-revealed scalar aggregates per Lawless 1987.
.nb_fullreg_nd_score <- function(theta, n_obs,
                                  datasources, dealer_ci, server_list,
                                  server_names, y_server, nl,
                                  ci_os, ci_nl, transport_pks,
                                  session_id,
                                  .dsAgg, .sendBlob, verbose = FALSE,
                                  yt_from_y_share = FALSE) {
  theta <- as.numeric(theta)
  if (!is.finite(theta) || theta <= 0)
    return(list(score = NA_real_, deriv = NA_real_, n = n_obs))

  # === Step 1: mu_share via .ring127_exp_round_keyed_extended ===
  .ring127_exp_round_keyed_extended(
    in_key = "k2_nb_eta_share_fp",
    out_key = "k2_nb_mu_share_fp",
    n = n_obs,
    datasources = datasources, dealer_ci = dealer_ci,
    server_list = server_list, server_names = server_names,
    y_server = y_server, nl = nl,
    transport_pks = transport_pks, session_id = session_id,
    .dsAgg = .dsAgg, .sendBlob = .sendBlob)

  # === Step 2: (mu + theta)_share -- party-0 (label) absorbs scalar theta ===
  theta_fp_b64 <- .to_b64url(dsVert:::.callMpcTool("k2-float-to-fp", list(
    values = array(as.numeric(theta), dim = 1L),
    frac_bits = 50L, ring = "ring127"))$fp_data)
  for (server in server_list) {
    ci <- which(server_names == server)
    is_p0 <- (server == y_server)
    .dsAgg(datasources[ci], call(name = "k2Ring127AffineCombineDS",
      a_key = "k2_nb_mu_share_fp", b_key = NULL,
      sign_a = 1, sign_b = 0,
      public_const_fp = if (is_p0) theta_fp_b64 else NULL,
      is_party0 = is_p0,
      output_key = "k2_nb_mupt_share_fp",
      n = as.numeric(n_obs), session_id = session_id))
  }

  # === Step 3: log(mu + theta)_share via NR-LOG (Pugh 2004 Sec.3) ===
  # Wide-Chebyshev seed on [0.1, 1000] degree 60 + 5 NR iters of
  # y_{n+1} = y_n + x * exp(-y_n) - 1 (quadratic convergence,
  # Goldschmidt 1964). Drives rel error to ~7.8e-27 << ULP 2^-50,
  # eliminating the plaintext-rescale single-scale precision floor
  # that bottlenecked the previous extended-log path at MARGINAL 11.6x.
  .ring127_log_round_keyed_nr(
    in_key = "k2_nb_mupt_share_fp",
    out_key = "k2_nb_log_mupt_share_fp",
    n = n_obs,
    datasources = datasources, dealer_ci = dealer_ci,
    server_list = server_list, server_names = server_names,
    y_server = y_server, nl = nl,
    transport_pks = transport_pks, session_id = session_id,
    .dsAgg = .dsAgg, .sendBlob = .sendBlob)

  # === Step 4: 1/(theta + mu)_share via .ring127_recip_round_keyed ===
  .ring127_recip_round_keyed(
    in_key = "k2_nb_mupt_share_fp",
    out_key = "k2_nb_recip_mupt_share_fp",
    n = n_obs,
    datasources = datasources, dealer_ci = dealer_ci,
    server_list = server_list, server_names = server_names,
    y_server = y_server, nl = nl,
    transport_pks = transport_pks, session_id = session_id,
    .dsAgg = .dsAgg, .sendBlob = .sendBlob)

  # === Step 5: (y + theta) shares ===
  if (isTRUE(yt_from_y_share)) {
    # K>=3 setup already shares y between the two selected DCF parties via
    # k2ShareInputDS(ring=127). Add the public scalar theta to party 0's y
    # share; the peer share is unchanged. No patient-level y or y+theta
    # vector is transported in plaintext.
    theta_fp_b64_yt <- .to_b64url(dsVert:::.callMpcTool("k2-float-to-fp", list(
      values = array(as.numeric(theta), dim = 1L),
      frac_bits = 50L, ring = "ring127"))$fp_data)
    for (server in server_list) {
      ci <- which(server_names == server)
      is_p0 <- (server == y_server)
      .dsAgg(datasources[ci], call(name = "k2Ring127AffineCombineDS",
        a_key = "k2_y_share_fp", b_key = NULL,
        sign_a = 1, sign_b = 0,
        public_const_fp = if (is_p0) theta_fp_b64_yt else NULL,
        is_party0 = is_p0,
        output_key = "k2_nb_yt_share_fp",
        n = as.numeric(n_obs), session_id = session_id))
    }
  } else {
    # K=2 legacy setup does not keep a y share; label re-shares y+theta as a
    # uniform Ring127 additive split to the non-label DCF party.
    yt_share_r <- .dsAgg(datasources[ci_os],
      call(name = "dsvertNBYThetaShareDS",
           theta = as.numeric(theta), target_pk = transport_pks[[nl]],
           session_id = session_id))
    if (is.list(yt_share_r) && length(yt_share_r) == 1L)
      yt_share_r <- yt_share_r[[1L]]
    yt_blob_slot <- "nb_yt_share_blob"
    .sendBlob(yt_share_r$sealed, yt_blob_slot, ci_nl)
    .dsAgg(datasources[ci_nl],
      call(name = "dsvertNBYThetaShareReceiveDS",
           peer_yt_share_blob_key = yt_blob_slot, session_id = session_id))
  }

  # === Step 6: (y + theta) * 1/(theta + mu) share via Beaver vecmul ===
  .ring127_vecmul(
    x_key = "k2_nb_yt_share_fp",
    y_key = "k2_nb_recip_mupt_share_fp",
    output_key = "k2_nb_ypt_over_tmu_share_fp",
    n = n_obs,
    datasources = datasources, dealer_ci = dealer_ci,
    server_list = server_list, server_names = server_names,
    y_server = y_server, nl = nl,
    transport_pks = transport_pks, session_id = session_id,
    .dsAgg = .dsAgg, .sendBlob = .sendBlob)

  # === Step 7: 1/(theta + mu)^2 share via Beaver vecmul ===
  .ring127_vecmul(
    x_key = "k2_nb_recip_mupt_share_fp",
    y_key = "k2_nb_recip_mupt_share_fp",
    output_key = "k2_nb_recip2_mupt_share_fp",
    n = n_obs,
    datasources = datasources, dealer_ci = dealer_ci,
    server_list = server_list, server_names = server_names,
    y_server = y_server, nl = nl,
    transport_pks = transport_pks, session_id = session_id,
    .dsAgg = .dsAgg, .sendBlob = .sendBlob)

  # === Step 8: (y + theta) * 1/(theta + mu)^2 share via Beaver vecmul ===
  .ring127_vecmul(
    x_key = "k2_nb_yt_share_fp",
    y_key = "k2_nb_recip2_mupt_share_fp",
    output_key = "k2_nb_ypt_over_tmu2_share_fp",
    n = n_obs,
    datasources = datasources, dealer_ci = dealer_ci,
    server_list = server_list, server_names = server_names,
    y_server = y_server, nl = nl,
    transport_pks = transport_pks, session_id = session_id,
    .dsAgg = .dsAgg, .sendBlob = .sendBlob)

  # === Step 9: Reveal four scalar aggregates via cooperative open ===
  rs <- function(key) .nb_fullreg_nd_reveal_sum(
    key, datasources, server_list, server_names, y_server, nl,
    session_id, .dsAgg)
  sum_log_mupt      <- rs("k2_nb_log_mupt_share_fp")
  sum_inv_tmu       <- rs("k2_nb_recip_mupt_share_fp")
  sum_ypt_over_tmu  <- rs("k2_nb_ypt_over_tmu_share_fp")
  sum_ypt_over_tmu2 <- rs("k2_nb_ypt_over_tmu2_share_fp")

  # === Step 10: Sumpsi(y+theta), Sumpsi_1(y+theta) plaintext at label ===
  psi_r <- .dsAgg(datasources[ci_os],
    call(name = "dsvertNBPsiAggregateDS", theta = theta, session_id = session_id))
  if (is.list(psi_r) && length(psi_r) == 1L) psi_r <- psi_r[[1L]]
  sum_psi <- as.numeric(psi_r$sum_psi)
  sum_tri <- as.numeric(psi_r$sum_tri)
  n_ret   <- as.integer(psi_r$n)

  # === Score + Hessian per Lawless 1987 ===
  score_val <- sum_psi - n_ret * digamma(theta) +
               n_ret * log(theta) - sum_log_mupt +
               n_ret - sum_ypt_over_tmu
  deriv_val <- sum_tri - n_ret * trigamma(theta) +
               n_ret / theta - 2 * sum_inv_tmu + sum_ypt_over_tmu2

  if (isTRUE(verbose))
    message(sprintf(
      "[NBFullRegND] theta=%.4f  score=%+.4e  deriv=%+.4e  n=%d  Sumlog(mu+theta)=%.3f  Sum1/(theta+mu)=%.3f  Sum(y+theta)/(theta+mu)=%.3f  Sum(y+theta)/(theta+mu)^2=%.3f  Sumpsi=%.3f  Sumpsi_1=%.3f",
      theta, score_val, deriv_val, n_ret,
      sum_log_mupt, sum_inv_tmu, sum_ypt_over_tmu,
      sum_ypt_over_tmu2, sum_psi, sum_tri))

  list(score = score_val, deriv = deriv_val, n = n_ret,
       diagnostics = list(
         sum_log_mupt = sum_log_mupt,
         sum_inv_tmu = sum_inv_tmu,
         sum_ypt_over_tmu = sum_ypt_over_tmu,
         sum_ypt_over_tmu2 = sum_ypt_over_tmu2,
         sum_psi = sum_psi, sum_tri = sum_tri))
}
