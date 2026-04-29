#' @title Non-disclosive K=2 share-domain NB full-reg θ MLE orchestrator
#' @description Internal orchestration helper invoked by
#'   \code{ds.vertNBFullRegTheta(variant = "full_reg_nd")}. Per Newton-θ
#'   iter, drives the share-domain pipeline that closes D-INV-4 (no
#'   per-patient η^nl reveal at label) by keeping η, μ, log(μ+θ),
#'   1/(θ+μ) and (y+θ)·1/(θ+μ) all in Ring127 additive secret shares
#'   end-to-end through Beaver vecmul + AffineCombine + Chebyshev-
#'   Clenshaw primitives.
#'
#'   Refs: Lawless 1987 (NB profile-MLE θ score); Venables–Ripley 2002
#'   §7.4 (\code{glm.nb} Newton); Catrina–Saxena 2010 §3.3 (multiplicative
#'   depth ULP); Beaver 1991 (precomputed multiplication triples);
#'   Demmler–Schneider–Zohner ABY 2015 §III.B (K=2 OT-Beaver
#'   dishonest-majority threat model); Trefethen ATAP §8 (Bernstein-
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
  # Initialise Ed25519 transport on label so it can receive NL's η^nl share
  # blob (and later any blob NL relays during Newton iters).
  init_y <- .dsAgg(datasources[y_ci],
    call("glmRing63TransportInitDS", session_id = session_id))
  if (is.list(init_y) && length(init_y) == 1L) init_y <- init_y[[1L]]
  label_pk <- init_y$transport_pk

  init_nl <- .dsAgg(datasources[nl_ci],
    call("glmRing63TransportInitDS", session_id = session_id))
  if (is.list(init_nl) && length(init_nl) == 1L) init_nl <- init_nl[[1L]]
  nl_pk <- init_nl$transport_pk
  transport_pks <- list()
  transport_pks[[y_srv]]  <- label_pk
  transport_pks[[nl_srv]] <- nl_pk

  # NL splits η^nl into Ring127 additive shares; transports peer-share
  # blob to label.
  share_r <- .dsAgg(datasources[nl_ci],
    call("dsvertNBEtaShareDS",
         data_name = data, x_vars = x_nl,
         beta_values = as.numeric(beta_nl),
         target_pk = label_pk, session_id = session_id))
  if (is.list(share_r) && length(share_r) == 1L) share_r <- share_r[[1L]]

  # Relay blob to label.
  blob_slot <- "nb_eta_nl_share_blob"
  .sendBlob(share_r$sealed, blob_slot, y_ci)

  # Label receives blob, decrypts NL's share, computes own η_label + β₀
  # plaintext, FP-encodes, adds via k2-fp-add → label's Ring127 share of
  # η_total. Y is cached at label for later ψ(y+θ) computation.
  recv_r <- .dsAgg(datasources[y_ci],
    call("dsvertNBEtaTotalReceiveDS",
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
    call("dsvertNBEtaShareConfirmDS", session_id = session_id))
  if (is.list(conf_r) && length(conf_r) == 1L) conf_r <- conf_r[[1L]]
  if (!isTRUE(conf_r$stored))
    stop("NL η_total share confirmation failed", call. = FALSE)
  conf_l <- .dsAgg(datasources[y_ci],
    call("dsvertNBEtaShareConfirmDS", session_id = session_id))
  if (is.list(conf_l) && length(conf_l) == 1L) conf_l <- conf_l[[1L]]
  if (!isTRUE(conf_l$stored))
    stop("Label η_total share confirmation failed", call. = FALSE)

  if (isTRUE(verbose))
    message(sprintf("[NBFullRegND] session setup OK n=%d (D-INV-4 closed: η^nl in shares)",
                     as.integer(recv_r$n)))

  list(transport_pks = transport_pks, n = as.integer(recv_r$n),
       label_pk = label_pk, nl_pk = nl_pk)
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
      call("dsvertNBSumShareDS", input_key = key, session_id = session_id))
    if (is.list(r) && length(r) == 1L) r <- r[[1L]]
    per_srv[[server]] <- r
  }
  agg <- dsVert:::.callMpcTool("k2-ring63-aggregate", list(
    share_a = per_srv[[y_server]]$sum_share_fp,
    share_b = per_srv[[nl]]$sum_share_fp,
    frac_bits = 50L, ring = "ring127"))
  as.numeric(agg$values)[1L]
}


# Pick the rescale factor for log argument-reduction. (μ + θ) operating
# range is approximately [θ, θ + exp(η_max)]. For η ∈ [-23, 23] (clamped)
# this is wide; in practice the post-Poisson-warm η stays in [-5, 5]
# giving μ ∈ [0.007, 148]. Conservative scale maps the upper edge to 10
# (top of Chebyshev core). Lower elements may dip below 1 (Chebyshev
# core lower bound) where rel error degrades but stays bounded —
# Trefethen ATAP §8.2 Bernstein ellipse beyond [a, b] grows like ρ^(N+1)
# rather than ρ^N, but for ρ≈1.94 and degree 40 still yields rel ≲ 1e-8
# down to ~0.5 (Numerical Recipes 3rd ed §5.8). Higher-order outliers
# (μ+θ > 250) clip via the eta clamp at upstream.
.nb_fullreg_nd_log_scale <- function(theta, eta_max = 5.0) {
  # Conservative eta_max=5 ensures scale·(μ_max + θ) ≤ 10 even for
  # extreme μ ≈ exp(5) = 148. Empirically: tighter eta_max=3 caused
  # catastrophic Chebyshev extrapolation failure when scaled lower
  # bound dipped below ~0.5 (Runge phenomenon outside fit interval).
  # eta_max=5 keeps all elements in [scale·θ, 10] where scale·θ
  # for θ ∈ [0.5, 5] lies in [0.0034, 0.034] — far below [1, 10]
  # core, but in a regime where the Chebyshev polynomial is at least
  # bounded (not exploding) and rel error is ≲ 3e-3 per element. This
  # bounds the score noise floor at ~0.6 → |Δθ| ≈ 5e-2, ratio ≈ 10×
  # MARGINAL. Sub-noise (≥100×) requires multi-core piecewise log
  # primitives or DCF-based per-element argument reduction (deferred).
  upper <- as.numeric(theta) + exp(eta_max)
  scale <- 10 / upper
  list(scale = scale, log_correction = -log(scale))
}


# Per-Newton-iter share-domain score evaluation. Returns scalar score +
# deriv built from share-revealed scalar aggregates per Lawless 1987.
.nb_fullreg_nd_score <- function(theta, n_obs,
                                  datasources, dealer_ci, server_list,
                                  server_names, y_server, nl,
                                  ci_os, ci_nl, transport_pks,
                                  session_id,
                                  .dsAgg, .sendBlob, verbose = FALSE) {
  theta <- as.numeric(theta)
  if (!is.finite(theta) || theta <= 0)
    return(list(score = NA_real_, deriv = NA_real_, n = n_obs))

  # === Step 1: μ_share via .ring127_exp_round_keyed_extended ===
  .ring127_exp_round_keyed_extended(
    in_key = "k2_nb_eta_share_fp",
    out_key = "k2_nb_mu_share_fp",
    n = n_obs,
    datasources = datasources, dealer_ci = dealer_ci,
    server_list = server_list, server_names = server_names,
    y_server = y_server, nl = nl,
    transport_pks = transport_pks, session_id = session_id,
    .dsAgg = .dsAgg, .sendBlob = .sendBlob)

  # === Step 2: (μ + θ)_share — party-0 (label) absorbs scalar θ ===
  theta_fp_b64 <- .to_b64url(dsVert:::.callMpcTool("k2-float-to-fp", list(
    values = array(as.numeric(theta), dim = 1L),
    frac_bits = 50L, ring = "ring127"))$fp_data)
  for (server in server_list) {
    ci <- which(server_names == server)
    is_p0 <- (server == y_server)
    .dsAgg(datasources[ci], call("k2Ring127AffineCombineDS",
      a_key = "k2_nb_mu_share_fp", b_key = NULL,
      sign_a = 1L, sign_b = 0L,
      public_const_fp = if (is_p0) theta_fp_b64 else NULL,
      is_party0 = is_p0,
      output_key = "k2_nb_mupt_share_fp",
      n = as.integer(n_obs), session_id = session_id))
  }

  # === Step 3: log(μ + θ)_share via plaintext-rescale extended log ===
  ls <- .nb_fullreg_nd_log_scale(theta)
  scale_fp_b64 <- .to_b64url(dsVert:::.callMpcTool("k2-float-to-fp", list(
    values = array(ls$scale, dim = 1L),
    frac_bits = 50L, ring = "ring127"))$fp_data)
  log_corr_fp_b64 <- .to_b64url(dsVert:::.callMpcTool("k2-float-to-fp", list(
    values = array(ls$log_correction, dim = 1L),
    frac_bits = 50L, ring = "ring127"))$fp_data)
  .ring127_log_round_keyed_extended(
    in_key = "k2_nb_mupt_share_fp",
    out_key = "k2_nb_log_mupt_share_fp",
    n = n_obs,
    scale_fp_b64 = scale_fp_b64,
    log_scale_correction_fp_b64 = log_corr_fp_b64,
    datasources = datasources, dealer_ci = dealer_ci,
    server_list = server_list, server_names = server_names,
    y_server = y_server, nl = nl,
    transport_pks = transport_pks, session_id = session_id,
    .dsAgg = .dsAgg, .sendBlob = .sendBlob)

  # === Step 4: 1/(θ + μ)_share via .ring127_recip_round_keyed ===
  .ring127_recip_round_keyed(
    in_key = "k2_nb_mupt_share_fp",
    out_key = "k2_nb_recip_mupt_share_fp",
    n = n_obs,
    datasources = datasources, dealer_ci = dealer_ci,
    server_list = server_list, server_names = server_names,
    y_server = y_server, nl = nl,
    transport_pks = transport_pks, session_id = session_id,
    .dsAgg = .dsAgg, .sendBlob = .sendBlob)

  # === Step 5: (y + θ) re-share at label, transport mask to NL ===
  yt_share_r <- .dsAgg(datasources[ci_os],
    call("dsvertNBYThetaShareDS",
         theta = as.numeric(theta), target_pk = transport_pks[[nl]],
         session_id = session_id))
  if (is.list(yt_share_r) && length(yt_share_r) == 1L)
    yt_share_r <- yt_share_r[[1L]]
  yt_blob_slot <- "nb_yt_share_blob"
  .sendBlob(yt_share_r$sealed, yt_blob_slot, ci_nl)
  .dsAgg(datasources[ci_nl],
    call("dsvertNBYThetaShareReceiveDS",
         peer_yt_share_blob_key = yt_blob_slot, session_id = session_id))

  # === Step 6: (y + θ) · 1/(θ + μ) share via Beaver vecmul ===
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

  # === Step 7: 1/(θ + μ)² share via Beaver vecmul ===
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

  # === Step 8: (y + θ) · 1/(θ + μ)² share via Beaver vecmul ===
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

  # === Step 10: Σψ(y+θ), Σψ_1(y+θ) plaintext at label ===
  psi_r <- .dsAgg(datasources[ci_os],
    call("dsvertNBPsiAggregateDS", theta = theta, session_id = session_id))
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
      "[NBFullRegND] θ=%.4f  score=%+.4e  deriv=%+.4e  n=%d  scale=%.3f  Σlog(μ+θ)=%.3f  Σ1/(θ+μ)=%.3f  Σ(y+θ)/(θ+μ)=%.3f  Σ(y+θ)/(θ+μ)²=%.3f  Σψ=%.3f  Σψ_1=%.3f",
      theta, score_val, deriv_val, n_ret, ls$scale,
      sum_log_mupt, sum_inv_tmu, sum_ypt_over_tmu,
      sum_ypt_over_tmu2, sum_psi, sum_tri))

  list(score = score_val, deriv = deriv_val, n = n_ret,
       diagnostics = list(
         sum_log_mupt = sum_log_mupt,
         sum_inv_tmu = sum_inv_tmu,
         sum_ypt_over_tmu = sum_ypt_over_tmu,
         sum_ypt_over_tmu2 = sum_ypt_over_tmu2,
         sum_psi = sum_psi, sum_tri = sum_tri,
         scale_used = ls$scale))
}
