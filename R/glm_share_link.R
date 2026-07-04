#' @title Share-domain (reveal-free) GLM link evaluation
#' @description Produces a secret-shared \code{mu = link^{-1}(eta)} directly on
#'   the additive shares, with NO reveal of any per-observation value — replacing
#'   the DCF wide-spline link (which relayed a masked \code{eta}, the F1/F1b
#'   leak). It reuses the exact share-domain Chebyshev primitives that the
#'   round-2-confirmed leak-free Newton family already ships
#'   (\code{.ring127_exp_round_keyed_extended}, \code{.ring127_recip_round_keyed}),
#'   mirroring \code{ds.vertNBFullRegThetaShare.R:.nb_fullreg_nd_score}.
#'
#'   Reads \code{ss$<in_eta_key>} and writes \code{ss$<out_mu_key>} — the same
#'   contract as \code{k2IdentityLinkDS} — so the downstream gradient / deviance
#'   / SE path is unchanged. It relays \code{dcf_masked} ZERO times.
#'
#'   Ring: the exp127/recip127 primitives are Ring127/frac50, so the incoming
#'   \code{eta} share MUST be Ring127-encoded (eta = X.beta is computed LOCALLY
#'   on each server from plaintext X and public beta, so its fixed-point encoding
#'   is a local choice — k2ComputeEtaShareDS must emit Ring127 for non-Gaussian
#'   families). The Gaussian identity link stays Ring63 and does not use this.
#'
#'   Sigmoid domain: for \code{eta in [-8,8]}, \code{denom = 1+exp(-eta)} lies in
#'   [1.0003, 2982] within the recip domain [1,3000]. Near-separable fits with
#'   |eta|>8 need a reveal-free share-domain clamp on eta first (TODO: via a
#'   secure comparison whose result stays shared — NOT a reveal).
#'
#' @keywords internal
.glm_share_link <- function(family, n, datasources, dealer_ci, server_list,
                            server_names, y_server, nl, transport_pks, session_id,
                            .dsAgg, .sendBlob,
                            in_eta_key = "k2_eta_share_fp",
                            out_mu_key = "secure_mu_share") {
  common <- list(
    n = n, datasources = datasources, dealer_ci = dealer_ci,
    server_list = server_list, server_names = server_names,
    y_server = y_server, nl = nl, transport_pks = transport_pks,
    session_id = session_id, .dsAgg = .dsAgg, .sendBlob = .sendBlob)

  if (identical(family, "poisson")) {
    # mu = exp(eta)  — one share-domain exp chain, no reveal.
    do.call(.ring127_exp_round_keyed_extended,
            c(list(in_key = in_eta_key, out_key = out_mu_key), common))
    return(invisible(out_mu_key))
  }

  # binomial: mu = sigmoid(eta) = 1 / (1 + exp(-eta)), all on shares.
  one_fp_b64 <- .to_b64url(dsVert:::.callMpcTool("k2-float-to-fp", list(
    values = array(1.0, dim = 1L), frac_bits = 50L, ring = "ring127"))$fp_data)

  # (a) neg_eta = -eta   (local affine on each server's share; no MPC round)
  for (server in server_list) {
    ci <- which(server_names == server)
    .dsAgg(datasources[ci], call(name = "k2Ring127AffineCombineDS",
      a_key = in_eta_key, b_key = NULL, sign_a = -1, sign_b = 0,
      public_const_fp = NULL, is_party0 = (server == y_server),
      output_key = "glm_neg_eta_share_fp",
      n = as.numeric(n), session_id = session_id))
  }
  # (b) exp(-eta)
  do.call(.ring127_exp_round_keyed_extended,
          c(list(in_key = "glm_neg_eta_share_fp", out_key = "glm_exp_neg_eta_fp"), common))
  # (c) denom = 1 + exp(-eta)   (party-0 = label server absorbs the public +1)
  for (server in server_list) {
    ci <- which(server_names == server)
    is_p0 <- (server == y_server)
    .dsAgg(datasources[ci], call(name = "k2Ring127AffineCombineDS",
      a_key = "glm_exp_neg_eta_fp", b_key = NULL, sign_a = 1, sign_b = 0,
      public_const_fp = if (is_p0) one_fp_b64 else NULL, is_party0 = is_p0,
      output_key = "glm_denom_share_fp",
      n = as.numeric(n), session_id = session_id))
  }
  # (d) mu = 1 / denom
  do.call(.ring127_recip_round_keyed,
          c(list(in_key = "glm_denom_share_fp", out_key = out_mu_key), common))
  invisible(out_mu_key)
}
