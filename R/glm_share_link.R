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

  # binomial: mu = sigmoid(eta) via the DIRECT share-domain Chebyshev sigmoid
  # (degree-29 Clenshaw = 29 Beaver rounds) instead of exp+recip (~85 rounds).
  # ~2.9x fewer rounds, reveal-free + dealer-free, Go-verified max abs 6.46e-6
  # (>> the ~1e-4 a GLM coefficient needs). Domain contract: eta in [-8, 8];
  # c_0 carries the +0.5 sigmoid baseline.
  do.call(.ring127_sigmoid_round_keyed,
          c(list(in_key = in_eta_key, out_key = out_mu_key), common))
  invisible(out_mu_key)
}
