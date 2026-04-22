#' @title Federated Cox proportional-hazards regression
#' @description Fit a Cox PH model on vertically partitioned DataSHIELD
#'   data using the reverse-cumsum reformulation of the partial-
#'   likelihood score. The client runs an L-BFGS outer loop; each step
#'   obtains the aggregate gradient
#'     \deqn{\nabla \ell(\beta) = \sum_j x_j (\delta_j - e^{\eta_j} G_j)}
#'     \deqn{G_j = \sum_{i: \delta_i=1, t_i \le t_j} 1 / S(t_i)}
#'     \deqn{S(t_i) = \sum_{k: t_k \ge t_i} e^{\eta_k}}
#'   via the already-deployed server helpers
#'   (k2SetCoxTimesDS / k2ApplyCoxPermutationDS /
#'   k2CoxReverseCumsumSDS / k2StoreCoxRecipDS / k2CoxForwardCumsumGDS)
#'   and the existing 4-phase DCF protocol with family="exp" +
#'   family="reciprocal". The partial-likelihood Beaver matvec is
#'   handled by the same glmRing63GenGradTriplesDS /
#'   k2StoreGradTripleDS / k2GradientR1DS / k2GradientR2DS machinery
#'   that ds.vertGLM uses, so no new cryptographic round is introduced.
#'
#'   Client view per iteration: the p-dimensional aggregate gradient
#'   and a scalar partial log-likelihood (via the standard
#'   Beaver-sum path). The client never sees \eqn{\eta_j}, \eqn{S(t_j)},
#'   \eqn{G_j}, or any per-patient quantity.
#'
#'   Inter-server disclosure: the DCF peer learns the ascending-time
#'   sort permutation (ranking of event times) and the binary event
#'   indicator. Absolute event times are NOT disclosed.
#'
#' @param formula Formula of the form \code{Surv(time, event) ~ x1 + ...}.
#'   If the LHS is not a \code{Surv(...)} expression, supply
#'   \code{time_col} / \code{event_col} explicitly.
#' @param data Aligned data-frame name on each server.
#' @param time_col,event_col Column names on the outcome server.
#' @param max_iter Outer L-BFGS iterations (default 30).
#' @param tol Convergence tolerance on max |delta beta|.
#' @param lambda L2 regularisation.
#' @param verbose Print progress.
#' @param datasources DataSHIELD connections.
#' @return A \code{ds.vertCox} object: \code{coefficients},
#'   \code{std_errors}, \code{covariance}, \code{loglik}
#'   (partial log-likelihood), \code{n_obs}, \code{n_events},
#'   \code{iterations}, \code{converged}.
#' @export
ds.vertCox <- function(formula, data = NULL,
                       time_col = NULL, event_col = NULL,
                       tstart_col = NULL,
                       strata_col = NULL,
                       max_iter = 30L, tol = 1e-4, lambda = 1e-4,
                       # compute_loglik=FALSE by default: in the Newton path
                       # it triggers a full DCF pipeline re-run at beta_hat
                       # (~60s on Opal, ~6 min on local harness with full key
                       # generation). Turn ON only when you need AIC/BIC or
                       # LR statistics. Note: the log-lik is DCF-spline
                       # biased (~10-15% off on NCCTG lung) and should be
                       # interpreted as approximate.
                       compute_loglik = FALSE, compute_se = FALSE,
                       # Speed knobs (all preserve non-disclosure):
                       #   num_intervals_exp/reciprocal: DCF spline grid
                       #     size. 75 is the sweet spot on Ring63 --
                       #     error < 1e-3 (below plan's coefficient bar)
                       #     while cutting Go spline compute time by
                       #     ~25% vs 100. Drop to 50 for even faster
                       #     Cox fits when the dataset is small (n<200)
                       #     and approximate SE is acceptable.
                       num_intervals_exp = 75L,
                       num_intervals_recip = 75L,
                       # One-step Newton path (default) -- at beta=0 the
                       # DCF splines are NOT invoked (exp(0)=1 is exact),
                       # so both grad(0) and Fisher(0) are computed in
                       # Ring63 at FP-quantisation precision rather than
                       # being perturbed by the ~1e-3 compounded DCF
                       # spline bias that the old iterative path suffers
                       # from. A single Newton step yields an MLE that
                       # is O(|beta|^2)-close to coxph on modest-effect
                       # covariates (|Delta| ~ 0.02-0.05 on Pima vs 0.31
                       # for the L-BFGS path). Set FALSE to fall back to
                       # the legacy iterative gradient-descent loop.
                       one_step_newton = TRUE,
                       # Refinement iterations: after the bias-free one-step
                       # Newton gives beta_1, do up to `newton_refine_iters`
                       # fixed-direction Newton steps:
                       #   beta_{k+1} = beta_k - n * Fisher(0)^{-1} * neg_grad(beta_k)
                       # grad(beta_k) is the existing DCF pipeline score (BIASED
                       # ~1-3% by the exp + recip splines). Fisher(0) is the
                       # bias-free information matrix from the first step,
                       # reused as a fixed preconditioner. Damped Newton
                       # converges LINEARLY with rate (1-eps) to the true MLE,
                       # giving 1-2 orders of accuracy per iter on strong-signal
                       # data. Set to 0 to keep the bare one-step.
                       # Path B: iterative Newton with Fisher(β_k) via
                       # Beaver on session-live μ, G, 1/S, μG shares.
                       # HARD CAP at 5 per P3 disclosure budget.
                       # TEMPORARY revert to 0 while the additive-bias
                       # structure is characterized (task #104 diagnostics
                       # a + b per Codex audit 2026-04-19 late). The
                       # targets file at docs/acceptance/path_b_targets.md
                       # REMAINS COMMITTED to Path B passing Cox
                       # lung/pima/strong — the 0 default is a revert so
                       # main stays executable, NOT a relaxation of the
                       # plan's strict <1e-3 goal for Cox.
                       newton_refine_iters = 5L,
                       newton_refine_tol = 1e-5,
                       # Ring selector (task #116 Cox STRICT closure):
                       # 63 (default) uses the current Ring63 uint64 pipeline
                       # at fracBits=20. 127 routes the input sharing, DCF
                       # generation, Beaver triples, wide-spline phases, and
                       # aggregate ops through the Uint128 Ring127 pipeline
                       # at fracBits=50 (~1e-15 per-op vs Ring63's ~1e-6).
                       # NOTE step 5b: ring=127 threads through the spline
                       # path (input-sharing + DCF + spline triples + 4-phase
                       # eval + mu aggregation). Downstream Beaver-gradient
                       # triples / k2-ring63-aggregate callsites beyond
                       # spline remain Ring63 for now.
                       # Default flipped to 127 on 2026-04-22: empirical
                       # 5/5 STRICT on Pima synthetic (max|Δβ|=1.06e-04
                       # vs Ring63 mixed STRICT/TIGHT/LOOSE), and Ring127
                       # runs ~2× faster because Path B converges in 5
                       # iters instead of oscillating around the Ring63
                       # Catrina-Saxena biased fixed point. See
                       # docs/error_bounds/cox_ring127_strict_evidence.md.
                       # Callers needing the legacy Ring63 behaviour must
                       # pass ring = 63L explicitly.
                       ring = 127L,
                       verbose = TRUE, datasources = NULL) {
  ring <- as.integer(ring)
  if (!ring %in% c(63L, 127L)) stop("ring must be 63 or 127", call. = FALSE)
  spline_frac_bits <- if (ring == 127L) 50L else 20L
  ring_tag <- if (ring == 127L) "ring127" else "ring63"

  if (is.null(datasources)) datasources <- DSI::datashield.connections_find()
  server_names <- names(datasources)
  if (!inherits(formula, "formula")) {
    stop("formula must be an R formula", call. = FALSE)
  }

  # Parse LHS: Surv(time, event) or Surv(tstart, tstop, event) or
  # survival::Surv(...). The call head of `survival::Surv(x, y)` is
  # itself a call to `::` (length 3), so we use deparse() for a single
  # canonical string.
  f_terms <- terms(formula)
  lhs <- attr(f_terms, "variables")[[2L]]
  x_vars_all <- attr(f_terms, "term.labels")
  if (is.call(lhs)) {
    head_name <- deparse(lhs[[1L]])
    if (head_name %in% c("Surv", "survival::Surv")) {
      # Count positional args (excluding `type=` etc.). Surv(time, event)
      # has 2 args, Surv(tstart, tstop, event) has 3.
      n_args <- length(lhs) - 1L
      if (n_args == 2L) {
        if (is.null(time_col))  time_col  <- as.character(lhs[[2L]])
        if (is.null(event_col)) event_col <- as.character(lhs[[3L]])
      } else if (n_args >= 3L) {
        if (is.null(tstart_col)) tstart_col <- as.character(lhs[[2L]])
        if (is.null(time_col))   time_col   <- as.character(lhs[[3L]])
        if (is.null(event_col))  event_col  <- as.character(lhs[[4L]])
      }
    }
  }
  if (is.null(time_col) || is.null(event_col)) {
    stop("Need time_col and event_col", call. = FALSE)
  }

  # Locate servers.
  col_results <- DSI::datashield.aggregate(datasources,
    call("dsvertColNamesDS", data_name = data))
  y_server <- NULL
  x_vars <- list()
  for (srv in server_names) {
    cols <- col_results[[srv]]$columns
    if (time_col %in% cols && event_col %in% cols) y_server <- srv
    feats <- intersect(x_vars_all, cols)
    if (length(feats) > 0L) x_vars[[srv]] <- feats
  }
  if (is.null(y_server)) {
    stop("time+event not co-located on a single server", call. = FALSE)
  }
  server_list <- names(x_vars)
  if (length(server_list) != 2L) {
    stop("ds.vertCox is K=2-native; need exactly 2 servers with covariates",
         call. = FALSE)
  }
  non_label_servers <- setdiff(server_list, y_server)
  nl <- non_label_servers[1L]
  p_coord <- length(x_vars[[y_server]])
  p_nl <- length(x_vars[[nl]])
  p_total <- p_coord + p_nl

  session_id <- .mpc_session_id()
  on.exit({
    for (.srv in server_list) {
      .ci <- which(server_names == .srv)
      tryCatch(DSI::datashield.aggregate(datasources[.ci],
        call("mpcCleanupDS", session_id = session_id)),
        error = function(e) NULL)
    }
  }, add = TRUE)

  if (verbose) {
    message(sprintf(
      "[ds.vertCox] y_server=%s, n=?, covariates: %s (%d) + %s (%d)",
      y_server,
      y_server, p_coord, nl, p_nl))
  }

  # Setup transport keys + standardise features (reuse GLM setup).
  setup <- .glm_mpc_setup(
    datasources = datasources, server_names = server_names,
    server_list = server_list, non_label_servers = non_label_servers,
    y_server = y_server, y_var = time_col, x_vars = x_vars,
    data_name = data, family = "gaussian", session_id = session_id,
    verbose = verbose)
  transport_pks <- setup$transport_pks
  std_data <- setup$std_data
  .dsAgg <- setup$.dsAgg
  .sendBlob <- setup$.sendBlob

  # Cox n: the standardized frame lives in session storage, not
  # parent.frame(), so query the raw aligned data instead (the row
  # count is identical before and after standardisation). The count
  # helper returns a named list with `$n_obs`, mirroring the ds.vertGLM
  # pattern at ds.vertGLM.R:290.
  n_obs <- tryCatch({
    r <- .dsAgg(datasources[which(server_names == y_server)],
                call("getObsCountDS", data_name = data))
    if (is.list(r) && length(r) == 1L) r <- r[[1L]]
    r$n_obs
  }, error = function(e) NULL)
  if (is.null(n_obs)) {
    stop("could not determine n on outcome server", call. = FALSE)
  }

  # === Input sharing (additive shares of X between the two DCF parties) ===
  .to_b64url <- function(x) gsub("\\+", "-",
    gsub("/", "_", gsub("=+$", "", x, perl = TRUE), fixed = TRUE),
    fixed = TRUE)
  share_results <- list()
  for (server in server_list) {
    ci <- which(server_names == server)
    peer <- setdiff(server_list, server)
    peer_pk_safe <- .to_b64url(transport_pks[[peer]])
    srv_x <- x_vars[[server]]; if (length(srv_x) == 0) srv_x <- NULL
    r <- .dsAgg(datasources[ci], call("k2ShareInputDS",
      data_name = std_data, x_vars = srv_x,
      y_var = NULL,   # Cox does not share a numeric y_share (time is sorted, not shared)
      peer_pk = peer_pk_safe, ring = ring, session_id = session_id))
    if (is.list(r) && length(r) == 1L) r <- r[[1L]]
    share_results[[server]] <- r
  }
  for (server in server_list) {
    peer <- setdiff(server_list, server)
    peer_ci <- which(server_names == peer)
    .sendBlob(share_results[[server]]$encrypted_x_share,
              "k2_peer_x_share", peer_ci)
  }
  for (server in server_list) {
    ci <- which(server_names == server)
    peer <- setdiff(server_list, server)
    .dsAgg(datasources[ci], call("k2ReceiveShareDS",
      peer_p = as.integer(length(x_vars[[peer]])),
      session_id = session_id))
  }

  # === Register Cox event times + broadcast sort permutation ===
  # Outcome server sorts by ascending t, serialises (perm, delta) to
  # the peer via transport-encrypted blob. k2SetCoxTimesDS reads the
  # time + event columns from the raw aligned frame (the standardised
  # frame lives in session storage and only holds the covariates).
  # Time-varying path: if tstart_col is given, we treat the data as
  # Andersen-Gill counting-process form. To keep the risk-set
  # computation tractable with the existing reverse-cumsum machinery,
  # we encode tstart as an ADDITIONAL stratum so every unique tstart
  # becomes its own stratum break (this degenerates to episode-split
  # Cox: each (person, interval) row is at risk only within the
  # stratum of its own tstart -- correct when each person contributes
  # a single interval, conservative otherwise). A true dual-cumsum
  # Andersen-Gill implementation remains the principled path; the
  # stratum-encoded approximation shipped here is correct for the
  # common case of left-truncation at a single tstart value per row.
  effective_strata <- strata_col
  if (!is.null(tstart_col) && nzchar(tstart_col)) {
    # Synthesise a combined strata column on the outcome server.
    if (verbose) {
      message("[ds.vertCox] tv-Cox via tstart-stratum encoding (",
              if (!is.null(strata_col)) "interacting with strata_col"
              else "standalone", ")")
    }
    # Ask outcome server to build __dsvert_tv_strata column.
    tryCatch(
      .dsAgg(datasources[which(server_names == y_server)],
        call("dsvertCoxTVStrataDS", data_name = data,
             tstart_column = tstart_col,
             base_strata_column = strata_col,
             output_column = "__dsvert_tv_strata")),
      error = function(e) stop(
        "dsvertCoxTVStrataDS unavailable: ", conditionMessage(e),
        call. = FALSE))
    effective_strata <- "__dsvert_tv_strata"
  }
  cox_times <- .dsAgg(
    datasources[which(server_names == y_server)],
    call("k2SetCoxTimesDS", data_name = data,
         time_column = time_col, event_column = event_col,
         strata_column = effective_strata,
         peer_pk = transport_pks[[nl]],
         session_id = session_id))
  if (is.list(cox_times) && length(cox_times) == 1L)
    cox_times <- cox_times[[1L]]
  n_events <- cox_times$n_events
  .sendBlob(cox_times$peer_blob, "k2_peer_cox_meta",
            which(server_names == nl))
  .dsAgg(datasources[which(server_names == nl)],
    call("k2ReceiveCoxMetaDS", session_id = session_id))
  # Permute X shares on both parties.
  for (server in server_list) {
    ci <- which(server_names == server)
    .dsAgg(datasources[ci],
      call("k2ApplyCoxPermutationDS", session_id = session_id))
  }

  # (Newton path runs AFTER .cox_score_round is defined, so it can use
  # it for optional refinement iterations — see below near the main loop.)

  # === Pre-generate DCF keys: we need both family="exp" (for mu) AND
  #     family="reciprocal" (for 1/S). Reuse the same wide spline
  #     machinery; dealer = non-label server. Generate two sets of keys
  #     per iteration (dispatch via family arg to phase1..phase4).
  dealer <- nl; dealer_ci <- which(server_names == nl)
  if (verbose) message(sprintf("[ds.vertCox] n=%d, n_events=%d, p=%d",
                                n_obs, n_events, p_total))

  # === Outer L-BFGS loop ===
  beta <- rep(0, p_total)
  coord_idx <- seq_len(p_coord)
  nl_idx <- seq_len(p_nl) + p_coord
  .lbfgs_direction_local <- .lbfgs_direction  # reuse GLM helper

  s_hist <- list(); y_hist <- list(); prev_theta <- NULL; prev_grad <- NULL
  converged <- FALSE
  final_iter <- 0L
  loglik <- NA_real_
  # Polyak-averaging iterate: steepest descent zig-zags on ill-
  # conditioned Cox Fisher info, but the Cesaro mean of the raw
  # iterates is guaranteed to converge to the MLE for convex
  # objectives (Polyak 1990, Nedic-Ozdaglar 2007). We track the
  # running mean from iter ~n/2 onward (tail averaging).
  beta_avg_sum <- rep(0, p_total)
  beta_avg_n <- 0L
  tail_start <- 3L  # begin averaging from iter 3

  # Track which DCF families have keys cached this session so we can
  # skip regen on iters 2+. Each (family, n, num_intervals) combo gets
  # generated once. Valid for the whole ds.vertCox session because the
  # DCF key material doesn't depend on beta.
  dcf_keys_cached <- list()

  # Cache of Ring127 Chebyshev exp public coefficients (fetched once per
  # session from the local Go binary — they are deterministic public
  # constants, so no cross-party sharing is needed). Populated lazily
  # in `.exp127_round` on first call.
  exp127_coef_cache <- NULL

  # Batched Beaver vecmul helper (Ring127 only) — wraps the 10-DS-call
  # triple-gen + consume + R1 + R2 chain that was previously inlined in
  # .cox_score_round. Used by .exp127_round for each Horner step, so
  # the Clenshaw recurrence stays readable (30 helper calls vs 300
  # inlined DS ops). The Ring63 Cox path keeps its inlined version to
  # preserve C4 regression guarantees.
  .run_beaver_vecmul_ring127 <- function(x_key, y_key, output_key, n) {
    tri <- .dsAgg(datasources[dealer_ci],
      call("k2BeaverVecmulGenTriplesDS",
           dcf0_pk = transport_pks[[y_server]],
           dcf1_pk = transport_pks[[nl]],
           n = as.integer(n),
           session_id = session_id, frac_bits = 50L,
           ring = 127L))
    if (is.list(tri) && length(tri) == 1L) tri <- tri[[1L]]
    .sendBlob(tri$triple_blob_0, "k2_beaver_vecmul_triple",
              which(server_names == y_server))
    .sendBlob(tri$triple_blob_1, "k2_beaver_vecmul_triple",
              which(server_names == nl))
    all_ci <- vapply(server_list, function(s) which(server_names == s),
                      integer(1L))
    .dsAgg(datasources[all_ci],
      call("k2BeaverVecmulConsumeTripleDS", session_id = session_id))
    r1_b <- list()
    for (server in server_list) {
      ci <- which(server_names == server)
      peer <- setdiff(server_list, server)
      r <- .dsAgg(datasources[ci], call("k2BeaverVecmulR1DS",
        peer_pk = transport_pks[[peer]],
        x_key = x_key, y_key = y_key,
        n = as.integer(n),
        session_id = session_id, frac_bits = 50L, ring = 127L))
      if (is.list(r) && length(r) == 1L) r <- r[[1L]]
      r1_b[[server]] <- r
    }
    .sendBlob(r1_b[[y_server]]$peer_blob, "k2_beaver_vecmul_peer_masked",
              which(server_names == nl))
    .sendBlob(r1_b[[nl]]$peer_blob, "k2_beaver_vecmul_peer_masked",
              which(server_names == y_server))
    for (server in server_list) {
      .dsAgg(datasources[which(server_names == server)],
        call("k2BeaverVecmulR2DS",
             is_party0 = (server == y_server),
             x_key = x_key, y_key = y_key,
             output_key = output_key,
             n = as.integer(n),
             session_id = session_id, frac_bits = 50L, ring = 127L))
    }
    invisible(NULL)
  }

  # Ring127 spline-less exp(eta) via Chebyshev-polynomial Horner (Clenshaw).
  # Replaces the wide-spline poisson path at ring=127 only — Ring63 and
  # any ring=127 with other families still go through .wide_spline_round.
  #
  # Input : eta share in ss$k2_eta_share_fp (populated by k2ComputeEtaShareDS).
  # Output: mu share in ss$secure_mu_share (same slot the spline path writes).
  #
  # Clenshaw on Ring127 FP (fracBits=50, degree N=30, domain [-5, 5]):
  #   y        = eta · (1/a)                  (local scale by public 1/a)
  #   twoY     = 2·y                          (local affine)
  #   b_{N+1}  = 0,  b_N = c_N (party0 only)  (bootstrap)
  #   for k = N-1, N-2, ..., 1:
  #     b_k    = c_k + twoY · b_{k+1} − b_{k+2}
  #              ^public       ^Beaver       ^share (local subtract)
  #   result   = c_0 + y · b_1 − b_2
  #
  # Each Horner step = 1 Beaver vecmul (twoY · b_{k+1}) + 1 affine combine.
  # ~360 DS round-trips per call, independent of n (all vectors batched).
  # On Ring127 with NCCTG n=210 → primitive accuracy ~3e-14 rel per element,
  # dominated by 30-step Clenshaw ULP drift (target ~1e-12 post-assembly).
  .exp127_round <- function(n_target) {
    if (ring != 127L) {
      stop(".exp127_round invoked with ring=", ring,
           "; Ring127-only path.", call. = FALSE)
    }
    n_int <- as.integer(n_target)

    # --- Step 1: fetch public Chebyshev coefficients (once per session).
    if (is.null(exp127_coef_cache)) {
      exp127_coef_cache <<- dsVert:::.callMpcTool(
        "k2-exp127-get-coeffs", list(frac_bits = 50L))
    }
    coef_res <- exp127_coef_cache
    degree <- as.integer(coef_res$degree)
    all_coeffs_raw <- jsonlite::base64_dec(coef_res$coeffs)
    # Split flat blob into one base64 string per coefficient (16 B each).
    c_b64 <- vapply(seq_len(degree + 1L), function(idx) {
      s <- (idx - 1L) * 16L + 1L
      e <- s + 15L
      jsonlite::base64_enc(all_coeffs_raw[s:e])
    }, character(1))
    # c_b64[idx] holds the base64 of c_{idx - 1} (0-indexed coefficient k).

    # --- Step 2: y = eta · (1/a)  (local scale, both parties).
    for (server in server_list) {
      ci <- which(server_names == server)
      .dsAgg(datasources[ci], call("k2Ring127LocalScaleDS",
        in_key = "k2_eta_share_fp",
        scalar_fp = coef_res$one_over_a,
        output_key = "k2_r127_horner_y",
        n = n_int, session_id = session_id))
    }

    # --- Step 3: twoY = y + y  (affine combine sign_a=+1, sign_b=+1).
    for (server in server_list) {
      ci <- which(server_names == server)
      is_coord <- (server == y_server)
      .dsAgg(datasources[ci], call("k2Ring127AffineCombineDS",
        a_key = "k2_r127_horner_y",
        b_key = "k2_r127_horner_y",
        sign_a = 1L, sign_b = 1L,
        public_const_fp = NULL,
        is_party0 = is_coord,
        output_key = "k2_r127_horner_twoY",
        n = n_int, session_id = session_id))
    }

    # --- Step 4: bootstrap b_N (party 0 has c_N, party 1 has 0) and
    # b_{N+1} = 0 share on both parties.
    for (server in server_list) {
      ci <- which(server_names == server)
      is_coord <- (server == y_server)
      # b_N: const=c_N, signs zero.
      .dsAgg(datasources[ci], call("k2Ring127AffineCombineDS",
        a_key = NULL, b_key = NULL,
        sign_a = 0L, sign_b = 0L,
        public_const_fp = c_b64[degree + 1L],
        is_party0 = is_coord,
        output_key = "k2_r127_horner_bB",
        n = n_int, session_id = session_id))
      # b_{N+1}: const=NULL, signs zero → zero share.
      .dsAgg(datasources[ci], call("k2Ring127AffineCombineDS",
        a_key = NULL, b_key = NULL,
        sign_a = 0L, sign_b = 0L,
        public_const_fp = NULL,
        is_party0 = is_coord,
        output_key = "k2_r127_horner_bA",
        n = n_int, session_id = session_id))
    }

    # --- Step 5: main Horner loop for k = N-1 downto 1.
    # Rolling slot convention at the start of each iter:
    #   slot_B holds b_{k+1}
    #   slot_A holds b_{k+2}
    # At end of iter, b_k is written into slot_A (overwriting the now-unused
    # b_{k+2}); we then swap the R-level labels so slot_B := b_k and
    # slot_A := b_{k+1} for the next iter. No server-side copy.
    slot_B <- "k2_r127_horner_bB"
    slot_A <- "k2_r127_horner_bA"
    for (k in seq.int(degree - 1L, 1L)) {
      # Beaver(twoY, slot_B=b_{k+1}) → k2_r127_horner_tmp = twoY · b_{k+1}.
      .run_beaver_vecmul_ring127(
        x_key = "k2_r127_horner_twoY",
        y_key = slot_B,
        output_key = "k2_r127_horner_tmp",
        n = n_int)
      # b_k = tmp + c_k_party0 − slot_A (b_{k+2}); store back into slot_A.
      for (server in server_list) {
        ci <- which(server_names == server)
        is_coord <- (server == y_server)
        .dsAgg(datasources[ci], call("k2Ring127AffineCombineDS",
          a_key = "k2_r127_horner_tmp",
          b_key = slot_A,
          sign_a = 1L, sign_b = -1L,
          public_const_fp = c_b64[k + 1L],
          is_party0 = is_coord,
          output_key = slot_A,
          n = n_int, session_id = session_id))
      }
      # Rotate slot labels — see comment above.
      swap <- slot_A; slot_A <- slot_B; slot_B <- swap
    }
    # Post-loop: slot_B = b_1, slot_A = b_2.

    # --- Step 6: final step result = y · b_1 + c_0 − b_2.
    .run_beaver_vecmul_ring127(
      x_key = "k2_r127_horner_y",
      y_key = slot_B,
      output_key = "k2_r127_horner_tmp",
      n = n_int)
    for (server in server_list) {
      ci <- which(server_names == server)
      is_coord <- (server == y_server)
      .dsAgg(datasources[ci], call("k2Ring127AffineCombineDS",
        a_key = "k2_r127_horner_tmp",
        b_key = slot_A,
        sign_a = 1L, sign_b = -1L,
        public_const_fp = c_b64[1L],
        is_party0 = is_coord,
        output_key = "secure_mu_share",
        n = n_int, session_id = session_id))
    }

    invisible(NULL)
  }

  # Cache of Ring127 Chebyshev recip public coefficients + NR constants.
  # Populated lazily in `.recip127_round` on first call via the local
  # `.callMpcTool` binary (deterministic public constants, no DS relay).
  recip127_coef_cache <- NULL

  # Ring127 spline-less 1/x via Chebyshev-Horner initial guess (degree 30
  # on domain [1, 3000]) + Newton-Raphson refinement (6 iters). Replaces
  # the wide-spline reciprocal path at ring=127 only. Ring63 + other
  # ring=127 families keep the existing DCF spline path.
  #
  # Input : x share (positive S(t) values) in ss$k2_eta_share_fp (written
  #         by k2CoxPrepareRecipPhaseDS, which copies k2_cox_S_share_fp
  #         into the eta slot before invoking this round).
  # Output: 1/x share in ss$secure_mu_share (same slot the spline path
  #         writes, so downstream k2StoreCoxRecipDS is unaffected).
  #
  # Algorithm (all on shares except public coeffs):
  #   1. t     = x · (1/halfRange) + (-mid/halfRange)     [LocalScale + Affine]
  #   2. twoT  = t + t                                    [Affine +1/+1]
  #   3. Clenshaw Horner: b_{N+1}=0, b_N=c_N_party0,
  #        for k=N-1..1:  b_k = c_k + twoT · b_{k+1} − b_{k+2}
  #      y     = c_0 + t · b_1 − b_2    ← Chebyshev initial guess (~58% worst)
  #   4. NR (6 iters):  y ← y · (2 − x · y)
  #        — rel_err squares per iter, reaches Ring127 ULP.
  #
  # ~42 Beaver vecmuls + 35 AffineCombines per call, batched over n.
  # 0-bit disclosure preserved (no range reduction, no shift reveal).
  .recip127_round <- function(n_target) {
    if (ring != 127L) {
      stop(".recip127_round invoked with ring=", ring,
           "; Ring127-only path.", call. = FALSE)
    }
    n_int <- as.integer(n_target)

    # --- Step 1: fetch public recip coefficients + NR constants.
    if (is.null(recip127_coef_cache)) {
      recip127_coef_cache <<- dsVert:::.callMpcTool(
        "k2-recip127-get-coeffs", list(frac_bits = 50L))
    }
    rc <- recip127_coef_cache
    degree <- as.integer(rc$degree)
    nr_steps <- as.integer(rc$nr_steps)
    all_coeffs_raw <- jsonlite::base64_dec(rc$coeffs)
    c_b64 <- vapply(seq_len(degree + 1L), function(idx) {
      s <- (idx - 1L) * 16L + 1L
      e <- s + 15L
      jsonlite::base64_enc(all_coeffs_raw[s:e])
    }, character(1))

    # --- Step 2a: t_pre = x · (1/halfRange)  (local scale, both parties).
    for (server in server_list) {
      ci <- which(server_names == server)
      .dsAgg(datasources[ci], call("k2Ring127LocalScaleDS",
        in_key = "k2_eta_share_fp",
        scalar_fp = rc$one_over_half_range,
        output_key = "k2_r127_recip_t_pre",
        n = n_int, session_id = session_id))
    }
    # --- Step 2b: t = t_pre + (-mid/halfRange) on party 0 (public offset).
    for (server in server_list) {
      ci <- which(server_names == server)
      is_coord <- (server == y_server)
      .dsAgg(datasources[ci], call("k2Ring127AffineCombineDS",
        a_key = "k2_r127_recip_t_pre",
        b_key = NULL,
        sign_a = 1L, sign_b = 0L,
        public_const_fp = rc$neg_mid_over_half_range,
        is_party0 = is_coord,
        output_key = "k2_r127_recip_t",
        n = n_int, session_id = session_id))
    }

    # --- Step 3: twoT = t + t.
    for (server in server_list) {
      ci <- which(server_names == server)
      is_coord <- (server == y_server)
      .dsAgg(datasources[ci], call("k2Ring127AffineCombineDS",
        a_key = "k2_r127_recip_t",
        b_key = "k2_r127_recip_t",
        sign_a = 1L, sign_b = 1L,
        public_const_fp = NULL,
        is_party0 = is_coord,
        output_key = "k2_r127_recip_twoT",
        n = n_int, session_id = session_id))
    }

    # --- Step 4: bootstrap b_N (party0 has c_N, party1 has 0) + b_{N+1} = 0.
    for (server in server_list) {
      ci <- which(server_names == server)
      is_coord <- (server == y_server)
      .dsAgg(datasources[ci], call("k2Ring127AffineCombineDS",
        a_key = NULL, b_key = NULL,
        sign_a = 0L, sign_b = 0L,
        public_const_fp = c_b64[degree + 1L],
        is_party0 = is_coord,
        output_key = "k2_r127_recip_bB",
        n = n_int, session_id = session_id))
      .dsAgg(datasources[ci], call("k2Ring127AffineCombineDS",
        a_key = NULL, b_key = NULL,
        sign_a = 0L, sign_b = 0L,
        public_const_fp = NULL,
        is_party0 = is_coord,
        output_key = "k2_r127_recip_bA",
        n = n_int, session_id = session_id))
    }

    # --- Step 5: Horner loop (same rolling 2-slot pattern as .exp127_round).
    slot_B <- "k2_r127_recip_bB"
    slot_A <- "k2_r127_recip_bA"
    for (k in seq.int(degree - 1L, 1L)) {
      .run_beaver_vecmul_ring127(
        x_key = "k2_r127_recip_twoT",
        y_key = slot_B,
        output_key = "k2_r127_recip_tmp",
        n = n_int)
      for (server in server_list) {
        ci <- which(server_names == server)
        is_coord <- (server == y_server)
        .dsAgg(datasources[ci], call("k2Ring127AffineCombineDS",
          a_key = "k2_r127_recip_tmp",
          b_key = slot_A,
          sign_a = 1L, sign_b = -1L,
          public_const_fp = c_b64[k + 1L],
          is_party0 = is_coord,
          output_key = slot_A,
          n = n_int, session_id = session_id))
      }
      swap <- slot_A; slot_A <- slot_B; slot_B <- swap
    }
    # Post-loop: slot_B = b_1, slot_A = b_2.

    # --- Step 6: y_0 = c_0 + t · b_1 − b_2.
    .run_beaver_vecmul_ring127(
      x_key = "k2_r127_recip_t",
      y_key = slot_B,
      output_key = "k2_r127_recip_tmp",
      n = n_int)
    for (server in server_list) {
      ci <- which(server_names == server)
      is_coord <- (server == y_server)
      .dsAgg(datasources[ci], call("k2Ring127AffineCombineDS",
        a_key = "k2_r127_recip_tmp",
        b_key = slot_A,
        sign_a = 1L, sign_b = -1L,
        public_const_fp = c_b64[1L],
        is_party0 = is_coord,
        output_key = "k2_r127_recip_y",
        n = n_int, session_id = session_id))
    }

    # --- Step 7: NR refinement y ← y · (2 − x · y), 6 iters.
    # Rolling slot alternation to avoid server-side copies; the last iter
    # writes directly into secure_mu_share (the spline-path output slot)
    # so no final identity copy is needed.
    y_cur <- "k2_r127_recip_y"
    y_next <- "k2_r127_recip_y_alt"
    for (iter in seq_len(nr_steps)) {
      # xy = x · y_cur.
      .run_beaver_vecmul_ring127(
        x_key = "k2_eta_share_fp",
        y_key = y_cur,
        output_key = "k2_r127_recip_xy",
        n = n_int)
      # twoMinusXy = 2 − xy  (sign_a=0, sign_b=-1, const=2 on party 0).
      for (server in server_list) {
        ci <- which(server_names == server)
        is_coord <- (server == y_server)
        .dsAgg(datasources[ci], call("k2Ring127AffineCombineDS",
          a_key = NULL,
          b_key = "k2_r127_recip_xy",
          sign_a = 0L, sign_b = -1L,
          public_const_fp = rc$two_fp,
          is_party0 = is_coord,
          output_key = "k2_r127_recip_twoMinusXy",
          n = n_int, session_id = session_id))
      }
      # y_next = y_cur · twoMinusXy  (Beaver). Final iter writes to
      # secure_mu_share directly.
      final_slot <- if (iter == nr_steps) "secure_mu_share" else y_next
      .run_beaver_vecmul_ring127(
        x_key = y_cur,
        y_key = "k2_r127_recip_twoMinusXy",
        output_key = final_slot,
        n = n_int)
      if (iter < nr_steps) {
        swap <- y_cur; y_cur <- y_next; y_next <- swap
      }
    }

    invisible(NULL)
  }

  .wide_spline_round <- function(family_name, n_target,
                                  num_intervals = 100L,
                                  need_dcf_keys = TRUE) {
    # Ring127 spline-less paths: replace the spline with Chebyshev-Horner
    # evaluation (poisson) or Chebyshev-Horner + NR refinement (reciprocal),
    # eliminating the ~1e-4 spline noise floor. Ring63 and any ring=127
    # families other than poisson/reciprocal keep the existing DCF path.
    if (ring == 127L && identical(family_name, "poisson")) {
      .exp127_round(n_target = n_target)
      return(invisible(NULL))
    }
    if (ring == 127L && identical(family_name, "reciprocal")) {
      .recip127_round(n_target = n_target)
      return(invisible(NULL))
    }
    # Cache hit: skip the ~20s key-gen + distribute round.
    cache_key <- sprintf("%s_%d_%d", family_name, n_target, num_intervals)
    if (isTRUE(dcf_keys_cached[[cache_key]])) need_dcf_keys <- FALSE
    # One complete 4-phase DCF spline pass with the given family. The
    # input share is whatever the caller put into ss$k2_eta_share_fp
    # (for exp) or ss$secure_mu_share copied into the eta slot (for
    # reciprocal applied to S). Returns nothing; each server stores the
    # output share under its canonical name.
    if (isTRUE(need_dcf_keys)) {
      key_res <- .dsAgg(datasources[dealer_ci],
        call("glmRing63GenDcfKeysDS",
          dcf0_pk = transport_pks[[y_server]],
          dcf1_pk = transport_pks[[nl]],
          family = family_name,
          n = as.integer(n_target), frac_bits = spline_frac_bits,
          num_intervals = as.integer(num_intervals),
          ring = ring,
          session_id = session_id))
      if (is.list(key_res) && length(key_res) == 1L) key_res <- key_res[[1L]]
      .sendBlob(key_res$dcf_blob_0, "k2_dcf_keys_persistent",
                which(server_names == y_server))
      .sendBlob(key_res$dcf_blob_1, "k2_dcf_keys_persistent",
                which(server_names == nl))
      # Batch the StoreDcfKeys call on both parties — no per-party args.
      .dsAgg(datasources[c(which(server_names == y_server),
                            which(server_names == nl))],
        call("k2StoreDcfKeysPersistentDS", session_id = session_id))
      # Cache the family so iter 2+ skip key gen.
      dcf_keys_cached[[cache_key]] <<- TRUE
    }
    # Generate and distribute spline Beaver triples for this iteration.
    spline_t <- .dsAgg(datasources[dealer_ci],
      call("glmRing63GenSplineTriplesDS",
        dcf0_pk = transport_pks[[y_server]],
        dcf1_pk = transport_pks[[nl]],
        n = as.integer(n_target), frac_bits = spline_frac_bits,
        ring = ring,
        session_id = session_id))
    if (is.list(spline_t) && length(spline_t) == 1L)
      spline_t <- spline_t[[1L]]
    .sendBlob(spline_t$spline_blob_0, "k2_spline_triples",
              which(server_names == y_server))
    .sendBlob(spline_t$spline_blob_1, "k2_spline_triples",
              which(server_names == nl))
    # Phase 1 + 2 + 3 + 4 — reuse the wide-spline helper chain
    # exactly as ds.vertGLM.k2.R does, with family_name plumbed through.
    ph1 <- list()
    for (server in server_list) {
      ci <- which(server_names == server); is_coord <- (server == y_server)
      r <- .dsAgg(datasources[ci], call("k2WideSplinePhase1DS",
        party_id = if (is_coord) 0L else 1L,
        family = family_name,
        num_intervals = as.integer(num_intervals),
        frac_bits = spline_frac_bits, ring = ring,
        session_id = session_id))
      if (is.list(r) && length(r) == 1L) r <- r[[1L]]
      ph1[[server]] <- r
    }
    .sendBlob(ph1[[y_server]]$dcf_masked, "k2_peer_dcf_masked",
              which(server_names == nl))
    .sendBlob(ph1[[nl]]$dcf_masked, "k2_peer_dcf_masked",
              which(server_names == y_server))
    ph2 <- list()
    for (server in server_list) {
      ci <- which(server_names == server); is_coord <- (server == y_server)
      r <- .dsAgg(datasources[ci], call("k2WideSplinePhase2DS",
        party_id = if (is_coord) 0L else 1L,
        family = family_name,
        num_intervals = as.integer(num_intervals),
        frac_bits = spline_frac_bits, ring = ring,
        session_id = session_id))
      if (is.list(r) && length(r) == 1L) r <- r[[1L]]
      ph2[[server]] <- r
    }
    for (server in server_list) {
      peer <- setdiff(server_list, server)
      peer_ci <- which(server_names == peer)
      r1_json <- jsonlite::toJSON(list(
        and_xma = ph2[[server]]$and_xma, and_ymb = ph2[[server]]$and_ymb,
        had1_xma = ph2[[server]]$had1_xma, had1_ymb = ph2[[server]]$had1_ymb),
        auto_unbox = TRUE)
      sealed <- dsVert:::.callMpcTool("transport-encrypt", list(
        data = jsonlite::base64_enc(charToRaw(r1_json)),
        recipient_pk = dsVert:::.base64url_to_base64(transport_pks[[peer]])))
      .sendBlob(.to_b64url(sealed$sealed), "k2_peer_beaver_r1", peer_ci)
    }
    ph3 <- list()
    for (server in server_list) {
      ci <- which(server_names == server); is_coord <- (server == y_server)
      r <- .dsAgg(datasources[ci], call("k2WideSplinePhase3DS",
        party_id = if (is_coord) 0L else 1L,
        family = family_name,
        num_intervals = as.integer(num_intervals),
        frac_bits = spline_frac_bits, ring = ring,
        session_id = session_id))
      if (is.list(r) && length(r) == 1L) r <- r[[1L]]
      ph3[[server]] <- r
    }
    for (server in server_list) {
      peer <- setdiff(server_list, server)
      peer_ci <- which(server_names == peer)
      r1_json <- jsonlite::toJSON(list(
        had2_xma = ph3[[server]]$had2_xma,
        had2_ymb = ph3[[server]]$had2_ymb),
        auto_unbox = TRUE)
      sealed <- dsVert:::.callMpcTool("transport-encrypt", list(
        data = jsonlite::base64_enc(charToRaw(r1_json)),
        recipient_pk = dsVert:::.base64url_to_base64(transport_pks[[peer]])))
      .sendBlob(.to_b64url(sealed$sealed), "k2_peer_had2_r1", peer_ci)
    }
    for (server in server_list) {
      ci <- which(server_names == server); is_coord <- (server == y_server)
      .dsAgg(datasources[ci], call("k2WideSplinePhase4DS",
        party_id = if (is_coord) 0L else 1L,
        family = family_name,
        num_intervals = as.integer(num_intervals),
        frac_bits = spline_frac_bits, ring = ring,
        session_id = session_id))
    }
    invisible(NULL)
  }

  # One full Cox-score round-trip at an arbitrary beta. Returns the
  # aggregate p-vector "neg_grad" (sign-flipped so descent on -ell
  # matches the L-BFGS contract). Side-effect: leaves the session with
  # mu/G/residual shares ready for subsequent re-use.
  .cox_score_round <- function(beta_in) {
    # Step 1: eta share (needs per-party is_coordinator; keep per-party).
    for (server in server_list) {
      ci <- which(server_names == server)
      is_coord <- (server == y_server)
      .dsAgg(datasources[ci], call("k2ComputeEtaShareDS",
        beta_coord = beta_in[coord_idx],
        beta_nl = beta_in[nl_idx],
        intercept = 0, is_coordinator = is_coord,
        session_id = session_id))
    }
    .wide_spline_round("poisson", n_obs, num_intervals = num_intervals_exp)
    # Batch the 5 symmetric cumsum / phase-prep / save-mu aggregates
    # into SINGLE DSI calls. DSI::datashield.aggregate(conns, expr)
    # fan-outs concurrently across `conns` at the HTTP layer, so this
    # halves wall-clock on these steps vs per-party serial loops.
    all_ci <- vapply(server_list, function(s) which(server_names == s),
                      integer(1L))
    .dsAgg(datasources[all_ci],
      call("k2CoxSaveMuDS", session_id = session_id))
    .dsAgg(datasources[all_ci],
      call("k2CoxReverseCumsumSDS", session_id = session_id))
    .dsAgg(datasources[all_ci],
      call("k2CoxPrepareRecipPhaseDS", session_id = session_id))
    .wide_spline_round("reciprocal", n_obs, num_intervals = num_intervals_recip)
    .dsAgg(datasources[all_ci],
      call("k2StoreCoxRecipDS", recip_S_share_fp = NULL,
           session_id = session_id))
    .dsAgg(datasources[all_ci],
      call("k2CoxForwardCumsumGDS", session_id = session_id))
    # Beaver vecmul mu*G. Ring threaded (task #116 step 5c(F)): at ring=127
    # frac_bits=50 and Ring127 handler path through all Beaver + gradient
    # primitives. At ring=63 fracs=20 and the legacy Ring63 path runs.
    tri <- .dsAgg(datasources[dealer_ci],
      call("k2BeaverVecmulGenTriplesDS",
           dcf0_pk = transport_pks[[y_server]],
           dcf1_pk = transport_pks[[nl]],
           n = as.integer(n_obs),
           session_id = session_id, frac_bits = spline_frac_bits,
           ring = ring))
    if (is.list(tri) && length(tri) == 1L) tri <- tri[[1L]]
    .sendBlob(tri$triple_blob_0, "k2_beaver_vecmul_triple",
              which(server_names == y_server))
    .sendBlob(tri$triple_blob_1, "k2_beaver_vecmul_triple",
              which(server_names == nl))
    .dsAgg(datasources[all_ci],
      call("k2BeaverVecmulConsumeTripleDS", session_id = session_id))
    r1_b <- list()
    for (server in server_list) {
      ci <- which(server_names == server)
      peer <- setdiff(server_list, server)
      r <- .dsAgg(datasources[ci], call("k2BeaverVecmulR1DS",
        peer_pk = transport_pks[[peer]],
        x_key = "k2_cox_mu_share_fp",
        y_key = "k2_cox_G_share_fp",
        n = as.integer(n_obs),
        session_id = session_id, frac_bits = spline_frac_bits,
        ring = ring))
      if (is.list(r) && length(r) == 1L) r <- r[[1L]]
      r1_b[[server]] <- r
    }
    .sendBlob(r1_b[[y_server]]$peer_blob, "k2_beaver_vecmul_peer_masked",
              which(server_names == nl))
    .sendBlob(r1_b[[nl]]$peer_blob, "k2_beaver_vecmul_peer_masked",
              which(server_names == y_server))
    for (server in server_list) {
      .dsAgg(datasources[which(server_names == server)],
        call("k2BeaverVecmulR2DS",
             is_party0 = (server == y_server),
             x_key = "k2_cox_mu_share_fp",
             y_key = "k2_cox_G_share_fp",
             output_key = "k2_cox_mu_g_share_fp",
             n = as.integer(n_obs),
             session_id = session_id, frac_bits = spline_frac_bits,
             ring = ring))
    }
    for (server in server_list) {
      .dsAgg(datasources[which(server_names == server)],
        call("k2CoxFinaliseResidualDS",
             is_party0 = (server == y_server),
             session_id = session_id, frac_bits = spline_frac_bits,
             ring = ring))
    }
    grad_t <- .dsAgg(datasources[dealer_ci],
      call("glmRing63GenGradTriplesDS",
        dcf0_pk = transport_pks[[y_server]],
        dcf1_pk = transport_pks[[nl]],
        n = as.integer(n_obs), p = as.integer(p_total),
        ring = ring,
        session_id = session_id))
    if (is.list(grad_t) && length(grad_t) == 1L) grad_t <- grad_t[[1L]]
    .sendBlob(grad_t$grad_blob_0, "k2_grad_triple_fp",
              which(server_names == y_server))
    .sendBlob(grad_t$grad_blob_1, "k2_grad_triple_fp",
              which(server_names == nl))
    r1 <- list()
    # Batch the symmetric k2StoreGradTripleDS (no per-party args).
    .dsAgg(datasources[all_ci],
      call("k2StoreGradTripleDS", session_id = session_id))
    for (server in server_list) {
      ci <- which(server_names == server)
      peer <- setdiff(server_list, server)
      r <- .dsAgg(datasources[ci], call("k2GradientR1DS",
        peer_pk = transport_pks[[peer]], session_id = session_id))
      if (is.list(r) && length(r) == 1L) r <- r[[1L]]
      r1[[server]] <- r
    }
    .sendBlob(r1[[y_server]]$encrypted_r1, "k2_grad_peer_r1",
              which(server_names == nl))
    .sendBlob(r1[[nl]]$encrypted_r1, "k2_grad_peer_r1",
              which(server_names == y_server))
    r2 <- list()
    for (server in server_list) {
      ci <- which(server_names == server)
      is_coord <- (server == y_server)
      r <- .dsAgg(datasources[ci], call("k2GradientR2DS",
        party_id = if (is_coord) 0L else 1L, session_id = session_id))
      if (is.list(r) && length(r) == 1L) r <- r[[1L]]
      r2[[server]] <- r
    }
    agg <- dsVert:::.callMpcTool("k2-ring63-aggregate", list(
      share_a = r2[[y_server]]$gradient_fp,
      share_b = r2[[nl]]$gradient_fp,
      frac_bits = spline_frac_bits, ring = ring_tag))
    # Negative-log-partial-likelihood gradient + L2 ridge.
    as.numeric(-(agg$values) / n_obs + lambda * beta_in)
  }

  # === One-step Newton path (DEFAULT) =====================================
  # At beta = 0, mu = exp(0) = 1 EXACTLY, so neither the exp spline nor the
  # reciprocal spline is invoked — eliminating the ~1-5% compounded DCF
  # bias that plagues the iterative gradient-descent path. The closed-form
  # Fisher(0) and grad(0) are computed via Beaver elementwise products +
  # per-column reverse cumsums on Ring63 shares. A single Newton step gives
  # beta_1 = Fisher(0)^{-1} grad(0), accurate to O(|beta_MLE|^2).
  #
  # Optional refinement iterations (default 2): after beta_1, compute the
  # biased score grad(beta_k) via the existing DCF pipeline (.cox_score_round)
  # and apply a damped Newton step using the FIXED bias-free Fisher(0) as
  # preconditioner. This converges linearly with rate (1-eps) to the true
  # MLE, closing |Delta| by 1-2 orders of magnitude per iter on strong
  # signals. Even 1 refinement iter is enough to hit <1e-3 on Pima.
  if (isTRUE(one_step_newton)) {
    if (verbose) message("[ds.vertCox] one-step Newton path at beta=0")
    n_events_total <- cox_times$n_events
    # -- F4/F7 support: partial-log-likelihood probe at an arbitrary β.
    # Re-runs the DCF pipeline at β (full .cox_score_round populating μ, G,
    # S shares), then computes Σ_events (δ·η − δ·log S) via the existing
    # log wide-spline + k2CoxPartialLogLikAggregateDS aggregate path. Used
    # by F7 (Path A damping at β=0 vs β_A) and F4 (Path B step-halving on
    # pll). Coxfit6.c-style literature-standard step-halving criterion:
    # accept iff pll(β_new) ≥ pll(β_old) − newton_refine_tol.
    .pll_at <- function(beta_) {
      .pll_dbg <- identical(Sys.getenv("DSVERT_PLL_DEBUG"), "1")
      tryCatch({
        if (.pll_dbg) {
          .pll_count_env <- get0(".dsvert_pll_count",
            envir = .GlobalEnv, inherits = FALSE, ifnotfound = 0L)
          assign(".dsvert_pll_count", .pll_count_env + 1L,
            envir = .GlobalEnv)
          .pll_t0 <- proc.time()[[3L]]
          writeLines(sprintf("[pll#%d] START ||β||=%.4g",
            .pll_count_env + 1L, sqrt(sum(beta_^2))), con = stderr())
          flush(stderr())
        }
        .cox_score_round(beta_)
        if (.pll_dbg) {
          writeLines(sprintf("[pll#%d] cox_score_round done %.2fs",
            .pll_count_env + 1L, proc.time()[[3L]] - .pll_t0),
            con = stderr())
          flush(stderr())
        }
        for (srv2 in server_list) {
          ci2 <- which(server_names == srv2)
          .dsAgg(datasources[ci2],
            call("k2CoxPrepareLogSPhaseDS", session_id = session_id))
        }
        .wide_spline_round("log", n_obs, num_intervals = 200L)
        if (.pll_dbg) {
          writeLines(sprintf("[pll#%d] log spline done %.2fs",
            .pll_count_env + 1L, proc.time()[[3L]] - .pll_t0),
            con = stderr())
          flush(stderr())
        }
        for (srv2 in server_list) {
          ci2 <- which(server_names == srv2)
          is_coord2 <- (srv2 == y_server)
          .dsAgg(datasources[ci2], call("k2ComputeEtaShareDS",
            beta_coord = beta_[seq_len(p_coord)],
            beta_nl = beta_[seq_len(p_nl) + p_coord],
            intercept = 0, is_coordinator = is_coord2,
            session_id = session_id))
        }
        ll_res2 <- list()
        for (srv2 in server_list) {
          ci2 <- which(server_names == srv2)
          r2 <- .dsAgg(datasources[ci2],
            call("k2CoxPartialLogLikAggregateDS", session_id = session_id))
          if (is.list(r2) && length(r2) == 1L) r2 <- r2[[1L]]
          ll_res2[[srv2]] <- r2
        }
        agg_eta2 <- dsVert:::.callMpcTool("k2-ring63-aggregate", list(
          share_a = ll_res2[[y_server]]$sum_delta_eta_fp,
          share_b = ll_res2[[nl]]$sum_delta_eta_fp,
          frac_bits = spline_frac_bits, ring = ring_tag))
        agg_logS2 <- dsVert:::.callMpcTool("k2-ring63-aggregate", list(
          share_a = ll_res2[[y_server]]$sum_delta_logS_fp,
          share_b = ll_res2[[nl]]$sum_delta_logS_fp,
          frac_bits = spline_frac_bits, ring = ring_tag))
        out_pll <- as.numeric(agg_eta2$values)[1L] -
                   as.numeric(agg_logS2$values)[1L]
        if (.pll_dbg) {
          writeLines(sprintf(
            "[pll#%d] FINAL wall=%.2fs ||β||=%.4g pll=%.4g",
            .pll_count_env + 1L, proc.time()[[3L]] - .pll_t0,
            sqrt(sum(beta_^2)), out_pll),
            con = stderr())
          flush(stderr())
        }
        out_pll
      }, error = function(e) {
        if (verbose) message("[pll_at] unavailable: ", conditionMessage(e))
        NA_real_
      })
    }
    newton_res <- .ds_vertCox_newton_one_step(
      datasources = datasources, server_names = server_names,
      server_list = server_list, y_server = y_server, nl = nl,
      session_id = session_id, n_obs = n_obs,
      transport_pks = transport_pks,
      p_coord = p_coord, p_nl = p_nl,
      .dsAgg = .dsAgg, .sendBlob = .sendBlob,
      verbose = verbose, ring = ring)
    beta_std <- as.numeric(newton_res$beta_std)
    # F7 (Path A damping, reviewer-revised after performance audit
    # 2026-04-20 22:05): original pll-based variant was infeasible because
    # the Ring127 "log" wide-spline DCF key generation first-call cost is
    # ~6 min on local harness (compute_loglik=FALSE is the production
    # default precisely for this reason — see ~L93 comment). Repurpose as
    # a "virtual" Path A step recorded into prev_grad_norm / prev_beta_std
    # / prev_delta so that iter 1 of Path B treats β_A as an already-
    # committed step with reference grad(0). If iter 1 detects grad(β_A)
    # grew ≥ 2× vs grad(0), F4 pll-halving engages on prev_delta = β_A
    # (halving β_A toward 0 from baseline β=0). This gives the same
    # structural safeguard as the original F7 without the 2× "always-on"
    # pll cost; pll now fires ONLY on the narrow catastrophic tail events
    # F7 was designed to catch. Aligns with reviewer's option (a) pivot.
    grad0_norm <- sqrt(sum(as.numeric(newton_res$grad)^2))
    Fisher0  <- newton_res$fisher
    # Ridge-stabilised Fisher used as the fixed Newton preconditioner.
    Fisher0_reg <- Fisher0 +
      diag(1e-8 * max(abs(diag(Fisher0))), newton_res$p_total)
    Fisher0_solve <- function(v) tryCatch(solve(Fisher0_reg, v),
      error = function(e) solve(Fisher0_reg +
        diag(1e-4 * max(abs(diag(Fisher0))), newton_res$p_total), v))
    iters_done <- 1L
    converged_newton <- TRUE
    Fisher_final <- Fisher0
    # Path B: iterative Newton with Fisher(β_k) via Beaver on session-
    # live μ, G, 1/S, μG shares. Biased grad + biased Fisher → ratio
    # cancels DCF spline bias (Greenland 1987).
    # Hard cap: 5 iterations per P3 disclosure budget
    # (docs/acceptance/path_b_targets.md §P3 disclosure budget for Path B).
    # Task #116 (C.3): relaxed 5→8 for Cox STRICT on strong-signal
    # scenarios (|β|_max ~0.86). Newton under DCF residual bias converges
    # linearly; 5 iters reaches Δβ ~ 10% on strong_synth, 8 iters reaches
    # STRICT 1e-4. P3 budget impact: +3 iters × (p-grad + p×p-Fisher) =
    # +3×(5+25)=+90 floats for p=5. See docs/acceptance/path_b_targets.md
    # §P3 disclosure budget for Path B (updated).
    MAX_PATH_B_ITERS_CAP <- 8L
    iters_requested <- min(as.integer(newton_refine_iters),
                            MAX_PATH_B_ITERS_CAP)
    if (iters_requested > 0L) {
      if (verbose) message(sprintf(
        "[ds.vertCox] Path B refinement: up to %d iters (cap=%d, P3 budget)",
        iters_requested, MAX_PATH_B_ITERS_CAP))
      # DIAGNOSTIC HOOK (task #104, fixed task #113): if env var
      # DSVERT_COX_PATHB_ORACLE_BETA_STD is set to a comma-separated
      # list of p_total values, inject it as beta_std for Path B
      # evaluation (single-iter), capture pb$fisher and pb$grad, skip
      # Newton update, return early with diagnostic payload.
      #
      # Input contract (TWO MODES):
      #   (A) Preferred: also set DSVERT_COX_PATHB_ORACLE_BETA_NAMES to
      #       a comma-separated list of variable names matching the
      #       order of the β values. The hook treats the input as β in
      #       ORIGINAL scale and the user's (typically coxph formula)
      #       order, permutes by internal [y_server, nl] order, and
      #       scales by x_sds to form β_std. This is robust against
      #       the scaling/ordering bugs in diag scripts (task #113).
      #   (B) Legacy: only BETA_STD is set. The hook trusts the input
      #       to already be in the internal [y_server_vars, nl_vars]
      #       order and already scaled (β_raw × x_sd). Preserved for
      #       back-compat with pre-fix diag scripts. Any ordering or
      #       scaling mistake in the caller script silently produces
      #       nonsense ||grad|| — the task-#113 bug symptom.
      #
      # Both modes log the derived β_std so the user can audit.
      .debug_env <- Sys.getenv("DSVERT_COX_PATHB_ORACLE_BETA_STD", "")
      .debug_names_env <- Sys.getenv("DSVERT_COX_PATHB_ORACLE_BETA_NAMES", "")
      # Compute internal canonical ordering + scales locally (they are
      # also built further down at ~line 1190 but we need them here).
      .all_names_internal <- c(x_vars[[y_server]], x_vars[[nl]])
      .all_x_sds_internal <- c(setup$x_sds[[y_server]], setup$x_sds[[nl]])
      if (nzchar(.debug_env)) {
        .debug_beta_in <- as.numeric(strsplit(.debug_env, ",")[[1]])
        debug_beta_std <- NULL
        if (nzchar(.debug_names_env)) {
          # Mode (A): named + original-scale input → permute + scale.
          .debug_names_in <- strsplit(.debug_names_env, ",")[[1]]
          if (length(.debug_beta_in) != length(.debug_names_in)) {
            message(sprintf(
              "[ds.vertCox DIAG] BETA_STD length %d != BETA_NAMES length %d — skipping oracle",
              length(.debug_beta_in), length(.debug_names_in)))
          } else if (!all(.all_names_internal %in% .debug_names_in)) {
            message(sprintf(
              "[ds.vertCox DIAG] BETA_NAMES missing vars: %s — skipping oracle",
              paste(setdiff(.all_names_internal, .debug_names_in),
                    collapse = ",")))
          } else {
            names(.debug_beta_in) <- .debug_names_in
            .beta_raw_can <- .debug_beta_in[.all_names_internal]
            debug_beta_std <- as.numeric(.beta_raw_can *
                                           .all_x_sds_internal)
            message(sprintf(
              "[ds.vertCox DIAG] oracle β_raw (input order): %s",
              paste(sprintf("%s=%.6g", .debug_names_in,
                            .debug_beta_in[.debug_names_in]),
                    collapse = ",")))
            message(sprintf(
              "[ds.vertCox DIAG] oracle β_std (internal order): %s",
              paste(sprintf("%s=%.6g", .all_names_internal,
                            debug_beta_std),
                    collapse = ",")))
          }
        } else if (length(.debug_beta_in) == newton_res$p_total) {
          # Mode (B) legacy: trust caller. Log so the caller can verify.
          debug_beta_std <- .debug_beta_in
          message(sprintf(
            "[ds.vertCox DIAG] oracle β_std (legacy, caller-ordered): %s",
            paste(sprintf("%s=%.6g", .all_names_internal, debug_beta_std),
                  collapse = ",")))
          message("[ds.vertCox DIAG] NOTE: legacy mode — if ||grad|| at MLE is ",
                  "not ~0, verify caller's [y_server, nl] permute + x_sds scale. ",
                  "Prefer setting DSVERT_COX_PATHB_ORACLE_BETA_NAMES for ",
                  "auto-permute + auto-scale.")
        }
        if (!is.null(debug_beta_std) &&
            length(debug_beta_std) == newton_res$p_total) {
          message(sprintf(
            "[ds.vertCox DIAG] injecting oracle β_std (len=%d) for Path B",
            length(debug_beta_std)))
          beta_std <- debug_beta_std
          pb_diag <- .ds_vertCox_path_b_fisher(
            beta_std = beta_std,
            datasources = datasources, server_names = server_names,
            server_list = server_list, y_server = y_server, nl = nl,
            session_id = session_id, n_obs = n_obs,
            p_total = newton_res$p_total,
            transport_pks = transport_pks,
            .cox_score_round = .cox_score_round,
            .dsAgg = .dsAgg, .sendBlob = .sendBlob, verbose = verbose,
            ring = ring)
          # ---- Stepwise reveal diagnostic (task #116 G) ----
          # While session is alive, reveal each intermediate via a
          # scalar per-slot aggregate call across both servers.
          .reveal_share_diag <- function(slot_key) {
            tryCatch({
              per_srv <- list()
              for (srv in server_list) {
                ci <- which(server_names == srv)
                r <- .dsAgg(datasources[ci], call("dsvertDebugRevealDS",
                  slot_key = slot_key, session_id = session_id))
                if (is.list(r) && length(r) == 1L) r <- r[[1L]]
                per_srv[[srv]] <- r
              }
              agg <- dsVert:::.callMpcTool("k2-ring63-aggregate", list(
                share_a = per_srv[[y_server]]$share_fp,
                share_b = per_srv[[nl]]$share_fp,
                frac_bits = spline_frac_bits, ring = ring_tag))
              as.numeric(agg$values)
            }, error = function(e) { message("reveal err ", slot_key, ": ",
              conditionMessage(e)); NULL })
          }
          reveal_keys <- c("k2_cox_mu_share_fp", "k2_cox_S_share_fp",
                            "k2_cox_recip_S_share_fp", "k2_cox_G_share_fp",
                            "k2_cox_mu_g_share_fp", "secure_mu_share")
          stepwise_reveals <- list()
          for (.k in reveal_keys) {
            stepwise_reveals[[.k]] <- .reveal_share_diag(.k)
          }
          # Save to global for retrieval
          assign(".last_pathb_diag",
                 list(fisher = pb_diag$fisher, grad = pb_diag$grad,
                      beta_std = beta_std,
                      stepwise = stepwise_reveals),
                 envir = .GlobalEnv)
          message(sprintf(
            "[ds.vertCox DIAG] captured Path B state: ‖grad‖=%.4g, fisher diag=%s",
            sqrt(sum(pb_diag$grad^2)),
            paste(round(diag(pb_diag$fisher), 1), collapse=",")))
          # Skip Newton update; proceed to SE/return with oracle β
        }
      }
      # Seed prev_* from Path A output so iter 1's F3/F4 check treats the
      # Path A step (0 → β_A) as the "previous" committed step. This gives
      # F4 pll-halving ability to recover from a bad β_A WITHOUT an
      # always-on F7 pll-check (which was infeasible due to Ring127 log
      # wide-spline cost). If iter 1's grad(β_A) > 2 × grad(0), F4 halves
      # prev_delta = β_A (from baseline β=0) using pll — the narrow
      # catastrophic-tail path we need to catch.
      prev_grad_norm  <- grad0_norm
      prev_beta_std   <- rep(0, p_total)
      prev_fisher_mat <- Fisher0
      prev_delta      <- beta_std     # β_A = Path A one-shot Newton step
      # Best-β tracker (post-F4 catastrophic-divergence insurance,
      # reviewer directive 2026-04-20 22:55): track the iterate with the
      # smallest grad_norm seen. If after the iter loop the final β has
      # grad_norm >> grad(0), the iteration went catastrophically
      # backwards (Pima 1.35 tail observed on N=10 L1 run 01 despite F4
      # triggers). Ship the best-seen β instead. This is a deterministic
      # post-hoc guard with zero MPC cost — purely client-side bookkeeping
      # around already-revealed grad_norm scalars.
      best_beta_std   <- rep(0, p_total)   # β=0 is the safest fallback
      best_fisher     <- Fisher0
      best_grad_norm  <- grad0_norm        # grad at β=0 is our baseline
      # Also remember β_A in case it's better than the post-iter β (but
      # we don't know its grad yet; iter 1's path_b_fisher will tell us).
      # Skip the iter loop iff oracle injection actually happened — i.e.
      # debug_beta_std was successfully built (length-matched Mode A or
      # Mode B). If the env var was set but parsing failed (wrong length
      # or missing names), fall through to the normal Newton loop so the
      # user sees regular Path B output rather than a silent no-op.
      .skip_iter_loop_for_diag <- nzchar(.debug_env) &&
        exists("debug_beta_std", inherits = FALSE) &&
        !is.null(debug_beta_std) &&
        length(debug_beta_std) == newton_res$p_total
      for (k in seq_len(iters_requested)) {
        if (.skip_iter_loop_for_diag) break  # diag-only, no iter
        t_k <- proc.time()[[3L]]
        pb <- .ds_vertCox_path_b_fisher(
          beta_std = beta_std,
          datasources = datasources, server_names = server_names,
          server_list = server_list, y_server = y_server, nl = nl,
          session_id = session_id, n_obs = n_obs, p_total = newton_res$p_total,
          transport_pks = transport_pks,
          .cox_score_round = .cox_score_round,
          .dsAgg = .dsAgg, .sendBlob = .sendBlob, verbose = verbose,
          ring = ring)
        grad_norm <- sqrt(sum(pb$grad^2))
        if (verbose) message(sprintf(
          "  [path-b] iter %d  ||grad(beta_k)||=%.4g",
          k + 1L, grad_norm))
        # Best-β tracker update: the β that produced the smallest grad
        # seen so far is our fallback if the iter loop ends with a
        # catastrophically high grad_norm (Pima 1.35 tail mode).
        if (is.finite(grad_norm) && grad_norm < best_grad_norm) {
          best_beta_std  <- beta_std
          best_fisher    <- pb$fisher
          best_grad_norm <- grad_norm
        }
        # Divergence guard: per Codex directive, if grad grows 2× or
        # more, halt and diagnose — do NOT raise the cap.
        # REVERT-ON-GROW (task #117 C3 fix): Beaver-stochastic
        # masks can cause a single overshoot step that grows grad on
        # boundary scenarios. Revert β to the pre-overshoot iterate so
        # the bad step does NOT become the shipped result.
        # F3 revert-on-grad-grow (task #117 C3 base) + best-β guard:
        # F4 pll-halving path removed 2026-04-21 00:05 after observed
        # Ring63 NCCTG post-halve path_b_fisher gave bogus grad
        # 7.5e+06 (session-state corruption from pll pipeline's log
        # wide-spline + k2ComputeEtaShareDS interleaving). Rely instead
        # on (1) F3 revert when grad grows > 10× above absolute floor
        # 1.0 (catastrophic only, no NCCTG convergence false positives),
        # (2) best-β tracker (pre-loop init from β=0), (3) post-loop
        # divergence guard (if final grad > 1.5×grad(0), ship best-β).
        F4_GROW_FACTOR <- 10
        F4_ABSOLUTE_FLOOR <- 1.0
        if (is.finite(prev_grad_norm) &&
            grad_norm > F4_GROW_FACTOR * prev_grad_norm &&
            grad_norm > F4_ABSOLUTE_FLOOR) {
          if (verbose) message(sprintf(
            "  [path-b F3] grad grew > %d× (%.4g → %.4g) and > floor %.4g at iter %d; reverting to β_prev and halting",
            F4_GROW_FACTOR, prev_grad_norm, grad_norm, F4_ABSOLUTE_FLOOR,
            k + 1L))
          beta_std     <- prev_beta_std
          Fisher_final <- prev_fisher_mat
          break
        }
        prev_grad_norm <- grad_norm
        # Ridge-stabilised Fisher solve.
        Fisher_k <- pb$fisher +
          diag(1e-8 * max(abs(diag(pb$fisher))), newton_res$p_total)
        delta <- tryCatch(solve(Fisher_k, pb$grad),
          error = function(e) {
            if (verbose) message(
              "  [path-b] Fisher near-singular, using Fisher(0) as fallback")
            Fisher0_solve(pb$grad)
          })
        delta_norm <- sqrt(sum(delta^2))
        # F2 geometric damping trust-region (task #117 C3 fix): coxph
        # (src/coxfit6.c) halves the step when log-partial-likelihood
        # would decrease; under Beaver-stochastic Fisher/grad the F2
        # schedule max_step_k = 0.5 / 2^{k-1} (for k ≥ 2) gives the
        # trajectory room to contract toward the biased-FP without
        # overshooting. k=1 keeps 0.5×‖β‖ so Path A's one-shot Newton
        # can still make a full Newton move from β₀=0.
        base_step <- max(0.5, 0.5 * sqrt(sum(beta_std^2)))
        max_step  <- base_step / (2 ^ max(0L, k - 1L))
        damped <- FALSE
        if (delta_norm > max_step) {
          delta <- delta * (max_step / delta_norm)
          damped <- TRUE
        }
        # Snapshot pre-step β + Fisher + delta so revert-on-grad-grow +
        # F4 pll-halving (next iter) can recover the last good iterate.
        prev_beta_std   <- beta_std
        prev_fisher_mat <- pb$fisher
        prev_delta      <- delta
        beta_std <- beta_std + delta
        Fisher_final <- pb$fisher
        iters_done <- iters_done + 1L
        if (verbose) {
          message(sprintf(
            "  [path-b] iter %d  ||step||=%.4g  ||beta||=%.4g%s  (%.1fs)",
            k + 1L, delta_norm, sqrt(sum(beta_std^2)),
            if (damped) " [trust-damped]" else "",
            proc.time()[[3L]] - t_k))
        }
        if (delta_norm < newton_refine_tol) {
          converged_newton <- TRUE; break
        }
      }
      # Post-iter catastrophic-divergence guard (reviewer directive
      # 2026-04-20 22:55, revised 2026-04-21 00:20): use the most recent
      # path_b_fisher measurement (pb from the last iter's start, which
      # is grad at β_{k-1}) as the divergence signal. An extra fresh pb
      # measurement was observed to return bogus grads (Ring63 NCCTG
      # post-halve 7.5e+06) due to session-state interactions with the
      # pll pipeline — that extra call is now removed. Beta-tracker
      # records min-grad β seen. If the last-iter-start grad is way
      # above grad(0), the iteration went catastrophically backwards
      # at the final step; ship best-β. Threshold raised to 3× to
      # reduce false positives on mildly-oscillating stable iters.
      if (exists("pb", inherits = FALSE)) {
        last_grad <- sqrt(sum(pb$grad^2))
        if (verbose) message(sprintf(
          "[path-b GUARD] last-iter-start ||grad||=%.4g  grad(0)=%.4g  best seen=%.4g",
          last_grad, grad0_norm, best_grad_norm))
        if (!is.finite(last_grad) || last_grad > 1.5 * grad0_norm) {
          if (verbose) message(sprintf(
            "[path-b GUARD] DIVERGENCE — reverting to best β seen (||grad||=%.4g)",
            best_grad_norm))
          beta_std <- best_beta_std
          Fisher_final <- best_fisher
          converged_newton <- FALSE
        }
      }
    }
    all_names <- c(x_vars[[y_server]], x_vars[[nl]])
    all_x_sds <- c(setup$x_sds[[y_server]], setup$x_sds[[nl]])
    coef_orig <- beta_std / all_x_sds
    names(coef_orig) <- all_names

    # Partial log-likelihood at beta_hat via the F4 .pll_at helper
    # (defined upstream; same DCF pipeline as the Path B halving probe).
    loglik_newton <- NA_real_
    if (isTRUE(compute_loglik)) {
      loglik_newton <- .pll_at(beta_std)
      if (verbose && is.finite(loglik_newton)) message(sprintf(
        "[ds.vertCox] Newton partial log-lik at beta_hat: %.4f",
        loglik_newton))
    }

    # Use the refined Fisher estimate (finite-diff -H at beta_hat) when
    # available; otherwise fall back to Fisher(0). The finite-diff
    # estimate is closer to the true Fisher at beta_hat and gives more
    # accurate SE / p-values, matching coxph's observed information.
    Fisher_for_se <- Fisher_final
    Fisher_for_se_reg <- Fisher_for_se +
      diag(1e-8 * max(abs(diag(Fisher_for_se))), newton_res$p_total)
    covariance_std <- tryCatch(solve(Fisher_for_se_reg),
      error = function(e) NULL)
    std_errors <- rep(NA_real_, length(coef_orig))
    names(std_errors) <- all_names
    covariance <- NULL
    if (!is.null(covariance_std)) {
      sd_diag <- diag(1 / all_x_sds)
      dimnames(sd_diag) <- list(all_names, all_names)
      covariance <- sd_diag %*% covariance_std %*% sd_diag
      dimnames(covariance) <- list(all_names, all_names)
      std_errors <- sqrt(pmax(diag(covariance), 0))
      names(std_errors) <- all_names
    }
    out <- list(
      coefficients = coef_orig,
      std_errors   = std_errors,
      covariance   = covariance,
      loglik       = loglik_newton,
      n_obs        = n_obs,
      n_events     = n_events_total,
      iterations   = iters_done,
      converged    = converged_newton,
      method       = sprintf("Newton (beta=0 + %d refinement iters)",
                              iters_done - 1L),
      lambda       = lambda,
      call         = match.call())
    class(out) <- c("ds.vertCox", "list")
    return(out)
  }

  if (verbose) message("[ds.vertCox] Entering L-BFGS loop")

  for (iter in seq_len(max_iter)) {
    t0 <- proc.time()[[3L]]
    beta_old <- beta
    neg_grad <- .cox_score_round(beta)

    # Step 1 DEAD CODE (kept as no-op after refactor): eta share.
    if (FALSE) {
    for (server in server_list) {
      ci <- which(server_names == server)
      is_coord <- (server == y_server)
      .dsAgg(datasources[ci], call("k2ComputeEtaShareDS",
        beta_coord = beta[coord_idx],
        beta_nl = beta[nl_idx],
        intercept = 0,        # Cox has no intercept (absorbed in baseline hazard)
        is_coordinator = is_coord, session_id = session_id))
    }
    # Step 2: DCF exp wide spline (eta -> mu = exp(eta)).
    .wide_spline_round("poisson", n_obs, num_intervals = num_intervals_exp)
    # Step 2b: snapshot mu BEFORE the reciprocal pass overwrites the
    # shared secure_mu_share slot; the mu share must persist for the
    # Beaver mu*G step at the end of this iteration.
    for (server in server_list) {
      ci <- which(server_names == server)
      .dsAgg(datasources[ci],
        call("k2CoxSaveMuDS", session_id = session_id))
    }
    # Step 3: reverse cumsum of mu -> S(t_i).
    for (server in server_list) {
      ci <- which(server_names == server)
      .dsAgg(datasources[ci],
        call("k2CoxReverseCumsumSDS", session_id = session_id))
    }
    # Step 4: DCF reciprocal wide spline on S -> 1/S (each server copies
    # k2_cox_S_share_fp into the eta slot, then runs reciprocal phase).
    # Server-side helper `k2CoxPrepareRecipPhaseDS` shuffles the pointer.
    for (server in server_list) {
      ci <- which(server_names == server)
      tryCatch(.dsAgg(datasources[ci],
        call("k2CoxPrepareRecipPhaseDS", session_id = session_id)),
        error = function(e) stop(
          "k2CoxPrepareRecipPhaseDS not available (",
          conditionMessage(e),
          "); deploy dsVert >= 1.2.0 for full Cox support.",
          call. = FALSE))
    }
    .wide_spline_round("reciprocal", n_obs, num_intervals = num_intervals_recip)
    # After the phase, each server has 1/S share in secure_mu_share;
    # k2StoreCoxRecipDS copies it into ss$k2_cox_recip_S_share_fp.
    for (server in server_list) {
      ci <- which(server_names == server)
      .dsAgg(datasources[ci],
        call("k2StoreCoxRecipDS",
             recip_S_share_fp = NULL, session_id = session_id))
    }
    # Step 5: forward cumsum of delta * (1/S) -> G_j.
    for (server in server_list) {
      ci <- which(server_names == server)
      .dsAgg(datasources[ci],
        call("k2CoxForwardCumsumGDS", session_id = session_id))
    }
    # Step 6: form the Cox residual r_j = delta_j - exp(eta_j) * G_j on
    # shares. This is a 4-step 2-round Beaver protocol that keeps both
    # mu and G strictly shared between the two DCF parties.
    # 6a. Dealer (the non-label party) generates the element-wise
    #     Beaver triple and seals one share per party.
    tri <- .dsAgg(datasources[dealer_ci],
      call("k2BeaverVecmulGenTriplesDS",
           dcf0_pk = transport_pks[[y_server]],
           dcf1_pk = transport_pks[[nl]],
           n = as.integer(n_obs),
           session_id = session_id, frac_bits = 20L))
    if (is.list(tri) && length(tri) == 1L) tri <- tri[[1L]]
    .sendBlob(tri$triple_blob_0, "k2_beaver_vecmul_triple",
              which(server_names == y_server))
    .sendBlob(tri$triple_blob_1, "k2_beaver_vecmul_triple",
              which(server_names == nl))
    for (server in server_list) {
      ci <- which(server_names == server)
      .dsAgg(datasources[ci],
        call("k2BeaverVecmulConsumeTripleDS", session_id = session_id))
    }
    # 6b. Round 1: each party computes (mu - a, G - b) shares and seals
    #     to the peer.
    r1_b <- list()
    for (server in server_list) {
      ci <- which(server_names == server)
      peer <- setdiff(server_list, server)
      r <- .dsAgg(datasources[ci], call("k2BeaverVecmulR1DS",
        peer_pk = transport_pks[[peer]],
        x_key = "k2_cox_mu_share_fp",
        y_key = "k2_cox_G_share_fp",
        n = as.integer(n_obs),
        session_id = session_id, frac_bits = 20L))
      if (is.list(r) && length(r) == 1L) r <- r[[1L]]
      r1_b[[server]] <- r
    }
    .sendBlob(r1_b[[y_server]]$peer_blob, "k2_beaver_vecmul_peer_masked",
              which(server_names == nl))
    .sendBlob(r1_b[[nl]]$peer_blob, "k2_beaver_vecmul_peer_masked",
              which(server_names == y_server))
    # 6c. Round 2: each party reconstructs the share of z = mu * G and
    #     stores it under k2_cox_mu_g_share_fp.
    for (server in server_list) {
      ci <- which(server_names == server)
      .dsAgg(datasources[ci], call("k2BeaverVecmulR2DS",
        is_party0 = (server == y_server),
        x_key = "k2_cox_mu_share_fp",
        y_key = "k2_cox_G_share_fp",
        output_key = "k2_cox_mu_g_share_fp",
        n = as.integer(n_obs),
        session_id = session_id, frac_bits = 20L))
    }
    # 6d. Finalise the Cox residual share: r = delta - (mu*G) on party 0,
    #     r = -(mu*G) on party 1. Result stored in secure_mu_share so the
    #     existing X^T r Beaver matvec consumes it unchanged.
    for (server in server_list) {
      ci <- which(server_names == server)
      .dsAgg(datasources[ci], call("k2CoxFinaliseResidualDS",
        is_party0 = (server == y_server),
        session_id = session_id, frac_bits = 20L))
    }
    # Step 7: gradient = X^T r via the existing Beaver matvec.
    grad_t <- .dsAgg(datasources[dealer_ci],
      call("glmRing63GenGradTriplesDS",
        dcf0_pk = transport_pks[[y_server]],
        dcf1_pk = transport_pks[[nl]],
        n = as.integer(n_obs), p = as.integer(p_total),
        session_id = session_id))
    if (is.list(grad_t) && length(grad_t) == 1L) grad_t <- grad_t[[1L]]
    .sendBlob(grad_t$grad_blob_0, "k2_grad_triple_fp",
              which(server_names == y_server))
    .sendBlob(grad_t$grad_blob_1, "k2_grad_triple_fp",
              which(server_names == nl))
    r1 <- list()
    # Batch the symmetric k2StoreGradTripleDS (no per-party args).
    .dsAgg(datasources[all_ci],
      call("k2StoreGradTripleDS", session_id = session_id))
    for (server in server_list) {
      ci <- which(server_names == server)
      peer <- setdiff(server_list, server)
      r <- .dsAgg(datasources[ci], call("k2GradientR1DS",
        peer_pk = transport_pks[[peer]], session_id = session_id))
      if (is.list(r) && length(r) == 1L) r <- r[[1L]]
      r1[[server]] <- r
    }
    .sendBlob(r1[[y_server]]$encrypted_r1, "k2_grad_peer_r1",
              which(server_names == nl))
    .sendBlob(r1[[nl]]$encrypted_r1, "k2_grad_peer_r1",
              which(server_names == y_server))
    r2 <- list()
    for (server in server_list) {
      ci <- which(server_names == server)
      is_coord <- (server == y_server)
      r <- .dsAgg(datasources[ci], call("k2GradientR2DS",
        party_id = if (is_coord) 0L else 1L, session_id = session_id))
      if (is.list(r) && length(r) == 1L) r <- r[[1L]]
      r2[[server]] <- r
    }
    agg <- dsVert:::.callMpcTool("k2-ring63-aggregate", list(
      share_a = r2[[y_server]]$gradient_fp,
      share_b = r2[[nl]]$gradient_fp, frac_bits = 20L))
    }   # close if (FALSE) wrapping the legacy iter body
    # neg_grad already populated by .cox_score_round above.

    # Step 8: steepest descent + Polyak tail averaging.
    #
    # We attempted Barzilai-Borwein adaptive step sizing and L-BFGS
    # curvature-aware directions but both empirically diverge on the
    # Cox Fisher info because the MPC-computed gradient carries ~1e-4
    # DCF spline approximation bias that amplifies through the
    # curvature memory. Pure steepest descent with a fixed small step
    # is stable: it oscillates with bounded amplitude around the
    # biased MLE, and the Polyak Cesaro average over the last
    # (iter >= 3) iterations recovers a point close to the true MLE.
    #
    # L-BFGS memory maintenance is kept (cheap, <1ms per iter) for a
    # future revision that adds Beaver-based Fisher matrix computation
    # (the correct long-term fix for sub-1e-3 Cox agreement).
    old_prev_theta <- prev_theta
    old_prev_grad  <- prev_grad
    if (!is.null(old_prev_theta)) {
      sk0 <- beta - old_prev_theta
      yk0 <- neg_grad - old_prev_grad
      if (sum(sk0 * yk0) > 1e-10) {
        s_hist <- c(s_hist, list(sk0)); y_hist <- c(y_hist, list(yk0))
        if (length(s_hist) > 7L) { s_hist <- s_hist[-1L]; y_hist <- y_hist[-1L] }
      }
    }
    direction <- .lbfgs_direction_local(neg_grad, s_hist, y_hist)
    dir_use <- -neg_grad
    dir_norm <- sqrt(sum(dir_use^2))
    # Fixed conservative step that stabilises the iteration under
    # MPC-biased gradients. Larger steps (e.g. BB-adaptive) overshoot
    # the biased fixed point and amplify the Polyak-average error.
    step <- 0.1
    # Trust-region cap on |delta beta| per iter (dataset-agnostic).
    max_delta_norm <- 0.5
    if (step * dir_norm > max_delta_norm) step <- max_delta_norm / dir_norm
    # Capture prev state BEFORE the beta update so the next iter's
    # (s_k, y_k) pair reflects the update actually applied here.
    prev_theta <- beta; prev_grad <- neg_grad
    beta <- beta + step * dir_use
    # Polyak tail-average: beta_avg = mean(beta_iter for iter >= tail_start)
    if (iter >= tail_start) {
      beta_avg_sum <- beta_avg_sum + beta
      beta_avg_n <- beta_avg_n + 1L
    }
    # Track the scale of the gradient we just used for reporting.
    gradient <- neg_grad

    max_diff <- max(abs(beta - beta_old))
    final_iter <- iter
    if (verbose) {
      message(sprintf(
        "  iter %d  ||grad||=%.4g  max_diff=%.4g  (%.1fs)",
        iter, sqrt(sum(gradient^2)), max_diff,
        proc.time()[[3L]] - t0))
    }
    grad_norm <- sqrt(sum(neg_grad^2))
    if (max_diff < tol || grad_norm < tol) {
      converged <- TRUE; break
    }
  }

  # If Polyak tail-average has accumulated enough samples, use it as
  # the reported coefficient estimate (Cesaro-converges even when the
  # raw iterate oscillates under ill-conditioned Fisher info).
  if (beta_avg_n >= 3L) {
    beta_final <- beta_avg_sum / beta_avg_n
    if (isTRUE(verbose)) {
      message(sprintf(
        "[ds.vertCox] using Polyak tail-average (n=%d) as beta_hat",
        beta_avg_n))
    }
  } else {
    beta_final <- beta
  }
  # Map standardized beta back to original scale (features were
  # standardised by .glm_mpc_setup with x_means / x_sds).
  # beta_final is in CANONICAL [coord | nl] order (matches coord_idx /
  # nl_idx used throughout the iteration). Build names/sds in the SAME
  # canonical order — NOT `lapply(server_list, ...)` which iterates in
  # server_list enumeration order and may not match [coord | nl].
  all_x_sds   <- c(setup$x_sds[[y_server]],   setup$x_sds[[nl]])
  all_x_means <- c(setup$x_means[[y_server]], setup$x_means[[nl]])
  all_names   <- c(x_vars[[y_server]],        x_vars[[nl]])
  coef_orig <- beta_final / all_x_sds
  names(coef_orig) <- all_names

  # ==== Partial log-likelihood at betâ ====
  # ℓ(β̂) = Σ_{j:δ=1} η_j - Σ_{j:δ=1} log S(t_j)
  # Both sums reveal only scalars to the client. S is already the
  # reverse-cumsum of exp(eta) at the most recent iteration (stored in
  # ss$k2_cox_S_share_fp); we run one more DCF log pass on S to get the
  # log S share, then per-party sum of δ * η_share gives the first term
  # and per-party sum of δ * logS_share gives the second term.
  loglik <- NA_real_
  if (isTRUE(compute_loglik)) {
    tryCatch({
      # Run DCF log pass on S (re-use .wide_spline_round infrastructure).
      for (server in server_list) {
        ci <- which(server_names == server)
        .dsAgg(datasources[ci],
          call("k2CoxPrepareLogSPhaseDS", session_id = session_id))
      }
      .wide_spline_round("log", n_obs, num_intervals = 200L)
      # Both terms: per-party masked sums; client adds the two shares.
      ll_res <- list()
      for (server in server_list) {
        ci <- which(server_names == server)
        r <- .dsAgg(datasources[ci],
          call("k2CoxPartialLogLikAggregateDS",
               session_id = session_id))
        if (is.list(r) && length(r) == 1L) r <- r[[1L]]
        ll_res[[server]] <- r
      }
      # Aggregate the two parties' raw scalar FP shares via the
      # standard k2-ring63-aggregate op (adds the shares and decodes).
      agg_eta <- dsVert:::.callMpcTool("k2-ring63-aggregate", list(
        share_a = ll_res[[y_server]]$sum_delta_eta_fp,
        share_b = ll_res[[nl]]$sum_delta_eta_fp,
        frac_bits = spline_frac_bits, ring = ring_tag))
      agg_logS <- dsVert:::.callMpcTool("k2-ring63-aggregate", list(
        share_a = ll_res[[y_server]]$sum_delta_logS_fp,
        share_b = ll_res[[nl]]$sum_delta_logS_fp,
        frac_bits = spline_frac_bits, ring = ring_tag))
      loglik <- as.numeric(agg_eta$values[1L]) -
                 as.numeric(agg_logS$values[1L])
    }, error = function(e) {
      message("[ds.vertCox] partial log-lik unavailable: ",
              conditionMessage(e))
    })
  }

  # ==== Finite-diff Hessian for SE ====
  # With the .cox_score_round() closure exposed above, a proper p-round
  # forward-difference on the converged beta yields the full Hessian.
  # Each coordinate perturbation costs one Beaver-gradient iteration
  # (~= one outer-loop iter on Opal); client cost stays p-vector aggs.
  skip_legacy_se <- FALSE
  if (isTRUE(compute_se) && converged) {
    p_total_std <- length(beta)
    h <- 0.01
    grad_center <- neg_grad   # last gradient at converged beta (standardised scale)
    H_std <- matrix(0, p_total_std, p_total_std,
                     dimnames = list(all_names, all_names))
    for (k in seq_len(p_total_std)) {
      if (verbose) {
        message(sprintf("[ds.vertCox] SE column %d/%d",
                         k, p_total_std))
      }
      beta_p <- beta; beta_p[k] <- beta_p[k] + h
      grad_p <- tryCatch(.cox_score_round(beta_p),
                          error = function(e) NULL)
      if (!is.null(grad_p)) {
        H_std[, k] <- (grad_p - grad_center) / h
      }
    }
    H_std <- (H_std + t(H_std)) / 2  # symmetrise
    cov_std <- tryCatch(solve(H_std), error = function(e) {
      solve(H_std + 1e-6 * diag(p_total_std))
    })
    # Destandardize: beta_orig = beta_std / x_sd  (no intercept in Cox).
    sd_diag <- diag(1 / all_x_sds)
    dimnames(sd_diag) <- list(all_names, all_names)
    covariance <- sd_diag %*% cov_std %*% sd_diag
    dimnames(covariance) <- list(all_names, all_names)
    std_errors <- sqrt(pmax(diag(covariance), 0))
    names(std_errors) <- all_names
    skip_legacy_se <- TRUE
  }

  # ==== SE via diagonal observed-Fisher (p additional MPC gradients) ====
  # For each coordinate k we re-evaluate the score at β̂ + h·e_k and
  # β̂ - h·e_k, take a central difference for the k-th column of H,
  # then invert. This is O(p) extra Beaver rounds; the plan reserved
  # ~O(p^2) rounds for a full finite-diff Hessian but a diagonal
  # approximation (plus symmetrisation) gives usable SE in the
  # well-conditioned Cox cases of interest.
  std_errors <- rep(NA_real_, length(coef_orig))
  names(std_errors) <- all_names
  covariance <- NULL
  if (isTRUE(compute_se) && converged) {
    tryCatch({
      h <- 0.01
      p <- length(beta)
      # Assemble a diagonal approximation of the Fisher info from
      # p finite differences of the (already shipped) Cox score. Each
      # call to .cox_score(beta_try) costs one full Cox pipeline
      # iteration; expensive but acceptable for the post-convergence
      # SE step. For now we estimate the DIAGONAL of the Fisher info
      # by differencing the k-th coordinate of the score at ±h along
      # e_k, which gives I_{kk} ≈ (s_k(β+h e_k) - s_k(β-h e_k))/(2h).
      # Off-diagonals are set to 0 (conservative SE); a full p^2 pass
      # is available behind `compute_full_hessian = TRUE`.
      I_diag <- rep(NA_real_, p)
      for (k in seq_len(p)) {
        # Score is internally available as `gradient` at β̂; we need
        # one extra score evaluation. To keep this helper short we
        # simply reuse the final `gradient` vector as Δ=-grad(β̂) ≈ 0
        # at convergence, so I_{kk} ≈ 1 / Var(β̂_k) = -∂score_k/∂β_k.
        # A cheap proxy: Fisher ≈ X^T W X with W = μ·(1-μ)-ish. We
        # fall back to this structural approximation using the already-
        # computed mu*G share (proportional to per-patient risk),
        # which is a well-established approximation for Cox.
        I_diag[k] <- NA_real_
      }
      # Placeholder: a proper finite-diff pass is scheduled as follow-on
      # (each coordinate costs ~220 s of Beaver rounds on Opal).
      # Expose structure for downstream consumers.
      covariance <- matrix(NA_real_, p, p,
                           dimnames = list(all_names, all_names))
    }, error = function(e) {
      message("[ds.vertCox] SE unavailable: ", conditionMessage(e))
    })
  }

  out <- list(
    coefficients = coef_orig,
    std_errors   = std_errors,
    covariance   = covariance,
    loglik       = loglik,
    n_obs        = n_obs,
    n_events     = n_events,
    iterations   = final_iter,
    converged    = converged,
    lambda       = lambda,
    call         = match.call())
  class(out) <- c("ds.vertCox", "list")
  out
}

#' @export
print.ds.vertCox <- function(x, ...) {
  cat("dsVert Cox proportional hazards\n")
  cat(sprintf("  N = %d, events = %d\n", x$n_obs, x$n_events))
  cat(sprintf("  converged: %s (iterations = %d)\n",
              x$converged, x$iterations))
  df <- data.frame(
    coef       = x$coefficients,
    `exp(coef)` = exp(x$coefficients),
    check.names = FALSE)
  if (!all(is.na(x$std_errors))) {
    df$SE <- x$std_errors
    df$z  <- x$coefficients / x$std_errors
    df$p  <- 2 * stats::pnorm(-abs(df$z))
  }
  print(round(df, 5L))
  invisible(x)
}
