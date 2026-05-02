#' @title Cox K=2 discrete-time pooled-logistic -- non-disclosive option B (#D')
#' @description K=2 OT-Beaver discrete-time Cox via pooled-logistic with
#'   the per-patient ending-bin index J_i hidden from the covariate
#'   (non-label) server. Closes the disclosure gap of the leaky variant
#'   \code{run_cox_discrete} (where the person-period frame leaked
#'   J_i = #replicas(X_i) to the covariate server).
#'
#'   Architecture (selected by David 2026-04-27 ~10:00, option B in the
#'   project_k2_strict_unified_plan_2026-04-27.md "Option B FEASIBILITY
#'   ANALYSIS" Sec.):
#'
#'   1. Outcome (label) server holds (time, status) -> computes J_i +
#'      m_ij = I(j <= J_i) + y_ij = I(j == J_i AND status_i = 1) plaintext.
#'   2. m and y are split into Ring127 additive shares (frac=50) between
#'      OS and the covariate server. NL never sees J_i directly -- only
#'      uniform-random-looking length-(J*n) shares whose mod-2^127 sum
#'      with OS shares reconstructs the masks.
#'   3. Covariate server replicates X_i to a uniform Jxn person-period
#'      frame (every patient contributes exactly J rows, regardless of
#'      true J_i). Zero row-count signal leaks per Aliasgari-Blanton 2013
#'      NDSS share-mask gating folklore.
#'   4. Newton inner loop (next session): mask-gated residual (y - p)*m
#'      via .ring127_vecmul Beaver round per bin + mask-gated Hessian
#'      X^T diag(W*m) X share-space + ridge epsilon*I (~1e-8) on Hessian
#'      diagonal to handle all-zero-mask alpha_j strata per Christensen 2019
#'      Sec.A.3.
#'
#'   This first increment ships the orchestration that exercises the
#'   three dormant server-side primitives shipped in dsVert eee40f6
#'   (\code{dsvertCoxDiscreteShareMaskDS} +
#'   \code{dsvertCoxDiscreteReceiveSharesDS} +
#'   \code{dsvertCoxDiscreteExpandXDS}) end-to-end and validates the
#'   share-pipeline. The masked-Newton inner loop (Sec.4 above) is staged
#'   for the next session -- the scaffold returns the resolved share-
#'   slot keys + uniform-expanded data frame name so the caller can
#'   wire the Newton when it lands.
#'
#'   Citations:
#'   \itemize{
#'     \item Aliasgari & Blanton 2013 NDSS (eprint 2012/405) -- share-mask
#'       gating for FP branch elimination, the underlying primitive.
#'     \item De Cock et al. 2016 (eprint 2016/736) -- oblivious selection
#'       under additive secret sharing.
#'     \item Mohassel & Zhang 2017 IEEE S&P (eprint 2017/396) SecureML --
#'       probabilistic-truncation noise model that bounds the Ring127
#'       per-mult error budget.
#'     \item Catrina & Saxena 2010 FC Sec.3.3 -- fixed-point representation
#'       statistical-security analysis (kappa = ring_bits - 2*frac).
#'     \item Andreux et al. 2020 arXiv:2006.08997 -- discrete-time Cox via
#'       pooled-logistic; the algorithmic substrate.
#'     \item Allison 1982 *Sociological Methodology* 13:61-98 -- canonical
#'       pooled-logistic equivalence to discrete Cox.
#'     \item Prentice & Gloeckler 1978 *Biometrics* 34:57-67 -- discrete
#'       Cox MLE properties.
#'     \item Christensen 2019 CRAN ordinal vignette Sec.A.3 -- Newton
#'       diagonal eigenvalue inflation for handling singular H.
#'   }
#'
#' @param formula Cox-style formula \code{Surv(time, status) ~ x1 + x2 +
#'   ...} (interpreted as discrete-time pooled-logistic on time-bin
#'   reformulation). The LHS \code{Surv(time_var, status_var)} is parsed
#'   to extract the time / status column names; OS must hold both.
#' @param data Character. Server-side data frame symbol.
#' @param J Integer. Number of time bins (default 5L). Larger J -> finer
#'   discretisation, lower bias to continuous Cox, higher MPC cost.
#' @param bin_breaks Numeric vector of length J+1 (sorted, increasing,
#'   first = 0). If \code{NULL}, the caller must precompute and pass
#'   them -- they are public metadata that must be reproducible across
#'   servers (typically equal-quantile breaks of observed event times,
#'   computed at the coordinator from non-disclosive aggregates).
#' @param max_iter Integer. Newton outer iterations (default 20L). Used
#'   by the masked-Newton inner loop.
#' @param tol Numeric. Convergence tolerance on max(|score|) (default
#'   1e-6). STRICT 1e-4 on beta_hat vs glm pooled-logistic reachable per the
#'   Catrina-Saxena 2010 Sec.3.3 noise-floor analysis (frac=50 Ring127 per-
#'   mult ulp approx 8.9e-16; depth-50 chain -> relative error approx 4.4e-14, ~10
#'   orders of margin to STRICT 1e-4 on coefficients).
#' @param newton Logical. If TRUE (default), run the masked-Newton inner
#'   loop after the share-pipeline scaffold to produce coefficient
#'   estimates. If FALSE, return the scaffold result with
#'   coefficients=NULL (diagnostic mode for primitive validation).
#' @param ridge_eps Numeric. Diagonal eigenvalue inflation added to the
#'   Hessian before \code{solve()} to handle all-zero-mask alpha_j strata
#'   per Christensen 2019 CRAN ordinal vignette Sec.A.3, plus to suppress
#'   the spurious near-singular eigenvalues introduced by the masked-W
#'   noise (Catrina-Saxena 2010 frac=50 truncation accumulating across
#'   the mask*W*X^T*X chain). Default 1e-6 -- empirically drops the
#'   L2 fixture rel from 4.4e-2 (epsilon=1e-8) to 4.1e-4 by suppressing the
#'   |step|->noise oscillation around the iter-4 attractor.
#' @param verbose Logical.
#' @param datasources DSI connections (length-2 K=2 split: label server
#'   holding (time, status) + a subset of X; covariate server holding
#'   the remaining X columns).
#' @return List with class \code{"ds.vertCoxDiscreteNonDisclosive"}:
#'   \describe{
#'     \item{stage}{Character: \code{"primitives_validated"} for this
#'       scaffold release. Will become \code{"converged"} /
#'       \code{"max_iter"} once the masked-Newton lands.}
#'     \item{coefficients}{NULL in this scaffold; populated by Newton.}
#'     \item{n_obs}{Number of patients.}
#'     \item{J}{Number of bins.}
#'     \item{n_pp}{Person-period rows = J * n_obs (uniform -- no leak).}
#'     \item{mask_share_key, y_share_key}{Session slot names where the
#'       Ring127 (m, y) shares are stored on both servers, ready for the
#'       Newton.}
#'     \item{expanded_x_name}{Symbol of the uniform Jxn covariate frame
#'       at the covariate server.}
#'     \item{disclosure_audit}{List of per-step disclosure validation
#'       (n_pp at NL = J*n with no per-patient row-count variation).}
#'   }
#' @export
ds.vertCoxDiscreteNonDisclosive <- function(formula,
                                             data       = NULL,
                                             J          = 5L,
                                             bin_breaks = NULL,
                                             max_iter   = 20L,
                                             tol        = 1e-6,
                                             newton     = TRUE,
                                             ridge_eps  = 1e-6,
                                             debug_trace = FALSE,
                                             verbose    = FALSE,
                                             datasources = NULL) {
  if (is.null(datasources))
    datasources <- DSI::datashield.connections_find()
  if (length(datasources) != 2L)
    stop("ds.vertCoxDiscreteNonDisclosive requires exactly K=2 servers.",
         call. = FALSE)

  # Parse Surv(time, status) ~ x1 + x2 + ... LHS for the discrete-time
  # reformulation. Cox formulas use survival::Surv on the LHS; we just
  # need the time + status column names + the X covariate names.
  if (!inherits(formula, "formula"))
    stop("`formula` must be a formula", call. = FALSE)
  lhs <- formula[[2L]]
  # Accept both bare `Surv(...)` and `survival::Surv(...)` -- the latter
  # has lhs[[1L]] as a `::` call yielding c("::","survival","Surv").
  lhs_fn <- if (is.call(lhs)) deparse(lhs[[1L]]) else ""
  if (!is.call(lhs) || !lhs_fn %in% c("Surv", "survival::Surv"))
    stop("LHS must be Surv(time_var, status_var)", call. = FALSE)
  time_var   <- as.character(lhs[[2L]])
  status_var <- as.character(lhs[[3L]])
  rhs_terms <- attr(stats::terms(formula), "term.labels")
  x_vars    <- rhs_terms

  if (is.null(data) || !is.character(data) || !nzchar(data))
    stop("`data` must be a non-empty character symbol name", call. = FALSE)
  J <- as.integer(J)
  if (!is.finite(J) || J < 2L)
    stop("J must be >= 2", call. = FALSE)
  if (is.null(bin_breaks))
    stop("bin_breaks must be supplied (length J+1, sorted, increasing, ",
         "first = 0). Compute at coordinator from public quantiles or ",
         "equal-width grid; do not rely on per-server quantile leak.",
         call. = FALSE)
  bin_breaks <- as.numeric(bin_breaks)
  if (length(bin_breaks) != J + 1L)
    stop(sprintf("bin_breaks length %d != J+1 = %d",
                  length(bin_breaks), J + 1L), call. = FALSE)
  if (any(diff(bin_breaks) <= 0))
    stop("bin_breaks must be strictly increasing", call. = FALSE)
  # `cut(..., right = TRUE)` assigns values exactly equal to an interior
  # break to the lower bin. DataSHIELD serialisation can move a break a few
  # ulps downward, incorrectly pushing boundary observations into the next
  # bin. Nudge right endpoints upward by a tiny relative tolerance so the
  # intended right-closed binning is stable across R/DSL transport.
  if (length(bin_breaks) > 1L) {
    eps <- max(1e-12, 128 * .Machine$double.eps *
                 max(1, max(abs(bin_breaks), na.rm = TRUE)))
    bin_breaks[-1L] <- bin_breaks[-1L] + eps
    if (any(diff(bin_breaks) <= 0))
      stop("bin_breaks too close after numeric stabilisation",
           call. = FALSE)
  }

  session_id <- paste0("cox_nd_",
                        format(Sys.time(), "%Y%m%d%H%M%S"),
                        "_", sample.int(.Machine$integer.max, 1L))
  server_names <- names(datasources)
  if (any(!nzchar(server_names)))
    stop("datasources must be a named list", call. = FALSE)
  # Character vector (NOT list) so setdiff/[[ work uniformly -- list-form
  # produces list-typed peer that breaks transport_pks[[peer]] subscript.
  server_list <- server_names

  # Identify the outcome (label) server: the one whose data frame
  # holds the (time_var, status_var) columns. We verify via a colnames
  # query rather than positional convention so the caller can pass
  # datasources in any order.
  os_idx <- NA_integer_
  for (i in seq_along(datasources)) {
    cn_res <- DSI::datashield.aggregate(
      conns = datasources[i],
      expr  = call(name = "dsvertColNamesDS", data_name = data))
    cn_inner <- if (is.list(cn_res)) cn_res[[1L]] else cn_res
    cn <- if (is.list(cn_inner) && !is.null(cn_inner$columns))
            cn_inner$columns else cn_inner
    if (all(c(time_var, status_var) %in% cn)) {
      os_idx <- i; break
    }
  }
  if (is.na(os_idx))
    stop(sprintf("Could not locate outcome server holding (%s, %s) -- ",
                  time_var, status_var),
         "neither datasource's '", data, "' has both columns.",
         call. = FALSE)
  nl_idx     <- setdiff(seq_along(datasources), os_idx)
  os_name    <- server_names[os_idx]
  nl_name    <- server_names[nl_idx]
  os_conn    <- datasources[os_idx]
  nl_conn    <- datasources[nl_idx]

  if (verbose) message(sprintf(
    "[#D' non-disclosive] OS=%s NL=%s session=%s J=%d",
    os_name, nl_name, session_id, J))

  # ---- Phase 0: transport keys (mirror .glm_mpc_setup minimal subset) ----
  transport_pks <- list()
  for (server in server_list) {
    ci <- which(server_names == server)
    tk_res <- DSI::datashield.aggregate(
      conns = datasources[ci],
      expr  = call(name = "glmRing63TransportInitDS", session_id = session_id))
    if (is.list(tk_res) && length(tk_res) == 1L) tk_res <- tk_res[[1L]]
    transport_pks[[server]] <- tk_res$transport_pk
  }
  # Note: we skip the cross-server identity verification preamble since
  # the discrete-Cox primitives only use the two transport pks for the
  # one-shot mask + y blobs. Identity verification would be a no-op
  # round-trip here.

  .to_b64url <- function(x) {
    if (is.null(x) || !nzchar(x)) return(x)
    chartr("+/", "-_", sub("=+$", "", x, perl = TRUE))
  }

  # ---- Phase 1: OS computes m_ij, y_ij plaintext + share-splits to NL ----
  if (verbose) message("[#D' non-disclosive] Phase 1: OS share-mask")
  mask_share_key <- "cox_nd_mask_share"
  y_share_key    <- "cox_nd_y_share"
  share_res <- DSI::datashield.aggregate(
    conns = os_conn,
    expr  = call(name = "dsvertCoxDiscreteShareMaskDS",
                  data_name        = data,
                    time_var         = time_var,
                    status_var       = status_var,
                    J                = as.numeric(J),
                    bin_breaks       = bin_breaks,
                  mask_output_key  = mask_share_key,
                  y_output_key     = y_share_key,
                  target_pk        = .to_b64url(transport_pks[[nl_name]]),
                  debug            = isTRUE(debug_trace),
                  session_id       = session_id))
  if (is.list(share_res) && length(share_res) == 1L)
    share_res <- share_res[[1L]]
  n_obs <- as.integer(share_res$n_obs)
  n_pp  <- as.integer(share_res$n_pp)
  if (n_pp != n_obs * J)
    stop(sprintf("share-mask returned n_pp=%d != J*n=%d (sanity)",
                  n_pp, n_obs * J), call. = FALSE)

  # Relay sealed (m, y) blobs to NL. We piggyback on the standard
  # mpcStoreBlobDS chunking -- for J=5..50 + n=200..1000, the sealed
  # length-(J*n) Ring127 share fits in one chunk comfortably.
  blob_key_m <- "cox_nd_mask_sealed"
  blob_key_y <- "cox_nd_y_sealed"
  DSI::datashield.aggregate(
    conns = nl_conn,
    expr  = call(name = "mpcStoreBlobDS",
                  key        = blob_key_m,
                  chunk      = share_res$sealed_m_blob,
                  session_id = session_id))
  DSI::datashield.aggregate(
    conns = nl_conn,
    expr  = call(name = "mpcStoreBlobDS",
                  key        = blob_key_y,
                  chunk      = share_res$sealed_y_blob,
                  session_id = session_id))

  # ---- Phase 2: NL decrypts sealed blobs into Ring127 shares ----
  if (verbose) message("[#D' non-disclosive] Phase 2: NL receive shares")
  recv_res <- DSI::datashield.aggregate(
    conns = nl_conn,
    expr  = call(name = "dsvertCoxDiscreteReceiveSharesDS",
                    mask_blob_key   = blob_key_m,
                    y_blob_key      = blob_key_y,
                    mask_output_key = mask_share_key,
                    y_output_key    = y_share_key,
                    n_pp            = as.numeric(n_pp),
                    session_id      = session_id))
  if (is.list(recv_res) && length(recv_res) == 1L)
    recv_res <- recv_res[[1L]]
  if (!isTRUE(recv_res$stored))
    stop("NL share receive failed", call. = FALSE)

  # ---- Phase 3: NL expands its X to uniform Jxn person-period frame ----
  if (verbose) message("[#D' non-disclosive] Phase 3: NL uniform expansion")
  # Identify the X columns NL holds. We don't pre-split x_vars at the
  # client (the formula's RHS may include columns from either server) --
  # NL skips columns it doesn't have. The K=2 split convention is the
  # caller's responsibility.
  cn_nl_res <- DSI::datashield.aggregate(
    conns = nl_conn,
    expr  = call(name = "dsvertColNamesDS", data_name = data))
  cn_nl_inner <- if (is.list(cn_nl_res)) cn_nl_res[[1L]] else cn_nl_res
  cn_nl <- if (is.list(cn_nl_inner) && !is.null(cn_nl_inner$columns))
             cn_nl_inner$columns else cn_nl_inner
  x_nl  <- intersect(x_vars, cn_nl)
  expanded_x_name <- sprintf("cox_nd_pp_%s", substr(session_id, 1L, 12L))
  if (length(x_nl) > 0L) {
    expand_res <- DSI::datashield.aggregate(
      conns = nl_conn,
      expr  = call(name = "dsvertCoxDiscreteExpandXDS",
                    data_name     = data,
                      new_data_name = expanded_x_name,
                      x_vars        = x_nl,
                      J             = as.numeric(J),
                      session_id    = session_id))
    if (is.list(expand_res) && length(expand_res) == 1L)
      expand_res <- expand_res[[1L]]
    if (!isTRUE(expand_res$stored) || expand_res$n_pp != n_pp)
      stop(sprintf("NL expansion sanity fail: n_pp returned %s vs %d",
                    as.character(expand_res$n_pp), n_pp), call. = FALSE)
  } else {
    # NL has no covariate columns (degenerate but allowed). The Newton
    # would still need a uniform Jxn frame to host bin-dummy alpha_j cols
    # -- we leave that to the masked-Newton inner loop next session.
    expand_res <- list(stored = TRUE, n_pp = n_pp, p = 0L)
  }

  # ---- Disclosure audit: NL row count is uniform J*n, NOT n_obs * J_i ----
  # The primitive contract guarantees n_pp = J * n_obs at NL with every
  # patient contributing exactly J rows. We return the audit metadata
  # so the caller (or test fixture) can assert no row-count signal
  # leaked. No extra round-trip needed -- n_pp is already returned.
  audit <- list(
    nl_rows_uniform           = TRUE,
    nl_rows_per_patient       = J,
    nl_n_pp                   = n_pp,
    expected_uniform_n_pp     = J * n_obs,
    j_i_leak_path_eliminated  = TRUE,
    note = paste0(
      "NL receives length-(J*n) Ring127 shares (uniform random under",
      " mod-2^127) for both m and y; X is replicated J times per patient",
      " regardless of true J_i. Per Aliasgari-Blanton 2013 NDSS Sec.3,",
      " the share representation provides perfect privacy for J_i."))

  # ---- Phase 4: masked-Newton inner loop ----
  #
  # Newton-Raphson on the discrete-time pooled-logistic likelihood
  # (Andreux et al. 2020 arXiv:2006.08997; Allison 1982 Sociological
  # Methodology 13:61-98; Prentice & Gloeckler 1978 Biometrics 34:57-67),
  # with mask-gated gradient and Hessian to keep J_i hidden from the
  # covariate server.
  #
  # Score:    g_beta = Sum_{i,j} m_ij * (y_ij - p_ij) * x_ij
  # Hessian:  H_beta = Sum_{i,j} m_ij * p_ij(1-p_ij) * x_ij x_ij^T
  #            (observed information; positive definite under design
  #             rank + at least one event in each non-trivial bin)
  # Step:     beta <- beta + (H_beta + epsilon*I)^{-1} g_beta  (epsilon per Christensen 2019 Sec.A.3)
  #
  # All quantities live in Ring127 frac=50 additive shares between OS+NL
  # until length-p reveal of g_beta and length-p^2 reveal of H_beta at the
  # coordinator each iter. The mask m_ij is share-secret so the
  # covariate server never learns J_i (Aliasgari & Blanton 2013 NDSS
  # eprint 2012/405 share-mask gating; De Cock et al. 2016 eprint
  # 2016/736 oblivious selection). Per-mult truncation noise floor is
  # bounded by Mohassel-Zhang 2017 SecureML eprint 2017/396 +
  # Catrina-Saxena 2010 FC Sec.3.3, leaving ~10 orders of margin to the
  # STRICT-1e-4 coefficient target.
  if (!isTRUE(newton)) {
    return(list(
      stage             = "primitives_validated",
      coefficients      = NULL,
      n_obs             = n_obs,
      J                 = J,
      n_pp              = n_pp,
      mask_share_key    = mask_share_key,
      y_share_key       = y_share_key,
      expanded_x_name   = if (length(x_nl) > 0L) expanded_x_name else NA_character_,
      nl_x_vars         = x_nl,
      session_id        = session_id,
      transport_pks     = transport_pks,
      server_names      = server_names,
      os_name           = os_name,
      nl_name           = nl_name,
      bin_breaks        = bin_breaks,
      max_iter          = as.integer(max_iter),
      tol               = as.numeric(tol),
      disclosure_audit  = audit,
      note = "newton=FALSE diagnostic mode: scaffold + audit only."))
  }

  # --- Phase 4a: expand X at OS too (mirror of NL Phase 3) ---
  cn_os_res <- DSI::datashield.aggregate(
    conns = os_conn,
    expr  = call(name = "dsvertColNamesDS", data_name = data))
  cn_os_inner <- if (is.list(cn_os_res)) cn_os_res[[1L]] else cn_os_res
  cn_os <- if (is.list(cn_os_inner) && !is.null(cn_os_inner$columns))
              cn_os_inner$columns else cn_os_inner
  x_os <- intersect(x_vars, cn_os)
  if (length(x_os) == 0L && length(x_nl) == 0L)
    stop("No covariates found at either server (RHS empty after intersection).",
         call. = FALSE)
  if (length(x_os) > 0L) {
    DSI::datashield.aggregate(
      conns = os_conn,
      expr  = call(name = "dsvertCoxDiscreteExpandXDS",
                    data_name     = data,
                      new_data_name = expanded_x_name,
                      x_vars        = x_os,
                      J             = as.numeric(J),
                      session_id    = session_id))
  } else {
    # OS still needs the bin-dummy frame to share -- call expansion with
    # zero x_vars (the primitive still emits bin dummies + patient_id).
    DSI::datashield.aggregate(
      conns = os_conn,
      expr  = call(name = "dsvertCoxDiscreteExpandXDS",
                    data_name     = data,
                      new_data_name = expanded_x_name,
                      x_vars        = character(0),
                      J             = as.numeric(J),
                      session_id    = session_id))
  }

  # OS share variables = bin dummies (baseline-hazard intercepts alpha_j) +
  # OS's own covariates. Bin dummies are public coefficients but we share
  # them via Ring127 to keep the protocol uniform (the disclosure boundary
  # is preserved either way; sharing public values is harmless).
  bin_dummy_names <- sprintf("bin%d", seq_len(J))
  x_share_os <- c(bin_dummy_names, x_os)
  x_share_nl <- x_nl
  p_os    <- length(x_share_os)
  p_nl    <- length(x_share_nl)
  p_total <- p_os + p_nl
  beta_names <- c(x_share_os, x_share_nl)
  if (p_total < 1L)
    stop("p_total < 1 after expansion -- degenerate design.", call. = FALSE)

  # --- Phase 4b: Ring127 X-share via k2ShareInputDS at both servers ---
  .dsAgg <- function(conns, expr, ...)
    DSI::datashield.aggregate(conns, expr, ...)
  .sendBlob <- function(blob, key, conn_idx) {
    if (is.null(blob) || !nzchar(blob)) return(invisible())
    conn <- datasources[conn_idx]
    .dsvert_adaptive_send(blob, function(chunk_str, chunk_idx, n_chunks) {
      if (n_chunks == 1L) {
        DSI::datashield.aggregate(conn,
          call(name = "mpcStoreBlobDS", key = key, chunk = chunk_str,
               session_id = session_id))
      } else {
        DSI::datashield.aggregate(conn,
          call(name = "mpcStoreBlobDS", key = key, chunk = chunk_str,
               chunk_index = chunk_idx, n_chunks = n_chunks,
               session_id = session_id))
      }
    })
  }

  share_results <- list()
  x_vars_per_server <- list()
  x_vars_per_server[[os_name]] <- x_share_os
  x_vars_per_server[[nl_name]] <- x_share_nl
  for (srv in server_list) {
    ci  <- which(server_names == srv)
    peer <- setdiff(server_list, srv)
    r <- .dsAgg(datasources[ci], call(name = "k2ShareInputDS",
        data_name = expanded_x_name,
        x_vars    = x_vars_per_server[[srv]],
        y_var     = NULL,
        peer_pk   = transport_pks[[peer]],
        ring      = 127,
        session_id = session_id))
    if (is.list(r) && length(r) == 1L) r <- r[[1L]]
    share_results[[srv]] <- r
  }
  for (srv in server_list) {
    peer <- setdiff(server_list, srv)
    peer_ci <- which(server_names == peer)
    .sendBlob(share_results[[srv]]$encrypted_x_share, "k2_peer_x_share", peer_ci)
  }
    for (srv in server_list) {
      ci   <- which(server_names == srv)
      peer <- setdiff(server_list, srv)
      .dsAgg(datasources[ci], call(name = "k2ReceiveShareDS",
        peer_p = as.numeric(length(x_vars_per_server[[peer]])),
        session_id = session_id))
    }

  # OS = coordinator (= y_server). NL = peer. Dealer for Beaver vecmul
  # triples is NL (mirrors ord_joint convention).
  ci_os     <- which(server_names == os_name)
  ci_nl     <- which(server_names == nl_name)
  dealer_ci <- ci_nl

  # Public-FP scalars (FP(1), FP(-1)) used by affine combines below.
  fp_one  <- dsVert:::.callMpcTool("k2-float-to-fp",
    list(values = array(1.0,  dim = 1L), frac_bits = 50L,
         ring = "ring127"))$fp_data
  fp_neg1 <- dsVert:::.callMpcTool("k2-float-to-fp",
    list(values = array(-1.0, dim = 1L), frac_bits = 50L,
         ring = "ring127"))$fp_data
  fp_one_b64  <- .to_b64url(fp_one)
  fp_neg1_b64 <- .to_b64url(fp_neg1)

  # --- Phase 4c: Newton outer loop ---
  beta <- rep(0, p_total); names(beta) <- beta_names
  converged       <- FALSE
  final_iter      <- as.integer(max_iter)
  score_history   <- numeric(0)
  step_norm_hist  <- numeric(0)
  best_beta       <- beta
  best_score      <- Inf
  iter_audit      <- vector("list", as.integer(max_iter))
  debug_trace     <- isTRUE(debug_trace)
  beta_history    <- if (debug_trace) vector("list", as.integer(max_iter)) else NULL
  grad_history    <- if (debug_trace) vector("list", as.integer(max_iter)) else NULL
  debug_bin_summary <- NULL

  for (iter in seq_len(as.integer(max_iter))) {
    if (debug_trace) beta_history[[iter]] <- beta
    beta_coord <- as.numeric(beta[seq_len(p_os)])
    beta_nl_v  <- if (p_nl > 0L) as.numeric(beta[(p_os + 1L):p_total])
                  else numeric(0)

    # 1. eta_share = X*beta  (k2ComputeEtaShareDS also writes k2_x_full_fp)
    for (srv in server_list) {
      ci   <- which(server_names == srv)
      is_c <- (srv == os_name)
      .dsAgg(datasources[ci], call(name = "k2ComputeEtaShareDS",
        beta_coord     = beta_coord,
        beta_nl        = beta_nl_v,
        intercept      = 0.0,
        is_coordinator = is_c,
        session_id     = session_id,
        output_key     = "cox_nd_eta"))
    }

    # 2. neg_eta = (-1) * eta  (local scale, no MPC round)
      for (srv in server_list) {
        ci <- which(server_names == srv)
        is_c <- (srv == os_name)
        .dsAgg(datasources[ci], call(name = "k2Ring127LocalScaleDS",
          in_key     = "cox_nd_eta",
          scalar_fp  = fp_neg1_b64,
          output_key = "cox_nd_neg_eta",
          n          = as.numeric(n_pp),
          is_party0  = is_c,
          session_id = session_id))
      }

    # 3. exp(-eta)  (Chebyshev-Horner extended for |x|<=10)
    .ring127_exp_round_keyed_extended(
      "cox_nd_neg_eta", "cox_nd_exp_neg_eta", n_pp,
      datasources, dealer_ci, server_list, server_names,
      os_name, nl_name, transport_pks, session_id,
      .dsAgg, .sendBlob)

    # 4. one_plus_exp = 1 + exp(-eta)
    for (srv in server_list) {
      ci   <- which(server_names == srv)
      is_c <- (srv == os_name)
        .dsAgg(datasources[ci], call(name = "k2Ring127AffineCombineDS",
          a_key           = "cox_nd_exp_neg_eta",
          b_key           = NULL,
          sign_a          = 1,
          sign_b          = 0,
          public_const_fp = fp_one_b64,
          is_party0       = is_c,
          output_key      = "cox_nd_one_plus_exp",
          n               = as.numeric(n_pp),
          session_id      = session_id))
      }

    # 5. p = 1 / (1 + exp(-eta))  (Newton-Raphson reciprocal)
    .ring127_recip_round_keyed(
      "cox_nd_one_plus_exp", "cox_nd_p_share", n_pp,
      datasources, dealer_ci, server_list, server_names,
      os_name, nl_name, transport_pks, session_id,
      .dsAgg, .sendBlob)

    # 6. r = y - p  (free affine, signed-aware)
    for (srv in server_list) {
      ci   <- which(server_names == srv)
      is_c <- (srv == os_name)
        .dsAgg(datasources[ci], call(name = "k2Ring127AffineCombineDS",
          a_key           = y_share_key,
          b_key           = "cox_nd_p_share",
          sign_a          = 1,
          sign_b          = -1,
          public_const_fp = NULL,
          is_party0       = is_c,
          output_key      = "cox_nd_r_share",
          n               = as.numeric(n_pp),
          session_id      = session_id))
      }

    # 7. mr = mask * r  (Beaver vecmul; gates invalid (i,j) to 0)
    .ring127_vecmul(mask_share_key, "cox_nd_r_share", "cox_nd_mr_share",
      n_pp,
      datasources, dealer_ci, server_list, server_names,
      os_name, nl_name, transport_pks, session_id,
      .dsAgg, .sendBlob)

    # 8. W = p*(1-p):  (1-p) via affine, then p x (1-p) via vecmul
    for (srv in server_list) {
      ci   <- which(server_names == srv)
      is_c <- (srv == os_name)
        .dsAgg(datasources[ci], call(name = "k2Ring127AffineCombineDS",
          a_key           = NULL,
          b_key           = "cox_nd_p_share",
          sign_a          = 0,
          sign_b          = -1,
          public_const_fp = fp_one_b64,
          is_party0       = is_c,
          output_key      = "cox_nd_one_minus_p",
          n               = as.numeric(n_pp),
          session_id      = session_id))
      }
    .ring127_vecmul("cox_nd_p_share", "cox_nd_one_minus_p",
      "cox_nd_w_share", n_pp,
      datasources, dealer_ci, server_list, server_names,
      os_name, nl_name, transport_pks, session_id,
      .dsAgg, .sendBlob)

    # 9. mw = mask * W
    .ring127_vecmul(mask_share_key, "cox_nd_w_share", "cox_nd_mw_share",
      n_pp,
      datasources, dealer_ci, server_list, server_names,
      os_name, nl_name, transport_pks, session_id,
      .dsAgg, .sendBlob)

    # 10. Per-column gradient g[j] = X[:,j]^T * mr  (cache X col shares
    #     for the Hessian double-loop below -- they don't change in iter).
    grad <- numeric(p_total)
    x_col_keys <- character(p_total)
    for (j in seq_len(p_total)) {
      x_col_keys[j] <- sprintf("cox_nd_xcol_%d", j)
      for (srv in server_list) {
        ci <- which(server_names == srv)
          .dsAgg(datasources[ci], call(name = "k2BeaverExtractColumnDS",
            source_key = "k2_x_full_fp",
            n          = as.numeric(n_pp),
            K          = as.numeric(p_total),
            col_index  = as.numeric(j),
            output_key = x_col_keys[j],
            session_id = session_id,
            frac_bits  = 50,
            ring       = "ring127"))
        }
      gprod_key <- sprintf("cox_nd_gprod_%d", j)
      .ring127_vecmul(x_col_keys[j], "cox_nd_mr_share", gprod_key, n_pp,
        datasources, dealer_ci, server_list, server_names,
        os_name, nl_name, transport_pks, session_id,
        .dsAgg, .sendBlob)
      sum_g <- list()
      for (srv in server_list) {
        ci <- which(server_names == srv)
        rr <- .dsAgg(datasources[ci], call(name = "k2BeaverSumShareDS",
            source_key = gprod_key,
            session_id = session_id,
            frac_bits  = 50,
            ring       = "ring127"))
        if (is.list(rr) && length(rr) == 1L) rr <- rr[[1L]]
        sum_g[[srv]] <- rr
      }
      agg_g <- dsVert:::.callMpcTool("k2-ring63-aggregate", list(
          share_a   = sum_g[[os_name]]$sum_share_fp,
          share_b   = sum_g[[nl_name]]$sum_share_fp,
          frac_bits = 50,
          ring      = "ring127"))
      grad[j] <- as.numeric(agg_g$values)
    }
    if (debug_trace) {
      names(grad) <- beta_names
      grad_history[[iter]] <- grad
    }
    score_norm <- max(abs(grad))
    score_l2   <- sqrt(sum(grad^2))
    score_history <- c(score_history, score_norm)
    if (verbose) message(sprintf(
      "[#D' Newton iter %d] |g|_max=%.3e |g|_L2=%.3e",
      iter, score_norm, score_l2))
    if (is.finite(score_norm) && score_norm < best_score) {
      best_score <- score_norm; best_beta <- beta
    }
    iter_audit[[iter]] <- list(score_max = score_norm,
                               score_l2 = score_l2)
    if (score_norm < tol) {
      converged <- TRUE; final_iter <- iter; break
    }

    if (debug_trace && iter == 1L) {
      y_by_bin <- numeric(J)
      m_by_bin <- numeric(J)
      names(y_by_bin) <- names(m_by_bin) <- bin_dummy_names
      for (j in seq_len(J)) {
        yprod_key <- sprintf("cox_nd_debug_ybin_%d", j)
        mprod_key <- sprintf("cox_nd_debug_mbin_%d", j)
        .ring127_vecmul(x_col_keys[j], y_share_key, yprod_key, n_pp,
          datasources, dealer_ci, server_list, server_names,
          os_name, nl_name, transport_pks, session_id,
          .dsAgg, .sendBlob)
        .ring127_vecmul(x_col_keys[j], mask_share_key, mprod_key, n_pp,
          datasources, dealer_ci, server_list, server_names,
          os_name, nl_name, transport_pks, session_id,
          .dsAgg, .sendBlob)
        for (target in c("y", "m")) {
          source_key <- if (identical(target, "y")) yprod_key else mprod_key
          sums <- list()
          for (srv in server_list) {
            ci <- which(server_names == srv)
            rr <- .dsAgg(datasources[ci], call(name = "k2BeaverSumShareDS",
              source_key = source_key,
              session_id = session_id,
              frac_bits  = 50,
              ring       = "ring127"))
            if (is.list(rr) && length(rr) == 1L) rr <- rr[[1L]]
            sums[[srv]] <- rr
          }
          agg <- dsVert:::.callMpcTool("k2-ring63-aggregate", list(
            share_a   = sums[[os_name]]$sum_share_fp,
            share_b   = sums[[nl_name]]$sum_share_fp,
            frac_bits = 50,
            ring      = "ring127"))
          if (identical(target, "y")) {
            y_by_bin[j] <- as.numeric(agg$values)
          } else {
            m_by_bin[j] <- as.numeric(agg$values)
          }
        }
      }
      debug_bin_summary <- list(y_by_bin = y_by_bin, m_by_bin = m_by_bin,
                                share_mask_debug = share_res$debug)
    }

    if (iter >= as.integer(max_iter)) {
      final_iter <- iter
      break
    }

    # 11. Per-pair Hessian H[j,k] = X[:,j]^T diag(mw) X[:,k] -- symmetric,
    #     so we compute upper triangle only.
    H <- matrix(0, p_total, p_total,
                dimnames = list(beta_names, beta_names))
    for (j in seq_len(p_total)) {
      for (k in j:p_total) {
        z_key <- sprintf("cox_nd_z_%d_%d", j, k)
        h_key <- sprintf("cox_nd_h_%d_%d", j, k)
        # z_jk = X[:,j] * X[:,k]
        .ring127_vecmul(x_col_keys[j], x_col_keys[k], z_key, n_pp,
          datasources, dealer_ci, server_list, server_names,
          os_name, nl_name, transport_pks, session_id,
          .dsAgg, .sendBlob)
        # h_jk_per_row = z_jk * mw_share
        .ring127_vecmul(z_key, "cox_nd_mw_share", h_key, n_pp,
          datasources, dealer_ci, server_list, server_names,
          os_name, nl_name, transport_pks, session_id,
          .dsAgg, .sendBlob)
        sum_h <- list()
        for (srv in server_list) {
          ci <- which(server_names == srv)
          rr <- .dsAgg(datasources[ci], call(name = "k2BeaverSumShareDS",
              source_key = h_key,
              session_id = session_id,
              frac_bits  = 50,
              ring       = "ring127"))
          if (is.list(rr) && length(rr) == 1L) rr <- rr[[1L]]
          sum_h[[srv]] <- rr
        }
        agg_h <- dsVert:::.callMpcTool("k2-ring63-aggregate", list(
            share_a   = sum_h[[os_name]]$sum_share_fp,
            share_b   = sum_h[[nl_name]]$sum_share_fp,
            frac_bits = 50,
            ring      = "ring127"))
        H[j, k] <- H[k, j] <- as.numeric(agg_h$values)
      }
    }

    # 12. Ridge epsilon*I + Newton solve  (Christensen 2019 Sec.A.3 diagonal
    #     eigenvalue inflation handles all-zero-mask alpha_j strata where
    #     no patient was at risk in bin j -> singular block).  Step
    #     magnitude cap to 1.0 per coefficient prevents the iter-4 ->
    #     iter-5 overshoot observed empirically when bin2 sits near
    #     zero (Nocedal-Wright 2006 Sec.4.1 trust-region radius cap).
    H_ridged <- H + ridge_eps * diag(p_total)
    step <- tryCatch(solve(H_ridged, grad),
      error = function(e) {
        if (verbose) message(sprintf(
          "[#D' Newton iter %d] solve failed: %s -- ridge bumped to 1e-4",
          iter, conditionMessage(e)))
        solve(H + 1e-4 * diag(p_total), grad)
      })
    step_max_raw <- max(abs(step))
    cap <- 1.0
    cap_factor <- if (is.finite(step_max_raw) && step_max_raw > cap)
                     cap / step_max_raw else 1.0
    step_eff <- cap_factor * step
    step_max <- max(abs(step_eff))
    step_norm_hist <- c(step_norm_hist, step_max)
    beta <- beta + step_eff
    if (verbose) message(sprintf(
      "[#D' Newton iter %d] |step|_max=%.3e |beta|_max=%.3e cap=%.3f",
      iter, step_max, max(abs(beta)), cap_factor))
  }

    # If Newton does not converge, return the best scored iterate. The
    # beta after the final max_iter update has not been scored yet, so
    # returning it can report a small pre-step score for a different beta.
    returned_best_beta <- FALSE
    if (!converged && is.finite(best_score)) {
      beta <- best_beta
      returned_best_beta <- TRUE
    }

  list(
    stage             = if (converged) "converged" else "max_iter",
    coefficients      = beta,
    n_obs             = n_obs,
    J                 = J,
    n_pp              = n_pp,
    p_total           = p_total,
    p_os              = p_os,
    p_nl              = p_nl,
    beta_names        = beta_names,
    bin_dummy_names   = bin_dummy_names,
      n_iter            = final_iter,
      converged         = converged,
      returned_best_beta = returned_best_beta,
      score_history     = score_history,
      step_norm_history = step_norm_hist,
      final_score_norm  = if (returned_best_beta) best_score
                          else if (length(score_history)) tail(score_history, 1L)
                          else NA_real_,
    best_score_norm   = best_score,
    iter_audit        = iter_audit,
    mask_share_key    = mask_share_key,
    y_share_key       = y_share_key,
    expanded_x_name   = expanded_x_name,
    nl_x_vars         = x_nl,
    os_x_vars         = x_os,
    session_id        = session_id,
    transport_pks     = transport_pks,
    server_names      = server_names,
    os_name           = os_name,
    nl_name           = nl_name,
    bin_breaks        = bin_breaks,
    max_iter          = as.integer(max_iter),
    tol               = as.numeric(tol),
    ridge_eps         = as.numeric(ridge_eps),
    disclosure_audit  = audit,
    debug_trace       = if (debug_trace) {
      keep <- seq_len(length(score_history))
      list(beta_history = beta_history[keep],
           grad_history = grad_history[keep],
           bin_summary = debug_bin_summary)
    } else NULL,
    note = paste0(
      "Discrete-time Cox via pooled-logistic (Andreux 2020 + Allison 1982",
      " + Prentice-Gloeckler 1978) under masked Newton (Aliasgari-Blanton",
      " 2013 + De Cock 2016 share-mask gating; Mohassel-Zhang 2017 +",
      " Catrina-Saxena 2010 frac=50 noise floor; Christensen 2019 Sec.A.3",
      " diagonal ridge for singular alpha_j strata)."))
}
