#' Internal Ring127 Chebyshev-Horner + Beaver-vecmul helpers reusable
#' across joint-Newton orchestrations (Multinom, Ordinal). Mirror of
#' \code{.exp127_round} / \code{.run_beaver_vecmul_ring127} from
#' \code{ds.vertCox.R} but parameterised on \code{input_key} /
#' \code{output_key} so callers can maintain K-1 parallel share slots.
#'
#' @keywords internal
#' @noRd

.ring127_exp_coef_cache <- new.env(parent = emptyenv())
.ring127_recip_coef_cache <- new.env(parent = emptyenv())
.ring127_log_coef_cache <- new.env(parent = emptyenv())
.ring127_log_wide_coef_cache <- new.env(parent = emptyenv())
.ring127_half_fp_cache <- new.env(parent = emptyenv())

# Fresh, non-secret namespace for one high-level Ring127 invocation. Exact-GC
# destinations are one-shot by design; a CSPRNG tag prevents a later call in
# the same session from colliding with an earlier intermediate.
.ring127_invocation_tag <- function() {
  sub("^op_", "", .dsvert_exact_gc_new_context()$operation_id)
}

# FP(0.5) for share-side argument reduction exp(x) = exp(x/2)^2.
# Cached once per session to avoid repeated Go round-trips.
.ring127_get_half_fp <- function() {
  if (is.null(.ring127_half_fp_cache$v)) {
    r <- dsVert:::.callMpcTool("k2-float-to-fp",
      list(values = array(0.5, dim = 1L), frac_bits = 50, ring = "ring127"))
    .ring127_half_fp_cache$v <- r$fp_data
  }
  .ring127_half_fp_cache$v
}

# Standard base64 -> base64url (Opal/Rock DSL parser chokes on `=`, `+`, `/`
# inside double-quoted string literals -- documented in
# dsVert/R/mpcUtils.R:175 alongside the inverse `.base64url_to_base64`
# helper). Mirror of the local closure used inside ds.vertCox.R line 231;
# extracted here so joint-Newton orchestrations share the canonical form.
.to_b64url <- function(x) {
  if (is.null(x) || !nzchar(x)) return(x)
  chartr("+/", "-_", sub("=+$", "", x, perl = TRUE))
}

# Run one independent local Ring127 phase at every computation peer in a
# single DSI aggregate cycle. The caller still invokes this helper once per
# dependency edge, so operations that consume a key written by an earlier
# phase are never coalesced into the same request.
.ring127_local_fanout <- function(datasources, server_list, server_names,
                                  y_server, .dsAgg, make_call) {
  conns <- match(server_list, server_names)
  if (anyNA(conns)) {
    stop("Ring127 computation peer is not connected", call. = FALSE)
  }
  expressions <- stats::setNames(lapply(server_list, function(server) {
    make_call(server, identical(server, y_server))
  }), server_list)
  .dsAgg(datasources[conns], expressions)
  invisible(NULL)
}

# Exact-GC helpers accept a low-level DSI aggregate callback and add their own
# strict result contract.  Ring127 callers pass an already-strict two-argument
# phase wrapper, so adapt the calling convention without forwarding duplicate
# async/error controls into that wrapper.
.ring127_aggregate_adapter <- function(.dsAgg) {
  force(.dsAgg)
  function(conns, expr, ...) .dsAgg(conns, expr)
}

# Establish the exact-GC transport once for one explicit high-level primitive.
# Readiness is propagated only through that call tree: there is deliberately no
# package-global cache whose lifetime could cross a session or peer binding.
.ring127_transport_once <- function(transport_ready, datasources, server_list,
                                    server_names, session_id, .dsAgg) {
  if (isTRUE(transport_ready)) return(TRUE)
  servers <- match(server_list, server_names)
  if (length(servers) != 2L || anyNA(servers) || anyDuplicated(servers)) {
    stop("Ring127 exact transport requires exactly two connected peers",
         call. = FALSE)
  }
  .dsvert_setup_exact_gc_transport(
    datasources = datasources, server_names = server_names,
    servers = servers, session_id = session_id,
    .aggregate = .ring127_aggregate_adapter(.dsAgg))
  TRUE
}

# Ring127 Beaver vecmul on shares under arbitrary session keys.
# Mirrors dsVertClient:::.run_beaver_vecmul_ring127 (which is a closure
# inside ds.vertCox.R and not callable externally).
.ring127_vecmul <- function(x_key, y_key, output_key, n,
                            datasources, dealer_ci, server_list,
                            server_names, y_server, nl, transport_pks,
                            session_id, .dsAgg, .sendBlob,
                            transport_ready = FALSE,
                            destination_fresh = FALSE) {
  n_int <- as.integer(n)
  all_ci <- vapply(server_list, function(s) which(server_names == s),
                    integer(1L))
  exact_output_key <- if (isTRUE(destination_fresh)) {
    output_key
  } else {
    paste0("__r127_mul_", .ring127_invocation_tag())
  }
  .dsvert_exact_gc_vecmul_run(
    datasources = datasources, server_names = server_names,
    servers = all_ci, session_id = session_id, total_n = n_int,
    x_key = x_key, y_key = y_key, output_key = exact_output_key,
    transport_ready = transport_ready,
    .aggregate = .ring127_aggregate_adapter(.dsAgg))
  if (!isTRUE(destination_fresh)) {
    .ring127_local_fanout(
      datasources, server_list, server_names, y_server, .dsAgg,
      function(server, is_coord) call(
        name = "k2Ring127AffineCombineDS",
        a_key = exact_output_key, b_key = NULL,
        sign_a = 1, sign_b = 0, public_const_fp = NULL,
        is_party0 = is_coord, output_key = output_key,
        n = as.numeric(n_int), session_id = session_id))
  }
  invisible(NULL)
}

# Exact multiplication by a public fixed-point scalar. The scalar is installed
# as additive shares (entirely on one peer, zero on the other), then routed
# through the same exact-GC truncation adapter as a private vecmul. This avoids
# truncating random additive shares independently.
.ring127_exact_public_scale <- function(
    in_key, scalar_fp, output_key, n,
    datasources, dealer_ci, server_list, server_names, y_server, nl,
    transport_pks, session_id, .dsAgg, .sendBlob,
    transport_ready = FALSE,
    destination_fresh = FALSE) {
  n_int <- as.integer(n)
  tag <- .ring127_invocation_tag()
  constant_key <- paste0("__r127_exact_public_", tag)
  .ring127_local_fanout(
    datasources, server_list, server_names, y_server, .dsAgg,
    function(server, is_coord) call(
      name = "k2Ring127AffineCombineDS",
      a_key = NULL, b_key = NULL, sign_a = 0, sign_b = 0,
      public_const_fp = scalar_fp,
      is_party0 = is_coord,
      output_key = constant_key, n = as.numeric(n_int),
      session_id = session_id))
  .ring127_vecmul(
    in_key, constant_key, output_key, n_int,
    datasources, dealer_ci, server_list, server_names, y_server, nl,
    transport_pks, session_id, .dsAgg, .sendBlob,
    transport_ready = transport_ready,
    destination_fresh = destination_fresh)
  invisible(output_key)
}

# Clenshaw exp Horner on Ring127 shares with keyed IO. Takes the
# Chebyshev coefficients from the Go binary once per session and runs
# `degree + 1` Beaver vecmul rounds against the input share.
.ring127_exp_round_keyed <- function(in_key, out_key, n,
                                     datasources, dealer_ci, server_list,
                                     server_names, y_server, nl,
                                     transport_pks, session_id,
                                     .dsAgg, .sendBlob,
                                     transport_ready = FALSE) {
  n_int <- as.integer(n)
  transport_ready <- .ring127_transport_once(
    transport_ready, datasources, server_list, server_names,
    session_id, .dsAgg)
  if (is.null(.ring127_exp_coef_cache$coef_res)) {
    .ring127_exp_coef_cache$coef_res <- dsVert:::.callMpcTool(
      "k2-exp127-get-coeffs", list(frac_bits = 50))
  }
  coef_res <- .ring127_exp_coef_cache$coef_res
  # Opal DSL "==" parser fix -- strip base64 padding client-side;
  # DS functions re-pad before decoding.
  coef_res_one_over_a <- .to_b64url(coef_res$one_over_a)
  degree <- as.integer(coef_res$degree)
  all_coeffs_raw <- jsonlite::base64_dec(coef_res$coeffs)
  c_b64 <- vapply(seq_len(degree + 1L), function(idx) {
    s <- (idx - 1L) * 16L + 1L; e <- s + 15L
    .to_b64url(jsonlite::base64_enc(all_coeffs_raw[s:e]))
  }, character(1))

  tag <- .ring127_invocation_tag()
  tmp_y    <- paste0("__r127_ey_",    tag)
  tmp_twoY <- paste0("__r127_etwoY_", tag)
  tmp_bB   <- paste0("__r127_ebB_",   tag)
  tmp_bA   <- paste0("__r127_ebA_",   tag)

  .ring127_exact_public_scale(
    in_key, coef_res_one_over_a, tmp_y, n_int,
    datasources, dealer_ci, server_list, server_names, y_server, nl,
    transport_pks, session_id, .dsAgg, .sendBlob,
    transport_ready = transport_ready, destination_fresh = TRUE)
  .ring127_local_fanout(
    datasources, server_list, server_names, y_server, .dsAgg,
    function(server, is_coord) call(name = "k2Ring127AffineCombineDS",
      a_key = tmp_y, b_key = tmp_y, sign_a = 1, sign_b = 1,
      public_const_fp = NULL, is_party0 = is_coord,
      output_key = tmp_twoY, n = as.numeric(n_int), session_id = session_id))
  .ring127_local_fanout(
    datasources, server_list, server_names, y_server, .dsAgg,
    function(server, is_coord) call(name = "k2Ring127AffineCombineDS",
      a_key = NULL, b_key = NULL, sign_a = 0, sign_b = 0,
      public_const_fp = c_b64[degree + 1L], is_party0 = is_coord,
      output_key = tmp_bB, n = as.numeric(n_int), session_id = session_id))
  .ring127_local_fanout(
    datasources, server_list, server_names, y_server, .dsAgg,
    function(server, is_coord) call(name = "k2Ring127AffineCombineDS",
      a_key = NULL, b_key = NULL, sign_a = 0, sign_b = 0,
      public_const_fp = NULL, is_party0 = is_coord,
      output_key = tmp_bA, n = as.numeric(n_int), session_id = session_id))

  slot_B <- tmp_bB; slot_A <- tmp_bA
  product_stage <- 0L
  for (k in seq.int(degree - 1L, 1L)) {
    product_stage <- product_stage + 1L
    product_key <- sprintf("__r127_emul%02d_%s", product_stage, tag)
    .ring127_vecmul(tmp_twoY, slot_B, product_key, n_int,
                    datasources, dealer_ci, server_list, server_names,
                    y_server, nl, transport_pks, session_id,
                    .dsAgg, .sendBlob, transport_ready = transport_ready,
                    destination_fresh = TRUE)
    .ring127_local_fanout(
      datasources, server_list, server_names, y_server, .dsAgg,
      function(server, is_coord) call(name = "k2Ring127AffineCombineDS",
        a_key = product_key, b_key = slot_A, sign_a = 1, sign_b = -1,
        public_const_fp = c_b64[k + 1L], is_party0 = is_coord,
        output_key = slot_A, n = as.numeric(n_int), session_id = session_id))
    swap <- slot_A; slot_A <- slot_B; slot_B <- swap
  }
  product_stage <- product_stage + 1L
  product_key <- sprintf("__r127_emul%02d_%s", product_stage, tag)
  .ring127_vecmul(tmp_y, slot_B, product_key, n_int,
                  datasources, dealer_ci, server_list, server_names,
                  y_server, nl, transport_pks, session_id,
                  .dsAgg, .sendBlob, transport_ready = transport_ready,
                  destination_fresh = TRUE)
  .ring127_local_fanout(
    datasources, server_list, server_names, y_server, .dsAgg,
    function(server, is_coord) call(name = "k2Ring127AffineCombineDS",
      a_key = product_key, b_key = slot_A, sign_a = 1, sign_b = -1,
      public_const_fp = c_b64[1L], is_party0 = is_coord,
      output_key = out_key, n = as.numeric(n_int), session_id = session_id))
  invisible(NULL)
}

# Coefficient cache for the direct-sigmoid Chebyshev primitive.
.ring127_sigmoid_coef_cache <- new.env(parent = emptyenv())

# Ring127 direct share-domain sigmoid(x) = 1/(1+exp(-x)) via a single
# Chebyshev Clenshaw pass on the shares.
# Structural clone of .ring127_exp_round_keyed: reveal-free (relays
# dcf_masked ZERO times), dealer-free (only .ring127_vecmul = IKNP), and
# env-agnostic (is_party0 derived per server from `server == y_server`).
# Replaces exp127+recip127 (~85 Beaver rounds) with `degree` (=29) rounds
# for the GLM logistic link. Domain contract: eta in [-8, 8]; c_0 carries
# the +0.5 sigmoid baseline (no extra affine term). Does NOT touch the
# high-accuracy exp127/recip127 used by the NB/ordinal/multinom/Cox family.
.ring127_sigmoid_round_keyed <- function(in_key, out_key, n,
                                         datasources, dealer_ci, server_list,
                                         server_names, y_server, nl,
                                         transport_pks, session_id,
                                         .dsAgg, .sendBlob,
                                         transport_ready = FALSE) {
  n_int <- as.integer(n)
  transport_ready <- .ring127_transport_once(
    transport_ready, datasources, server_list, server_names,
    session_id, .dsAgg)
  if (is.null(.ring127_sigmoid_coef_cache$coef_res)) {
    .ring127_sigmoid_coef_cache$coef_res <- dsVert:::.callMpcTool(
      "k2-sigmoid127-get-coeffs", list(frac_bits = 50))
  }
  coef_res <- .ring127_sigmoid_coef_cache$coef_res
  coef_res_one_over_a <- .to_b64url(coef_res$one_over_a)
  degree <- as.integer(coef_res$degree)
  all_coeffs_raw <- jsonlite::base64_dec(coef_res$coeffs)
  c_b64 <- vapply(seq_len(degree + 1L), function(idx) {
    s <- (idx - 1L) * 16L + 1L; e <- s + 15L
    .to_b64url(jsonlite::base64_enc(all_coeffs_raw[s:e]))
  }, character(1))

  tag <- .ring127_invocation_tag()
  tmp_y    <- paste0("__r127_sgy_",    tag)
  tmp_twoY <- paste0("__r127_sgtwoY_", tag)
  tmp_bB   <- paste0("__r127_sgbB_",   tag)
  tmp_bA   <- paste0("__r127_sgbA_",   tag)

  .ring127_exact_public_scale(
    in_key, coef_res_one_over_a, tmp_y, n_int,
    datasources, dealer_ci, server_list, server_names, y_server, nl,
    transport_pks, session_id, .dsAgg, .sendBlob,
    transport_ready = transport_ready, destination_fresh = TRUE)
  .ring127_local_fanout(
    datasources, server_list, server_names, y_server, .dsAgg,
    function(server, is_coord) call(name = "k2Ring127AffineCombineDS",
      a_key = tmp_y, b_key = tmp_y, sign_a = 1, sign_b = 1,
      public_const_fp = NULL, is_party0 = is_coord,
      output_key = tmp_twoY, n = as.numeric(n_int), session_id = session_id))
  .ring127_local_fanout(
    datasources, server_list, server_names, y_server, .dsAgg,
    function(server, is_coord) call(name = "k2Ring127AffineCombineDS",
      a_key = NULL, b_key = NULL, sign_a = 0, sign_b = 0,
      public_const_fp = c_b64[degree + 1L], is_party0 = is_coord,
      output_key = tmp_bB, n = as.numeric(n_int), session_id = session_id))
  .ring127_local_fanout(
    datasources, server_list, server_names, y_server, .dsAgg,
    function(server, is_coord) call(name = "k2Ring127AffineCombineDS",
      a_key = NULL, b_key = NULL, sign_a = 0, sign_b = 0,
      public_const_fp = NULL, is_party0 = is_coord,
      output_key = tmp_bA, n = as.numeric(n_int), session_id = session_id))

  slot_B <- tmp_bB; slot_A <- tmp_bA
  product_stage <- 0L
  for (k in seq.int(degree - 1L, 1L)) {
    product_stage <- product_stage + 1L
    product_key <- sprintf("__r127_sgmul%02d_%s", product_stage, tag)
    .ring127_vecmul(tmp_twoY, slot_B, product_key, n_int,
                    datasources, dealer_ci, server_list, server_names,
                    y_server, nl, transport_pks, session_id,
                    .dsAgg, .sendBlob, transport_ready = transport_ready,
                    destination_fresh = TRUE)
    .ring127_local_fanout(
      datasources, server_list, server_names, y_server, .dsAgg,
      function(server, is_coord) call(name = "k2Ring127AffineCombineDS",
        a_key = product_key, b_key = slot_A, sign_a = 1, sign_b = -1,
        public_const_fp = c_b64[k + 1L], is_party0 = is_coord,
        output_key = slot_A, n = as.numeric(n_int), session_id = session_id))
    swap <- slot_A; slot_A <- slot_B; slot_B <- swap
  }
  product_stage <- product_stage + 1L
  product_key <- sprintf("__r127_sgmul%02d_%s", product_stage, tag)
  .ring127_vecmul(tmp_y, slot_B, product_key, n_int,
                  datasources, dealer_ci, server_list, server_names,
                  y_server, nl, transport_pks, session_id,
                  .dsAgg, .sendBlob, transport_ready = transport_ready,
                  destination_fresh = TRUE)
  .ring127_local_fanout(
    datasources, server_list, server_names, y_server, .dsAgg,
    function(server, is_coord) call(name = "k2Ring127AffineCombineDS",
      a_key = product_key, b_key = slot_A, sign_a = 1, sign_b = -1,
      public_const_fp = c_b64[1L], is_party0 = is_coord,
      output_key = out_key, n = as.numeric(n_int), session_id = session_id))
  invisible(NULL)
}

# Coefficient cache for the direct-softplus Chebyshev primitive.
.ring127_softplus_coef_cache <- new.env(parent = emptyenv())

.ring127_glm_softplus_manifest_mul <- function(
    stage, invocation_id, output_key, n, datasources, server_list, server_names,
    session_id, .dsAgg) {
  .dsvert_block_retired_remote_route("legacy_glm")
  stage <- as.integer(stage)
  if (length(stage) != 1L || is.na(stage) || stage < 0L || stage > 36L) {
    stop("Invalid producer-bound GLM softplus stage", call. = FALSE)
  }
  selected <- match(server_list, server_names)
  if (length(selected) != 2L || anyNA(selected) || anyDuplicated(selected)) {
    stop("GLM softplus requires exactly two connected computation peers",
         call. = FALSE)
  }
  selected_names <- server_names[selected]
  conns <- datasources[selected]
  names(conns) <- selected_names
  requests <- stats::setNames(lapply(selected_names, function(server) {
    call(name = "exactGCGLMSoftplusPrepareDS", stage = stage,
         invocation_id = invocation_id,
         session_id = session_id)
  }), selected_names)
  prepared <- .dsvert_fanout_by_site(
    conns, requests,
    operation = "GLM producer-bound softplus multiplication preparation",
    .aggregate = .dsAgg)
  manifests <- stats::setNames(lapply(selected_names, function(server) {
    value <- prepared[[server]]
    if (!is.list(value) || !is.character(value$manifest_handle) ||
        length(value$manifest_handle) != 1L) {
      stop("A GLM peer did not mint its softplus multiplication manifest.",
           call. = FALSE)
    }
    value
  }), selected_names)
  .dsvert_exact_gc_vecmul_run(
    datasources = datasources, server_names = server_names,
    servers = selected, session_id = session_id, total_n = as.integer(n),
    output_key = output_key, input_manifests = manifests,
    transport_ready = TRUE, .aggregate = .dsAgg)
  invisible(output_key)
}

# Ring127 direct share-domain softplus(x) = log(1+exp(x)) via a single
# Chebyshev Clenshaw pass on the shares (binomial-deviance link). Structural
# clone of .ring127_sigmoid_round_keyed: reveal-free (dcf_masked 0×), dealer-
# free (IKNP), env-agnostic. Domain contract: eta in [-8, 8]. degree=36 rounds.
.ring127_softplus_round_keyed <- function(in_key, out_key, n,
                                          datasources, dealer_ci, server_list,
                                          server_names, y_server, nl,
                                          transport_pks, session_id,
                                          .dsAgg, .sendBlob,
                                          transport_ready = FALSE,
                                          producer_bound_glm = FALSE) {
  n_int <- as.integer(n)
  transport_ready <- .ring127_transport_once(
    transport_ready, datasources, server_list, server_names,
    session_id, .dsAgg)
  if (is.null(.ring127_softplus_coef_cache$coef_res)) {
    .ring127_softplus_coef_cache$coef_res <- dsVert:::.callMpcTool(
      "k2-softplus127-get-coeffs", list(frac_bits = 50))
  }
  coef_res <- .ring127_softplus_coef_cache$coef_res
  coef_res_one_over_a <- .to_b64url(coef_res$one_over_a)
  degree <- as.integer(coef_res$degree)
  all_coeffs_raw <- jsonlite::base64_dec(coef_res$coeffs)
  c_b64 <- vapply(seq_len(degree + 1L), function(idx) {
    s <- (idx - 1L) * 16L + 1L; e <- s + 15L
    .to_b64url(jsonlite::base64_enc(all_coeffs_raw[s:e]))
  }, character(1))

  tag <- .ring127_invocation_tag()
  tmp_y    <- paste0("__r127_spy_",    tag)
  tmp_twoY <- paste0("__r127_sptwoY_", tag)
  tmp_bB   <- paste0("__r127_spbB_",   tag)
  tmp_bA   <- paste0("__r127_spbA_",   tag)

  if (isTRUE(producer_bound_glm)) {
    if (!identical(in_key, "k2_eta_share_fp") ||
        !identical(out_key, "softplus_share_fp") || degree != 36L) {
      stop("Producer-bound GLM softplus requires its fixed Ring127 contract",
           call. = FALSE)
    }
    constant_key <- paste0("__r127_spconst_", tag)
    .ring127_local_fanout(
      datasources, server_list, server_names, y_server, .dsAgg,
      function(server, is_coord) call(
        name = "k2Ring127AffineCombineDS",
        a_key = NULL, b_key = NULL, sign_a = 0, sign_b = 0,
        public_const_fp = coef_res_one_over_a,
        is_party0 = is_coord, output_key = constant_key,
        n = as.numeric(n_int), session_id = session_id))
    .ring127_glm_softplus_manifest_mul(
      0L, tag, tmp_y, n_int, datasources, server_list, server_names,
      session_id, .dsAgg)
  } else {
    .ring127_exact_public_scale(
      in_key, coef_res_one_over_a, tmp_y, n_int,
      datasources, dealer_ci, server_list, server_names, y_server, nl,
      transport_pks, session_id, .dsAgg, .sendBlob,
      transport_ready = transport_ready, destination_fresh = TRUE)
  }
  .ring127_local_fanout(
    datasources, server_list, server_names, y_server, .dsAgg,
    function(server, is_coord) call(name = "k2Ring127AffineCombineDS",
      a_key = tmp_y, b_key = tmp_y, sign_a = 1, sign_b = 1,
      public_const_fp = NULL, is_party0 = is_coord,
      output_key = tmp_twoY, n = as.numeric(n_int), session_id = session_id))
  .ring127_local_fanout(
    datasources, server_list, server_names, y_server, .dsAgg,
    function(server, is_coord) call(name = "k2Ring127AffineCombineDS",
      a_key = NULL, b_key = NULL, sign_a = 0, sign_b = 0,
      public_const_fp = c_b64[degree + 1L], is_party0 = is_coord,
      output_key = tmp_bB, n = as.numeric(n_int), session_id = session_id))
  .ring127_local_fanout(
    datasources, server_list, server_names, y_server, .dsAgg,
    function(server, is_coord) call(name = "k2Ring127AffineCombineDS",
      a_key = NULL, b_key = NULL, sign_a = 0, sign_b = 0,
      public_const_fp = NULL, is_party0 = is_coord,
      output_key = tmp_bA, n = as.numeric(n_int), session_id = session_id))

  slot_B <- tmp_bB; slot_A <- tmp_bA
  product_stage <- 0L
  for (k in seq.int(degree - 1L, 1L)) {
    product_stage <- product_stage + 1L
    product_key <- sprintf("__r127_spmul%02d_%s", product_stage, tag)
    if (isTRUE(producer_bound_glm)) {
      .ring127_glm_softplus_manifest_mul(
        product_stage, tag, product_key, n_int, datasources, server_list,
        server_names, session_id, .dsAgg)
    } else {
      .ring127_vecmul(tmp_twoY, slot_B, product_key, n_int,
                      datasources, dealer_ci, server_list, server_names,
                      y_server, nl, transport_pks, session_id,
                      .dsAgg, .sendBlob, transport_ready = transport_ready,
                      destination_fresh = TRUE)
    }
    .ring127_local_fanout(
      datasources, server_list, server_names, y_server, .dsAgg,
      function(server, is_coord) call(name = "k2Ring127AffineCombineDS",
        a_key = product_key, b_key = slot_A, sign_a = 1, sign_b = -1,
        public_const_fp = c_b64[k + 1L], is_party0 = is_coord,
        output_key = slot_A, n = as.numeric(n_int), session_id = session_id))
    swap <- slot_A; slot_A <- slot_B; slot_B <- swap
  }
  final_product <- sprintf("__r127_spmul36_%s", tag)
  if (isTRUE(producer_bound_glm)) {
    .ring127_glm_softplus_manifest_mul(
      36L, tag, final_product, n_int, datasources, server_list, server_names,
      session_id, .dsAgg)
  } else {
    .ring127_vecmul(tmp_y, slot_B, final_product, n_int,
                    datasources, dealer_ci, server_list, server_names,
                    y_server, nl, transport_pks, session_id,
                    .dsAgg, .sendBlob, transport_ready = transport_ready,
                    destination_fresh = TRUE)
  }
  .ring127_local_fanout(
    datasources, server_list, server_names, y_server, .dsAgg,
    function(server, is_coord) call(name = "k2Ring127AffineCombineDS",
      a_key = final_product, b_key = slot_A, sign_a = 1, sign_b = -1,
      public_const_fp = c_b64[1L], is_party0 = is_coord,
      output_key = out_key, n = as.numeric(n_int), session_id = session_id))
  invisible(NULL)
}

# Ring127 spline-less 1/x via Chebyshev-Horner initial guess + 6 NR
# iterations. Mirror of ds.vertCox.R's .recip127_round but keyed I/O.
.ring127_recip_round_keyed <- function(in_key, out_key, n,
                                       datasources, dealer_ci, server_list,
                                       server_names, y_server, nl,
                                       transport_pks, session_id,
                                       .dsAgg, .sendBlob,
                                       transport_ready = FALSE) {
  n_int <- as.integer(n)
  transport_ready <- .ring127_transport_once(
    transport_ready, datasources, server_list, server_names,
    session_id, .dsAgg)
  if (is.null(.ring127_recip_coef_cache$coef_res)) {
    .ring127_recip_coef_cache$coef_res <- dsVert:::.callMpcTool(
      "k2-recip127-get-coeffs", list(frac_bits = 50))
  }
  rc <- .ring127_recip_coef_cache$coef_res
  rc_one_over_half_range <- .to_b64url(rc$one_over_half_range)
  rc_neg_mid_over_half_range <- .to_b64url(rc$neg_mid_over_half_range)
  rc_two_fp <- .to_b64url(rc$two_fp)
  degree <- as.integer(rc$degree)
  nr_steps <- as.integer(rc$nr_steps)
  all_coeffs_raw <- jsonlite::base64_dec(rc$coeffs)
  c_b64 <- vapply(seq_len(degree + 1L), function(idx) {
    s <- (idx - 1L) * 16L + 1L; e <- s + 15L
    .to_b64url(jsonlite::base64_enc(all_coeffs_raw[s:e]))
  }, character(1))

  tag <- .ring127_invocation_tag()
  t_pre <- paste0("__r127_rp_tpre_", tag)
  t_key <- paste0("__r127_rp_t_",    tag)
  twoT  <- paste0("__r127_rp_twoT_", tag)
  bB    <- paste0("__r127_rp_bB_",   tag)
  bA    <- paste0("__r127_rp_bA_",   tag)
  y_seed <- paste0("__r127_rp_yseed_", tag)

  .ring127_exact_public_scale(
    in_key, rc_one_over_half_range, t_pre, n_int,
    datasources, dealer_ci, server_list, server_names, y_server, nl,
    transport_pks, session_id, .dsAgg, .sendBlob,
    transport_ready = transport_ready, destination_fresh = TRUE)
  .ring127_local_fanout(
    datasources, server_list, server_names, y_server, .dsAgg,
    function(server, is_coord) call(name = "k2Ring127AffineCombineDS",
      a_key = t_pre, b_key = NULL, sign_a = 1, sign_b = 0,
      public_const_fp = rc_neg_mid_over_half_range,
      is_party0 = is_coord, output_key = t_key,
      n = as.numeric(n_int), session_id = session_id))
  .ring127_local_fanout(
    datasources, server_list, server_names, y_server, .dsAgg,
    function(server, is_coord) call(name = "k2Ring127AffineCombineDS",
      a_key = t_key, b_key = t_key, sign_a = 1, sign_b = 1,
      public_const_fp = NULL, is_party0 = is_coord,
      output_key = twoT, n = as.numeric(n_int), session_id = session_id))
  .ring127_local_fanout(
    datasources, server_list, server_names, y_server, .dsAgg,
    function(server, is_coord) call(name = "k2Ring127AffineCombineDS",
      a_key = NULL, b_key = NULL, sign_a = 0, sign_b = 0,
      public_const_fp = c_b64[degree + 1L], is_party0 = is_coord,
      output_key = bB, n = as.numeric(n_int), session_id = session_id))
  .ring127_local_fanout(
    datasources, server_list, server_names, y_server, .dsAgg,
    function(server, is_coord) call(name = "k2Ring127AffineCombineDS",
      a_key = NULL, b_key = NULL, sign_a = 0, sign_b = 0,
      public_const_fp = NULL, is_party0 = is_coord,
      output_key = bA, n = as.numeric(n_int), session_id = session_id))

  slot_B <- bB; slot_A <- bA
  product_stage <- 0L
  for (k in seq.int(degree - 1L, 1L)) {
    product_stage <- product_stage + 1L
    product_key <- sprintf("__r127_rpmul%02d_%s", product_stage, tag)
    .ring127_vecmul(twoT, slot_B, product_key, n_int,
                    datasources, dealer_ci, server_list, server_names,
                    y_server, nl, transport_pks, session_id,
                    .dsAgg, .sendBlob, transport_ready = transport_ready,
                    destination_fresh = TRUE)
    .ring127_local_fanout(
      datasources, server_list, server_names, y_server, .dsAgg,
      function(server, is_coord) call(name = "k2Ring127AffineCombineDS",
        a_key = product_key, b_key = slot_A, sign_a = 1, sign_b = -1,
        public_const_fp = c_b64[k + 1L], is_party0 = is_coord,
        output_key = slot_A, n = as.numeric(n_int), session_id = session_id))
    swap <- slot_A; slot_A <- slot_B; slot_B <- swap
  }
  product_stage <- product_stage + 1L
  product_key <- sprintf("__r127_rpmul%02d_%s", product_stage, tag)
  .ring127_vecmul(t_key, slot_B, product_key, n_int,
                  datasources, dealer_ci, server_list, server_names,
                  y_server, nl, transport_pks, session_id,
                  .dsAgg, .sendBlob, transport_ready = transport_ready,
                  destination_fresh = TRUE)
  .ring127_local_fanout(
    datasources, server_list, server_names, y_server, .dsAgg,
    function(server, is_coord) call(name = "k2Ring127AffineCombineDS",
      a_key = product_key, b_key = slot_A, sign_a = 1, sign_b = -1,
      public_const_fp = c_b64[1L], is_party0 = is_coord,
      output_key = y_seed, n = as.numeric(n_int), session_id = session_id))

  y_cur <- y_seed
  for (iter in seq_len(nr_steps)) {
    xy <- sprintf("__r127_rp_xy%02d_%s", iter, tag)
    tmXY <- sprintf("__r127_rp_twoMinusXY%02d_%s", iter, tag)
    .ring127_vecmul(in_key, y_cur, xy, n_int,
                    datasources, dealer_ci, server_list, server_names,
                    y_server, nl, transport_pks, session_id,
                    .dsAgg, .sendBlob, transport_ready = transport_ready,
                    destination_fresh = TRUE)
    .ring127_local_fanout(
      datasources, server_list, server_names, y_server, .dsAgg,
      function(server, is_coord) call(name = "k2Ring127AffineCombineDS",
        a_key = NULL, b_key = xy, sign_a = 0, sign_b = -1,
        public_const_fp = rc_two_fp, is_party0 = is_coord,
        output_key = tmXY, n = as.numeric(n_int), session_id = session_id))
    next_slot <- sprintf("__r127_rp_nr%02d_%s", iter, tag)
    .ring127_vecmul(y_cur, tmXY, next_slot, n_int,
                    datasources, dealer_ci, server_list, server_names,
                    y_server, nl, transport_pks, session_id,
                    .dsAgg, .sendBlob, transport_ready = transport_ready,
                    destination_fresh = TRUE)
    y_cur <- next_slot
  }
  .ring127_local_fanout(
    datasources, server_list, server_names, y_server, .dsAgg,
    function(server, is_coord) call(name = "k2Ring127AffineCombineDS",
      a_key = y_cur, b_key = NULL, sign_a = 1, sign_b = 0,
      public_const_fp = NULL, is_party0 = is_coord,
      output_key = out_key, n = as.numeric(n_int), session_id = session_id))
  invisible(NULL)
}

# Argument-reduced exp for share-space: evaluates exp(x) for |x| <= ~10
# by scaling input by 1/2, running the core [-5,5] Chebyshev, and
# squaring the result via one Beaver vecmul.
#
# Mirrors Ring127ExpPlaintextExtended (dsVert/inst/dsvert-mpc/k2_exp127.go
# lines 87-119) but operates on secret-shared inputs. Preserves interior
# accuracy (rel <= 1e-12 per Trefethen ATAP Sec.8 bound) by keeping the
# Chebyshev core at a=5 rather than degrading coefficients at a=10.
#
# Use this in place of .ring127_exp_round_keyed when eta may leave the
# [-5, 5] band during Newton iterations (softmax / proportional-odds
# sigmoid) -- fixes bug #9 (max_step=NaN from exp overflow/wraparound at
# |eta|>5 on Ring127 fixed-point).
#
# Cost: one extra exact public-scale multiplication plus the final square.
.ring127_exp_round_keyed_extended <- function(in_key, out_key, n,
                                              datasources, dealer_ci, server_list,
                                              server_names, y_server, nl,
                                              transport_pks, session_id,
                                              .dsAgg, .sendBlob,
                                              transport_ready = FALSE) {
  n_int <- as.integer(n)
  transport_ready <- .ring127_transport_once(
    transport_ready, datasources, server_list, server_names,
    session_id, .dsAgg)
  half_fp_b64 <- .to_b64url(.ring127_get_half_fp())
  tag <- .ring127_invocation_tag()
  half_key <- paste0("__r127_exthalf_", tag)
  exp_half_key <- paste0("__r127_extexp_", tag)
  square_key <- paste0("__r127_extsquare_", tag)

  # Step 1: scale by 0.5 with exact joint truncation.
  .ring127_exact_public_scale(
    in_key, half_fp_b64, half_key, n_int,
    datasources, dealer_ci, server_list, server_names, y_server, nl,
    transport_pks, session_id, .dsAgg, .sendBlob,
    transport_ready = transport_ready, destination_fresh = TRUE)
  # Step 2: evaluate exp on half-input (interior Chebyshev [-5,5]).
  .ring127_exp_round_keyed(half_key, exp_half_key, n,
                           datasources, dealer_ci, server_list,
                           server_names, y_server, nl,
                           transport_pks, session_id,
                           .dsAgg, .sendBlob,
                           transport_ready = transport_ready)
  # Step 3: square via Beaver vecmul -> exp(x) = exp(x/2)^2.
  .ring127_vecmul(exp_half_key, exp_half_key, square_key, n_int,
                  datasources, dealer_ci, server_list, server_names,
                  y_server, nl, transport_pks, session_id,
                  .dsAgg, .sendBlob, transport_ready = transport_ready,
                  destination_fresh = TRUE)
  .ring127_local_fanout(
    datasources, server_list, server_names, y_server, .dsAgg,
    function(server, is_coord) call(name = "k2Ring127AffineCombineDS",
      a_key = square_key, b_key = NULL, sign_a = 1, sign_b = 0,
      public_const_fp = NULL, is_party0 = is_coord,
      output_key = out_key, n = as.numeric(n_int), session_id = session_id))
  invisible(NULL)
}

# Clenshaw log-shift Horner on Ring127 shares with keyed IO. Evaluates
# log(x) for x in the public Chebyshev core domain [Ring127LogShiftMin,
# Ring127LogShiftMax] = [1, 10] at fracBits=50 (rel <~ 1e-12 per
# Trefethen ATAP Sec.8 + Bernstein ellipse rhoapprox1.94 at degree 40).
#
# Mirrors .ring127_recip_round_keyed exactly through the Clenshaw stage --
# the affine "scale + offset" maps x in [1, 10] onto t in [-1, 1], then a
# 41-step Clenshaw recurrence assembles log(x) via Beaver vecmul +
# AffineCombine. No NR refinement needed: Chebyshev-only evaluation
# already achieves the Catrina-Saxena 2010 Sec.3.3 ULP floor at fracBits=50.
#
# Disclosure note (per K=2 OT-Beaver dishonest-majority threat model):
# coefficients + affine constants are public deterministic values from
# the Go init() -- distributing them leaks nothing. Beaver vecmul +
# AffineCombine pipeline is the same threat-model footing as exp/recip.
#
# Caller responsibility: ensure share-domain input encodes a value in
# [1, 10]. For NB full-regression theta MLE, mu+theta may exceed this range and
# the orchestrator must perform argument reduction first (e.g. divide
# by a known plaintext rescale factor and add the corresponding
# log-correction post-hoc, leveraging log(c*x) = log(c) + log(x)).
.ring127_log_round_keyed <- function(in_key, out_key, n,
                                     datasources, dealer_ci, server_list,
                                     server_names, y_server, nl,
                                     transport_pks, session_id,
                                     .dsAgg, .sendBlob,
                                     transport_ready = FALSE) {
  n_int <- as.integer(n)
  transport_ready <- .ring127_transport_once(
    transport_ready, datasources, server_list, server_names,
    session_id, .dsAgg)
  if (is.null(.ring127_log_coef_cache$coef_res)) {
    .ring127_log_coef_cache$coef_res <- dsVert:::.callMpcTool(
      "k2-log-shift-coeffs", list(frac_bits = 50))
  }
  rc <- .ring127_log_coef_cache$coef_res
  rc_one_over_half_range <- .to_b64url(rc$one_over_half_range)
  rc_neg_mid_over_half_range <- .to_b64url(rc$neg_mid_over_half_range)
  degree <- as.integer(rc$degree)
  all_coeffs_raw <- jsonlite::base64_dec(rc$coeffs)
  c_b64 <- vapply(seq_len(degree + 1L), function(idx) {
    s <- (idx - 1L) * 16L + 1L; e <- s + 15L
    .to_b64url(jsonlite::base64_enc(all_coeffs_raw[s:e]))
  }, character(1))

  tag <- .ring127_invocation_tag()
  t_pre <- paste0("__r127_lp_tpre_", tag)
  t_key <- paste0("__r127_lp_t_",    tag)
  twoT  <- paste0("__r127_lp_twoT_", tag)
  bB    <- paste0("__r127_lp_bB_",   tag)
  bA    <- paste0("__r127_lp_bA_",   tag)

  .ring127_exact_public_scale(
    in_key, rc_one_over_half_range, t_pre, n_int,
    datasources, dealer_ci, server_list, server_names, y_server, nl,
    transport_pks, session_id, .dsAgg, .sendBlob,
    transport_ready = transport_ready, destination_fresh = TRUE)
  .ring127_local_fanout(
    datasources, server_list, server_names, y_server, .dsAgg,
    function(server, is_coord) call(name = "k2Ring127AffineCombineDS",
      a_key = t_pre, b_key = NULL, sign_a = 1, sign_b = 0,
      public_const_fp = rc_neg_mid_over_half_range,
      is_party0 = is_coord, output_key = t_key,
      n = as.numeric(n_int), session_id = session_id))
  .ring127_local_fanout(
    datasources, server_list, server_names, y_server, .dsAgg,
    function(server, is_coord) call(name = "k2Ring127AffineCombineDS",
      a_key = t_key, b_key = t_key, sign_a = 1, sign_b = 1,
      public_const_fp = NULL, is_party0 = is_coord,
      output_key = twoT, n = as.numeric(n_int), session_id = session_id))
  .ring127_local_fanout(
    datasources, server_list, server_names, y_server, .dsAgg,
    function(server, is_coord) call(name = "k2Ring127AffineCombineDS",
      a_key = NULL, b_key = NULL, sign_a = 0, sign_b = 0,
      public_const_fp = c_b64[degree + 1L], is_party0 = is_coord,
      output_key = bB, n = as.numeric(n_int), session_id = session_id))
  .ring127_local_fanout(
    datasources, server_list, server_names, y_server, .dsAgg,
    function(server, is_coord) call(name = "k2Ring127AffineCombineDS",
      a_key = NULL, b_key = NULL, sign_a = 0, sign_b = 0,
      public_const_fp = NULL, is_party0 = is_coord,
      output_key = bA, n = as.numeric(n_int), session_id = session_id))

  slot_B <- bB; slot_A <- bA
  product_stage <- 0L
  for (k in seq.int(degree - 1L, 1L)) {
    product_stage <- product_stage + 1L
    product_key <- sprintf("__r127_lpmul%02d_%s", product_stage, tag)
    .ring127_vecmul(twoT, slot_B, product_key, n_int,
                    datasources, dealer_ci, server_list, server_names,
                    y_server, nl, transport_pks, session_id,
                    .dsAgg, .sendBlob, transport_ready = transport_ready,
                    destination_fresh = TRUE)
    .ring127_local_fanout(
      datasources, server_list, server_names, y_server, .dsAgg,
      function(server, is_coord) call(name = "k2Ring127AffineCombineDS",
        a_key = product_key, b_key = slot_A, sign_a = 1, sign_b = -1,
        public_const_fp = c_b64[k + 1L], is_party0 = is_coord,
        output_key = slot_A, n = as.numeric(n_int), session_id = session_id))
    swap <- slot_A; slot_A <- slot_B; slot_B <- swap
  }
  product_stage <- product_stage + 1L
  product_key <- sprintf("__r127_lpmul%02d_%s", product_stage, tag)
  .ring127_vecmul(t_key, slot_B, product_key, n_int,
                  datasources, dealer_ci, server_list, server_names,
                  y_server, nl, transport_pks, session_id,
                  .dsAgg, .sendBlob, transport_ready = transport_ready,
                  destination_fresh = TRUE)
  .ring127_local_fanout(
    datasources, server_list, server_names, y_server, .dsAgg,
    function(server, is_coord) call(name = "k2Ring127AffineCombineDS",
      a_key = product_key, b_key = slot_A, sign_a = 1, sign_b = -1,
      public_const_fp = c_b64[1L], is_party0 = is_coord,
      output_key = out_key, n = as.numeric(n_int), session_id = session_id))
  invisible(NULL)
}

# Argument-reduced log for share-space: evaluates log(x) for x in a wider
# operating range by rescaling the shared input by a public plaintext
# scale factor that maps x into the [1, 10] Chebyshev core, then
# correcting via log(c*x) = log(c) + log(x) post-hoc as a public
# constant added by party-0 in a final AffineCombine.
#
# The rescale factor `scale_fp_b64` is supplied by the caller as a public
# Ring127 FP value; the caller is responsible for selecting it so that
# scale*x in [1, 10] across the operating range. A common pattern for
# NB full-regression theta MLE: mu + theta in [theta, ~theta + exp(eta_max)]; if theta >= 1 and
# eta in [-5, 5] gives mu <= ~150, scale = 1/15 maps mu+theta into roughly
# [theta/15, 10]. The caller computes log_scale_correction = log(1/scale)
# at plaintext (it is public -- `scale` is a public const) and passes
# its FP encoding via `log_scale_correction_fp_b64` as the additive
# party-0 const in the final affine step.
#
# Disclosure note: scale_fp and log_scale_correction_fp are public
# plaintext constants; supplying them does not leak share state.
.ring127_log_round_keyed_extended <- function(in_key, out_key, n,
                                              scale_fp_b64,
                                              log_scale_correction_fp_b64,
                                              datasources, dealer_ci, server_list,
                                              server_names, y_server, nl,
                                              transport_pks, session_id,
                                              .dsAgg, .sendBlob,
                                              transport_ready = FALSE) {
  n_int <- as.integer(n)
  transport_ready <- .ring127_transport_once(
    transport_ready, datasources, server_list, server_names,
    session_id, .dsAgg)
  tag <- .ring127_invocation_tag()
  scaled_key <- paste0("__r127_logscaled_", tag)
  log_scaled_key <- paste0("__r127_logvalue_", tag)

  # Step 1: rescale with exact joint truncation. The caller guarantees the
  # scaled value falls in [1, 10].
  .ring127_exact_public_scale(
    in_key, scale_fp_b64, scaled_key, n_int,
    datasources, dealer_ci, server_list, server_names, y_server, nl,
    transport_pks, session_id, .dsAgg, .sendBlob,
    transport_ready = transport_ready, destination_fresh = TRUE)

  # Step 2: evaluate log on the rescaled share (Chebyshev core [1, 10]).
  .ring127_log_round_keyed(scaled_key, log_scaled_key, n,
                           datasources, dealer_ci, server_list,
                           server_names, y_server, nl,
                           transport_pks, session_id,
                           .dsAgg, .sendBlob,
                           transport_ready = transport_ready)

  # Step 3: add the public log-correction log(1/scale) so that the
  # final share encodes log(x) = log(scale*x) + log(1/scale). One
  # AffineCombine round (party-0 absorbs the constant; party-1's
  # contribution is identity).
  .ring127_local_fanout(
    datasources, server_list, server_names, y_server, .dsAgg,
    function(server, is_coord) call(name = "k2Ring127AffineCombineDS",
      a_key = log_scaled_key, b_key = NULL, sign_a = 1, sign_b = 0,
      public_const_fp = log_scale_correction_fp_b64,
      is_party0 = is_coord, output_key = out_key,
      n = as.numeric(n_int), session_id = session_id))
  invisible(NULL)
}


# Wide-Chebyshev seed + Newton-Raphson refinement log on Ring127 shares.
# Operates on inputs in [0.1, 1000] without per-element argument
# reduction, achieving ULP precision via quadratic NR convergence on
# the wide-Chebyshev seed.
#
# Architecture (Pugh 2004 PhD Sec.3 NR-on-Chebyshev for elementary
# functions; Goldschmidt 1964 NR division pattern; mirror of
# .ring127_recip_round_keyed which uses identical wide-Chebyshev +
# NR-refine pattern for 1/x):
#
#   y_0 = wide_chebyshev(x)           [60-deg Clenshaw on [0.1, 1000];
#                                       initial rel approx 30%]
#   y_{n+1} = y_n + x * exp(-y_n) - 1  [NR on f(y) = exp(y) - x;
#                                       quadratic convergence epsilon_{n+1}
#                                       approx -epsilon_n^2/2; Trefethen & Bau Sec.16]
#
# 5 NR iterations drive epsilon from ~0.30 -> 0.045 -> 0.001 -> 5e-7 -> 1.25e-13
# -> 7.8e-27, well below ULP 2^-50 approx 8.9e-16. Each NR iter costs:
#   - 1 .ring127_exp_round_keyed_extended (approx30 vecmul Clenshaw + 1
#     squaring vecmul approx 31 Beaver rounds)
#   - 1 .ring127_vecmul (1 Beaver round) for x * exp(-y_n)
#   - 1 k2Ring127AffineCombineDS round (no MPC) for y_n + result - 1
#   - 1 local-affine for -y_n
# Per-call total: 60 (init Cheb) + 5 x 32 = 220 Beaver rounds. About
# 5x the [1, 10] core single-call cost; trades runtime for the ULP
# precision floor across the full [0.1, 1000] operating range.
#
# Disclosure footing: identical to existing exp/recip/log primitives --
# all coefficients + affine constants are public deterministic values
# from the Go init(); share-side evaluation runs over the K=2 OT-Beaver
# dishonest-majority threat model (Demmler-Schneider-Zohner ABY 2015
# Sec.III.B) preserving D-INV-1..5.
#
# Caller responsibility: ensure the share-encoded input falls within
# [0.1, 1000]. For NB full-regression theta MLE with Poisson-warm
# eta in [-5, 5] -> mu in [0.0067, 148] and theta in [0.5, 5], (mu + theta) in
# [0.51, 153] is well within domain. Tighter eta clamps are encouraged
# at the orchestrator level for safety.
.ring127_log_round_keyed_nr <- function(in_key, out_key, n,
                                         datasources, dealer_ci, server_list,
                                         server_names, y_server, nl,
                                         transport_pks, session_id,
                                         .dsAgg, .sendBlob,
                                         nr_iters = 5L,
                                         transport_ready = FALSE) {
  n_int <- as.integer(n)
  nr_iters <- as.integer(nr_iters)
  transport_ready <- .ring127_transport_once(
    transport_ready, datasources, server_list, server_names,
    session_id, .dsAgg)

  # Step 0: fetch wide-Chebyshev coefficients (cache once per session).
  if (is.null(.ring127_log_wide_coef_cache$coef_res)) {
    .ring127_log_wide_coef_cache$coef_res <- dsVert:::.callMpcTool(
      "k2-log-shift-coeffs-wide", list(frac_bits = 50))
  }
  rc <- .ring127_log_wide_coef_cache$coef_res
  rc_one_over_half_range <- .to_b64url(rc$one_over_half_range)
  rc_neg_mid_over_half_range <- .to_b64url(rc$neg_mid_over_half_range)
  degree <- as.integer(rc$degree)
  all_coeffs_raw <- jsonlite::base64_dec(rc$coeffs)
  c_b64 <- vapply(seq_len(degree + 1L), function(idx) {
    s <- (idx - 1L) * 16L + 1L; e <- s + 15L
    .to_b64url(jsonlite::base64_enc(all_coeffs_raw[s:e]))
  }, character(1))

  tag <- .ring127_invocation_tag()
  t_pre  <- paste0("__r127_lpnr_tpre_", tag)
  t_key  <- paste0("__r127_lpnr_t_",    tag)
  twoT   <- paste0("__r127_lpnr_twoT_", tag)
  bB     <- paste0("__r127_lpnr_bB_",   tag)
  bA     <- paste0("__r127_lpnr_bA_",   tag)
  y_seed <- paste0("__r127_lpnr_yseed_", tag)

  # Step 1: wide-Cheb seed via "scale + affine + Clenshaw" pipeline.
  .ring127_exact_public_scale(
    in_key, rc_one_over_half_range, t_pre, n_int,
    datasources, dealer_ci, server_list, server_names, y_server, nl,
    transport_pks, session_id, .dsAgg, .sendBlob,
    transport_ready = transport_ready, destination_fresh = TRUE)
  .ring127_local_fanout(
    datasources, server_list, server_names, y_server, .dsAgg,
    function(server, is_coord) call(name = "k2Ring127AffineCombineDS",
      a_key = t_pre, b_key = NULL, sign_a = 1, sign_b = 0,
      public_const_fp = rc_neg_mid_over_half_range,
      is_party0 = is_coord, output_key = t_key,
      n = as.numeric(n_int), session_id = session_id))
  .ring127_local_fanout(
    datasources, server_list, server_names, y_server, .dsAgg,
    function(server, is_coord) call(name = "k2Ring127AffineCombineDS",
      a_key = t_key, b_key = t_key, sign_a = 1, sign_b = 1,
      public_const_fp = NULL, is_party0 = is_coord,
      output_key = twoT, n = as.numeric(n_int), session_id = session_id))
  .ring127_local_fanout(
    datasources, server_list, server_names, y_server, .dsAgg,
    function(server, is_coord) call(name = "k2Ring127AffineCombineDS",
      a_key = NULL, b_key = NULL, sign_a = 0, sign_b = 0,
      public_const_fp = c_b64[degree + 1L], is_party0 = is_coord,
      output_key = bB, n = as.numeric(n_int), session_id = session_id))
  .ring127_local_fanout(
    datasources, server_list, server_names, y_server, .dsAgg,
    function(server, is_coord) call(name = "k2Ring127AffineCombineDS",
      a_key = NULL, b_key = NULL, sign_a = 0, sign_b = 0,
      public_const_fp = NULL, is_party0 = is_coord,
      output_key = bA, n = as.numeric(n_int), session_id = session_id))

  slot_B <- bB; slot_A <- bA
  product_stage <- 0L
  for (k in seq.int(degree - 1L, 1L)) {
    product_stage <- product_stage + 1L
    product_key <- sprintf("__r127_lpnrmul%02d_%s", product_stage, tag)
    .ring127_vecmul(twoT, slot_B, product_key, n_int,
                    datasources, dealer_ci, server_list, server_names,
                    y_server, nl, transport_pks, session_id,
                    .dsAgg, .sendBlob, transport_ready = transport_ready,
                    destination_fresh = TRUE)
    .ring127_local_fanout(
      datasources, server_list, server_names, y_server, .dsAgg,
      function(server, is_coord) call(name = "k2Ring127AffineCombineDS",
        a_key = product_key, b_key = slot_A, sign_a = 1, sign_b = -1,
        public_const_fp = c_b64[k + 1L], is_party0 = is_coord,
        output_key = slot_A, n = as.numeric(n_int), session_id = session_id))
    swap <- slot_A; slot_A <- slot_B; slot_B <- swap
  }
  product_stage <- product_stage + 1L
  product_key <- sprintf("__r127_lpnrmul%02d_%s", product_stage, tag)
  .ring127_vecmul(t_key, slot_B, product_key, n_int,
                  datasources, dealer_ci, server_list, server_names,
                  y_server, nl, transport_pks, session_id,
                  .dsAgg, .sendBlob, transport_ready = transport_ready,
                  destination_fresh = TRUE)
  .ring127_local_fanout(
    datasources, server_list, server_names, y_server, .dsAgg,
    function(server, is_coord) call(name = "k2Ring127AffineCombineDS",
      a_key = product_key, b_key = slot_A, sign_a = 1, sign_b = -1,
      public_const_fp = c_b64[1L], is_party0 = is_coord,
      output_key = y_seed, n = as.numeric(n_int), session_id = session_id))

  # Step 2: 5 NR refinement iterations on shares.
  # y_{n+1} = y_n + x * exp(-y_n) - 1
  neg_one_fp_b64 <- .to_b64url(dsVert:::.callMpcTool("k2-float-to-fp", list(
    values = array(-1.0, dim = 1L), frac_bits = 50,
    ring = "ring127"))$fp_data)

  cur <- y_seed
  for (iter in seq_len(nr_iters)) {
    is_last <- (iter == nr_iters)
    negY <- sprintf("__r127_lpnr_negY%02d_%s", iter, tag)
    expN <- sprintf("__r127_lpnr_expNegY%02d_%s", iter, tag)
    xExp <- sprintf("__r127_lpnr_xExp%02d_%s", iter, tag)
    next_slot <- if (is_last) out_key else
      sprintf("__r127_lpnr_y%02d_%s", iter, tag)
    # Step 2a: -y_n share via local affine (sign_a = -1).
    .ring127_local_fanout(
      datasources, server_list, server_names, y_server, .dsAgg,
      function(server, is_coord) call(name = "k2Ring127AffineCombineDS",
        a_key = cur, b_key = NULL, sign_a = -1, sign_b = 0,
        public_const_fp = NULL, is_party0 = is_coord,
        output_key = negY, n = as.numeric(n_int), session_id = session_id))
    # Step 2b: exp(-y_n) share via existing primitive.
    .ring127_exp_round_keyed_extended(negY, expN, n_int,
                                       datasources, dealer_ci,
                                       server_list, server_names,
                                       y_server, nl, transport_pks,
                                       session_id,
                                       .dsAgg, .sendBlob,
                                       transport_ready = transport_ready)
    # Step 2c: x * exp(-y_n) share via Beaver vecmul.
    .ring127_vecmul(in_key, expN, xExp, n_int,
                    datasources, dealer_ci, server_list, server_names,
                    y_server, nl, transport_pks, session_id,
                    .dsAgg, .sendBlob, transport_ready = transport_ready,
                    destination_fresh = TRUE)
    # Step 2d: y_{n+1} = y_n + xExp + (-1) via affine-combine.
    # Note: a_key=cur, b_key=xExp, sign_a=1, sign_b=1, public_const=-1
    # at party-0 only. Uses canonical 3-term affine combine.
    .ring127_local_fanout(
      datasources, server_list, server_names, y_server, .dsAgg,
      function(server, is_coord) call(name = "k2Ring127AffineCombineDS",
        a_key = cur, b_key = xExp, sign_a = 1, sign_b = 1,
        public_const_fp = if (is_coord) neg_one_fp_b64 else NULL,
        is_party0 = is_coord,
        output_key = next_slot, n = as.numeric(n_int), session_id = session_id))
    cur <- next_slot
  }

  invisible(NULL)
}
