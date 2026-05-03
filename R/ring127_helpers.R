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

# Ring127 Beaver vecmul on shares under arbitrary session keys.
# Mirrors dsVertClient:::.run_beaver_vecmul_ring127 (which is a closure
# inside ds.vertCox.R and not callable externally).
.ring127_vecmul <- function(x_key, y_key, output_key, n,
                            datasources, dealer_ci, server_list,
                            server_names, y_server, nl, transport_pks,
                            session_id, .dsAgg, .sendBlob) {
  n_int <- as.integer(n)
  tri <- .dsAgg(datasources[dealer_ci],
    call(name = "k2BeaverVecmulGenTriplesDS",
         dcf0_pk = transport_pks[[y_server]],
         dcf1_pk = transport_pks[[nl]],
         n = as.numeric(n_int), session_id = session_id,
         frac_bits = 50, ring = 127))
  if (is.list(tri) && length(tri) == 1L) tri <- tri[[1L]]
  .sendBlob(tri$triple_blob_0, "k2_beaver_vecmul_triple",
            which(server_names == y_server))
  .sendBlob(tri$triple_blob_1, "k2_beaver_vecmul_triple",
            which(server_names == nl))
  all_ci <- vapply(server_list, function(s) which(server_names == s),
                    integer(1L))
  .dsAgg(datasources[all_ci],
         call(name = "k2BeaverVecmulConsumeTripleDS", session_id = session_id))
  r1_b <- list()
  for (server in server_list) {
    ci <- which(server_names == server)
    peer <- setdiff(server_list, server)
    r <- .dsAgg(datasources[ci], call(name = "k2BeaverVecmulR1DS",
      peer_pk = transport_pks[[peer]],
      x_key = x_key, y_key = y_key,
      n = as.numeric(n_int), session_id = session_id,
      frac_bits = 50, ring = 127))
    if (is.list(r) && length(r) == 1L) r <- r[[1L]]
    r1_b[[server]] <- r
  }
  .sendBlob(r1_b[[y_server]]$peer_blob, "k2_beaver_vecmul_peer_masked",
            which(server_names == nl))
  .sendBlob(r1_b[[nl]]$peer_blob, "k2_beaver_vecmul_peer_masked",
            which(server_names == y_server))
  for (server in server_list) {
    ci <- which(server_names == server)
    is_coord <- (server == y_server)
    .dsAgg(datasources[ci], call(name = "k2BeaverVecmulR2DS",
      is_party0 = is_coord, x_key = x_key, y_key = y_key,
      output_key = output_key,
      n = as.numeric(n_int), session_id = session_id,
      frac_bits = 50, ring = 127))
  }
  invisible(NULL)
}

# Clenshaw exp Horner on Ring127 shares with keyed IO. Takes the
# Chebyshev coefficients from the Go binary once per session and runs
# `degree + 1` Beaver vecmul rounds against the input share.
.ring127_exp_round_keyed <- function(in_key, out_key, n,
                                     datasources, dealer_ci, server_list,
                                     server_names, y_server, nl,
                                     transport_pks, session_id,
                                     .dsAgg, .sendBlob) {
  n_int <- as.integer(n)
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

  tag <- make.names(out_key, unique = TRUE)
  tmp_y    <- paste0("__r127_ey_",    tag)
  tmp_twoY <- paste0("__r127_etwoY_", tag)
  tmp_bB   <- paste0("__r127_ebB_",   tag)
  tmp_bA   <- paste0("__r127_ebA_",   tag)
  tmp_res  <- paste0("__r127_etmp_",  tag)

  for (server in server_list) {
    ci <- which(server_names == server)
    is_coord <- (server == y_server)
    .dsAgg(datasources[ci], call(name = "k2Ring127LocalScaleDS",
      in_key = in_key, scalar_fp = coef_res_one_over_a,
      output_key = tmp_y, n = as.numeric(n_int), session_id = session_id,
      is_party0 = is_coord))
    .dsAgg(datasources[ci], call(name = "k2Ring127AffineCombineDS",
      a_key = tmp_y, b_key = tmp_y, sign_a = 1, sign_b = 1,
      public_const_fp = NULL, is_party0 = is_coord,
      output_key = tmp_twoY, n = as.numeric(n_int), session_id = session_id))
    .dsAgg(datasources[ci], call(name = "k2Ring127AffineCombineDS",
      a_key = NULL, b_key = NULL, sign_a = 0, sign_b = 0,
      public_const_fp = c_b64[degree + 1L], is_party0 = is_coord,
      output_key = tmp_bB, n = as.numeric(n_int), session_id = session_id))
    .dsAgg(datasources[ci], call(name = "k2Ring127AffineCombineDS",
      a_key = NULL, b_key = NULL, sign_a = 0, sign_b = 0,
      public_const_fp = NULL, is_party0 = is_coord,
      output_key = tmp_bA, n = as.numeric(n_int), session_id = session_id))
  }

  slot_B <- tmp_bB; slot_A <- tmp_bA
  for (k in seq.int(degree - 1L, 1L)) {
    .ring127_vecmul(tmp_twoY, slot_B, tmp_res, n_int,
                    datasources, dealer_ci, server_list, server_names,
                    y_server, nl, transport_pks, session_id,
                    .dsAgg, .sendBlob)
    for (server in server_list) {
      ci <- which(server_names == server)
      is_coord <- (server == y_server)
      .dsAgg(datasources[ci], call(name = "k2Ring127AffineCombineDS",
        a_key = tmp_res, b_key = slot_A, sign_a = 1, sign_b = -1,
        public_const_fp = c_b64[k + 1L], is_party0 = is_coord,
        output_key = slot_A, n = as.numeric(n_int), session_id = session_id))
    }
    swap <- slot_A; slot_A <- slot_B; slot_B <- swap
  }
  .ring127_vecmul(tmp_y, slot_B, tmp_res, n_int,
                  datasources, dealer_ci, server_list, server_names,
                  y_server, nl, transport_pks, session_id,
                  .dsAgg, .sendBlob)
  for (server in server_list) {
    ci <- which(server_names == server)
    is_coord <- (server == y_server)
    .dsAgg(datasources[ci], call(name = "k2Ring127AffineCombineDS",
      a_key = tmp_res, b_key = slot_A, sign_a = 1, sign_b = -1,
      public_const_fp = c_b64[1L], is_party0 = is_coord,
      output_key = out_key, n = as.numeric(n_int), session_id = session_id))
  }
  invisible(NULL)
}

# Ring127 spline-less 1/x via Chebyshev-Horner initial guess + 6 NR
# iterations. Mirror of ds.vertCox.R's .recip127_round but keyed I/O.
.ring127_recip_round_keyed <- function(in_key, out_key, n,
                                       datasources, dealer_ci, server_list,
                                       server_names, y_server, nl,
                                       transport_pks, session_id,
                                       .dsAgg, .sendBlob) {
  n_int <- as.integer(n)
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

  tag <- make.names(out_key, unique = TRUE)
  t_pre <- paste0("__r127_rp_tpre_", tag)
  t_key <- paste0("__r127_rp_t_",    tag)
  twoT  <- paste0("__r127_rp_twoT_", tag)
  bB    <- paste0("__r127_rp_bB_",   tag)
  bA    <- paste0("__r127_rp_bA_",   tag)
  tmp   <- paste0("__r127_rp_tmp_",  tag)
  y_cur <- paste0("__r127_rp_y_",    tag)
  y_alt <- paste0("__r127_rp_yalt_", tag)
  xy    <- paste0("__r127_rp_xy_",   tag)
  tmXY  <- paste0("__r127_rp_twoMinusXY_", tag)

  for (server in server_list) {
    ci <- which(server_names == server)
    is_coord <- (server == y_server)
    .dsAgg(datasources[ci], call(name = "k2Ring127LocalScaleDS",
      in_key = in_key, scalar_fp = rc_one_over_half_range,
      output_key = t_pre, n = as.numeric(n_int), session_id = session_id,
      is_party0 = is_coord))
    .dsAgg(datasources[ci], call(name = "k2Ring127AffineCombineDS",
      a_key = t_pre, b_key = NULL, sign_a = 1, sign_b = 0,
      public_const_fp = rc_neg_mid_over_half_range,
      is_party0 = is_coord, output_key = t_key,
      n = as.numeric(n_int), session_id = session_id))
    .dsAgg(datasources[ci], call(name = "k2Ring127AffineCombineDS",
      a_key = t_key, b_key = t_key, sign_a = 1, sign_b = 1,
      public_const_fp = NULL, is_party0 = is_coord,
      output_key = twoT, n = as.numeric(n_int), session_id = session_id))
    .dsAgg(datasources[ci], call(name = "k2Ring127AffineCombineDS",
      a_key = NULL, b_key = NULL, sign_a = 0, sign_b = 0,
      public_const_fp = c_b64[degree + 1L], is_party0 = is_coord,
      output_key = bB, n = as.numeric(n_int), session_id = session_id))
    .dsAgg(datasources[ci], call(name = "k2Ring127AffineCombineDS",
      a_key = NULL, b_key = NULL, sign_a = 0, sign_b = 0,
      public_const_fp = NULL, is_party0 = is_coord,
      output_key = bA, n = as.numeric(n_int), session_id = session_id))
  }

  slot_B <- bB; slot_A <- bA
  for (k in seq.int(degree - 1L, 1L)) {
    .ring127_vecmul(twoT, slot_B, tmp, n_int,
                    datasources, dealer_ci, server_list, server_names,
                    y_server, nl, transport_pks, session_id,
                    .dsAgg, .sendBlob)
    for (server in server_list) {
      ci <- which(server_names == server)
      is_coord <- (server == y_server)
      .dsAgg(datasources[ci], call(name = "k2Ring127AffineCombineDS",
        a_key = tmp, b_key = slot_A, sign_a = 1, sign_b = -1,
        public_const_fp = c_b64[k + 1L], is_party0 = is_coord,
        output_key = slot_A, n = as.numeric(n_int), session_id = session_id))
    }
    swap <- slot_A; slot_A <- slot_B; slot_B <- swap
  }
  .ring127_vecmul(t_key, slot_B, tmp, n_int,
                  datasources, dealer_ci, server_list, server_names,
                  y_server, nl, transport_pks, session_id,
                  .dsAgg, .sendBlob)
  for (server in server_list) {
    ci <- which(server_names == server)
    is_coord <- (server == y_server)
    .dsAgg(datasources[ci], call(name = "k2Ring127AffineCombineDS",
      a_key = tmp, b_key = slot_A, sign_a = 1, sign_b = -1,
      public_const_fp = c_b64[1L], is_party0 = is_coord,
      output_key = y_cur, n = as.numeric(n_int), session_id = session_id))
  }

  for (iter in seq_len(nr_steps)) {
    .ring127_vecmul(in_key, y_cur, xy, n_int,
                    datasources, dealer_ci, server_list, server_names,
                    y_server, nl, transport_pks, session_id,
                    .dsAgg, .sendBlob)
    for (server in server_list) {
      ci <- which(server_names == server)
      is_coord <- (server == y_server)
      .dsAgg(datasources[ci], call(name = "k2Ring127AffineCombineDS",
        a_key = NULL, b_key = xy, sign_a = 0, sign_b = -1,
        public_const_fp = rc_two_fp, is_party0 = is_coord,
        output_key = tmXY, n = as.numeric(n_int), session_id = session_id))
    }
    final_slot <- if (iter == nr_steps) out_key else y_alt
    .ring127_vecmul(y_cur, tmXY, final_slot, n_int,
                    datasources, dealer_ci, server_list, server_names,
                    y_server, nl, transport_pks, session_id,
                    .dsAgg, .sendBlob)
    if (iter < nr_steps) {
      swap <- y_cur; y_cur <- y_alt; y_alt <- swap
    }
  }
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
# Cost: 1 extra k2Ring127LocalScaleDS + 1 extra k2BeaverVecmul (approx 3%
# overhead vs the degree-30 Clenshaw core, which uses 30 vecmuls).
.ring127_exp_round_keyed_extended <- function(in_key, out_key, n,
                                              datasources, dealer_ci, server_list,
                                              server_names, y_server, nl,
                                              transport_pks, session_id,
                                              .dsAgg, .sendBlob) {
  n_int <- as.integer(n)
  half_fp_b64 <- .to_b64url(.ring127_get_half_fp())
  half_key <- paste0(in_key, "_ext_half")
  exp_half_key <- paste0(out_key, "_ext_exphalf")

  # Step 1: scale input by 0.5 locally on each server (linear, no MPC round).
  for (server in server_list) {
    ci <- which(server_names == server)
    is_coord <- (server == y_server)
    .dsAgg(datasources[ci], call(name = "k2Ring127LocalScaleDS",
      in_key = in_key, scalar_fp = half_fp_b64,
      output_key = half_key, n = as.numeric(n_int), session_id = session_id,
      is_party0 = is_coord))
  }
  # Step 2: evaluate exp on half-input (interior Chebyshev [-5,5]).
  .ring127_exp_round_keyed(half_key, exp_half_key, n,
                           datasources, dealer_ci, server_list,
                           server_names, y_server, nl,
                           transport_pks, session_id,
                           .dsAgg, .sendBlob)
  # Step 3: square via Beaver vecmul -> exp(x) = exp(x/2)^2.
  .ring127_vecmul(exp_half_key, exp_half_key, out_key, n_int,
                  datasources, dealer_ci, server_list, server_names,
                  y_server, nl, transport_pks, session_id,
                  .dsAgg, .sendBlob)
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
                                     .dsAgg, .sendBlob) {
  n_int <- as.integer(n)
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

  tag <- make.names(out_key, unique = TRUE)
  t_pre <- paste0("__r127_lp_tpre_", tag)
  t_key <- paste0("__r127_lp_t_",    tag)
  twoT  <- paste0("__r127_lp_twoT_", tag)
  bB    <- paste0("__r127_lp_bB_",   tag)
  bA    <- paste0("__r127_lp_bA_",   tag)
  tmp   <- paste0("__r127_lp_tmp_",  tag)

  for (server in server_list) {
    ci <- which(server_names == server)
    is_coord <- (server == y_server)
    .dsAgg(datasources[ci], call(name = "k2Ring127LocalScaleDS",
      in_key = in_key, scalar_fp = rc_one_over_half_range,
      output_key = t_pre, n = as.numeric(n_int), session_id = session_id,
      is_party0 = is_coord))
    .dsAgg(datasources[ci], call(name = "k2Ring127AffineCombineDS",
      a_key = t_pre, b_key = NULL, sign_a = 1, sign_b = 0,
      public_const_fp = rc_neg_mid_over_half_range,
      is_party0 = is_coord, output_key = t_key,
      n = as.numeric(n_int), session_id = session_id))
    .dsAgg(datasources[ci], call(name = "k2Ring127AffineCombineDS",
      a_key = t_key, b_key = t_key, sign_a = 1, sign_b = 1,
      public_const_fp = NULL, is_party0 = is_coord,
      output_key = twoT, n = as.numeric(n_int), session_id = session_id))
    .dsAgg(datasources[ci], call(name = "k2Ring127AffineCombineDS",
      a_key = NULL, b_key = NULL, sign_a = 0, sign_b = 0,
      public_const_fp = c_b64[degree + 1L], is_party0 = is_coord,
      output_key = bB, n = as.numeric(n_int), session_id = session_id))
    .dsAgg(datasources[ci], call(name = "k2Ring127AffineCombineDS",
      a_key = NULL, b_key = NULL, sign_a = 0, sign_b = 0,
      public_const_fp = NULL, is_party0 = is_coord,
      output_key = bA, n = as.numeric(n_int), session_id = session_id))
  }

  slot_B <- bB; slot_A <- bA
  for (k in seq.int(degree - 1L, 1L)) {
    .ring127_vecmul(twoT, slot_B, tmp, n_int,
                    datasources, dealer_ci, server_list, server_names,
                    y_server, nl, transport_pks, session_id,
                    .dsAgg, .sendBlob)
    for (server in server_list) {
      ci <- which(server_names == server)
      is_coord <- (server == y_server)
      .dsAgg(datasources[ci], call(name = "k2Ring127AffineCombineDS",
        a_key = tmp, b_key = slot_A, sign_a = 1, sign_b = -1,
        public_const_fp = c_b64[k + 1L], is_party0 = is_coord,
        output_key = slot_A, n = as.numeric(n_int), session_id = session_id))
    }
    swap <- slot_A; slot_A <- slot_B; slot_B <- swap
  }
  .ring127_vecmul(t_key, slot_B, tmp, n_int,
                  datasources, dealer_ci, server_list, server_names,
                  y_server, nl, transport_pks, session_id,
                  .dsAgg, .sendBlob)
  for (server in server_list) {
    ci <- which(server_names == server)
    is_coord <- (server == y_server)
    .dsAgg(datasources[ci], call(name = "k2Ring127AffineCombineDS",
      a_key = tmp, b_key = slot_A, sign_a = 1, sign_b = -1,
      public_const_fp = c_b64[1L], is_party0 = is_coord,
      output_key = out_key, n = as.numeric(n_int), session_id = session_id))
  }
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
                                              .dsAgg, .sendBlob) {
  n_int <- as.integer(n)
  scaled_key <- paste0(in_key, "_ext_scaled")
  log_scaled_key <- paste0(out_key, "_ext_logscaled")

  # Step 1: locally rescale input by `scale` (no MPC round). After this
  # step, scaled_share encodes scale * x, which the caller has guaranteed
  # falls in [1, 10].
  for (server in server_list) {
    ci <- which(server_names == server)
    is_coord <- (server == y_server)
    .dsAgg(datasources[ci], call(name = "k2Ring127LocalScaleDS",
      in_key = in_key, scalar_fp = scale_fp_b64,
      output_key = scaled_key, n = as.numeric(n_int), session_id = session_id,
      is_party0 = is_coord))
  }

  # Step 2: evaluate log on the rescaled share (Chebyshev core [1, 10]).
  .ring127_log_round_keyed(scaled_key, log_scaled_key, n,
                           datasources, dealer_ci, server_list,
                           server_names, y_server, nl,
                           transport_pks, session_id,
                           .dsAgg, .sendBlob)

  # Step 3: add the public log-correction log(1/scale) so that the
  # final share encodes log(x) = log(scale*x) + log(1/scale). One
  # AffineCombine round (party-0 absorbs the constant; party-1's
  # contribution is identity).
  for (server in server_list) {
    ci <- which(server_names == server)
    is_coord <- (server == y_server)
    .dsAgg(datasources[ci], call(name = "k2Ring127AffineCombineDS",
      a_key = log_scaled_key, b_key = NULL, sign_a = 1, sign_b = 0,
      public_const_fp = log_scale_correction_fp_b64,
      is_party0 = is_coord, output_key = out_key,
      n = as.numeric(n_int), session_id = session_id))
  }
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
                                         nr_iters = 5L) {
  n_int <- as.integer(n)
  nr_iters <- as.integer(nr_iters)

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

  tag <- make.names(out_key, unique = TRUE)
  t_pre  <- paste0("__r127_lpnr_tpre_", tag)
  t_key  <- paste0("__r127_lpnr_t_",    tag)
  twoT   <- paste0("__r127_lpnr_twoT_", tag)
  bB     <- paste0("__r127_lpnr_bB_",   tag)
  bA     <- paste0("__r127_lpnr_bA_",   tag)
  tmp    <- paste0("__r127_lpnr_tmp_",  tag)
  y_seed <- paste0("__r127_lpnr_yseed_", tag)
  y_cur  <- paste0("__r127_lpnr_ycur_",  tag)
  y_alt  <- paste0("__r127_lpnr_yalt_",  tag)
  negY   <- paste0("__r127_lpnr_negY_",  tag)
  expN   <- paste0("__r127_lpnr_expNegY_", tag)
  xExp   <- paste0("__r127_lpnr_xExp_",  tag)

  # Step 1: wide-Cheb seed via "scale + affine + Clenshaw" pipeline.
  for (server in server_list) {
    ci <- which(server_names == server)
    is_coord <- (server == y_server)
    .dsAgg(datasources[ci], call(name = "k2Ring127LocalScaleDS",
      in_key = in_key, scalar_fp = rc_one_over_half_range,
      output_key = t_pre, n = as.numeric(n_int), session_id = session_id,
      is_party0 = is_coord))
    .dsAgg(datasources[ci], call(name = "k2Ring127AffineCombineDS",
      a_key = t_pre, b_key = NULL, sign_a = 1, sign_b = 0,
      public_const_fp = rc_neg_mid_over_half_range,
      is_party0 = is_coord, output_key = t_key,
      n = as.numeric(n_int), session_id = session_id))
    .dsAgg(datasources[ci], call(name = "k2Ring127AffineCombineDS",
      a_key = t_key, b_key = t_key, sign_a = 1, sign_b = 1,
      public_const_fp = NULL, is_party0 = is_coord,
      output_key = twoT, n = as.numeric(n_int), session_id = session_id))
    .dsAgg(datasources[ci], call(name = "k2Ring127AffineCombineDS",
      a_key = NULL, b_key = NULL, sign_a = 0, sign_b = 0,
      public_const_fp = c_b64[degree + 1L], is_party0 = is_coord,
      output_key = bB, n = as.numeric(n_int), session_id = session_id))
    .dsAgg(datasources[ci], call(name = "k2Ring127AffineCombineDS",
      a_key = NULL, b_key = NULL, sign_a = 0, sign_b = 0,
      public_const_fp = NULL, is_party0 = is_coord,
      output_key = bA, n = as.numeric(n_int), session_id = session_id))
  }

  slot_B <- bB; slot_A <- bA
  for (k in seq.int(degree - 1L, 1L)) {
    .ring127_vecmul(twoT, slot_B, tmp, n_int,
                    datasources, dealer_ci, server_list, server_names,
                    y_server, nl, transport_pks, session_id,
                    .dsAgg, .sendBlob)
    for (server in server_list) {
      ci <- which(server_names == server)
      is_coord <- (server == y_server)
      .dsAgg(datasources[ci], call(name = "k2Ring127AffineCombineDS",
        a_key = tmp, b_key = slot_A, sign_a = 1, sign_b = -1,
        public_const_fp = c_b64[k + 1L], is_party0 = is_coord,
        output_key = slot_A, n = as.numeric(n_int), session_id = session_id))
    }
    swap <- slot_A; slot_A <- slot_B; slot_B <- swap
  }
  .ring127_vecmul(t_key, slot_B, tmp, n_int,
                  datasources, dealer_ci, server_list, server_names,
                  y_server, nl, transport_pks, session_id,
                  .dsAgg, .sendBlob)
  for (server in server_list) {
    ci <- which(server_names == server)
    is_coord <- (server == y_server)
    .dsAgg(datasources[ci], call(name = "k2Ring127AffineCombineDS",
      a_key = tmp, b_key = slot_A, sign_a = 1, sign_b = -1,
      public_const_fp = c_b64[1L], is_party0 = is_coord,
      output_key = y_seed, n = as.numeric(n_int), session_id = session_id))
  }

  # Step 2: 5 NR refinement iterations on shares.
  # y_{n+1} = y_n + x * exp(-y_n) - 1
  one_fp_b64 <- .to_b64url(dsVert:::.callMpcTool("k2-float-to-fp", list(
    values = array(1.0, dim = 1L), frac_bits = 50,
    ring = "ring127"))$fp_data)
  neg_one_fp_b64 <- .to_b64url(dsVert:::.callMpcTool("k2-float-to-fp", list(
    values = array(-1.0, dim = 1L), frac_bits = 50,
    ring = "ring127"))$fp_data)

  cur <- y_seed
  for (iter in seq_len(nr_iters)) {
    is_last <- (iter == nr_iters)
    next_slot <- if (is_last) out_key else (if (cur == y_cur) y_alt else y_cur)
    # Step 2a: -y_n share via local affine (sign_a = -1).
    for (server in server_list) {
      ci <- which(server_names == server)
      is_coord <- (server == y_server)
      .dsAgg(datasources[ci], call(name = "k2Ring127AffineCombineDS",
        a_key = cur, b_key = NULL, sign_a = -1, sign_b = 0,
        public_const_fp = NULL, is_party0 = is_coord,
        output_key = negY, n = as.numeric(n_int), session_id = session_id))
    }
    # Step 2b: exp(-y_n) share via existing primitive.
    .ring127_exp_round_keyed_extended(negY, expN, n_int,
                                       datasources, dealer_ci,
                                       server_list, server_names,
                                       y_server, nl, transport_pks,
                                       session_id,
                                       .dsAgg, .sendBlob)
    # Step 2c: x * exp(-y_n) share via Beaver vecmul.
    .ring127_vecmul(in_key, expN, xExp, n_int,
                    datasources, dealer_ci, server_list, server_names,
                    y_server, nl, transport_pks, session_id,
                    .dsAgg, .sendBlob)
    # Step 2d: y_{n+1} = y_n + xExp + (-1) via affine-combine.
    # Note: a_key=cur, b_key=xExp, sign_a=1, sign_b=1, public_const=-1
    # at party-0 only. Uses canonical 3-term affine combine.
    for (server in server_list) {
      ci <- which(server_names == server)
      is_coord <- (server == y_server)
      .dsAgg(datasources[ci], call(name = "k2Ring127AffineCombineDS",
        a_key = cur, b_key = xExp, sign_a = 1, sign_b = 1,
        public_const_fp = if (is_coord) neg_one_fp_b64 else NULL,
        is_party0 = is_coord,
        output_key = next_slot, n = as.numeric(n_int), session_id = session_id))
    }
    cur <- next_slot
  }

  invisible(NULL)
}
