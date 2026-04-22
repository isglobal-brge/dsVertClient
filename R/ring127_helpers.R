#' Internal Ring127 Chebyshev-Horner + Beaver-vecmul helpers reusable
#' across joint-Newton orchestrations (Multinom, Ordinal). Mirror of
#' \code{.exp127_round} / \code{.run_beaver_vecmul_ring127} from
#' \code{ds.vertCox.R} but parameterised on \code{input_key} /
#' \code{output_key} so callers can maintain K−1 parallel share slots.
#'
#' @keywords internal
#' @noRd

.ring127_exp_coef_cache <- new.env(parent = emptyenv())
.ring127_recip_coef_cache <- new.env(parent = emptyenv())

# Ring127 Beaver vecmul on shares under arbitrary session keys.
# Mirrors dsVertClient:::.run_beaver_vecmul_ring127 (which is a closure
# inside ds.vertCox.R and not callable externally).
.ring127_vecmul <- function(x_key, y_key, output_key, n,
                            datasources, dealer_ci, server_list,
                            server_names, y_server, nl, transport_pks,
                            session_id, .dsAgg, .sendBlob) {
  n_int <- as.integer(n)
  tri <- .dsAgg(datasources[dealer_ci],
    call("k2BeaverVecmulGenTriplesDS",
         dcf0_pk = transport_pks[[y_server]],
         dcf1_pk = transport_pks[[nl]],
         n = n_int, session_id = session_id,
         frac_bits = 50L, ring = 127L))
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
      n = n_int, session_id = session_id,
      frac_bits = 50L, ring = 127L))
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
    .dsAgg(datasources[ci], call("k2BeaverVecmulR2DS",
      is_party0 = is_coord, output_key = output_key,
      n = n_int, session_id = session_id,
      frac_bits = 50L, ring = 127L))
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
      "k2-exp127-get-coeffs", list(frac_bits = 50L))
  }
  coef_res <- .ring127_exp_coef_cache$coef_res
  degree <- as.integer(coef_res$degree)
  all_coeffs_raw <- jsonlite::base64_dec(coef_res$coeffs)
  c_b64 <- vapply(seq_len(degree + 1L), function(idx) {
    s <- (idx - 1L) * 16L + 1L; e <- s + 15L
    jsonlite::base64_enc(all_coeffs_raw[s:e])
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
    .dsAgg(datasources[ci], call("k2Ring127LocalScaleDS",
      in_key = in_key, scalar_fp = coef_res$one_over_a,
      output_key = tmp_y, n = n_int, session_id = session_id))
    .dsAgg(datasources[ci], call("k2Ring127AffineCombineDS",
      a_key = tmp_y, b_key = tmp_y, sign_a = 1L, sign_b = 1L,
      public_const_fp = NULL, is_party0 = is_coord,
      output_key = tmp_twoY, n = n_int, session_id = session_id))
    .dsAgg(datasources[ci], call("k2Ring127AffineCombineDS",
      a_key = NULL, b_key = NULL, sign_a = 0L, sign_b = 0L,
      public_const_fp = c_b64[degree + 1L], is_party0 = is_coord,
      output_key = tmp_bB, n = n_int, session_id = session_id))
    .dsAgg(datasources[ci], call("k2Ring127AffineCombineDS",
      a_key = NULL, b_key = NULL, sign_a = 0L, sign_b = 0L,
      public_const_fp = NULL, is_party0 = is_coord,
      output_key = tmp_bA, n = n_int, session_id = session_id))
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
      .dsAgg(datasources[ci], call("k2Ring127AffineCombineDS",
        a_key = tmp_res, b_key = slot_A, sign_a = 1L, sign_b = -1L,
        public_const_fp = c_b64[k + 1L], is_party0 = is_coord,
        output_key = slot_A, n = n_int, session_id = session_id))
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
    .dsAgg(datasources[ci], call("k2Ring127AffineCombineDS",
      a_key = tmp_res, b_key = slot_A, sign_a = 1L, sign_b = -1L,
      public_const_fp = c_b64[1L], is_party0 = is_coord,
      output_key = out_key, n = n_int, session_id = session_id))
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
      "k2-recip127-get-coeffs", list(frac_bits = 50L))
  }
  rc <- .ring127_recip_coef_cache$coef_res
  degree <- as.integer(rc$degree)
  nr_steps <- as.integer(rc$nr_steps)
  all_coeffs_raw <- jsonlite::base64_dec(rc$coeffs)
  c_b64 <- vapply(seq_len(degree + 1L), function(idx) {
    s <- (idx - 1L) * 16L + 1L; e <- s + 15L
    jsonlite::base64_enc(all_coeffs_raw[s:e])
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
    .dsAgg(datasources[ci], call("k2Ring127LocalScaleDS",
      in_key = in_key, scalar_fp = rc$one_over_half_range,
      output_key = t_pre, n = n_int, session_id = session_id))
    .dsAgg(datasources[ci], call("k2Ring127AffineCombineDS",
      a_key = t_pre, b_key = NULL, sign_a = 1L, sign_b = 0L,
      public_const_fp = rc$neg_mid_over_half_range,
      is_party0 = is_coord, output_key = t_key,
      n = n_int, session_id = session_id))
    .dsAgg(datasources[ci], call("k2Ring127AffineCombineDS",
      a_key = t_key, b_key = t_key, sign_a = 1L, sign_b = 1L,
      public_const_fp = NULL, is_party0 = is_coord,
      output_key = twoT, n = n_int, session_id = session_id))
    .dsAgg(datasources[ci], call("k2Ring127AffineCombineDS",
      a_key = NULL, b_key = NULL, sign_a = 0L, sign_b = 0L,
      public_const_fp = c_b64[degree + 1L], is_party0 = is_coord,
      output_key = bB, n = n_int, session_id = session_id))
    .dsAgg(datasources[ci], call("k2Ring127AffineCombineDS",
      a_key = NULL, b_key = NULL, sign_a = 0L, sign_b = 0L,
      public_const_fp = NULL, is_party0 = is_coord,
      output_key = bA, n = n_int, session_id = session_id))
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
      .dsAgg(datasources[ci], call("k2Ring127AffineCombineDS",
        a_key = tmp, b_key = slot_A, sign_a = 1L, sign_b = -1L,
        public_const_fp = c_b64[k + 1L], is_party0 = is_coord,
        output_key = slot_A, n = n_int, session_id = session_id))
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
    .dsAgg(datasources[ci], call("k2Ring127AffineCombineDS",
      a_key = tmp, b_key = slot_A, sign_a = 1L, sign_b = -1L,
      public_const_fp = c_b64[1L], is_party0 = is_coord,
      output_key = y_cur, n = n_int, session_id = session_id))
  }

  for (iter in seq_len(nr_steps)) {
    .ring127_vecmul(in_key, y_cur, xy, n_int,
                    datasources, dealer_ci, server_list, server_names,
                    y_server, nl, transport_pks, session_id,
                    .dsAgg, .sendBlob)
    for (server in server_list) {
      ci <- which(server_names == server)
      is_coord <- (server == y_server)
      .dsAgg(datasources[ci], call("k2Ring127AffineCombineDS",
        a_key = NULL, b_key = xy, sign_a = 0L, sign_b = -1L,
        public_const_fp = rc$two_fp, is_party0 = is_coord,
        output_key = tmXY, n = n_int, session_id = session_id))
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
