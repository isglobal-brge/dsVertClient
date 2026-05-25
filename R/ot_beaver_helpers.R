#' Client-side OT-Beaver preprocessing orchestration
#'
#' These helpers execute the two OT cross-term directions required to generate
#' Beaver triples without a trusted dealer. They install the resulting triple
#' shares into the same server-side session slots consumed by the existing
#' online Beaver R1/R2 methods.
#'
#' @keywords internal
#' @noRd

.beaver_preprocessing_mode <- function(kind, n, p, ring) {
  mode <- getOption("dsvert.beaver_preprocessing", "auto")
  mode <- match.arg(tolower(as.character(mode)[1L]),
                    c("auto", "dealer", "ot", "iknp", "direct_ot"))
  if (identical(mode, "ot")) mode <- "iknp"
  if (identical(mode, "auto")) {
    bit_ops <- as.numeric(n) * max(1, as.numeric(p)) * as.numeric(ring) * 2
    threshold <- getOption("dsvert.ot_direct_max_bit_ops", 512L)
    mode <- if (is.finite(bit_ops) && bit_ops <= threshold) {
      "direct_ot"
    } else {
      "iknp"
    }
  }
  mode
}

.dealer_prepare_vecmul <- function(datasources, party_conns, party_names,
                                   transport_pks, session_id, n, ring,
                                   .dsAgg, .sendBlob, dealer_conn = NULL) {
  if (is.null(dealer_conn)) dealer_conn <- party_conns[[2L]]
  dealer_party <- match(as.integer(dealer_conn), as.integer(party_conns),
                        nomatch = NA_integer_) - 1L
  if (is.na(dealer_party)) dealer_party <- NULL
  ring <- as.integer(ring)
  frac_bits <- if (ring == 127L) 50L else 20L
  tri <- .dsAgg(datasources[dealer_conn],
    call(name = "k2BeaverVecmulGenTriplesDS",
         dcf0_pk = transport_pks[[party_names[[1L]]]],
         dcf1_pk = transport_pks[[party_names[[2L]]]],
         n = as.numeric(as.integer(n)), session_id = session_id,
         frac_bits = frac_bits, ring = ring, dealer_party = dealer_party))
  if (is.list(tri) && length(tri) == 1L) tri <- tri[[1L]]
  if (!is.null(tri$triple_blob_0) && nzchar(tri$triple_blob_0)) {
    .sendBlob(tri$triple_blob_0, "k2_beaver_vecmul_triple",
              party_conns[[1L]])
  }
  if (!is.null(tri$triple_blob_1) && nzchar(tri$triple_blob_1)) {
    .sendBlob(tri$triple_blob_1, "k2_beaver_vecmul_triple",
              party_conns[[2L]])
  }
  invisible(NULL)
}

.dealer_prepare_grad <- function(datasources, party_conns, party_names,
                                 transport_pks, session_id, n, p, ring,
                                 .dsAgg, .sendBlob, dealer_conn = NULL,
                                 grad_triple_key = "k2_grad_triple_fp") {
  if (is.null(dealer_conn)) dealer_conn <- party_conns[[2L]]
  dealer_party <- match(as.integer(dealer_conn), as.integer(party_conns),
                        nomatch = NA_integer_) - 1L
  if (is.na(dealer_party)) dealer_party <- NULL
  grad_t <- .dsAgg(datasources[dealer_conn],
    call(name = "glmRing63GenGradTriplesDS",
         dcf0_pk = transport_pks[[party_names[[1L]]]],
         dcf1_pk = transport_pks[[party_names[[2L]]]],
         n = as.integer(n), p = as.integer(p),
         ring = as.integer(ring), session_id = session_id,
         dealer_party = dealer_party))
  if (is.list(grad_t) && length(grad_t) == 1L) grad_t <- grad_t[[1L]]
  if (!is.null(grad_t$grad_blob_0) && nzchar(grad_t$grad_blob_0)) {
    .sendBlob(grad_t$grad_blob_0, grad_triple_key, party_conns[[1L]])
  }
  if (!is.null(grad_t$grad_blob_1) && nzchar(grad_t$grad_blob_1)) {
    .sendBlob(grad_t$grad_blob_1, grad_triple_key, party_conns[[2L]])
  }
  invisible(NULL)
}

.dealer_prepare_spline <- function(datasources, party_conns, party_names,
                                   transport_pks, session_id, n, ring,
                                   .dsAgg, .sendBlob, dealer_conn = NULL) {
  if (is.null(dealer_conn)) dealer_conn <- party_conns[[2L]]
  dealer_party <- match(as.integer(dealer_conn), as.integer(party_conns),
                        nomatch = NA_integer_) - 1L
  if (is.na(dealer_party)) dealer_party <- NULL
  ring <- as.integer(ring)
  frac_bits <- if (ring == 127L) 50L else 20L
  spline_t <- .dsAgg(datasources[dealer_conn],
    call(name = "glmRing63GenSplineTriplesDS",
         dcf0_pk = transport_pks[[party_names[[1L]]]],
         dcf1_pk = transport_pks[[party_names[[2L]]]],
         n = as.integer(n), frac_bits = frac_bits,
         ring = ring, session_id = session_id, dealer_party = dealer_party))
  if (is.list(spline_t) && length(spline_t) == 1L) spline_t <- spline_t[[1L]]
  if (!is.null(spline_t$spline_blob_0) && nzchar(spline_t$spline_blob_0)) {
    .sendBlob(spline_t$spline_blob_0, "k2_spline_triples",
              party_conns[[1L]])
  }
  if (!is.null(spline_t$spline_blob_1) && nzchar(spline_t$spline_blob_1)) {
    .sendBlob(spline_t$spline_blob_1, "k2_spline_triples",
              party_conns[[2L]])
  }
  invisible(NULL)
}

.ot_beaver_prepare_vecmul <- function(datasources, party_conns, party_names,
                                      transport_pks, session_id, n, ring,
                                      .dsAgg, .sendBlob,
                                      beaver_key = NULL,
                                      dealer_conn = NULL) {
  mode <- .beaver_preprocessing_mode("vecmul", n, 1L, ring)
  if (identical(mode, "dealer")) {
    return(.dealer_prepare_vecmul(
      datasources = datasources, party_conns = party_conns,
      party_names = party_names, transport_pks = transport_pks,
      session_id = session_id, n = n, ring = ring, .dsAgg = .dsAgg,
      .sendBlob = .sendBlob, dealer_conn = dealer_conn))
  }
  prepare <- if (identical(mode, "iknp")) .iknp_beaver_prepare else
    .ot_beaver_prepare
  prepare(
    datasources = datasources,
    party_conns = party_conns,
    party_names = party_names,
    transport_pks = transport_pks,
    session_id = session_id,
    kind = "vecmul",
    n = n,
    p = 0L,
    ring = ring,
    .dsAgg = .dsAgg,
    .sendBlob = .sendBlob,
    beaver_key = beaver_key)
}

.ot_beaver_prepare_grad <- function(datasources, party_conns, party_names,
                                    transport_pks, session_id, n, p, ring,
                                    .dsAgg, .sendBlob,
                                    beaver_key = NULL,
                                    dealer_conn = NULL,
                                    grad_triple_key = "k2_grad_triple_fp") {
  if (!is.null(beaver_key) &&
      identical(grad_triple_key, "k2_grad_triple_fp")) {
    grad_triple_key <- beaver_key
  }
  mode <- .beaver_preprocessing_mode("grad", n, p, ring)
  if (identical(mode, "dealer")) {
    return(.dealer_prepare_grad(
      datasources = datasources, party_conns = party_conns,
      party_names = party_names, transport_pks = transport_pks,
      session_id = session_id, n = n, p = p, ring = ring, .dsAgg = .dsAgg,
      .sendBlob = .sendBlob, dealer_conn = dealer_conn,
      grad_triple_key = grad_triple_key))
  }
  prepare <- if (identical(mode, "iknp")) .iknp_beaver_prepare else
    .ot_beaver_prepare
  prepare(
    datasources = datasources,
    party_conns = party_conns,
    party_names = party_names,
    transport_pks = transport_pks,
    session_id = session_id,
    kind = "matvec",
    n = n,
    p = p,
    ring = ring,
    .dsAgg = .dsAgg,
    .sendBlob = .sendBlob,
    beaver_key = beaver_key)
}

.ot_beaver_prepare_spline <- function(datasources, party_conns, party_names,
                                      transport_pks, session_id, n, ring,
                                      .dsAgg, .sendBlob,
                                      beaver_key = NULL,
                                      dealer_conn = NULL) {
  mode <- .beaver_preprocessing_mode("spline", n, 3L, ring)
  if (identical(mode, "dealer")) {
    return(.dealer_prepare_spline(
      datasources = datasources, party_conns = party_conns,
      party_names = party_names, transport_pks = transport_pks,
      session_id = session_id, n = n, ring = ring, .dsAgg = .dsAgg,
      .sendBlob = .sendBlob, dealer_conn = dealer_conn))
  }
  if (is.null(beaver_key)) {
    beaver_key <- paste0("k2_ot_spline_",
                         format(Sys.time(), "%Y%m%d%H%M%OS3"),
                         "_", sample.int(.Machine$integer.max, 1L))
  }
  ops <- c(and = "spline_and", had1 = "spline_had1", had2 = "spline_had2")
  prepare <- if (identical(mode, "iknp")) .iknp_beaver_prepare else
    .ot_beaver_prepare
  for (op in names(ops)) {
    prepare(
      datasources = datasources,
      party_conns = party_conns,
      party_names = party_names,
      transport_pks = transport_pks,
      session_id = session_id,
      kind = "vecmul",
      n = n,
      p = 0L,
      ring = ring,
      .dsAgg = .dsAgg,
      .sendBlob = .sendBlob,
      beaver_key = paste0(beaver_key, "_", op),
      target = ops[[op]])
  }
  invisible(beaver_key)
}

.ot_beaver_prepare <- function(datasources, party_conns, party_names,
                               transport_pks, session_id, kind, n, p, ring,
                               .dsAgg, .sendBlob, beaver_key = NULL,
                               target = NULL) {
  if (length(party_conns) != 2L || length(party_names) != 2L) {
    stop("OT-Beaver preprocessing requires exactly two DCF parties",
         call. = FALSE)
  }
  ring <- as.integer(ring)
  if (!ring %in% c(63L, 127L)) stop("ring must be 63 or 127", call. = FALSE)
  n <- as.integer(n)
  p <- as.integer(p)
  if (is.null(beaver_key)) {
    beaver_key <- paste0("k2_ot_beaver_", kind, "_",
                         format(Sys.time(), "%Y%m%d%H%M%OS3"),
                         "_", sample.int(.Machine$integer.max, 1L))
  }
  beaver_key <- make.names(beaver_key)
  ot_n <- if (identical(kind, "matvec")) n * p else n

  for (ci in party_conns) {
    .dsAgg(datasources[ci], call(name = "k2OtBeaverSampleDS",
      kind = kind, n = as.integer(n), p = as.integer(p), ring = ring,
      beaver_key = beaver_key, session_id = session_id))
  }

  a_key <- paste0(beaver_key, "_a")
  b_key <- if (identical(kind, "matvec")) {
    paste0(beaver_key, "_b_expanded")
  } else {
    paste0(beaver_key, "_b")
  }
  cross_send_key <- paste0(beaver_key, "_cross_send")
  cross_receive_key <- paste0(beaver_key, "_cross_receive")

  .ot_beaver_direction(
    datasources = datasources,
    sender_ci = party_conns[[1L]],
    receiver_ci = party_conns[[2L]],
    sender_name = party_names[[1L]],
    receiver_name = party_names[[2L]],
    x_key = a_key,
    y_key = b_key,
    output_sender_key = cross_send_key,
    output_receiver_key = cross_receive_key,
    ot_key = paste0(beaver_key, "_dir12"),
    n = ot_n,
    ring = ring,
    session_id = session_id,
    .dsAgg = .dsAgg,
    .sendBlob = .sendBlob)

  .ot_beaver_direction(
    datasources = datasources,
    sender_ci = party_conns[[2L]],
    receiver_ci = party_conns[[1L]],
    sender_name = party_names[[2L]],
    receiver_name = party_names[[1L]],
    x_key = a_key,
    y_key = b_key,
    output_sender_key = cross_send_key,
    output_receiver_key = cross_receive_key,
    ot_key = paste0(beaver_key, "_dir21"),
    n = ot_n,
    ring = ring,
    session_id = session_id,
    .dsAgg = .dsAgg,
    .sendBlob = .sendBlob)

  if (is.null(target)) target <- if (identical(kind, "matvec")) "grad" else "vecmul"
  for (ci in party_conns) {
    .dsAgg(datasources[ci], call(name = "k2OtBeaverFinalizeDS",
      beaver_key = beaver_key,
      target = target,
      cross_send_key = cross_send_key,
      cross_receive_key = cross_receive_key,
      session_id = session_id))
  }
  invisible(beaver_key)
}

.iknp_beaver_prepare <- function(datasources, party_conns, party_names,
                                 transport_pks, session_id, kind, n, p, ring,
                                 .dsAgg, .sendBlob, beaver_key = NULL,
                                 target = NULL) {
  if (length(party_conns) != 2L || length(party_names) != 2L) {
    stop("IKNP OT-Beaver preprocessing requires exactly two DCF parties",
         call. = FALSE)
  }
  ring <- as.integer(ring)
  if (!ring %in% c(63L, 127L)) stop("ring must be 63 or 127", call. = FALSE)
  n <- as.integer(n)
  p <- as.integer(p)
  if (is.null(beaver_key)) {
    beaver_key <- paste0("k2_iknp_beaver_", kind, "_",
                         format(Sys.time(), "%Y%m%d%H%M%OS3"),
                         "_", sample.int(.Machine$integer.max, 1L))
  }
  beaver_key <- make.names(beaver_key)
  ot_n <- if (identical(kind, "matvec")) n * p else n

  for (ci in party_conns) {
    .dsAgg(datasources[ci], call(name = "k2OtBeaverSampleDS",
      kind = kind, n = as.integer(n), p = as.integer(p), ring = ring,
      beaver_key = beaver_key, session_id = session_id))
  }

  a_key <- paste0(beaver_key, "_a")
  b_key <- if (identical(kind, "matvec")) {
    paste0(beaver_key, "_b_expanded")
  } else {
    paste0(beaver_key, "_b")
  }
  cross_send_key <- paste0(beaver_key, "_cross_send")
  cross_receive_key <- paste0(beaver_key, "_cross_receive")

  .iknp_beaver_direction(
    datasources = datasources,
    sender_ci = party_conns[[1L]],
    receiver_ci = party_conns[[2L]],
    x_key = a_key,
    y_key = b_key,
    output_sender_key = cross_send_key,
    output_receiver_key = cross_receive_key,
    iknp_key = paste0(beaver_key, "_dir12"),
    n = ot_n,
    ring = ring,
    session_id = session_id,
    .dsAgg = .dsAgg,
    .sendBlob = .sendBlob)

  .iknp_beaver_direction(
    datasources = datasources,
    sender_ci = party_conns[[2L]],
    receiver_ci = party_conns[[1L]],
    x_key = a_key,
    y_key = b_key,
    output_sender_key = cross_send_key,
    output_receiver_key = cross_receive_key,
    iknp_key = paste0(beaver_key, "_dir21"),
    n = ot_n,
    ring = ring,
    session_id = session_id,
    .dsAgg = .dsAgg,
    .sendBlob = .sendBlob)

  if (is.null(target)) target <- if (identical(kind, "matvec")) "grad" else "vecmul"
  for (ci in party_conns) {
    .dsAgg(datasources[ci], call(name = "k2OtBeaverFinalizeDS",
      beaver_key = beaver_key,
      target = target,
      cross_send_key = cross_send_key,
      cross_receive_key = cross_receive_key,
      session_id = session_id))
  }
  invisible(beaver_key)
}

.ot_beaver_direction <- function(datasources, sender_ci, receiver_ci,
                                 sender_name, receiver_name,
                                 x_key, y_key,
                                 output_sender_key, output_receiver_key,
                                 ot_key, n, ring, session_id,
                                 .dsAgg, .sendBlob) {
  setup <- .dsAgg(datasources[sender_ci], call(name = "k2OtMulSenderSetupDS",
    ot_key = ot_key, session_id = session_id))
  if (is.list(setup) && length(setup) == 1L) setup <- setup[[1L]]

  choices <- .dsAgg(datasources[receiver_ci],
    call(name = "k2OtMulReceiverChoicesDS",
         public_setup = setup$public_setup,
         y_key = y_key,
         ot_key = ot_key,
         n = as.integer(n),
         ring = ring,
         session_id = session_id))
  if (is.list(choices) && length(choices) == 1L) choices <- choices[[1L]]

  points_key <- paste0(ot_key, "_points")
  .sendBlob(choices$points, points_key, sender_ci)
  enc <- .dsAgg(datasources[sender_ci],
    call(name = "k2OtMulSenderEncryptDS",
         points_blob_key = points_key,
         x_key = x_key,
         ot_key = ot_key,
         output_key = output_sender_key,
         n = as.integer(n),
         ring = ring,
         session_id = session_id))
  if (is.list(enc) && length(enc) == 1L) enc <- enc[[1L]]

  cts_key <- paste0(ot_key, "_ciphertexts")
  .sendBlob(enc$ciphertexts, cts_key, receiver_ci)
  .dsAgg(datasources[receiver_ci],
    call(name = "k2OtMulReceiverDecryptDS",
         ciphertexts_blob_key = cts_key,
         ot_key = ot_key,
         output_key = output_receiver_key,
         n = as.integer(n),
         ring = ring,
         session_id = session_id))
  invisible(NULL)
}

.iknp_beaver_direction <- function(datasources, sender_ci, receiver_ci,
                                   x_key, y_key,
                                   output_sender_key, output_receiver_key,
                                   iknp_key, n, ring, session_id,
                                   .dsAgg, .sendBlob) {
  setup <- .dsAgg(datasources[receiver_ci],
    call(name = "k2IknpBaseReceiverSetupDS",
         iknp_key = iknp_key,
         session_id = session_id))
  if (is.list(setup) && length(setup) == 1L) setup <- setup[[1L]]

  choices <- .dsAgg(datasources[sender_ci],
    call(name = "k2IknpBaseSenderChoicesDS",
         public_setup = setup$public_setup,
         iknp_key = iknp_key,
         session_id = session_id))
  if (is.list(choices) && length(choices) == 1L) choices <- choices[[1L]]

  points_key <- paste0(iknp_key, "_base_points")
  .sendBlob(choices$points, points_key, receiver_ci)
  base_ct <- .dsAgg(datasources[receiver_ci],
    call(name = "k2IknpBaseReceiverEncryptDS",
         points_blob_key = points_key,
         iknp_key = iknp_key,
         session_id = session_id))
  if (is.list(base_ct) && length(base_ct) == 1L) base_ct <- base_ct[[1L]]

  base_ct_key <- paste0(iknp_key, "_base_ciphertexts")
  .sendBlob(base_ct$ciphertexts, base_ct_key, sender_ci)
  .dsAgg(datasources[sender_ci],
    call(name = "k2IknpBaseSenderFinalizeDS",
         ciphertexts_blob_key = base_ct_key,
         iknp_key = iknp_key,
         session_id = session_id))

  ext <- .dsAgg(datasources[receiver_ci],
    call(name = "k2IknpReceiverExtendDS",
         y_key = y_key,
         iknp_key = iknp_key,
         n = as.integer(n),
         ring = ring,
         session_id = session_id))
  if (is.list(ext) && length(ext) == 1L) ext <- ext[[1L]]

  u_key <- paste0(iknp_key, "_u_matrix")
  .sendBlob(ext$u_matrix, u_key, sender_ci)
  enc <- .dsAgg(datasources[sender_ci],
    call(name = "k2IknpSenderEncryptDS",
         u_matrix_blob_key = u_key,
         x_key = x_key,
         iknp_key = iknp_key,
         output_key = output_sender_key,
         n = as.integer(n),
         ring = ring,
         session_id = session_id))
  if (is.list(enc) && length(enc) == 1L) enc <- enc[[1L]]

  cts_key <- paste0(iknp_key, "_ciphertexts")
  .sendBlob(enc$ciphertexts, cts_key, receiver_ci)
  .dsAgg(datasources[receiver_ci],
    call(name = "k2IknpReceiverDecryptDS",
         ciphertexts_blob_key = cts_key,
         iknp_key = iknp_key,
         output_key = output_receiver_key,
         n = as.integer(n),
         ring = ring,
         session_id = session_id))
  invisible(NULL)
}
