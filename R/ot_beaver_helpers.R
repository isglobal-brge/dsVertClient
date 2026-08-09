#' Client-side Beaver preprocessing orchestration
#'
#' These helpers negotiate the IKNP OT-extension backend required by every
#' participating server. Participating-party dealer preprocessing is retired:
#' it can expose peer operands through the online Beaver openings.
#'
#' @keywords internal
#' @noRd

.iknp_base_cache <- new.env(parent = emptyenv())
.beaver_policy_cache <- new.env(parent = emptyenv())

.normalise_beaver_mode <- function(mode, allow_auto = TRUE) {
  mode <- tolower(as.character(mode)[1L])
  if (!nzchar(mode)) mode <- "auto"
  if (identical(mode, "ot")) mode <- "iknp"
  valid <- "iknp"
  if (allow_auto) valid <- c("auto", valid)
  if (!mode %in% valid) {
    stop("Invalid dsvert.beaver_preprocessing value: ", mode,
         ". Supported values are 'auto', 'iknp' and 'ot' ",
         "('ot' is an alias for 'iknp').",
         call. = FALSE)
  }
  mode
}

.beaver_fetch_policies <- function(datasources, party_conns, .dsAgg,
                                   session_id = NULL) {
  if (is.null(datasources) || is.null(party_conns) || is.null(.dsAgg)) {
    return(NULL)
  }
  cache_key <- paste(session_id %||% "no_session",
                     paste(as.integer(party_conns), collapse = ","),
                     sep = "|")
  cached <- get0(cache_key, envir = .beaver_policy_cache, inherits = FALSE)
  if (!is.null(cached)) return(cached)

  policies <- lapply(party_conns, function(ci) {
    out <- tryCatch(
      .dsAgg(datasources[ci], call(name = "dsvertBeaverPolicyDS")),
      error = function(e) {
        if (inherits(e, c("dsvert_peer_not_recognized",
                          "dsvert_dsi_poisoned_session"))) {
          stop(e)
        }
        stop("Could not verify the server IKNP policy: ",
             conditionMessage(e), call. = FALSE)
      })
    if (is.list(out) && length(out) == 1L &&
        is.list(out[[1L]]) && !is.null(out[[1L]]$allowed)) {
      out <- out[[1L]]
    }
    out$allowed <- unique(vapply(out$allowed, .normalise_beaver_mode,
                                 character(1), allow_auto = FALSE))
    out$preferred <- .normalise_beaver_mode(out$preferred %||% "iknp",
                                            allow_auto = FALSE)
    out$minimum <- .normalise_beaver_mode(out$minimum %||% "iknp",
                                          allow_auto = FALSE)
    out$requires_iknp <- TRUE
    out
  })
  assign(cache_key, policies, envir = .beaver_policy_cache)
  policies
}

.beaver_preprocessing_mode <- function(kind, n, p, ring,
                                       datasources = NULL,
                                       party_conns = NULL,
                                       .dsAgg = NULL,
                                       session_id = NULL) {
  mode <- getOption("dsvert.beaver_preprocessing", "auto")
  mode <- .normalise_beaver_mode(mode)
  policies <- .beaver_fetch_policies(datasources, party_conns, .dsAgg,
                                     session_id = session_id)
  if (is.null(policies)) {
    if (identical(mode, "auto")) return("iknp")
    return(.normalise_beaver_mode(mode, allow_auto = FALSE))
  }

  allowed_all <- Reduce(intersect, lapply(policies, `[[`, "allowed"))
  if (!length(allowed_all)) {
    stop("No common Beaver preprocessing backend is allowed by all servers",
         call. = FALSE)
  }
  if (!"iknp" %in% allowed_all) {
    stop("At least one server does not allow IKNP preprocessing",
         call. = FALSE)
  }
  "iknp"
}

.ot_beaver_prepare_vecmul <- function(datasources, party_conns, party_names,
                                      transport_pks, session_id, n, ring,
                                      .dsAgg, .sendBlob,
                                      beaver_key = NULL,
                                      dealer_conn = NULL) {
  mode <- .beaver_preprocessing_mode(
    "vecmul", n, 1L, ring, datasources, party_conns, .dsAgg, session_id)
  iknp <- function() .iknp_beaver_prepare(
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
  # DEALER MODE REMOVED (F2/F17): a participating-party dealer reconstructs the
  # peer's operands from the online openings d=x-a, e=y-b. IKNP OT-extension is
  # now the SOLE, dealer-free Beaver backend — no party ever holds a full triple.
  iknp()
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
  mode <- .beaver_preprocessing_mode(
    "grad", n, p, ring, datasources, party_conns, .dsAgg, session_id)
  iknp <- function() .iknp_beaver_prepare(
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
  # DEALER MODE REMOVED (F2/F17): IKNP OT-extension is the sole dealer-free
  # Beaver backend — no party ever holds a full triple.
  iknp()
}

.ot_beaver_prepare_spline <- function(datasources, party_conns, party_names,
                                      transport_pks, session_id, n, ring,
                                      .dsAgg, .sendBlob,
                                      beaver_key = NULL,
                                      dealer_conn = NULL) {
  mode <- .beaver_preprocessing_mode(
    "spline", n, 3L, ring, datasources, party_conns, .dsAgg, session_id)
  iknp <- function() {
    if (is.null(beaver_key)) {
      beaver_key <- paste0("k2_ot_spline_", .dsvert_uuid4())
    }
    ops <- c(and = "spline_and", had1 = "spline_had1", had2 = "spline_had2")
    for (op in names(ops)) {
      .iknp_beaver_prepare(
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
  # DEALER MODE REMOVED (F2/F17): IKNP OT-extension is the sole dealer-free
  # Beaver backend — no party ever holds a full triple.
  iknp()
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
    beaver_key <- paste0("k2_iknp_beaver_", kind, "_", .dsvert_uuid4())
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
    sender_name = party_names[[1L]],
    receiver_name = party_names[[2L]],
    sender_pk = transport_pks[[party_names[[1L]]]],
    receiver_pk = transport_pks[[party_names[[2L]]]],
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
    sender_name = party_names[[2L]],
    receiver_name = party_names[[1L]],
    sender_pk = transport_pks[[party_names[[2L]]]],
    receiver_pk = transport_pks[[party_names[[1L]]]],
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

.iknp_beaver_direction <- function(datasources, sender_ci, receiver_ci,
                                   x_key, y_key,
                                   output_sender_key, output_receiver_key,
                                   iknp_key, n, ring, session_id,
                                   sender_name, receiver_name,
                                   sender_pk, receiver_pk,
                                   .dsAgg, .sendBlob) {
  base_key <- make.names(sprintf("k2_iknp_base_%s_%s_ring%d",
                                 sender_name,
                                 receiver_name,
                                 as.integer(ring)))
  cache_key <- paste(session_id, sender_name, receiver_name, ring, sep = "|")
  base_cached <- isTRUE(get0(cache_key, envir = .iknp_base_cache,
                             inherits = FALSE))
  if (!base_cached) {
    setup <- .dsAgg(datasources[receiver_ci],
      call(name = "k2IknpBaseReceiverSetupDS",
           iknp_key = base_key,
           session_id = session_id))
    if (is.list(setup) && length(setup) == 1L) setup <- setup[[1L]]

    choices <- .dsAgg(datasources[sender_ci],
      call(name = "k2IknpBaseSenderChoicesDS",
           public_setup = setup$public_setup,
           iknp_key = base_key,
           recipient_pk = receiver_pk,
           ring = ring,
           session_id = session_id))
    if (is.list(choices) && length(choices) == 1L) choices <- choices[[1L]]

    .sendBlob(choices$points, choices$points_transfer, receiver_ci)
    base_ct <- .dsAgg(datasources[receiver_ci],
      call(name = "k2IknpBaseReceiverEncryptDS",
           iknp_key = base_key,
           producer_name = sender_name,
           recipient_pk = sender_pk,
           ring = ring,
           session_id = session_id))
    if (is.list(base_ct) && length(base_ct) == 1L) base_ct <- base_ct[[1L]]

    .sendBlob(base_ct$ciphertexts, base_ct$ciphertexts_transfer, sender_ci)
    .dsAgg(datasources[sender_ci],
      call(name = "k2IknpBaseSenderFinalizeDS",
           iknp_key = base_key,
           producer_name = receiver_name,
           ring = ring,
           session_id = session_id))
    assign(cache_key, TRUE, envir = .iknp_base_cache)
  }

  ext <- .dsAgg(datasources[receiver_ci],
    call(name = "k2IknpReceiverExtendDS",
         y_key = y_key,
         iknp_key = iknp_key,
         base_key = base_key,
         recipient_pk = sender_pk,
         n = as.integer(n),
         ring = ring,
         session_id = session_id))
  if (is.list(ext) && length(ext) == 1L) ext <- ext[[1L]]

  .sendBlob(ext$u_matrix, ext$u_matrix_transfer, sender_ci)
  # Relay the mandatory KOS15 consistency-check opener with the U matrix (no
  # extra round). There is intentionally no analyst/custodian downgrade knob.
  if (is.null(ext$kos_check) || !nzchar(ext$kos_check)) {
    stop("IKNP KOS consistency-check opener missing", call. = FALSE)
  }
  enc <- .dsAgg(datasources[sender_ci],
    call(name = "k2IknpSenderEncryptDS",
         x_key = x_key,
         iknp_key = iknp_key,
         base_key = base_key,
         output_key = output_sender_key,
         n = as.integer(n),
         ring = ring,
         producer_name = receiver_name,
         recipient_pk = receiver_pk,
         kos_check = ext$kos_check,
         session_id = session_id))
  if (is.list(enc) && length(enc) == 1L) enc <- enc[[1L]]

  .sendBlob(enc$ciphertexts, enc$ciphertexts_transfer, receiver_ci)
  .dsAgg(datasources[receiver_ci],
    call(name = "k2IknpReceiverDecryptDS",
         iknp_key = iknp_key,
         output_key = output_receiver_key,
         n = as.integer(n),
         ring = ring,
         producer_name = sender_name,
         session_id = session_id))
  invisible(NULL)
}
