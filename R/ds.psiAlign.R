#' Check whether a data symbol is pinned-padded-PSI aligned
#'
#' Each server validates a persistent, secret-token-authenticated manifest
#' against the current identifier values and row order. The client accepts the
#' result only when every named server returns the same fixed-schema public
#' attestation. No row count, identifier commitment or intersection size is
#' released.
#'
#' @param newobj Character. Symbol to validate on every server.
#' @param datasources Named DataSHIELD connections. When `NULL`, the current
#'   connections returned by [DSI::datashield.connections_find()] are used.
#' @return A list with `aligned` and `n_common`. `n_common` is always
#'   `NA_integer_`; the alignment-status route never releases cardinality.
#' @export
ds.isPsiAligned <- function(newobj = "DA", datasources = NULL) {
  if (is.null(datasources)) datasources <- DSI::datashield.connections_find()
  status <- .psi_alignment_status(newobj, datasources)
  list(aligned = status$aligned, n_common = status$n_common)
}

.psi_alignment_status <- function(newobj, datasources,
                                  .aggregate = DSI::datashield.aggregate) {
  invalid <- function() list(
    aligned = FALSE, n_common = NA_integer_, manifests = NULL)
  if (!is.character(newobj) || length(newobj) != 1L || is.na(newobj) ||
      !nzchar(newobj) || !is.list(datasources) || length(datasources) < 2L ||
      !is.function(.aggregate)) return(invalid())
  datasource_names <- names(datasources)
  if (is.null(datasource_names) || anyNA(datasource_names) ||
      any(!nzchar(datasource_names)) || anyDuplicated(datasource_names)) {
    return(invalid())
  }
  attestations <- tryCatch(
    .dsvert_aggregate_strict(
      conns = datasources,
      expr = call(
        name = "psiPaddedAttestationDS", data_name = newobj,
        session_id = ""),
      operation = "pinned padded PSI alignment attestation",
      .aggregate = .aggregate),
    error = function(e) NULL)
  if (!is.list(attestations) || length(attestations) != length(datasources) ||
      is.null(names(attestations)) ||
      !setequal(names(attestations), datasource_names)) return(invalid())
  validated <- tryCatch(
    lapply(attestations, .dsvert_validate_psi_padded_attestation),
    error = function(e) NULL)
  if (is.null(validated) ||
      !all(vapply(validated[-1L], identical, logical(1L), validated[[1L]]))) {
    return(invalid())
  }
  list(aligned = TRUE, n_common = NA_integer_, manifests = validated)
}

#' Align vertically partitioned records with pinned, fixed-capacity PSI
#'
#' Runs the sole public dsVert alignment protocol. Every participating server
#' pads its private input to a server-owned power-of-two capacity bucket. The
#' relay sees only authenticated ciphertexts with fixed shapes for that public
#' bucket; it never receives identifiers, row maps, membership bits or an exact
#' cardinality. Peer identities, the complete contract, phase receipts and
#' every server-to-server envelope are cryptographically bound before use.
#'
#' @param data_name Character. Source data-frame symbol on every server.
#' @param id_col Character. Identifier-column name on every server.
#' @param newobj Character. Destination symbol for the aligned data frames.
#' @param ref_server Deprecated compatibility argument. It must be `NULL`:
#'   the reference is selected deterministically from authenticated peer
#'   identities and cannot be chosen by the relay.
#' @param verbose Logical. Emit fixed, data-independent progress messages.
#' @param datasources Named DataSHIELD connections. At least two are required.
#' @param na.action Compatibility argument. Its value is validated but the
#'   server-owned `dsvert.psi.padded_missing_policy` determines eligibility;
#'   the analyst cannot override that disclosure policy.
#'
#' @return Invisibly, the common public protocol attestation. It contains the
#'   public capacity bucket and authenticated contract identifiers, never an
#'   input count or intersection cardinality.
#'
#' @details
#' The global all-server membership vector is evaluated with the purpose-bound
#' exact GC/OT backend between the two contract compute peers. Results are
#' delivered as authenticated shares to the reference and materialized in one
#' deterministic order on every server. A malicious relay can delay, drop or
#' replay traffic and thereby cause denial of service, but cannot forge a peer,
#' alter a contract or silently substitute a protocol message. As with any
#' two-party secure computation, collusion of both designated compute peers is
#' outside the non-collusion confidentiality claim. The public capacity bucket
#' is an upper bound on local input shape.
#'
#' This function is privacy-preserving preprocessing, not itself a statistical
#' release and not a differential-privacy mechanism. Subsequent result methods
#' remain subject to their own disclosure and DP contracts.
#'
#' @examples
#' \dontrun{
#' aligned <- ds.psiAlign(
#'   "D", "patient_id", "D_aligned", datasources = connections)
#' ds.isPsiAligned("D_aligned", datasources = connections)
#' }
#'
#' @importFrom DSI datashield.aggregate datashield.assign.expr datashield.connections_find
#' @export
ds.psiAlign <- function(data_name, id_col, newobj = "D_aligned",
                        ref_server = NULL, verbose = TRUE,
                        datasources = NULL,
                        na.action = c("na.omit", "na.fail", "none")) {
  values <- list(data_name = data_name, id_col = id_col, newobj = newobj)
  labels <- names(values)
  for (index in seq_along(values)) {
    value <- values[[index]]
    if (!is.character(value) || length(value) != 1L || is.na(value) ||
        !nzchar(value)) {
      stop(labels[[index]], " must be a single non-empty character string",
           call. = FALSE)
    }
  }
  if (!is.null(ref_server)) {
    stop(
      "ref_server must be NULL: pinned padded PSI selects its reference ",
      "from authenticated peer identities.", call. = FALSE)
  }
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("verbose must be TRUE or FALSE", call. = FALSE)
  }
  tryCatch(match.arg(na.action), error = function(e) stop(
    "na.action must be one of 'na.omit', 'na.fail' or 'none'",
    call. = FALSE))
  if (is.null(datasources)) {
    datasources <- DSI::datashield.connections_find()
  }
  server_names <- names(datasources)
  if (!is.list(datasources) || length(datasources) < 2L ||
      is.null(server_names) || anyNA(server_names) ||
      any(!nzchar(server_names)) || anyDuplicated(server_names)) {
    stop(
      "Pinned padded PSI requires at least two uniquely named DataSHIELD connections.",
      call. = FALSE)
  }
  .dsvert_maybe_negotiate_dsi_chunk_size(datasources)
  on.exit(.dsvert_reset_chunk_size(), add = TRUE)
  if (verbose) message("Starting pinned, fixed-capacity PSI alignment.")
  attestation <- .dsvert_psi_padded_align(
    data_name = data_name, id_col = id_col, newobj = newobj,
    datasources = datasources)
  if (verbose) message("Pinned PSI alignment completed and attested.")
  invisible(attestation)
}
