# Public protocol negotiation and closed provenance epoch schemas.  No source
# lookup is permitted until every participating connection has passed preflight.
.DSVERT_CLIENT_PROVENANCE_PROTOCOL <- "dsvert-provenance-epochs-v1"

.dsvert_provenance_preflight <- function(
    datasources, .aggregate = DSI::datashield.aggregate) {
  expected <- list(
    protocol = .DSVERT_CLIENT_PROVENANCE_PROTOCOL,
    version = "1.4.1",
    families = c("psi", "count", "frequency", "synopsis"))
  values <- tryCatch(.dsvert_aggregate_strict(
    datasources, call(name = "dsvertProvenanceCapabilitiesDS"),
    operation = "provenance epoch protocol preflight", .aggregate = .aggregate),
    error = function(error) stop(
      "All peers require the dsVert 1.4.1 provenance epoch protocol.",
      call. = FALSE))
  if (!all(vapply(values, identical, logical(1L), expected))) {
    stop("All peers require the dsVert 1.4.1 provenance epoch protocol.",
         call. = FALSE)
  }
  invisible(TRUE)
}

.dsvert_provenance_epochs_validate <- function(value, identities = NULL) {
  invalid <- function() stop("Invalid complete authorization epoch vector.",
                             call. = FALSE)
  if (!is.list(value) || !is.null(names(value)) || !length(value) ||
      length(value) > 4096L) invalid()
  peers <- vapply(value, function(entry) {
    fields <- c("peer_identity", "authorization_epoch", "source_contract_id")
    if (!is.list(entry) || is.null(names(entry)) || anyDuplicated(names(entry)) ||
        !setequal(names(entry), fields)) invalid()
    identity <- tryCatch(.dsvert_dp_analysis_client_identity_pk(
      entry$peer_identity), error = function(error) NULL)
    if (is.null(identity) || !identical(identity, entry$peer_identity) ||
        !all(vapply(entry[c("authorization_epoch", "source_contract_id")],
          function(field) is.character(field) && length(field) == 1L &&
            !is.na(field) && grepl("^[0-9a-f]{64}$", field), logical(1L)))) {
      invalid()
    }
    identity
  }, character(1L))
  if (anyDuplicated(peers) ||
      !identical(peers, sort(peers, method = "radix"))) invalid()
  epochs <- vapply(value, `[[`, character(1L), "authorization_epoch")
  if (anyDuplicated(epochs)) invalid()
  if (!is.null(identities)) {
    expected <- sort(unname(vapply(identities,
      .dsvert_dp_analysis_client_identity_pk, character(1L))), method = "radix")
    if (anyDuplicated(expected) || !identical(peers, expected)) invalid()
  }
  value
}

.dsvert_provenance_snapshot_commitment <- function(
    provenance, config, peer_name, request) {
  epochs <- .dsvert_provenance_epochs_validate(
    provenance$authorization_epochs, unname(config$peer_pins))
  value <- list(
    protocol = .DSVERT_CLIENT_PROVENANCE_PROTOCOL,
    domain = config$domain, cohort_id = config$cohort_id,
    owner_identity_pk = unname(config$peer_pins[[peer_name]]),
    dataset_id = config$dataset_id, dataset_version = config$dataset_version,
    alignment_purpose = config$alignment_purpose,
    derivation_policy_id = provenance$derivation_policy_id,
    authorization_epochs = epochs, request = request)
  digest::digest(charToRaw(paste0(
    "dsVert/provenance/stable-snapshot/v1|",
    .dsvert_joint_dp_client_json(value))), algo = "sha256", serialize = FALSE)
}

# The caller first verifies the server release signatures and artifact binding.
# This receipt separates stable public semantics from per-run authentication.
.dsvert_provenance_stable_receipt <- function(artifact) {
  semantic <- .dsvert_joint_dp_client_canonical(artifact$semantic)
  expected <- digest::digest(charToRaw(paste0(
    "dsVert/analysis-artifact-key/v1|",
    .dsvert_joint_dp_client_json(semantic))), algo = "sha256", serialize = FALSE)
  if (!identical(artifact$artifact_key, expected)) {
    stop("Invalid provenance receipt artifact binding.", call. = FALSE)
  }
  list(version = "dsvert-provenance-stable-receipt-v1",
       artifact_key = expected, semantic = semantic)
}
