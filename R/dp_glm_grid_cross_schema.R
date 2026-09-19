# Client-side authority checks for the versioned cross-owner grid contract.
# These helpers inspect public signed objects only and never access a dataset.

.dsvert_dp_glm_grid_cross_transcript_stop <- function(error) {
  stop(structure(list(
    message = "[dsvert_dp_public_failure:v1] Protected capsule operation failed.",
    call = NULL, code = "dp_protected_operation_failed"),
    class = c("dsvert_dp_public_failure", "error", "condition")))
}

.dsvert_dp_glm_grid_cross_verify <- function(message, public_key, signature) {
  tryCatch({
    public <- .dsvert_joint_dp_client_b64url(
      public_key, 32L, "cross grid identity key")
    signature <- .dsvert_joint_dp_client_b64url(
      signature, 64L, "cross grid signature")
    isTRUE(openssl::ed25519_verify(
      message, signature, openssl::read_ed25519_pubkey(public)))
  }, error = function(error) FALSE)
}

.dsvert_dp_glm_grid_cross_schema_levels <- function(value) {
  fail <- function() stop("Invalid cross grid signed schema.", call. = FALSE)
  # Canonical JSON arrays decode as unnamed lists. Their signed labels remain
  # character scalars; reject objects and mixed/coerced label types.
  if (is.list(value) && is.null(names(value)) && length(value) &&
      all(vapply(value, function(label) is.character(label) && length(label) == 1L &&
        !is.na(label), logical(1L)))) value <- unlist(value, use.names = FALSE)
  if (!is.atomic(value) || !is.null(dim(value))) fail()
  declared <- attr(value, "class", exact = TRUE)
  factor_value <- identical(declared, "factor") ||
    identical(declared, c("ordered", "factor"))
  plain <- is.null(declared)
  levels <- if (factor_value) {
    labels <- attr(value, "levels", exact = TRUE)
    codes <- value
    attributes(codes) <- NULL
    if (typeof(codes) != "integer" || !is.character(labels) ||
        anyNA(labels) || anyDuplicated(labels) ||
        any(!is.na(codes) & (codes < 1L | codes > length(labels)))) fail()
    unname(labels[codes])
  } else if (plain && is.character(value)) {
    unname(value)
  } else if (plain && is.logical(value)) {
    ifelse(is.na(value), NA_character_, ifelse(value, "TRUE", "FALSE"))
  } else if (plain && typeof(value) == "integer") {
    result <- rep(NA_character_, length(value))
    present <- !is.na(value)
    result[present] <- sprintf("%d", value[present])
    result
  } else if (plain && typeof(value) == "double") {
    present <- !is.na(value)
    if (any(!is.finite(value[present])) ||
        any(value[present] != floor(value[present])) ||
        any(abs(value[present]) > 2^53 - 1)) fail()
    result <- rep(NA_character_, length(value))
    normalized <- value[present]
    normalized[normalized == 0] <- 0
    result[present] <- sprintf("%.0f", normalized)
    result
  } else fail()
  levels <- enc2utf8(levels)
  if (!length(levels) || anyNA(levels) || anyDuplicated(levels) ||
      any(!nzchar(trimws(levels)))) fail()
  sort(unname(levels), method = "radix")
}

.dsvert_dp_glm_grid_cross_schema_validate <- function(
    policy, logical_snapshot, schema_manifest, signature_verifier) {
  tryCatch({
  exact <- function(value, fields) {
    is.list(value) && !is.null(names(value)) && !anyNA(names(value)) &&
      !anyDuplicated(names(value)) && setequal(names(value), fields)
  }
  fail <- function() stop("Invalid cross grid signed schema.", call. = FALSE)
  identifier <- function(value) {
    is.character(value) && length(value) == 1L && !is.na(value) &&
      grepl("^[A-Za-z0-9][A-Za-z0-9._:-]{0,127}$", value)
  }
  named <- function(value) {
    is.list(value) && length(value) > 0L && !is.null(names(value)) &&
      !anyNA(names(value)) && !anyDuplicated(names(value)) &&
      all(nzchar(names(value)))
  }
  if (!exact(schema_manifest, c("version", "logical_snapshot",
      "peer_pinset_sha256", "datasets", "signatures")) ||
      !identical(schema_manifest$version,
                 .DSVERT_CLIENT_DP_CAPSULE_SCHEMA_VERSION) ||
      !is.function(signature_verifier)) fail()
  snapshot <- schema_manifest$logical_snapshot
  if (!exact(snapshot, c("logical_snapshot_id", "version",
      "alignment_protocol_version")) ||
      !identifier(snapshot$logical_snapshot_id) ||
      !identifier(snapshot$version) ||
      !is.numeric(snapshot$alignment_protocol_version) ||
      length(snapshot$alignment_protocol_version) != 1L ||
      !is.finite(snapshot$alignment_protocol_version) ||
      snapshot$alignment_protocol_version < 1 ||
      snapshot$alignment_protocol_version !=
        floor(snapshot$alignment_protocol_version) ||
      !identical(.dsvert_joint_dp_client_canonical(snapshot),
                 .dsvert_joint_dp_client_canonical(logical_snapshot))) fail()
  pins <- policy$peer_pinset
  if (!is.character(pins) || length(pins) < 2L ||
      is.null(names(pins)) || anyNA(pins) || anyNA(names(pins)) ||
      any(!nzchar(names(pins))) || anyDuplicated(names(pins)) ||
      anyDuplicated(unname(pins))) fail()
  for (pin in pins) .dsvert_joint_dp_client_b64url(
    pin, 32L, "cross grid pinned identity")
  pins <- pins[order(names(pins), method = "radix")]
  pin_hash <- .dsvert_dp_capsule_source_hash(as.list(pins))
  if (!identical(schema_manifest$peer_pinset_sha256, pin_hash) ||
      !identical(policy$peer_pinset_sha256, pin_hash)) fail()
  datasets <- schema_manifest$datasets
  if (!named(datasets)) fail()
  seen_references <- seen_physical <- alignment_groups <- character()
  for (data_name in names(datasets)) {
    dataset <- datasets[[data_name]]
    if (!identifier(data_name) || !exact(dataset, c("dataset_id",
        "dataset_version", "schema_version", "alignment_group",
        "patient_keys", "columns")) ||
        !all(vapply(dataset[c("dataset_id", "dataset_version",
          "schema_version", "alignment_group")], identifier, logical(1L))) ||
        !named(dataset$patient_keys) ||
        !all(names(dataset$patient_keys) %in% names(pins)) ||
        !all(vapply(dataset$patient_keys, identifier, logical(1L))) ||
        !named(dataset$columns)) fail()
    owners <- character()
    for (column_name in names(dataset$columns)) {
      column <- dataset$columns[[column_name]]
      parts <- strsplit(column_name, "$", fixed = TRUE)[[1L]]
      if (!length(parts) %in% 1:2 ||
          !identical(paste(parts, collapse = "$"), column_name) ||
          !all(vapply(parts, identifier, logical(1L))) ||
          !is.list(column) || !identifier(column$owner_peer) ||
          !column$owner_peer %in% names(pins) ||
          (length(parts) == 2L &&
           !identical(parts[[1L]], column$owner_peer))) fail()
      if (identical(column$kind, "numeric")) {
        if (!exact(column, c("kind", "owner_peer", "lower", "upper")) ||
            !is.numeric(column$lower) || length(column$lower) != 1L ||
            !is.numeric(column$upper) || length(column$upper) != 1L ||
            !is.finite(column$lower) || !is.finite(column$upper) ||
            column$lower >= column$upper ||
            !is.finite((column$upper - column$lower)^2)) fail()
      } else if (identical(column$kind, "categorical")) {
        if (!exact(column, c("kind", "owner_peer", "levels"))) fail()
        column$levels <- .dsvert_dp_glm_grid_cross_schema_levels(column$levels)
        dataset$columns[[column_name]] <- column
      } else fail()
      physical <- paste(column$owner_peer, data_name,
                        tail(parts, 1L), sep = "\r")
      if (column_name %in% seen_references || physical %in% seen_physical)
        fail()
      seen_references <- c(seen_references, column_name)
      seen_physical <- c(seen_physical, physical)
      owners <- c(owners, column$owner_peer)
    }
    if (!setequal(names(dataset$patient_keys), unique(owners))) fail()
    alignment_groups <- c(alignment_groups, dataset$alignment_group)
    datasets[[data_name]] <- dataset
  }
  if (length(unique(alignment_groups)) != 1L) fail()
  unsigned <- .dsvert_joint_dp_client_canonical(list(
    version = .DSVERT_CLIENT_DP_CAPSULE_SCHEMA_VERSION,
    logical_snapshot = snapshot, peer_pinset_sha256 = pin_hash,
    datasets = datasets))
  signatures <- schema_manifest$signatures
  if (!named(signatures) || !setequal(names(signatures), names(pins)) ||
      !all(vapply(signatures, function(value) {
        is.character(value) && length(value) == 1L && !is.na(value) &&
          grepl("^[A-Za-z0-9_-]{86}$", value)
      }, logical(1L)))) fail()
  message <- charToRaw(paste0(.DSVERT_CLIENT_DP_CAPSULE_SCHEMA_SIGNATURE_DOMAIN,
                             .dsvert_joint_dp_client_json(unsigned)))
  for (peer in names(pins)) {
    if (!isTRUE(tryCatch(signature_verifier(message, unname(pins[[peer]]),
          signatures[[peer]]), error = function(error) FALSE))) fail()
  }
  list(unsigned = unsigned, sha256 = .dsvert_dp_capsule_source_hash(unsigned),
       signers = names(pins), signatures =
         .dsvert_joint_dp_client_canonical(signatures[names(pins)]))
  }, error = .dsvert_dp_glm_grid_cross_transcript_stop)
}
