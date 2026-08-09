# Public, machine-readable handling for a regenerated peer identity.  Only
# logical names and SHA-256 fingerprints of public Ed25519 keys cross this
# boundary; private keys and seeds never do.

.DSVERT_CLIENT_PEER_NOT_RECOGNIZED_VERSION <-
  "DSVERT_PEER_NOT_RECOGNIZED_V1"

.dsvert_client_identity_fingerprint <- function(identity_pk) {
  raw <- .dsvert_joint_dp_client_b64url(
    identity_pk, 32L, "peer identity public key")
  digest::digest(raw, algo = "sha256", serialize = FALSE)
}

.dsvert_client_peer_not_recognized_condition <- function(
    peer_name, observed_fingerprint_sha256,
    expected_fingerprint_sha256 = "unconfigured") {
  valid_peer <- is.character(peer_name) && length(peer_name) == 1L &&
    !is.na(peer_name) &&
    grepl("^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$", peer_name)
  valid_fingerprint <- function(value, allow_unconfigured = FALSE) {
    is.character(value) && length(value) == 1L && !is.na(value) &&
      (grepl("^[0-9a-f]{64}$", value) ||
         (allow_unconfigured && value %in% c("unconfigured", "inconsistent")))
  }
  if (!valid_peer ||
      !valid_fingerprint(observed_fingerprint_sha256) ||
      !valid_fingerprint(expected_fingerprint_sha256, TRUE)) {
    stop("Invalid peer-not-recognized error contract", call. = FALSE)
  }
  mismatch <- switch(
    expected_fingerprint_sha256,
    unconfigured = paste0(
      "logical peer '", peer_name,
      "' has no name-bound Ed25519 pin. Observed SHA-256 ",
      observed_fingerprint_sha256),
    inconsistent = paste0(
      "logical peer '", peer_name,
      "' has inconsistent name-bound Ed25519 pins across the participating ",
      "servers. Observed SHA-256 ", observed_fingerprint_sha256),
    paste0(
      "logical peer '", peer_name,
      "' does not match its name-bound Ed25519 pin. Expected SHA-256 ",
      expected_fingerprint_sha256, "; observed ",
      observed_fingerprint_sha256))
  message <- paste0(
    "peer_not_recognized: ", mismatch, ". Server administrator: from a ",
    "trusted administrative dsVertClient session, call ds.getIdentityPks(), ",
    "then verify the observed fingerprint out of band directly at the ",
    "affected site's console. On each other participating server that pins '",
    peer_name, "', persist the new name-bound dsvert.trusted_peers (or ",
    "dsvert.trusted_peer_", peer_name, ") value; the affected server must not ",
    "pin its own identity. Restart or reconnect the affected DataSHIELD ",
    "services and sessions, then retry. Never approve a replacement key ",
    "solely from the analyst/relay.")
  structure(list(
    message = message, call = NULL, code = "peer_not_recognized",
    protocol = .DSVERT_CLIENT_PEER_NOT_RECOGNIZED_VERSION,
    peer_name = peer_name,
    expected_fingerprint_sha256 = expected_fingerprint_sha256,
    observed_fingerprint_sha256 = observed_fingerprint_sha256,
    admin_action =
      "verify_out_of_band_then_update_name_bound_pin_and_restart"),
    class = c("dsvert_peer_not_recognized", "error", "condition"))
}

.dsvert_client_stop_peer_not_recognized <- function(...) {
  stop(.dsvert_client_peer_not_recognized_condition(...))
}

.dsvert_client_parse_peer_not_recognized <- function(message) {
  if (!is.character(message) || !length(message) || is.na(message[[1L]])) {
    return(NULL)
  }
  message <- substr(as.character(message[[1L]]), 1L, 8192L)
  pattern <- paste0(
    .DSVERT_CLIENT_PEER_NOT_RECOGNIZED_VERSION,
    "\\|peer=([A-Za-z0-9][A-Za-z0-9._-]{0,127})",
    "\\|expected=(unconfigured|inconsistent|[0-9a-f]{64})",
    "\\|observed=([0-9a-f]{64})(?![0-9a-f])")
  match <- regexec(pattern, message, perl = TRUE)
  fields <- regmatches(message, match)[[1L]]
  if (length(fields) != 4L) return(NULL)
  .dsvert_client_peer_not_recognized_condition(
    fields[[2L]], fields[[4L]], fields[[3L]])
}
