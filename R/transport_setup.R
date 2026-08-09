# Establish the identity-verified peer transport set on every participating
# server. Each server publishes its transport public key (and Ed25519 identity
# signature) via glmRing63TransportInitDS; the sorted set is stored back on each
# server via mpcStoreTransportKeysDS, which verifies identities and records the
# peers in ss$peer_transport_pks. Sealing primitives then pin recipients to that
# set, so a share can never be sealed to an analyst-supplied key.
#
# Returns the named list of transport public keys (by server name).

#' @keywords internal
.dsvert_setup_peer_transport <- function(datasources, server_names, servers,
                                         session_id,
                                         .aggregate = DSI::datashield.aggregate) {
  .dsvert_maybe_negotiate_dsi_chunk_size(
    datasources[servers], .aggregate)
  pks <- list()
  identity_info <- list()
  init_results <- .dsvert_aggregate_strict(
    datasources[servers],
    call(name = "glmRing63TransportInitDS", session_id = session_id),
    operation = "peer transport initialization", .aggregate = .aggregate)
  for (srv in servers) {
    r <- init_results[[srv]]
    pks[[srv]] <- r$transport_pk
    if (!is.null(r$identity_pk)) {
      identity_info[[srv]] <- list(identity_pk = r$identity_pk,
                                   signature = r$signature)
    }
  }
  pk_sorted <- pks[sort(names(pks))]
  id_sorted <- if (length(identity_info) > 0) identity_info[sort(names(identity_info))] else NULL
  to_b64url <- function(x) gsub("\\+", "-", gsub("/", "_", gsub("=+$", "", x, perl = TRUE), fixed = TRUE), fixed = TRUE)
  json_b64 <- function(x) to_b64url(gsub("\n", "", jsonlite::base64_enc(charToRaw(jsonlite::toJSON(x, auto_unbox = TRUE))), fixed = TRUE))
  pk_b64 <- json_b64(pk_sorted)
  id_b64 <- if (!is.null(id_sorted)) json_b64(id_sorted) else ""
  .dsvert_aggregate_strict(
    datasources[servers],
    call(name = "mpcStoreTransportKeysDS",
         transport_keys_b64 = pk_b64, identity_info_b64 = id_b64,
         session_id = session_id),
    operation = "peer transport pinning", result_contract = "logical_true",
    .aggregate = .aggregate)
  pks
}
