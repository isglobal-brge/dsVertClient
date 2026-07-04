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
                                         session_id) {
  pks <- list()
  identity_info <- list()
  for (srv in servers) {
    ci <- which(server_names == srv)
    r <- DSI::datashield.aggregate(datasources[ci],
      call(name = "glmRing63TransportInitDS", session_id = session_id))
    if (is.list(r) && length(r) == 1L) r <- r[[1L]]
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
  for (srv in servers) {
    ci <- which(server_names == srv)
    DSI::datashield.aggregate(datasources[ci],
      call(name = "mpcStoreTransportKeysDS",
           transport_keys_b64 = pk_b64, identity_info_b64 = id_b64,
           session_id = session_id))
  }
  pks
}
