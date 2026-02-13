#' @title ECDH-PSI Record Alignment (Blind Relay)
#' @description Privacy-preserving record alignment using Elliptic Curve
#'   Diffie-Hellman Private Set Intersection (ECDH-PSI) with blind-relay
#'   transport encryption. Aligns data frames across vertically partitioned
#'   DataSHIELD servers so that rows correspond to the same individuals.
#'   The client never sees raw EC points — only opaque encrypted blobs.
#'
#' @param data_name Character string. Name of the data frame on each server.
#' @param id_col Character string. Name of the identifier column.
#' @param newobj Character string. Name for the aligned data frame on servers.
#'   Default is \code{"D_aligned"}.
#' @param ref_server Character string or NULL. Name of the reference server.
#'   If NULL (default), the first connection is used.
#' @param datasources DataSHIELD connection object or list of connections.
#'   If NULL, uses all available connections.
#'
#' @return Invisibly returns a list with alignment statistics for each server:
#'   \itemize{
#'     \item \code{n_matched}: Number of records matched
#'     \item \code{n_total}: Number of records on that server
#'   }
#'
#' @details
#' This function performs privacy-preserving record alignment in a single call,
#' using ECDH-PSI with blind-relay transport encryption.
#'
#' \subsection{Protocol overview}{
#' ECDH-PSI exploits the commutativity of elliptic curve scalar multiplication:
#' \eqn{\alpha \cdot (\beta \cdot H(id)) = \beta \cdot (\alpha \cdot H(id))}.
#'
#' All EC point exchanges are encrypted server-to-server (X25519 + AES-256-GCM
#' ECIES). The client acts as a blind relay, seeing only opaque blobs.
#'
#' \enumerate{
#'   \item \strong{Phase 0}: Each server generates an X25519 transport keypair.
#'     Public keys are exchanged via the client.
#'   \item \strong{Phase 1}: The reference server masks IDs with scalar
#'     \eqn{\alpha}. Points are stored server-side (not returned to client).
#'   \item For each target server:
#'     \itemize{
#'       \item The reference encrypts masked points under the target's PK.
#'       \item The target decrypts, generates scalar \eqn{\beta}, double-masks
#'         ref points (stores locally), masks own IDs, encrypts them under
#'         the ref's PK.
#'       \item The reference decrypts, double-masks with \eqn{\alpha},
#'         encrypts result under target's PK.
#'       \item The target decrypts, matches double-masked sets, aligns data.
#'     }
#'   \item A multi-server intersection ensures only records present on ALL
#'     servers are retained.
#' }
#' }
#'
#' \subsection{Security (DDH assumption on P-256, malicious-client model)}{
#' \itemize{
#'   \item The client sees only opaque encrypted blobs — not EC points.
#'   \item Each server's scalar never leaves the server.
#'   \item PSI firewall: phase ordering + one-shot semantics prevent OPRF
#'     oracle attacks.
#' }
#' }
#'
#' @seealso \code{\link{ds.vertCor}}, \code{\link{ds.vertGLM}} for analysis
#'   functions that operate on aligned data.
#'
#' @examples
#' \dontrun{
#' # Align records across all servers using PSI
#' ds.psiAlign("D", "patient_id", "D_aligned", datasources = connections)
#'
#' # Now "D_aligned" on all servers has matching, ordered observations
#' }
#'
#' @importFrom DSI datashield.aggregate datashield.assign datashield.connections_find
#' @export
ds.psiAlign <- function(data_name, id_col, newobj = "D_aligned",
                         ref_server = NULL, datasources = NULL) {
  # Validate inputs
  if (!is.character(data_name) || length(data_name) != 1) {
    stop("data_name must be a single character string", call. = FALSE)
  }
  if (!is.character(id_col) || length(id_col) != 1) {
    stop("id_col must be a single character string", call. = FALSE)
  }
  if (!is.character(newobj) || length(newobj) != 1) {
    stop("newobj must be a single character string", call. = FALSE)
  }

  # Get datasources
  if (is.null(datasources)) {
    datasources <- DSI::datashield.connections_find()
  }
  if (length(datasources) == 0) {
    stop("No DataSHIELD connections found", call. = FALSE)
  }

  server_names <- names(datasources)
  n_servers <- length(datasources)

  # Determine reference server
  if (is.null(ref_server)) {
    ref_server <- server_names[1]
  }
  if (!ref_server %in% server_names) {
    stop("ref_server '", ref_server, "' not found in connections", call. = FALSE)
  }

  ref_conn <- datasources[ref_server]
  target_names <- setdiff(server_names, ref_server)

  chunk_size <- 10000  # DataSHIELD parser-safe chunk size

  # Store a large blob on a server via chunked mheStoreBlobDS calls.
  .storeLargeBlob <- function(key, data, conn) {
    n_chars <- nchar(data)
    n_chunks <- ceiling(n_chars / chunk_size)
    for (ch in seq_len(n_chunks)) {
      start <- (ch - 1) * chunk_size + 1
      end <- min(ch * chunk_size, n_chars)
      DSI::datashield.aggregate(
        conns = conn,
        expr = call("mheStoreBlobDS",
                    key = key,
                    chunk = substr(data, start, end),
                    chunk_index = as.integer(ch),
                    n_chunks = as.integer(n_chunks))
      )
    }
  }

  cat("=== ECDH-PSI Record Alignment (Blind Relay) ===\n")
  cat(sprintf("Reference: %s, Targets: %s\n\n", ref_server,
              paste(target_names, collapse = ", ")))

  # ==================================================================
  # Phase 0: PSI Transport Key Exchange
  # ==================================================================
  # Each server generates an X25519 transport keypair. The client
  # collects PKs and distributes them so servers can encrypt messages
  # for each other. The client never receives the secret keys.
  cat("[Phase 0] Transport key exchange...\n")
  psi_transport_pks <- list()
  for (name in server_names) {
    conn <- datasources[name]
    result <- DSI::datashield.aggregate(
      conns = conn,
      expr = call("psiInitDS")
    )
    psi_transport_pks[[name]] <- result[[1]]$transport_pk
  }

  # Distribute transport PKs to all servers
  for (name in server_names) {
    conn <- datasources[name]
    # Build key map: all server names + "ref" alias
    keys_for_server <- psi_transport_pks
    keys_for_server[["ref"]] <- psi_transport_pks[[ref_server]]
    DSI::datashield.aggregate(
      conns = conn,
      expr = call("psiStoreTransportKeysDS", keys_for_server)
    )
  }
  cat("  Transport keys exchanged.\n")

  # ==================================================================
  # Phase 1: Reference server masks its IDs
  # ==================================================================
  # The ref server generates a random P-256 scalar α and computes
  # α·H(id) for each ID. Points are stored server-side (NOT returned).
  cat("[Phase 1] Reference server masking IDs...\n")
  ref_result <- DSI::datashield.aggregate(
    conns = ref_conn,
    expr = call("psiMaskIdsDS", data_name, id_col)
  )
  ref_result <- ref_result[[1]]
  cat(sprintf("  %s: %d IDs masked (stored server-side)\n",
              ref_server, ref_result$n))

  # For each target server: Phases 2-7 (all via encrypted blobs)
  for (target_name in target_names) {
    target_conn <- datasources[target_name]

    # ==================================================================
    # Phase 2: Ref exports encrypted masked points for this target
    # ==================================================================
    cat(sprintf("[Phase 2] %s: exporting encrypted points for %s...\n",
                ref_server, target_name))
    export_result <- DSI::datashield.aggregate(
      conns = ref_conn,
      expr = call("psiExportMaskedDS", target_name)
    )
    encrypted_ref_blob <- export_result[[1]]$encrypted_blob

    # ==================================================================
    # Phase 3: Client relays encrypted blob to target.
    # Target decrypts, processes (masks own IDs, double-masks ref points),
    # returns encrypted own masked points under ref's PK.
    # ==================================================================
    cat(sprintf("[Phase 3] %s: processing (blind relay)...\n", target_name))
    .storeLargeBlob("ref_encrypted_blob", encrypted_ref_blob, target_conn)

    target_result <- DSI::datashield.aggregate(
      conns = target_conn,
      expr = call("psiProcessTargetDS", data_name, id_col,
                  from_storage = TRUE)
    )
    target_result <- target_result[[1]]
    cat(sprintf("  %s: %d IDs masked\n", target_name, target_result$n))
    encrypted_target_blob <- target_result$encrypted_blob

    # ==================================================================
    # Phases 4+5: Client relays target's encrypted blob to ref.
    # Ref decrypts, double-masks with α, re-encrypts under target's PK.
    # One-shot per target (PSI firewall enforced).
    # ==================================================================
    cat(sprintf("[Phase 4-5] %s: double-masking via %s (blind relay)...\n",
                target_name, ref_server))
    .storeLargeBlob("target_encrypted_blob", encrypted_target_blob, ref_conn)

    dm_result <- DSI::datashield.aggregate(
      conns = ref_conn,
      expr = call("psiDoubleMaskDS",
                  target_name = target_name,
                  from_storage = TRUE)
    )
    encrypted_dm_blob <- dm_result[[1]]$encrypted_blob

    # ==================================================================
    # Phases 6+7: Client relays encrypted double-masked blob to target.
    # Target decrypts, matches against stored ref double-masked points,
    # and reorders its data to match ref order.
    # ==================================================================
    cat(sprintf("[Phase 6-7] %s: matching and aligning (blind relay)...\n",
                target_name))
    .storeLargeBlob("dm_encrypted_blob", encrypted_dm_blob, target_conn)

    DSI::datashield.assign(
      conns = target_conn,
      symbol = newobj,
      value = call("psiMatchAndAlignDS", data_name,
                   from_storage = TRUE)
    )
  }

  # ==================================================================
  # Phase 7 for ref: Self-align (identity operation)
  # ==================================================================
  cat(sprintf("[Phase 7] %s: self-aligning...\n", ref_server))
  DSI::datashield.assign(
    conns = ref_conn,
    symbol = newobj,
    value = call("psiSelfAlignDS", data_name)
  )

  # ==================================================================
  # Phase 8: Multi-server intersection
  # ==================================================================
  # Integer indices are safe aggregate statistics (no EC points).
  cat("[Phase 8] Computing multi-server intersection...\n")
  all_indices <- list()
  for (name in server_names) {
    conn <- datasources[name]
    idx <- DSI::datashield.aggregate(
      conns = conn,
      expr = call("psiGetMatchedIndicesDS")
    )
    all_indices[[name]] <- idx[[1]]
  }

  # Intersect all index sets
  common_indices <- Reduce(intersect, all_indices)
  common_indices <- sort(as.integer(common_indices))

  cat(sprintf("  Common records: %d\n", length(common_indices)))

  # Broadcast common indices to all servers via blob storage and filter.
  for (name in server_names) {
    conn <- datasources[name]
    indices_blob <- paste(common_indices, collapse = ",")
    .storeLargeBlob("common_indices", indices_blob, conn)
    DSI::datashield.assign(
      conns = conn,
      symbol = newobj,
      value = call("psiFilterCommonDS", newobj, from_storage = TRUE)
    )
  }

  # Verify alignment
  stats <- list()
  for (name in server_names) {
    conn <- datasources[name]
    count <- DSI::datashield.aggregate(
      conns = conn,
      expr = call("getObsCountDS", newobj)
    )
    count <- count[[1]]
    n_server <- all_indices[[name]]
    stats[[name]] <- list(
      n_matched = count$n_obs,
      n_total = length(n_server)
    )
    message(sprintf(
      "Server '%s': %d of %d records matched (%.1f%%)",
      name, count$n_obs, length(n_server),
      100 * count$n_obs / max(length(n_server), 1)
    ))
  }

  cat("PSI alignment complete.\n")

  invisible(stats)
}
