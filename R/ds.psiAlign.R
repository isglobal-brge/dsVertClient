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
#' @param verbose Logical. If TRUE (default), print progress messages.
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
                         ref_server = NULL, verbose = TRUE,
                         datasources = NULL) {
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

  on.exit({
    .dsvert_reset_chunk_size()
  }, add = TRUE)

  # Determine reference server
  if (is.null(ref_server)) {
    ref_server <- server_names[1]
  }
  if (!ref_server %in% server_names) {
    stop("ref_server '", ref_server, "' not found in connections", call. = FALSE)
  }

  ref_conn <- datasources[ref_server]
  target_names <- setdiff(server_names, ref_server)

  # Generate session_id for PSI protocol state isolation
  session_id <- local({
    hex <- sample(c(0:9, letters[1:6]), 32, replace = TRUE)
    hex[13] <- "4"  # UUID v4
    hex[17] <- sample(c("8","9","a","b"), 1)  # variant 1
    paste0(
      paste(hex[1:8], collapse = ""), "-",
      paste(hex[9:12], collapse = ""), "-",
      paste(hex[13:16], collapse = ""), "-",
      paste(hex[17:20], collapse = ""), "-",
      paste(hex[21:32], collapse = "")
    )
  })

  # Store a large blob on a server with adaptive chunking and fallback
  .storeLargeBlob <- function(key, data, conn) {
    .dsvert_adaptive_send(data, function(chunk_str, chunk_idx, n_chunks) {
      if (n_chunks == 1L) {
        DSI::datashield.aggregate(
          conns = conn,
          expr = call("mpcStoreBlobDS", key = key, chunk = chunk_str,
                      session_id = session_id)
        )
      } else {
        DSI::datashield.aggregate(
          conns = conn,
          expr = call("mpcStoreBlobDS",
                      key = key,
                      chunk = chunk_str,
                      chunk_index = chunk_idx,
                      n_chunks = n_chunks,
                      session_id = session_id)
        )
      }
    })
  }

  if (verbose) message("=== ECDH-PSI Record Alignment (Blind Relay) ===")
  if (verbose) message(sprintf("Reference: %s, Targets: %s\n", ref_server,
              paste(target_names, collapse = ", ")))

  # ==================================================================
  # Phase 0: PSI Transport Key Exchange
  # ==================================================================
  # Each server generates an X25519 transport keypair. The client
  # collects PKs and distributes them so servers can encrypt messages
  # for each other. The client never receives the secret keys.
  if (verbose) message("[Phase 0] Transport key exchange...")
  # Initialize PSI on all servers in parallel
  init_results <- DSI::datashield.aggregate(
    conns = datasources,
    expr = call("psiInitDS", session_id = session_id)
  )
  psi_transport_pks <- list()
  psi_identity_info <- list()
  for (name in server_names) {
    psi_transport_pks[[name]] <- init_results[[name]]$transport_pk
    if (!is.null(init_results[[name]]$identity_pk))
      psi_identity_info[[name]] <- list(
        identity_pk = init_results[[name]]$identity_pk,
        signature   = init_results[[name]]$signature)
  }

  # Distribute transport PKs + identity info to all servers in parallel
  keys_for_server <- psi_transport_pks
  keys_for_server[["ref"]] <- psi_transport_pks[[ref_server]]
  id_for_server <- if (length(psi_identity_info) > 0) psi_identity_info else NULL
  if (!is.null(id_for_server))
    id_for_server[["ref"]] <- psi_identity_info[[ref_server]]
  # Store transport keys per-server (list args as base64url JSON to avoid Opal parser)
  .to_b64url <- function(x) gsub("\\+","-",gsub("/","_",gsub("=+$","",x,perl=TRUE),fixed=TRUE),fixed=TRUE)
  .json_to_b64url <- function(x) .to_b64url(gsub("\n","",jsonlite::base64_enc(charToRaw(jsonlite::toJSON(x, auto_unbox = TRUE))),fixed=TRUE))
  keys_b64 <- .json_to_b64url(keys_for_server)
  id_b64 <- if (!is.null(id_for_server)) .json_to_b64url(id_for_server) else ""
  for (name in server_names) {
    ci <- which(names(datasources) == name)
    DSI::datashield.aggregate(
      conns = datasources[ci],
      expr = call("psiStoreTransportKeysDS",
                    transport_keys_b64 = keys_b64,
                    identity_info_b64 = id_b64,
                    session_id = session_id))
  }
  if (verbose) message("  Transport keys exchanged.")

  # ==================================================================
  # Phase 1: Reference server masks its IDs
  # ==================================================================
  # The ref server generates a random P-256 scalar α and computes
  # α·H(id) for each ID. Points are stored server-side (NOT returned).
  if (verbose) message("[Phase 1] Reference server masking IDs...")
  ref_result <- DSI::datashield.aggregate(
    conns = ref_conn,
    expr = call("psiMaskIdsDS", data_name, id_col,
                  session_id = session_id)
  )
  ref_result <- ref_result[[1]]
  if (verbose) message(sprintf("  %s: %d IDs masked (stored server-side)",
              ref_server, ref_result$n))

  # Phase 2: Ref exports encrypted masked points for each target
  encrypted_ref_blobs <- list()
  for (target_name in target_names) {
    if (verbose) message(sprintf("[Phase 2] %s: exporting encrypted points for %s...",
                ref_server, target_name))
    export_result <- DSI::datashield.aggregate(
      conns = ref_conn,
      expr = call("psiExportMaskedDS", target_name, session_id = session_id)
    )
    encrypted_ref_blobs[[target_name]] <- export_result[[1]]$encrypted_blob
  }

  # Deliver blobs to all targets (sequential — different blobs per target)
  for (target_name in target_names) {
    .storeLargeBlob("ref_encrypted_blob", encrypted_ref_blobs[[target_name]],
                    datasources[target_name])
  }

  # Phase 3: ALL targets process in PARALLEL (independent operations)
  if (verbose) message(sprintf("[Phase 3] %s: processing (parallel)...",
              paste(target_names, collapse = ", ")))
  target_results <- DSI::datashield.aggregate(
    conns = datasources[target_names],
    expr = call("psiProcessTargetDS", data_name, id_col,
                from_storage = TRUE, session_id = session_id)
  )
  for (target_name in target_names) {
    tr <- target_results[[target_name]]
    if (verbose) message(sprintf("  %s: %d IDs masked", target_name, tr$n))
  }

  # Phases 4-7: Sequential per target (ref must process each individually)
  for (target_name in target_names) {
    encrypted_target_blob <- target_results[[target_name]]$encrypted_blob

    # Phase 4-5: Relay target blob to ref for double-masking
    if (verbose) message(sprintf("[Phase 4-5] %s: double-masking via %s (blind relay)...",
                target_name, ref_server))
    .storeLargeBlob("target_encrypted_blob", encrypted_target_blob, ref_conn)

    dm_result <- DSI::datashield.aggregate(
      conns = ref_conn,
      expr = call("psiDoubleMaskDS", target_name = target_name,
                  from_storage = TRUE, session_id = session_id)
    )
    encrypted_dm_blob <- dm_result[[1]]$encrypted_blob

    # Phase 6-7: Relay double-masked blob to target for matching
    if (verbose) message(sprintf("[Phase 6-7] %s: matching and aligning (blind relay)...",
                target_name))
    .storeLargeBlob("dm_encrypted_blob", encrypted_dm_blob,
                    datasources[target_name])

    DSI::datashield.assign(
      conns = datasources[target_name],
      symbol = newobj,
      value = call("psiMatchAndAlignDS", data_name,
                   from_storage = TRUE, session_id = session_id)
    )
  }

  # ==================================================================
  # Phase 7 for ref: Self-align (identity operation)
  # ==================================================================
  if (verbose) message(sprintf("[Phase 7] %s: self-aligning...", ref_server))
  DSI::datashield.assign(
    conns = ref_conn,
    symbol = newobj,
    value = call("psiSelfAlignDS", data_name,
                   session_id = session_id)
  )

  # ==================================================================
  # Phase 8: Multi-server intersection
  # ==================================================================
  # Integer indices are safe aggregate statistics (no EC points).
  if (verbose) message("[Phase 8] Computing multi-server intersection...")
  all_indices <- list()
  # Collect matched indices from all servers in parallel
  idx_results <- DSI::datashield.aggregate(
    conns = datasources,
    expr = call("psiGetMatchedIndicesDS", session_id = session_id)
  )
  for (name in server_names) all_indices[[name]] <- idx_results[[name]]

  # Intersect all index sets
  common_indices <- Reduce(intersect, all_indices)
  common_indices <- sort(as.integer(common_indices))

  if (verbose) message(sprintf("  Common records: %d", length(common_indices)))

  # Broadcast common indices to all servers via blob storage
  indices_blob <- paste(common_indices, collapse = ",")
  for (name in server_names) {
    .storeLargeBlob("common_indices", indices_blob, datasources[name])
  }
  # Filter in parallel
  DSI::datashield.assign(
    conns = datasources,
    symbol = newobj,
    value = call("psiFilterCommonDS", newobj, from_storage = TRUE,
                   session_id = session_id)
  )

  # Verify alignment in parallel
  count_results <- DSI::datashield.aggregate(
    conns = datasources,
    expr = call("getObsCountDS", newobj)
  )
  stats <- list()
  for (name in server_names) {
    count <- count_results[[name]]
    n_server <- all_indices[[name]]
    stats[[name]] <- list(
      n_matched = count$n_obs,
      n_total = length(n_server)
    )
    if (verbose) message(sprintf(
      "Server '%s': %d of %d records matched (%.1f%%)",
      name, count$n_obs, length(n_server),
      100 * count$n_obs / max(length(n_server), 1)
    ))
  }

  if (verbose) message("PSI alignment complete.")

  stats$n_common <- length(common_indices)
  invisible(stats)
}
