#' @title ECDH-PSI Record Alignment
#' @description Privacy-preserving record alignment using Elliptic Curve
#'   Diffie-Hellman Private Set Intersection (ECDH-PSI). Aligns data frames
#'   across vertically partitioned DataSHIELD servers so that rows correspond
#'   to the same individuals, without exposing identifiers.
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
#' using ECDH-PSI instead of SHA-256 hashing for stronger privacy guarantees.
#'
#' \subsection{Protocol overview}{
#' ECDH-PSI exploits the commutativity of elliptic curve scalar multiplication:
#' \eqn{\alpha \cdot (\beta \cdot H(id)) = \beta \cdot (\alpha \cdot H(id))}.
#'
#' \enumerate{
#'   \item The reference server hashes IDs to P-256 curve points and multiplies
#'     by a random scalar \eqn{\alpha}. Returns masked points to client.
#'   \item For each target server:
#'     \itemize{
#'       \item The target generates scalar \eqn{\beta}, double-masks ref points
#'         with \eqn{\beta} (stores locally), and masks own IDs with \eqn{\beta}.
#'       \item The reference double-masks target points with \eqn{\alpha}.
#'       \item The target matches double-masked sets to find the intersection,
#'         then reorders its data to match the reference order.
#'     }
#'   \item A multi-server intersection ensures only records present on ALL
#'     servers are retained.
#' }
#' }
#'
#' \subsection{Security (DDH assumption on P-256)}{
#' \itemize{
#'   \item The client sees only opaque elliptic curve points — not reversible
#'     to identifiers, not vulnerable to dictionary attacks.
#'   \item Each server's scalar never leaves the server.
#'   \item The DDH assumption prevents the client from linking single-masked
#'     points across servers.
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
  # This avoids R's C stack overflow for large call() expressions.
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

  cat("=== ECDH-PSI Record Alignment ===\n")
  cat(sprintf("Reference: %s, Targets: %s\n\n", ref_server,
              paste(target_names, collapse = ", ")))

  # ==================================================================
  # Phase 1: Reference server masks its IDs
  # ==================================================================
  # The ref server generates a random P-256 scalar α and computes
  # α·H(id) for each ID. These masked points are returned to the client.
  # The scalar α stays on the server (NEVER transmitted).
  cat("[Phase 1] Reference server masking IDs...\n")
  ref_result <- DSI::datashield.aggregate(
    conns = ref_conn,
    expr = call("psiMaskIdsDS", data_name, id_col)
  )
  ref_result <- ref_result[[1]]
  cat(sprintf("  %s: %d IDs masked\n", ref_server, ref_result$n))

  # For each target server: Phases 2-7
  for (target_name in target_names) {
    target_conn <- datasources[target_name]

    # ==================================================================
    # Phases 2+3: Client relays ref masked points to target server.
    # EC points are sent via chunked blob storage to prevent R's C stack
    # overflow with large datasets. Each point is ~44 bytes base64url;
    # at n=100K records the vector would be ~14MB inline — well above
    # R's call expression parser limit.
    # ==================================================================
    cat(sprintf("[Phase 2-3] %s: processing...\n", target_name))
    ref_points_blob <- paste(ref_result$masked_points, collapse = ",")
    .storeLargeBlob("ref_masked_points", ref_points_blob, target_conn)

    target_result <- DSI::datashield.aggregate(
      conns = target_conn,
      expr = call("psiProcessTargetDS", data_name, id_col,
                  from_storage = TRUE)
    )
    target_result <- target_result[[1]]
    cat(sprintf("  %s: %d IDs masked\n", target_name, target_result$n))

    # ==================================================================
    # Phases 4+5: Client relays target's masked IDs to the ref server.
    # Ref server double-masks them: α·(β·H(id)). Same blob storage
    # pattern for the n-length EC point vector.
    # ==================================================================
    cat(sprintf("[Phase 4-5] %s: double-masking via %s...\n",
                target_name, ref_server))
    target_points_blob <- paste(target_result$own_masked_points, collapse = ",")
    .storeLargeBlob("target_masked_points", target_points_blob, ref_conn)

    dm_result <- DSI::datashield.aggregate(
      conns = ref_conn,
      expr = call("psiDoubleMaskDS", from_storage = TRUE)
    )
    dm_result <- dm_result[[1]]

    # ==================================================================
    # Phases 6+7: Client relays double-masked own points to target.
    # Target matches them against stored double-masked ref points,
    # identifies common IDs, and reorders its data to match ref order.
    # ==================================================================
    cat(sprintf("[Phase 6-7] %s: matching and aligning...\n", target_name))
    dm_points_blob <- paste(dm_result$double_masked_points, collapse = ",")
    .storeLargeBlob("double_masked_points", dm_points_blob, target_conn)

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
  # After Phases 2-7, each target has aligned to the reference, but
  # different target servers may have matched different subsets of ref IDs
  # (e.g. server2 has patients 1-40, server3 has patients 10-50).
  # We intersect all matched index sets to keep only records present
  # on ALL servers. This ensures the final aligned data is consistent.
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
  # The index vector grows with n_common; at n=200K the integer vector
  # would be ~15MB inline, exceeding R's call expression limit.
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
