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
#' This function replaces \code{\link{ds.hashId}} + \code{\link{ds.alignRecords}}
#' with a single call that provides stronger privacy guarantees.
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
    # Target generates its own scalar β, double-masks ref points
    # (β·(α·H(id)) — stored locally for matching), and returns its own
    # masked IDs (β·H(id)) to the client.
    # ==================================================================
    cat(sprintf("[Phase 2-3] %s: processing...\n", target_name))
    target_result <- DSI::datashield.aggregate(
      conns = target_conn,
      expr = call("psiProcessTargetDS", data_name, id_col,
                  ref_result$masked_points)
    )
    target_result <- target_result[[1]]
    cat(sprintf("  %s: %d IDs masked\n", target_name, target_result$n))

    # ==================================================================
    # Phases 4+5: Client relays target's masked IDs to the ref server.
    # Ref server double-masks them: α·(β·H(id)). By ECDH commutativity,
    # α·β·H(id) = β·α·H(id), so the target can match these against its
    # stored double-masked ref points to find common IDs.
    # ==================================================================
    cat(sprintf("[Phase 4-5] %s: double-masking via %s...\n",
                target_name, ref_server))
    dm_result <- DSI::datashield.aggregate(
      conns = ref_conn,
      expr = call("psiDoubleMaskDS", target_result$own_masked_points)
    )
    dm_result <- dm_result[[1]]

    # ==================================================================
    # Phases 6+7: Client relays double-masked own points to target.
    # Target matches them against stored double-masked ref points,
    # identifies common IDs, and reorders its data to match ref order.
    # ==================================================================
    cat(sprintf("[Phase 6-7] %s: matching and aligning...\n", target_name))
    DSI::datashield.assign(
      conns = target_conn,
      symbol = newobj,
      value = call("psiMatchAndAlignDS", data_name,
                   dm_result$double_masked_points)
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

  # Broadcast common indices to all servers and filter
  for (name in server_names) {
    conn <- datasources[name]
    DSI::datashield.assign(
      conns = conn,
      symbol = newobj,
      value = call("psiFilterCommonDS", newobj, common_indices)
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
