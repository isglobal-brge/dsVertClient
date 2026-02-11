#' @title Align Records Across Servers
#' @description Client-side function that aligns records on DataSHIELD servers
#'   to match a reference set of hashed identifiers. This ensures observations
#'   are properly matched across vertically partitioned data.
#'
#' @param data_name Character string. Name of the data frame on each server.
#' @param id_col Character string. Name of the identifier column.
#' @param reference_hashes Character vector. Hashes from the reference server
#'   (obtained via \code{\link{ds.hashId}}).
#' @param newobj Character string. Name for the aligned data frame on servers.
#'   Default is "D_aligned".
#' @param algo Character string. Hash algorithm (must match reference).
#'   Default is "sha256".
#' @param datasources DataSHIELD connection object or list of connections.
#'   If NULL, uses all available connections.
#'
#' @return Invisibly returns a list with alignment statistics for each server:
#'   \itemize{
#'     \item \code{n_matched}: Number of records matched
#'     \item \code{n_reference}: Number of reference hashes
#'   }
#'
#' @details
#' This function performs record alignment for vertically partitioned data:
#'
#' \enumerate{
#'   \item Takes reference hashes (from \code{ds.hashId} on one server)
#'   \item For each target server:
#'     \itemize{
#'       \item Hashes local identifiers
#'       \item Matches against reference hashes
#'       \item Reorders data to match reference order
#'       \item Creates new aligned data frame
#'     }
#' }
#'
#' After alignment, all servers will have:
#' \itemize{
#'   \item Same number of observations
#'   \item Observations in the same order (by identifier)
#'   \item Only observations present in all partitions
#' }
#'
#' @seealso \code{\link{ds.hashId}} for obtaining reference hashes
#'
#' @examples
#' \dontrun{
#' # Step 1: Get reference hashes from first server
#' ref <- ds.hashId("D", "patient_id", datasource = conns$server1)
#'
#' # Step 2: Align all servers (including reference) to these hashes
#' ds.alignRecords("D", "patient_id", ref$hashes,
#'                 newobj = "D_aligned", datasources = conns)
#'
#' # Now "D_aligned" on all servers has matching, ordered observations
#' }
#'
#' @importFrom DSI datashield.assign datashield.aggregate datashield.connections_find
#' @export
ds.alignRecords <- function(data_name, id_col, reference_hashes,
                            newobj = "D_aligned", algo = "sha256",
                            datasources = NULL) {
  # Validate inputs
  if (!is.character(data_name) || length(data_name) != 1) {
    stop("data_name must be a single character string", call. = FALSE)
  }
  if (!is.character(id_col) || length(id_col) != 1) {
    stop("id_col must be a single character string", call. = FALSE)
  }
  if (!is.character(reference_hashes) || length(reference_hashes) == 0) {
    stop("reference_hashes must be a non-empty character vector", call. = FALSE)
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

  # For each server, call the align function
  stats <- list()

  for (i in seq_along(datasources)) {
    conn <- datasources[i]
    server_name <- names(datasources)[i]
    if (is.null(server_name)) server_name <- paste0("server", i)

    # Build the call expression using call()
    call_expr <- call("alignRecordsDS", data_name, id_col, reference_hashes, algo)

    # Execute assign on server
    DSI::datashield.assign(
      conns = conn,
      symbol = newobj,
      value = call_expr
    )

    # Get count to verify alignment
    count_call <- call("getObsCountDS", newobj)
    count_result <- DSI::datashield.aggregate(conn, count_call)

    if (is.list(count_result) && length(count_result) == 1) {
      count_result <- count_result[[1]]
    }

    stats[[server_name]] <- list(
      n_matched = count_result$n_obs,
      n_reference = length(reference_hashes)
    )

    message(sprintf(
      "Server '%s': %d of %d records matched (%.1f%%)",
      server_name,
      count_result$n_obs,
      length(reference_hashes),
      100 * count_result$n_obs / length(reference_hashes)
    ))
  }

  invisible(stats)
}
