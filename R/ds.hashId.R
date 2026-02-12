#' @title Hash Identifier Column
#' @description Client-side function that retrieves hashed identifiers from
#'   a DataSHIELD server. Used as the first step in record matching for
#'   vertically partitioned data.
#'
#' @param data_name Character string. Name of the data frame on the server.
#' @param id_col Character string. Name of the identifier column to hash.
#' @param algo Character string. Hash algorithm. Default is "sha256".
#' @param datasource A DataSHIELD connection object (DSConnection).
#'   If NULL, uses the first available connection.
#'
#' @return A list containing:
#'   \itemize{
#'     \item \code{hashes}: Character vector of hashed identifiers
#'     \item \code{n}: Number of observations
#'   }
#'
#' @details
#' This function is typically called on the "reference" server in a
#' vertical partitioning scenario. The returned hashes are then passed
#' to \code{\link{ds.alignRecords}} on other servers to align the data.
#'
#' @seealso \code{\link{ds.alignRecords}} for aligning records across servers
#'
#' @examples
#' \dontrun{
#' # Get hashes from reference server
#' ref_hashes <- ds.hashId("D", "patient_id", datasource = conns$server1)
#'
#' # Use hashes to align other servers
#' ds.alignRecords("D", "patient_id", ref_hashes$hashes,
#'                 datasources = conns[c("server2", "server3")])
#' }
#'
#' @importFrom DSI datashield.aggregate datashield.connections_find
#' @export
ds.hashId <- function(data_name, id_col, algo = "sha256", datasource = NULL) {
  .Deprecated("ds.psiAlign",
    msg = "ds.hashId() is deprecated. Use ds.psiAlign() for ECDH-PSI alignment (stronger privacy).")
  # Validate inputs
  if (!is.character(data_name) || length(data_name) != 1) {
    stop("data_name must be a single character string", call. = FALSE)
  }
  if (!is.character(id_col) || length(id_col) != 1) {
    stop("id_col must be a single character string", call. = FALSE)
  }

  # Get datasource
  if (is.null(datasource)) {
    datasource <- DSI::datashield.connections_find()
    if (length(datasource) == 0) {
      stop("No DataSHIELD connections found", call. = FALSE)
    }
    # Use only the first connection for reference hashes
    datasource <- datasource[1]
  }

  # Build the call using substitute to properly construct the expression
  call_expr <- call("hashIdDS", data_name, id_col, algo)

  # Execute on server
  result <- DSI::datashield.aggregate(
    conns = datasource,
    expr = call_expr
  )

  # Return the first (and should be only) result
  if (is.list(result) && length(result) == 1) {
    return(result[[1]])
  }
  return(result)
}
