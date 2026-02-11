#' @title Validate Identifier Format Across Servers
#' @description Client-side function that validates identifier format
#'   consistency across all DataSHIELD servers before record alignment.
#'
#' @param data_name Character string. Name of the data frame on each server.
#' @param id_col Character string. Name of the identifier column.
#' @param pattern Character string (optional). Regular expression pattern
#'   that IDs should match on all servers.
#' @param datasources DataSHIELD connection object(s). If NULL, uses all
#'   available connections.
#'
#' @return A list with class "ds.id.validation" containing:
#'   \itemize{
#'     \item \code{valid}: Logical, TRUE if formats are consistent across servers
#'     \item \code{servers}: Data frame with validation results per server
#'     \item \code{format_match}: Logical, TRUE if all servers have same format signature
#'     \item \code{pattern_match}: Logical, TRUE if all servers match the pattern
#'       (only if pattern provided)
#'     \item \code{warnings}: Character vector of any warnings detected
#'   }
#'
#' @details
#' This function should be called before \code{ds.alignRecords} to ensure
#' that identifier formats are consistent across all data partitions. It
#' helps catch common issues like:
#' \itemize{
#'   \item Different ID formats (e.g., "001" vs "1" vs "ID001")
#'   \item Missing or duplicate identifiers
#'   \item Type mismatches (numeric vs character)
#' }
#'
#' The function uses format signatures (hashed format characteristics) to
#' compare formats without revealing actual identifier values.
#'
#' @seealso \code{\link{ds.hashId}}, \code{\link{ds.alignRecords}}
#'
#' @examples
#' \dontrun{
#' # Validate ID format before alignment
#' validation <- ds.validateIdFormat("D", "patient_id", datasources = conns)
#' print(validation)
#'
#' if (validation$valid) {
#'   # Proceed with alignment
#'   ref_hashes <- ds.hashId("D", "patient_id", datasource = conns["server1"])
#'   ds.alignRecords("D", "patient_id", ref_hashes$hashes, "D_aligned", conns)
#' }
#'
#' # With pattern validation
#' validation <- ds.validateIdFormat("D", "patient_id",
#'                                   pattern = "^[A-Z]{2}[0-9]{6}$",
#'                                   datasources = conns)
#' }
#'
#' @importFrom DSI datashield.aggregate datashield.connections_find
#' @export
ds.validateIdFormat <- function(data_name, id_col, pattern = NULL,
                                 datasources = NULL) {
  # Validate inputs
  if (!is.character(data_name) || length(data_name) != 1) {
    stop("data_name must be a single character string", call. = FALSE)
  }
  if (!is.character(id_col) || length(id_col) != 1) {
    stop("id_col must be a single character string", call. = FALSE)
  }

  # Get datasources
  if (is.null(datasources)) {
    datasources <- DSI::datashield.connections_find()
  }

  if (length(datasources) == 0) {
    stop("No DataSHIELD connections found", call. = FALSE)
  }

  server_names <- names(datasources)

  # Build call
  if (is.null(pattern)) {
    call_expr <- call("validateIdFormatDS", data_name, id_col)
  } else {
    call_expr <- call("validateIdFormatDS", data_name, id_col, pattern)
  }

  # Execute on all servers
  results <- DSI::datashield.aggregate(conns = datasources, expr = call_expr)

  # Compile results
  servers_df <- data.frame(
    server = server_names,
    n_obs = sapply(results, function(x) x$n_obs),
    n_unique = sapply(results, function(x) x$n_unique),
    n_missing = sapply(results, function(x) x$n_missing),
    id_class = sapply(results, function(x) x$id_class),
    format_signature = sapply(results, function(x) x$format_signature),
    stringsAsFactors = FALSE
  )

  if (!is.null(pattern)) {
    servers_df$pct_match <- sapply(results, function(x) x$pct_match)
    servers_df$all_match <- sapply(results, function(x) x$all_match)
  }

  # Check format consistency
  unique_signatures <- unique(servers_df$format_signature)
  format_match <- length(unique_signatures) == 1 && !is.na(unique_signatures[1])

  # Check pattern match if provided
  if (!is.null(pattern)) {
    pattern_match <- all(sapply(results, function(x) x$all_match))
  } else {
    pattern_match <- NA
  }

  # Generate warnings
  warnings <- character()

  if (!format_match) {
    warnings <- c(warnings,
      "Format signatures differ across servers - IDs may have inconsistent formats")
  }

  if (any(servers_df$n_missing > 0)) {
    servers_with_missing <- servers_df$server[servers_df$n_missing > 0]
    warnings <- c(warnings,
      paste("Missing IDs detected on:", paste(servers_with_missing, collapse = ", ")))
  }

  if (any(servers_df$n_unique != servers_df$n_obs - servers_df$n_missing)) {
    warnings <- c(warnings,
      "Duplicate IDs detected on one or more servers")
  }

  unique_classes <- unique(servers_df$id_class)
  if (length(unique_classes) > 1) {
    warnings <- c(warnings,
      paste("ID column has different types across servers:",
            paste(unique_classes, collapse = ", ")))
  }

  if (!is.null(pattern) && !pattern_match) {
    servers_not_matching <- servers_df$server[!servers_df$all_match]
    warnings <- c(warnings,
      paste("Pattern not matched on:", paste(servers_not_matching, collapse = ", ")))
  }

  # Determine overall validity
  valid <- format_match && (is.na(pattern_match) || pattern_match)

  result <- list(
    valid = valid,
    servers = servers_df,
    format_match = format_match,
    pattern_match = pattern_match,
    warnings = warnings
  )

  class(result) <- c("ds.id.validation", "list")
  return(result)
}

#' @title Print Method for ds.id.validation Objects
#' @description Prints a summary of ID format validation results.
#' @param x A ds.id.validation object
#' @param ... Additional arguments (ignored)
#' @export
print.ds.id.validation <- function(x, ...) {
  cat("\nIdentifier Format Validation\n")
  cat("============================\n\n")

  cat("Overall Status:", ifelse(x$valid, "VALID", "ISSUES DETECTED"), "\n")
  cat("Format Consistency:", ifelse(x$format_match, "Yes", "No"), "\n")

  if (!is.na(x$pattern_match)) {
    cat("Pattern Match:", ifelse(x$pattern_match, "Yes", "No"), "\n")
  }

  cat("\nServer Details:\n")
  print(x$servers, row.names = FALSE)

  if (length(x$warnings) > 0) {
    cat("\nWarnings:\n")
    for (w in x$warnings) {
      cat(" -", w, "\n")
    }
  }

  invisible(x)
}
