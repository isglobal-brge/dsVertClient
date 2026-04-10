#' @title Query Server Identity Public Keys
#' @description Queries all connected DataSHIELD servers for their Ed25519
#'   identity public keys. Used by administrators to discover PKs for
#'   configuring the \code{dsvert.trusted_peers} option on each server.
#'
#' @param datasources DataSHIELD connections. If NULL, uses all available.
#' @return Named list: server_name -> identity_pk (base64url string).
#'
#' @examples
#' \dontrun{
#' pks <- ds.getIdentityPks(datasources = conns)
#' # Set trusted peers on each server via opalr:
#' # dsadmin.set_option(opal, "dsvert.trusted_peers", paste(pks, collapse=","))
#' }
#'
#' @importFrom DSI datashield.aggregate datashield.connections_find
#' @export
ds.getIdentityPks <- function(datasources = NULL) {
  if (is.null(datasources))
    datasources <- DSI::datashield.connections_find()
  results <- DSI::datashield.aggregate(
    conns = datasources,
    expr = call("dsvertIdentityPkDS")
  )
  pks <- list()
  for (name in names(results))
    pks[[name]] <- results[[name]]$identity_pk
  pks
}
