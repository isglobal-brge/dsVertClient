#' @title Query Server Identity Public Keys
#' @description Queries all connected DataSHIELD servers for their Ed25519
#'   identity public keys. Used by administrators to discover PKs for
#'   configuring the name-bound \code{dsvert.trusted_peers} option on each
#'   server. A new or changed key must be verified independently at the
#'   corresponding server console before updating any pin. The analyst/relay
#'   is not a trusted source for accepting a replacement identity.
#'
#' @param datasources DataSHIELD connections. If NULL, uses all available.
#' @return Named list: server_name -> identity_pk (base64url string).
#'
#' @examples
#' \dontrun{
#' pks <- ds.getIdentityPks(datasources = conns)
#' # After out-of-band verification, configure every server with a named map
#' # of the *other* logical peers. For example, on site_a:
#' site_a_pins <- unlist(pks[setdiff(names(pks), "site_a")])
#' names(site_a_pins) <- setdiff(names(pks), "site_a")
#' # Persist site_a_pins as the server option dsvert.trusted_peers, restart
#' # the affected DataSHIELD service, and repeat for every participating site.
#' }
#'
#' @importFrom DSI datashield.aggregate datashield.connections_find
#' @export
ds.getIdentityPks <- function(datasources = NULL) {
  if (is.null(datasources))
    datasources <- DSI::datashield.connections_find()
  results <- .dsvert_aggregate_strict(
    conns = datasources,
    expr = call(name = "dsvertIdentityPkDS"),
    operation = "server identity discovery"
  )
  pks <- lapply(names(results), function(server) {
    value <- results[[server]]
    if (!is.list(value) || !identical(names(value), "identity_pk")) {
      stop("Server '", server, "' returned an invalid public identity",
           call. = FALSE)
    }
    normalized <- .dsvert_dp_normalize_identity_pk(value$identity_pk)
    if (is.null(normalized)) {
      stop("Server '", server, "' returned an invalid public identity",
           call. = FALSE)
    }
    normalized
  })
  stats::setNames(pks, names(results))
}
