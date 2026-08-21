# Normalize a local routing alias for every named DataSHIELD connection.
#
# These names are not part of the common PSI semantics: each server resolves
# its own configured symbol and then signs the same logical contract.  Keeping
# the map complete and ordered prevents positional routing across sites.
.dsvert_site_character <- function(value, datasources, what) {
  server_names <- names(datasources)
  if (!is.list(datasources) || !length(datasources) ||
      is.null(server_names) || anyNA(server_names) ||
      any(!nzchar(server_names)) || anyDuplicated(server_names)) {
    stop(what, " requires uniquely named DataSHIELD connections", call. = FALSE)
  }
  scalar <- is.character(value) && length(value) == 1L &&
    is.null(names(value)) && !is.na(value) && nzchar(value)
  if (isTRUE(scalar)) {
    return(stats::setNames(rep.int(value, length(server_names)), server_names))
  }
  if (!(is.character(value) || is.list(value)) || is.null(names(value)) ||
      length(value) != length(server_names) || anyNA(names(value)) ||
      any(!nzchar(names(value))) || anyDuplicated(names(value)) ||
      !setequal(names(value), server_names)) {
    stop(what, " must be a complete named per-site map", call. = FALSE)
  }
  result <- vapply(value[server_names], function(site_value) {
    if (!is.character(site_value) || length(site_value) != 1L ||
        is.na(site_value) || !nzchar(site_value) || !is.null(names(site_value))) {
      stop(what, " must contain one non-empty character alias per site",
           call. = FALSE)
    }
    site_value
  }, character(1L))
  stats::setNames(unname(result), server_names)
}
