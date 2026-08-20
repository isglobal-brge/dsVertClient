.dsvert_new_federation <- function(symbol, sites, source_symbols,
                                    id_columns, attestation,
                                    public_schema = NULL) {
  if (is.null(public_schema)) {
    version <- 1L
    value <- list(
      version = version,
      symbol = symbol,
      sites = sites,
      source_symbols = source_symbols,
      id_columns = id_columns,
      attestation = attestation)
  } else {
    version <- 2L
    value <- list(
      version = version,
      symbol = symbol,
      sites = sites,
      source_symbols = source_symbols,
      id_columns = id_columns,
      public_schema = public_schema,
      attestation = attestation)
  }
  structure(
    value,
    class = c("ds.vertFederation", "list"))
}

.dsvert_federation_public_schema_validate <- function(
    schema, sites, id_columns) {
  expected <- c("server", "column", "kind", "role")
  invalid <- function() stop(
    "Invalid ds.vertFederation: malformed public_schema", call. = FALSE)
  if (!is.data.frame(schema) || !identical(names(schema), expected) ||
      nrow(schema) < length(sites) ||
      any(!vapply(schema, is.character, logical(1L))) ||
      anyNA(schema) ||
      any(!schema$server %in% sites) ||
      any(!vapply(schema$column, .dsvert_formula_identifier, logical(1L))) ||
      any(!schema$kind %in% c("identifier", "numeric", "categorical")) ||
      any(!schema$role %in% c("id", "data")) ||
      anyDuplicated(paste(schema$server, schema$column, sep = "\r")) ||
      !identical(unique(schema$server), sites)) {
    invalid()
  }
  for (site in sites) {
    local <- schema[schema$server == site, , drop = FALSE]
    identifier <- local$role == "id"
    if (sum(identifier) != 1L ||
        !identical(local$column[identifier], unname(id_columns[[site]])) ||
        any((local$kind == "identifier") != identifier)) {
      invalid()
    }
  }
  schema
}

.dsvert_federation_public_schema <- function(
    symbol, datasources, id_columns, attestation) {
  # This catalogue is routing metadata, never an authorization decision. Each
  # analysis still verifies its signed server-owned artifact and ownership.
  sites <- names(datasources)
  responses <- .dsvert_aggregate_strict(
    conns = datasources,
    expr = call(name = "dsvertColNamesDS", data_name = symbol),
    operation = "federation public-schema discovery")
  expected <- c(
    "version", "peer_name", "dataset_id", "dataset_version", "columns",
    "kinds", "roles", "data_access")
  if (!is.list(responses) || !identical(names(responses), sites)) {
    stop("Public schema responses do not cover the complete federation.",
         call. = FALSE)
  }
  rows <- lapply(sites, function(site) {
    response <- responses[[site]]
    if (!is.list(response) || !identical(names(response), expected)) {
      stop("A custodian returned an invalid public column catalogue.",
           call. = FALSE)
    }
    valid_columns <- is.character(response$columns) &&
      length(response$columns) >= 1L && !anyNA(response$columns) &&
      !anyDuplicated(response$columns) &&
      identical(response$columns,
                sort(response$columns, method = "radix")) &&
      all(vapply(response$columns, .dsvert_formula_identifier, logical(1L)))
    if (!identical(response$version, "dsvert-public-column-catalog-v1") ||
        !identical(response$peer_name, site) ||
        !identical(response$dataset_id, attestation$dataset_id) ||
        !identical(response$dataset_version, attestation$dataset_version) ||
        !identical(response$data_access, FALSE) || !valid_columns ||
        !is.character(response$kinds) ||
        !identical(names(response$kinds), response$columns) ||
        anyNA(response$kinds) ||
        any(!response$kinds %in% c(
          "identifier", "numeric", "categorical")) ||
        !is.character(response$roles) ||
        !identical(names(response$roles), response$columns) ||
        anyNA(response$roles) ||
        any(!response$roles %in% c("id", "data"))) {
      stop("A custodian returned an invalid public column catalogue.",
           call. = FALSE)
    }
    data.frame(
      server = rep(site, length(response$columns)),
      column = response$columns,
      kind = unname(response$kinds),
      role = unname(response$roles),
      stringsAsFactors = FALSE)
  })
  schema <- do.call(rbind, rows)
  rownames(schema) <- NULL
  .dsvert_federation_public_schema_validate(schema, sites, id_columns)
}

.dsvert_federation_formula_schema <- function(federation) {
  if (!inherits(federation, "ds.vertFederation") ||
      !identical(federation$version, 2L)) {
    stop(
      "This federation predates public schema discovery; align it again.",
      call. = FALSE)
  }
  federation$public_schema
}

.dsvert_federation_status <- function(federation, datasources) {
  invalid <- function(message) stop(
    "Invalid ds.vertFederation: ", message, call. = FALSE)
  version <- if (is.list(federation)) federation$version else NULL
  required <- if (identical(version, 1L)) {
    c("version", "symbol", "sites", "source_symbols", "id_columns",
      "attestation")
  } else if (identical(version, 2L)) {
    c("version", "symbol", "sites", "source_symbols", "id_columns",
      "public_schema", "attestation")
  } else {
    character()
  }
  if (!inherits(federation, "ds.vertFederation") ||
      !identical(names(federation), required) ||
      !is.character(federation$symbol) || length(federation$symbol) != 1L ||
      is.na(federation$symbol) || !nzchar(federation$symbol) ||
      !is.character(federation$sites) || length(federation$sites) < 2L ||
      anyNA(federation$sites) || any(!nzchar(federation$sites)) ||
      anyDuplicated(federation$sites)) {
    invalid("malformed object")
  }
  sites <- federation$sites
  for (field in c("source_symbols", "id_columns")) {
    value <- federation[[field]]
    if (!is.character(value) || !identical(names(value), sites) ||
        length(value) != length(sites) || anyNA(value) ||
        any(!nzchar(value))) {
      invalid(paste0("malformed ", field))
    }
  }
  if (identical(version, 2L)) {
    .dsvert_federation_public_schema_validate(
      federation$public_schema, sites, federation$id_columns)
  }
  attestation <- tryCatch(
    .dsvert_validate_psi_padded_attestation(federation$attestation),
    error = function(e) NULL)
  if (is.null(attestation) || attestation$peer_count != length(sites)) {
    invalid("malformed attestation")
  }
  datasource_sites <- names(datasources)
  if (!is.list(datasources) || is.null(datasource_sites) ||
      anyNA(datasource_sites) || any(!nzchar(datasource_sites)) ||
      anyDuplicated(datasource_sites) ||
      !setequal(datasource_sites, sites)) {
    stop(
      "A ds.vertFederation requires the same named datasource set used for alignment.",
      call. = FALSE)
  }
  datasources <- datasources[sites]
  status <- .psi_alignment_status(federation$symbol, datasources)
  current <- if (is.list(status$manifests) && length(status$manifests)) {
    status$manifests[[1L]]
  } else {
    NULL
  }
  if (!isTRUE(status$aligned) || !identical(status$n_common, NA_integer_) ||
      !identical(current, attestation)) {
    stop(
      "The ds.vertFederation symbol is no longer PSI-aligned to its original attestation.",
      call. = FALSE)
  }
  list(symbol = federation$symbol, datasources = datasources)
}

.dsvert_resolve_federation <- function(value, datasources) {
  if (!inherits(value, "ds.vertFederation")) {
    return(list(value = value, datasources = datasources))
  }
  resolved <- .dsvert_federation_status(value, datasources)
  list(value = resolved$symbol, datasources = resolved$datasources)
}

.dsvert_federation_argument <- function(value, datasources) {
  if (!inherits(value, "ds.vertFederation")) {
    return(list(value = value, datasources = datasources))
  }
  datasources <- .dsvert_dp_datasources(datasources)
  .dsvert_resolve_federation(value, datasources)
}

#' @export
print.ds.vertFederation <- function(x, ...) {
  cat("Vertically aligned dsVert federation\n")
  cat("  Symbol:", x$symbol, "\n")
  cat("  Sites:", paste(x$sites, collapse = ", "), "\n")
  if (identical(x$version, 2L) && is.data.frame(x$public_schema)) {
    cat("  Public columns:", nrow(x$public_schema), "\n")
  }
  cat("  PSI attested:", isTRUE(x$attestation$alignment_attested), "\n")
  invisible(x)
}
