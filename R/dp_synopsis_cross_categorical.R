# Canonical signed projection for one cross-owner categorical pair.  The full
# bootstrap catalog authorizes the projection but is deliberately absent from
# its public manifest identity.

.dsvert_dp_synopsis_pair_reference_v1 <- function(value, what) {
  if (!.dsvert_dp_is_string(value) ||
      nchar(value, type = "bytes") > 257L) {
    stop("Invalid ", what, call. = FALSE)
  }
  expression <- tryCatch(
    parse(text = value, keep.source = FALSE),
    error = function(error) NULL)
  if (is.null(expression) || length(expression) != 1L) {
    stop(what, " must be a plain column or an exact server$column reference",
         call. = FALSE)
  }
  parsed <- .dsvert_formula_reference(
    expression[[1L]], qualified = TRUE, what = what)
  list(
    reference = parsed$reference, owner_peer = parsed$server,
    column = parsed$column)
}

.dsvert_dp_synopsis_pair_federation_references_v1 <- function(
    federation, values) {
  if (!inherits(federation, "ds.vertFederation") ||
      !identical(federation$version, 2L)) {
    stop("A version-2 ds.vertFederation is required", call. = FALSE)
  }
  schema <- .dsvert_federation_public_schema_validate(
    .dsvert_federation_formula_schema(federation), federation$sites,
    federation$id_columns)
  if (!is.character(values) || length(values) != 2L || anyNA(values) ||
      any(!nzchar(values))) {
    stop("Two non-empty contingency variables are required", call. = FALSE)
  }
  resolve <- function(value) {
    reference <- .dsvert_dp_synopsis_pair_reference_v1(
      value, "A contingency variable")
    matches <- schema$column == reference$column
    if (!is.null(reference$owner_peer)) {
      if (!reference$owner_peer %in% schema$server) {
        stop("Unknown contingency server: ", reference$owner_peer,
             call. = FALSE)
      }
      matches <- matches & schema$server == reference$owner_peer
      if (!any(matches)) {
        if (reference$column %in% schema$column) {
          stop(
            "Column ", reference$column, " is not owned by server ",
            reference$owner_peer, call. = FALSE)
        }
        stop("Contingency variable is missing from public schema: ",
             reference$reference, call. = FALSE)
      }
    } else if (!any(matches)) {
      stop("Contingency variable is missing from public schema: ",
           reference$column, call. = FALSE)
    } else if (sum(matches) != 1L) {
      stop(
        "Contingency variable has multiple owners; qualify it as ",
        "server$column: ", reference$column, call. = FALSE)
    }
    row <- schema[which(matches), , drop = FALSE]
    if (!identical(row$role[[1L]], "data")) {
      stop("Identifier columns cannot be used in a contingency table: ",
           reference$reference, call. = FALSE)
    }
    if (!identical(row$kind[[1L]], "categorical")) {
      stop(
        "Contingency variable has incompatible public kind '",
        row$kind[[1L]], "': ", reference$reference,
        "; expected categorical", call. = FALSE)
    }
    list(
      reference = paste0(row$server[[1L]], "$", row$column[[1L]]),
      owner_peer = row$server[[1L]], column = row$column[[1L]])
  }
  resolved <- lapply(values, resolve)
  physical <- vapply(resolved, function(reference) paste(
    reference$owner_peer, reference$column, sep = "\r"), character(1L))
  if (anyDuplicated(physical)) {
    stop("row_var and col_var must identify two distinct variables",
         call. = FALSE)
  }
  resolved
}

.dsvert_dp_synopsis_cross_pair_request_v1 <- function(value) {
  fields <- c("version", "family", "dataset", "columns")
  if (!is.list(value) || !.dsvert_dp_has_exact_names(value, fields) ||
      !identical(
        value$version, .DSVERT_CLIENT_SYNOPSIS_CROSS_PAIR_REQUEST_VERSION) ||
      !identical(value$family, "categorical_pair")) {
    stop("Invalid cross-owner categorical Synopsis projection request",
         call. = FALSE)
  }
  columns <- value$columns
  if (is.list(columns) && is.null(names(columns)) &&
      all(vapply(columns, .dsvert_dp_is_string, logical(1L)))) {
    columns <- unname(unlist(columns, use.names = FALSE))
  }
  if (!is.character(columns) || length(columns) != 2L || anyNA(columns) ||
      !is.null(names(columns))) {
    stop("Invalid cross-owner categorical Synopsis projection columns",
         call. = FALSE)
  }
  references <- lapply(columns, .dsvert_dp_synopsis_pair_reference_v1,
                       what = "A cross-owner categorical Synopsis column")
  canonical <- vapply(references, `[[`, character(1L), "reference")
  if (anyDuplicated(canonical) ||
      (all(vapply(references, function(reference) {
        is.null(reference$owner_peer)
      }, logical(1L))) &&
       identical(references[[1L]]$column, references[[2L]]$column))) {
    stop("A cross-owner categorical Synopsis projection needs two columns",
         call. = FALSE)
  }
  columns <- sort(unname(canonical), method = "radix")
  list(
    version = .DSVERT_CLIENT_SYNOPSIS_CROSS_PAIR_REQUEST_VERSION,
    family = "categorical_pair",
    dataset = .dsvert_dp_synopsis_local_pair_id_v1(
      value$dataset, "cross-owner categorical Synopsis dataset"),
    columns = unname(columns))
}

.dsvert_dp_synopsis_cross_pair_selector_hash_v1 <- function(value) {
  .dsvert_dp_synopsis_client_hash(
    .DSVERT_CLIENT_SYNOPSIS_CROSS_PAIR_SELECTOR_DOMAIN, value)
}

.dsvert_dp_synopsis_cross_pair_public_id_v1 <- function(columns) {
  semantic <- lapply(columns, function(column) {
    column[c("owner_peer", "dataset", "column")]
  })
  keys <- vapply(semantic, function(column) paste(
    column$owner_peer, column$dataset, column$column, sep = "\r"),
    character(1L))
  semantic <- semantic[order(keys, method = "radix")]
  paste0("categorical_pair_", substr(.dsvert_vector_hash(list(
    version = "dsvert-synopsis-cross-categorical-semantic-id-v1",
    family = "categorical_pair", columns = semantic)), 1L, 48L))
}

.dsvert_dp_synopsis_cross_pair_spec_v1 <- function(schema, request) {
  vertical <- tryCatch(
    schema$workload_contract$value$vertical_cross,
    error = function(error) NULL)
  datasets <- if (is.list(schema) && is.list(schema$unsigned)) {
    schema$unsigned$datasets
  } else NULL
  owners <- if (is.list(datasets)) unique(unlist(lapply(
    datasets, function(dataset) names(dataset$patient_keys)),
    use.names = FALSE)) else character()
  if (!is.list(vertical) || !length(vertical) || !is.list(datasets)) {
    stop("The signed cross-owner categorical pair is missing or ambiguous",
         call. = FALSE)
  }
  normalized <- lapply(names(vertical), function(analysis_id) {
    entry <- vertical[[analysis_id]]
    raw <- if (is.list(entry)) entry$spec else NULL
    fields <- c(
      "version", "left_dataset", "right_dataset", "left", "right",
      "family")
    valid_entry <- .dsvert_dp_has_exact_names(
      entry, c("owner_peer", "spec")) &&
      .dsvert_dp_is_string(entry$owner_peer) &&
      entry$owner_peer %in% owners &&
      .dsvert_dp_has_exact_names(raw, fields) &&
      identical(raw$version, "v2") &&
      identical(raw$family, "categorical_pair")
    if (!isTRUE(valid_entry)) {
      stop("The signed cross-owner categorical workload is invalid",
           call. = FALSE)
    }
    left_dataset <- .dsvert_dp_synopsis_local_pair_id_v1(
      raw$left_dataset, "signed cross-owner categorical dataset")
    right_dataset <- .dsvert_dp_synopsis_local_pair_id_v1(
      raw$right_dataset, "signed cross-owner categorical dataset")
    left <- .dsvert_dp_synopsis_pair_reference_v1(
      raw$left, "A signed cross-owner categorical column")$reference
    right <- .dsvert_dp_synopsis_pair_reference_v1(
      raw$right, "A signed cross-owner categorical column")$reference
    left_column <- tryCatch(datasets[[left_dataset]]$columns[[left]],
                            error = function(error) NULL)
    right_column <- tryCatch(datasets[[right_dataset]]$columns[[right]],
                             error = function(error) NULL)
    valid_columns <- is.list(left_column) && is.list(right_column) &&
      identical(left_column$kind, "categorical") &&
      identical(right_column$kind, "categorical") &&
      !identical(left_column$owner_peer, right_column$owner_peer) &&
      identical(datasets[[left_dataset]]$alignment_group,
                datasets[[right_dataset]]$alignment_group) &&
      left_column$owner_peer %in%
        names(datasets[[left_dataset]]$patient_keys) &&
      right_column$owner_peer %in%
        names(datasets[[right_dataset]]$patient_keys)
    if (!isTRUE(valid_columns)) {
      stop("The signed cross-owner categorical workload is invalid",
           call. = FALSE)
    }
    left_physical <- .dsvert_dp_synopsis_pair_reference_v1(
      left, "A signed cross-owner categorical column")$column
    right_physical <- .dsvert_dp_synopsis_pair_reference_v1(
      right, "A signed cross-owner categorical column")$column
    references <- c(
      paste(left_column$owner_peer, left_dataset, left_physical, sep = "::"),
      paste(right_column$owner_peer, right_dataset, right_physical,
            sep = "::"))
    order_index <- order(references, method = "radix")
    side_columns <- c(left, right)[order_index]
    side_datasets <- c(left_dataset, right_dataset)[order_index]
    side_owners <- c(
      left_column$owner_peer, right_column$owner_peer)[order_index]
    list(
      analysis_id = analysis_id, owner_peer = entry$owner_peer,
      spec = list(
        version = "v2", family = "categorical_pair",
        left_dataset = side_datasets[[1L]], left = side_columns[[1L]],
        right_dataset = side_datasets[[2L]], right = side_columns[[2L]],
        left_owner = side_owners[[1L]], right_owner = side_owners[[2L]]))
  })
  requested <- lapply(request$columns, .dsvert_dp_synopsis_pair_reference_v1,
                      what = "A cross-owner categorical Synopsis column")
  schema_columns <- unlist(lapply(names(datasets), function(dataset_name) {
    dataset <- datasets[[dataset_name]]
    lapply(names(dataset$columns), function(reference) {
      column <- dataset$columns[[reference]]
      parsed <- .dsvert_dp_synopsis_pair_reference_v1(
        reference, "A signed categorical schema column")
      list(
        dataset = dataset_name, owner_peer = column$owner_peer,
        column = parsed$column)
    })
  }), recursive = FALSE)
  bare_owner_unique <- function(reference) {
    if (!is.null(reference$owner_peer)) return(TRUE)
    matches <- schema_columns[vapply(schema_columns, function(column) {
      identical(column$column, reference$column)
    }, logical(1L))]
    owners <- unique(vapply(
      matches, `[[`, character(1L), "owner_peer"))
    length(owners) == 1L
  }
  if (!all(vapply(requested, bare_owner_unique, logical(1L)))) {
    stop("A bare cross-owner categorical column must identify one signed owner",
         call. = FALSE)
  }
  matches <- vapply(normalized, function(candidate) {
    spec <- candidate$spec
    sides <- lapply(c("left", "right"), function(side) {
      reference <- spec[[side]]
      column <- datasets[[spec[[paste0(side, "_dataset")]]]]$
        columns[[reference]]
      list(
        dataset = spec[[paste0(side, "_dataset")]],
        owner_peer = column$owner_peer,
        column = .dsvert_dp_synopsis_pair_reference_v1(
          reference, "A signed categorical schema column")$column)
    })
    accepts <- function(reference, side) {
      identical(reference$column, side$column) &&
        (is.null(reference$owner_peer) ||
           identical(reference$owner_peer, side$owner_peer))
    }
    bindings <- c(
      accepts(requested[[1L]], sides[[1L]]) &&
        accepts(requested[[2L]], sides[[2L]]),
      accepts(requested[[1L]], sides[[2L]]) &&
        accepts(requested[[2L]], sides[[1L]]))
    request$dataset %in% vapply(
      sides, `[[`, character(1L), "dataset") && sum(bindings) == 1L
  }, logical(1L))
  selected <- normalized[matches]
  if (length(selected) != 1L) {
    stop("The signed cross-owner categorical pair is missing or ambiguous",
         call. = FALSE)
  }
  selected[[1L]]
}

.dsvert_dp_synopsis_cross_pair_project_components_v1 <- function(
    schema, selected, context) {
  spec <- selected$spec
  sides <- list(
    left = list(dataset = spec$left_dataset, reference = spec$left,
                owner_peer = spec$left_owner),
    right = list(dataset = spec$right_dataset, reference = spec$right,
                 owner_peer = spec$right_owner))
  datasets <- list()
  projected_references <- list()
  for (side_name in names(sides)) {
    side <- sides[[side_name]]
    source <- schema$unsigned$datasets[[side$dataset]]
    column <- if (is.list(source)) source$columns[[side$reference]] else NULL
    if (!is.list(source) || !is.list(column) ||
        !identical(column$kind, "categorical") ||
        !identical(column$owner_peer, side$owner_peer) ||
        !is.atomic(column$levels) || !length(column$levels) ||
        !side$owner_peer %in% names(source$patient_keys)) {
      stop("Invalid cross-owner categorical Synopsis projection source",
           call. = FALSE)
    }
    if (is.null(datasets[[side$dataset]])) {
      datasets[[side$dataset]] <- source[c(
        "dataset_id", "dataset_version", "schema_version",
        "alignment_group")]
      datasets[[side$dataset]]$patient_keys <- list()
      datasets[[side$dataset]]$columns <- list()
    }
    datasets[[side$dataset]]$patient_keys[[side$owner_peer]] <-
      source$patient_keys[[side$owner_peer]]
    column$levels <- .dsvert_dp_capsule_manifest_string_array(
      as.list(column$levels), "cross categorical projection levels")
    physical <- .dsvert_dp_synopsis_pair_reference_v1(
      side$reference, "A cross-owner categorical projection column")$column
    projected_reference <- paste0(side$owner_peer, "$", physical)
    projected_references[[side_name]] <- projected_reference
    datasets[[side$dataset]]$columns[[projected_reference]] <-
      .dsvert_joint_dp_client_canonical(column)
  }
  datasets <- datasets[order(names(datasets), method = "radix")]
  datasets <- lapply(datasets, function(dataset) {
    dataset$patient_keys <- dataset$patient_keys[
      order(names(dataset$patient_keys), method = "radix")]
    dataset$columns <- dataset$columns[
      order(names(dataset$columns), method = "radix")]
    dataset
  })
  semantic_columns <- lapply(names(sides), function(side_name) {
    side <- sides[[side_name]]
    list(
      owner_peer = side$owner_peer, dataset = side$dataset,
      column = .dsvert_dp_synopsis_pair_reference_v1(
        side$reference, "A cross-owner categorical projection column")$column)
  })
  public_analysis_id <-
    .dsvert_dp_synopsis_cross_pair_public_id_v1(semantic_columns)
  public_owner <- sort(vapply(
    semantic_columns, `[[`, character(1L), "owner_peer"),
    method = "radix")[[1L]]
  raw_spec <- list(
    version = "v2", left_dataset = spec$left_dataset,
    right_dataset = spec$right_dataset,
    left = projected_references$left,
    right = projected_references$right, family = "categorical_pair")
  workload <- .dsvert_joint_dp_client_canonical(list(
    version = .DSVERT_CLIENT_DP_CAPSULE_WORKLOAD_CONTRACT_VERSION,
    describe = list(), survival = list(), gaussian = list(),
    vertical_cross = stats::setNames(list(list(
      owner_peer = public_owner, spec = raw_spec)),
      public_analysis_id)))
  parent_snapshot <- schema$unsigned$logical_snapshot
  alignment <- tryCatch(
    as.integer(parent_snapshot$alignment_protocol_version),
    error = function(error) NA_integer_)
  if (!is.list(parent_snapshot) || !identical(
        parent_snapshot$logical_snapshot_id, context$policy$cohort_id) ||
      length(alignment) != 1L || is.na(alignment) || alignment < 1L ||
      !identical(as.numeric(alignment),
                 as.numeric(parent_snapshot$alignment_protocol_version))) {
    stop("Invalid cross-owner categorical Synopsis projection snapshot",
         call. = FALSE)
  }
  fingerprint <- .dsvert_dp_capsule_manifest_hash(list(
    protocol = "dsvert-biomedical-capsule-logical-snapshot-v1",
    domain = context$policy$domain,
    cohort_id = context$policy$cohort_id,
    peer_pinset_sha256 = context$policy$peer_pinset_sha256,
    alignment_protocol_version = alignment,
    datasets = datasets, workload_contract = workload))
  logical_snapshot <- .dsvert_joint_dp_client_canonical(list(
    logical_snapshot_id = context$policy$cohort_id,
    version = paste0("schema-v1-", fingerprint),
    alignment_protocol_version = alignment))
  projected_schema <- .dsvert_joint_dp_client_canonical(list(
    version = .DSVERT_CLIENT_DP_CAPSULE_SCHEMA_VERSION,
    logical_snapshot = logical_snapshot,
    peer_pinset_sha256 = context$policy$peer_pinset_sha256,
    datasets = datasets))
  scope <- .dsvert_joint_dp_client_canonical(list(
    mode = "catalog_v1", numeric_moments = character(),
    categorical_marginals = character(),
    strict_missing_categorical = character(),
    categorical_pairs = list(),
    correlations = list()))
  projected_policy <- context$policy
  projected_policy$primitive_scope <- scope
  schema_json <- .dsvert_joint_dp_client_json(projected_schema)
  workload_json <- .dsvert_joint_dp_client_json(workload)
  list(
    unsigned = projected_schema, json = schema_json,
    sha256 = digest::digest(
      schema_json, algo = "sha256", serialize = FALSE),
    logical_snapshot = logical_snapshot,
    workload_contract = list(
      value = workload, json = workload_json,
      sha256 = digest::digest(
        workload_json, algo = "sha256", serialize = FALSE)),
    primitive_scope = scope,
    policy_sha256 = .dsvert_vector_hash(projected_policy))
}

.dsvert_dp_synopsis_cross_pair_resolve_v1 <- function(
    request, schema, context) {
  request <- .dsvert_dp_synopsis_cross_pair_request_v1(request)
  selected <- .dsvert_dp_synopsis_cross_pair_spec_v1(schema, request)
  projection <- .dsvert_dp_synopsis_cross_pair_project_components_v1(
    schema, selected, context)
  spec <- selected$spec
  columns <- lapply(c("left", "right"), function(side) {
    dataset <- spec[[paste0(side, "_dataset")]]
    reference <- spec[[side]]
    column <- schema$unsigned$datasets[[dataset]]$columns[[reference]]
    levels <- .dsvert_dp_capsule_manifest_string_array(
      as.list(column$levels), "cross categorical levels")
    list(
      side = side, dataset = dataset,
      reference = paste0(
        column$owner_peer, "$",
        .dsvert_dp_synopsis_pair_reference_v1(
          reference, "A signed categorical schema column")$column),
      column = .dsvert_dp_synopsis_pair_reference_v1(
        reference, "A signed categorical schema column")$column,
      owner_peer = column$owner_peer,
      levels_sha256 = .dsvert_vector_hash(as.list(levels)))
  })
  parent <- list(
    schema_sha256 = projection$sha256,
    workload_contract_sha256 = projection$workload_contract$sha256,
    logical_snapshot = projection$logical_snapshot,
    policy_sha256 = projection$policy_sha256)
  unsigned <- .dsvert_joint_dp_client_canonical(list(
    version = .DSVERT_CLIENT_SYNOPSIS_CROSS_PAIR_SELECTOR_VERSION,
    family = "categorical_pair",
    analysis_id = .dsvert_dp_synopsis_cross_pair_public_id_v1(columns),
    columns = columns, parent = parent))
  .dsvert_joint_dp_client_canonical(c(unsigned, list(
    sha256 = .dsvert_dp_synopsis_cross_pair_selector_hash_v1(unsigned))))
}

.dsvert_dp_synopsis_cross_pair_project_v1 <- function(
    schema, selector, context) {
  physical <- tryCatch(vapply(
    selector$columns, `[[`, character(1L), "column"),
    error = function(error) character())
  datasets <- tryCatch(vapply(
    selector$columns, `[[`, character(1L), "dataset"),
    error = function(error) character())
  expected <- if (is.list(selector) && identical(
      selector$version,
      .DSVERT_CLIENT_SYNOPSIS_CROSS_PAIR_SELECTOR_VERSION) &&
      length(physical) == 2L && length(datasets) == 2L) tryCatch(
    .dsvert_dp_synopsis_cross_pair_resolve_v1(list(
      version = .DSVERT_CLIENT_SYNOPSIS_CROSS_PAIR_REQUEST_VERSION,
      family = "categorical_pair", dataset = datasets[[1L]],
      columns = vapply(
        selector$columns, `[[`, character(1L), "reference")),
      schema, context),
    error = function(error) NULL) else NULL
  if (is.null(expected) || !identical(
        .dsvert_joint_dp_client_json(selector),
        .dsvert_joint_dp_client_json(expected))) {
    stop("The cross-owner categorical Synopsis selector is detached from signed metadata",
         call. = FALSE)
  }
  selected <- .dsvert_dp_synopsis_cross_pair_spec_v1(
    schema, .dsvert_dp_synopsis_cross_pair_request_v1(list(
      version = .DSVERT_CLIENT_SYNOPSIS_CROSS_PAIR_REQUEST_VERSION,
      family = "categorical_pair", dataset = datasets[[1L]],
      columns = vapply(
        selector$columns, `[[`, character(1L), "reference"))))
  projection <- .dsvert_dp_synopsis_cross_pair_project_components_v1(
    schema, selected, context)
  projection$primitive_scope <- NULL
  projection$policy_sha256 <- NULL
  projection
}

.dsvert_dp_synopsis_projection_is_cross_v1 <- function(value) {
  is.list(value) && identical(
    value$version, .DSVERT_CLIENT_SYNOPSIS_CROSS_PAIR_SELECTOR_VERSION)
}

.dsvert_dp_synopsis_categorical_pair_request_v1 <- function(value) {
  fields <- c("version", "family", "dataset", "columns", "owner_peer")
  if (!is.list(value) || !.dsvert_dp_has_exact_names(value, fields) ||
      !identical(
        value$version,
        .DSVERT_CLIENT_SYNOPSIS_CATEGORICAL_PAIR_REQUEST_VERSION) ||
      !identical(value$family, "categorical_pair")) {
    stop("Invalid categorical Synopsis projection request", call. = FALSE)
  }
  normalized <- .dsvert_dp_synopsis_cross_pair_request_v1(list(
    version = .DSVERT_CLIENT_SYNOPSIS_CROSS_PAIR_REQUEST_VERSION,
    family = "categorical_pair", dataset = value$dataset,
    columns = value$columns))
  owner <- value$owner_peer
  if (!is.null(owner)) owner <- .dsvert_dp_synopsis_local_pair_id_v1(
    owner, "categorical Synopsis owner")
  list(
    version = .DSVERT_CLIENT_SYNOPSIS_CATEGORICAL_PAIR_REQUEST_VERSION,
    family = "categorical_pair", dataset = normalized$dataset,
    columns = normalized$columns, owner_peer = owner)
}

.dsvert_dp_synopsis_categorical_pair_resolve_v1 <- function(
    request, schema, context) {
  request <- .dsvert_dp_synopsis_categorical_pair_request_v1(request)
  candidates <- list()
  references <- lapply(
    request$columns, .dsvert_dp_synopsis_pair_reference_v1,
    what = "A categorical Synopsis column")
  datasets <- tryCatch(schema$unsigned$datasets,
                       error = function(error) NULL)
  signed_columns <- if (is.list(datasets)) unlist(lapply(
    names(datasets), function(dataset_name) {
      dataset <- datasets[[dataset_name]]
      lapply(names(dataset$columns), function(reference) {
        column <- dataset$columns[[reference]]
        list(
          owner_peer = column$owner_peer,
          column = .dsvert_dp_synopsis_pair_reference_v1(
            reference, "A signed categorical schema column")$column)
      })
    }), recursive = FALSE) else list()
  for (reference in references) {
    if (is.null(reference$owner_peer)) {
      matches <- signed_columns[vapply(signed_columns, function(column) {
        identical(column$column, reference$column)
      }, logical(1L))]
      owners <- unique(vapply(
        matches, `[[`, character(1L), "owner_peer"))
      if (length(owners) != 1L) {
        stop("A bare categorical Synopsis column is missing or ambiguous across signed owners",
             call. = FALSE)
      }
    }
  }
  qualified_owners <- unique(Filter(Negate(is.null), lapply(
    references, `[[`, "owner_peer")))
  local_owner <- request$owner_peer
  if (is.null(local_owner) && length(qualified_owners) == 1L) {
    local_owner <- qualified_owners[[1L]]
  }
  local_allowed <- length(qualified_owners) <= 1L && all(vapply(
    references, function(reference) {
      is.null(reference$owner_peer) || is.null(local_owner) ||
        identical(reference$owner_peer, local_owner)
    }, logical(1L)))
  local <- if (isTRUE(local_allowed)) tryCatch(
    .dsvert_dp_synopsis_local_pair_resolve_v1(list(
      version = .DSVERT_CLIENT_SYNOPSIS_LOCAL_PAIR_REQUEST_VERSION,
      family = "categorical_pairs", dataset = request$dataset,
      columns = vapply(references, `[[`, character(1L), "column"),
      owner_peer = local_owner), schema, context),
    error = function(error) NULL) else NULL
  if (!is.null(local)) candidates$local <- local
  if (is.null(request$owner_peer)) {
    cross <- tryCatch(.dsvert_dp_synopsis_cross_pair_resolve_v1(list(
      version = .DSVERT_CLIENT_SYNOPSIS_CROSS_PAIR_REQUEST_VERSION,
      family = "categorical_pair", dataset = request$dataset,
      columns = request$columns), schema, context),
      error = function(error) NULL)
    if (!is.null(cross)) candidates$cross <- cross
  }
  if (length(candidates) != 1L) {
    stop("The signed categorical Synopsis pair is missing or ambiguous",
         call. = FALSE)
  }
  candidates[[1L]]
}

.dsvert_dp_synopsis_projection_resolve_v1 <- function(
    request, schema, context) {
  if (is.list(request) && identical(
        request$version,
        .DSVERT_CLIENT_SYNOPSIS_CATEGORICAL_PAIR_REQUEST_VERSION)) {
    return(.dsvert_dp_synopsis_categorical_pair_resolve_v1(
      request, schema, context))
  }
  if (is.list(request) && identical(
        request$version,
        .DSVERT_CLIENT_SYNOPSIS_CROSS_PAIR_REQUEST_VERSION)) {
    return(.dsvert_dp_synopsis_cross_pair_resolve_v1(
      request, schema, context))
  }
  .dsvert_dp_synopsis_local_pair_resolve_v1(request, schema, context)
}

.dsvert_dp_synopsis_projection_request_v1 <- function(value) {
  if (is.list(value) && identical(
        value$version,
        .DSVERT_CLIENT_SYNOPSIS_CATEGORICAL_PAIR_REQUEST_VERSION)) {
    return(.dsvert_dp_synopsis_categorical_pair_request_v1(value))
  }
  if (is.list(value) && identical(
        value$version,
        .DSVERT_CLIENT_SYNOPSIS_CROSS_PAIR_REQUEST_VERSION)) {
    return(.dsvert_dp_synopsis_cross_pair_request_v1(value))
  }
  .dsvert_dp_synopsis_local_pair_request_v1(value)
}

.dsvert_dp_synopsis_projection_request_compatible_v1 <- function(
    selector, request) {
  request <- .dsvert_dp_synopsis_projection_request_v1(request)
  sides <- tryCatch(lapply(selector$columns, function(column) list(
    dataset = column$dataset %||% selector$dataset,
    owner_peer = column$owner_peer %||% selector$owner_peer,
    column = column$column)), error = function(error) list())
  references <- tryCatch(lapply(
    request$columns, .dsvert_dp_synopsis_pair_reference_v1,
    what = "A categorical Synopsis column"),
    error = function(error) list())
  accepts <- function(reference, side) {
    identical(reference$column, side$column) &&
      (is.null(reference$owner_peer) ||
         identical(reference$owner_peer, side$owner_peer))
  }
  bindings <- if (length(sides) == 2L && length(references) == 2L) c(
    accepts(references[[1L]], sides[[1L]]) &&
      accepts(references[[2L]], sides[[2L]]),
    accepts(references[[1L]], sides[[2L]]) &&
      accepts(references[[2L]], sides[[1L]])) else FALSE
  tuple_match <- sum(bindings) == 1L
  if (.dsvert_dp_synopsis_projection_is_cross_v1(selector)) {
    cross_request <- request$version %in% c(
      .DSVERT_CLIENT_SYNOPSIS_CROSS_PAIR_REQUEST_VERSION,
      .DSVERT_CLIENT_SYNOPSIS_CATEGORICAL_PAIR_REQUEST_VERSION)
    return(isTRUE(cross_request) && is.null(request$owner_peer) &&
      length(sides) == 2L &&
      request$dataset %in% vapply(
        sides, `[[`, character(1L), "dataset") && tuple_match)
  }
  local_request <- request$version %in% c(
    .DSVERT_CLIENT_SYNOPSIS_LOCAL_PAIR_REQUEST_VERSION,
    .DSVERT_CLIENT_SYNOPSIS_CATEGORICAL_PAIR_REQUEST_VERSION)
  is.list(selector) && length(sides) == 2L &&
    isTRUE(local_request) &&
    identical(selector$dataset, request$dataset) &&
    tuple_match &&
    (is.null(request$owner_peer) ||
       identical(selector$owner_peer, request$owner_peer))
}

.dsvert_dp_synopsis_projection_project_v1 <- function(
    schema, selector, context) {
  if (.dsvert_dp_synopsis_projection_is_cross_v1(selector)) {
    return(.dsvert_dp_synopsis_cross_pair_project_v1(
      schema, selector, context))
  }
  .dsvert_dp_synopsis_local_pair_project_v1(schema, selector, context)
}

.dsvert_dp_synopsis_cross_pair_manifest_v1 <- function(manifest, selector) {
  artifacts <- tryCatch(
    manifest$workload$families$categorical_pairs$cross_artifacts,
    error = function(error) NULL)
  vertical <- tryCatch(
    manifest$workload$vertical_crosses$configured_crosses,
    error = function(error) NULL)
  scope <- tryCatch(manifest$workload$primitive_scope,
                    error = function(error) NULL)
  explicit <- tryCatch(scope$selection$explicit_catalog,
                       error = function(error) NULL)
  artifact <- if (is.list(artifacts) && length(artifacts) == 1L &&
      identical(names(artifacts), selector$analysis_id)) {
    artifacts[[1L]]
  } else NULL
  empty <- function(value) {
    (is.list(value) || is.character(value)) && length(value) == 0L
  }
  expected_sides <- lapply(selector$columns, function(column) list(
    dataset = column$dataset, column = column$column,
    owner_peer = column$owner_peer))
  actual_sides <- if (is.list(artifact)) list(
    artifact$left[c("dataset", "column", "owner_peer")],
    artifact$right[c("dataset", "column", "owner_peer")]) else NULL
  valid <- identical(
      .dsvert_joint_dp_client_json(manifest$logical_snapshot),
      .dsvert_joint_dp_client_json(selector$parent$logical_snapshot)) &&
    is.list(artifact) && identical(
      artifact$version,
      .DSVERT_CLIENT_DP_CATEGORICAL_CROSS_ARTIFACT_VERSION) &&
    identical(actual_sides, expected_sides) &&
    is.list(vertical) && length(vertical) == 1L &&
    identical(names(vertical), selector$analysis_id) &&
    identical(scope$mode, "catalog_v1") && is.list(explicit) &&
    all(vapply(c(
      "numeric_moments", "categorical_marginals", "categorical_pairs",
      "correlations"), function(name) empty(explicit[[name]]), logical(1L)))
  if (!isTRUE(valid)) {
    stop("The signed cross-owner categorical Synopsis projection is invalid",
         call. = FALSE)
  }
  invisible(manifest)
}

.dsvert_dp_synopsis_supported_categorical_cross_v1 <- function(manifest) {
  workload <- if (is.list(manifest)) manifest$workload else NULL
  families <- if (is.list(workload)) workload$families else NULL
  categorical <- if (is.list(families)) {
    families$categorical_pairs
  } else NULL
  artifacts <- if (is.list(categorical)) {
    categorical$cross_artifacts
  } else NULL
  artifact <- if (is.list(artifacts) && length(artifacts) == 1L) {
    artifacts[[1L]]
  } else NULL
  vertical <- if (is.list(workload)) workload$vertical_crosses else NULL
  configured <- if (is.list(vertical)) vertical$configured_crosses else NULL
  configured_artifact <- if (is.list(configured) &&
      length(configured) == 1L) configured[[1L]] else NULL
  scope <- if (is.list(workload)) workload$primitive_scope else NULL
  explicit <- tryCatch(
    scope$selection$explicit_catalog, error = function(error) NULL)
  empty <- function(value) {
    (is.list(value) || is.character(value)) && length(value) == 0L
  }
  structured <- c(
    "numeric_moments", "numeric_pair_moments", "gaussian_models",
    "fixed_numeric_histograms")
  direct <- c(
    "correlation_artifacts", "describe_artifacts", "survival_artifacts")
  other_empty <- is.list(families) && all(vapply(
    structured, function(name) {
      candidate <- families[[name]]
      is.list(candidate) && is.list(candidate$artifacts) &&
        length(candidate$artifacts) == 0L
    }, logical(1L))) && all(vapply(direct, function(name) {
      candidate <- families[[name]]
      is.list(candidate) && length(candidate) == 0L
    }, logical(1L)))
  explicit_empty <- is.list(explicit) && all(vapply(c(
    "numeric_moments", "categorical_marginals", "categorical_pairs",
    "correlations"), function(name) empty(explicit[[name]]), logical(1L)))
  is.list(artifact) && identical(
    artifact$version,
    .DSVERT_CLIENT_DP_CATEGORICAL_CROSS_ARTIFACT_VERSION) &&
    identical(artifact$spec_version, "v2") &&
    identical(artifact$implementation_state,
              "cross_owner_exact_gc_materialized") &&
    is.list(configured_artifact) &&
    identical(configured_artifact$version, artifact$version) &&
    identical(configured_artifact$analysis_id, artifact$analysis_id) &&
    identical(vertical$implementation_state,
              "categorical_pair_exact_gc_materialized") &&
    identical(as.numeric(vertical$included_coordinate_count),
              as.numeric(artifact$coordinate_count)) &&
    identical(scope$mode, "catalog_v1") && isTRUE(explicit_empty) &&
    is.list(categorical$sets) && length(categorical$sets) == 0L &&
    isTRUE(other_empty)
}

.dsvert_dp_synopsis_categorical_cross_remote_context_v1 <- function(value) {
  fields <- c("manifest_sha256", "claim_set_json", "compilation_json")
  if (!is.list(value) || !.dsvert_dp_has_exact_names(value, fields) ||
      !.dsvert_vector_hex(value$manifest_sha256) ||
      !all(vapply(value[c("claim_set_json", "compilation_json")],
                  .dsvert_dp_is_string, logical(1L)))) {
    stop("Invalid categorical cross Synopsis remote context", call. = FALSE)
  }
  value
}

.dsvert_dp_synopsis_categorical_cross_bind_call_v1 <- function(
    context, analysis_id, session_id) {
  context <- .dsvert_dp_synopsis_categorical_cross_remote_context_v1(context)
  call(
    name = "dsvertDPSynopsisCategoricalCrossBindDS",
    manifest_sha256 = context$manifest_sha256,
    claim_set_json = context$claim_set_json,
    compilation_json = context$compilation_json,
    analysis_id = analysis_id, session_id = session_id)
}

.dsvert_dp_synopsis_categorical_cross_finalize_call_v1 <- function(
    context, analysis_id, session_id) {
  context <- .dsvert_dp_synopsis_categorical_cross_remote_context_v1(context)
  call(
    name = "dsvertDPSynopsisCategoricalCrossFinalizeDS",
    manifest_sha256 = context$manifest_sha256,
    claim_set_json = context$claim_set_json,
    compilation_json = context$compilation_json,
    analysis_id = analysis_id, session_id = session_id)
}

.dsvert_dp_synopsis_categorical_cross_source_receipt_v1 <- function(
    source_manifest, tickets, authorities, sources) {
  first <- tickets[[authorities[[1L]]]]$value
  valid <- is.list(source_manifest) && isTRUE(source_manifest$cross_enabled) &&
    is.list(tickets) && all(authorities %in% names(tickets)) &&
    is.list(first) && .dsvert_dp_capsule_source_hex(first$contract_hash) &&
    identical(unname(unlist(first$source_peers, use.names = FALSE)), sources)
  if (!isTRUE(valid)) {
    stop("Invalid categorical cross Synopsis source receipt", call. = FALSE)
  }
  list(
    version = .DSVERT_CLIENT_DP_CAPSULE_SOURCE_RECEIPT_VERSION,
    phase = "recipient_aggregates_committed",
    purpose = source_manifest$purpose,
    capsule_id = source_manifest$capsule_id,
    contract_hash = first$contract_hash,
    peer_pinset_sha256 = first$peer_pinset_sha256,
    source_peers = sources,
    designated_noise_peers = authorities,
    coordinate_count = as.numeric(first$coordinate_count),
    chunk_coordinates = as.numeric(first$chunk_coordinates),
    chunk_count = as.numeric(first$chunk_count), ring_bits = 128,
    record_encoding = "little_endian_unsigned_fixed_16_bytes",
    durable_replay = TRUE, operation_or_request_limit = FALSE,
    history_can_deny_operation = FALSE, payload_exposed = FALSE,
    sampler_handoff_ready = FALSE,
    release_coordinate_count = source_manifest$release_coordinate_count,
    release_coordinate_order_sha256 =
      source_manifest$release_coordinate_order_sha256,
    private_layout_sha256 = source_manifest$private_layout_sha256)
}
