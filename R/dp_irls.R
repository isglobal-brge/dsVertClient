# Thin B3 selector-only orchestration. All statistical choices come from the
# unanimous signed contract; the client independently verifies public updates.
.dsvert_dp_irls_client_id <- function(analysis_id) {
  if (!is.character(analysis_id) || length(analysis_id) != 1L ||
      is.na(analysis_id) || !grepl("^[A-Za-z0-9][A-Za-z0-9._:-]{0,127}$", analysis_id))
    stop("analysis_id must be one bounded ASCII identifier.", call. = FALSE)
  analysis_id
}

.dsvert_dp_irls_client_pinned <- function(datasources) {
  pins <- getOption("dsvert.dp_irls.trusted_pins")
  if (is.null(pins)) stop(
    "Custodian-verified dsvert.dp_irls.trusted_pins are required before IRLS startup.",
    call. = FALSE)
  pins <- .dsvert_dp_irls_client_pins(unlist(pins))
  if (!setequal(names(pins), names(datasources)))
    stop("IRLS datasource names do not match the trusted pinset.", call. = FALSE)
  pins
}

.dsvert_dp_irls_client_call <- function(endpoint, analysis_id, datasources,
                                       .aggregate = DSI::datashield.aggregate) {
  .dsvert_aggregate_strict(conns = datasources,
    expr = call(endpoint, analysis_id = analysis_id),
    operation = paste("DP IRLS", endpoint), .aggregate = .aggregate)
}

.dsvert_dp_irls_client_metadata <- function(analysis_id, datasources, pins,
                                           .aggregate = DSI::datashield.aggregate) {
  values <- .dsvert_dp_irls_client_call("dpIRLSMetadataDS", analysis_id,
                                      datasources, .aggregate)
  checked <- lapply(values, function(value) {
    if (!.dsvert_dp_has_exact_names(value, c("version", "analysis_id", "root",
                                             "signed_contract", "stages")) ||
        !identical(value$version, "dsvert-dp-irls-metadata-v1") ||
        !identical(value$analysis_id, analysis_id))
      stop("Invalid IRLS metadata response.", call. = FALSE)
    verified <- .dsvert_dp_irls_verify_chain(value$signed_contract, pins,
      value$stages, expected_root = value$root)
    list(metadata = value, verified = verified)
  })
  identity <- vapply(checked, function(x) .dsvert_dp_irls_client_state_json(
    x$metadata[c("analysis_id", "root", "signed_contract")]), character(1L))
  if (length(unique(identity)) != 1L)
    stop("Authorities returned different IRLS contracts.", call. = FALSE)
  # A interrupted last acknowledgment may leave one authority behind. Only a
  # common prefix is treated as public until publication confirms both peers.
  first <- checked[[1L]]
  common <- min(vapply(checked, function(x) length(x$metadata$stages), integer(1L)))
  for (i in seq_len(common)) {
    bytes <- vapply(checked, function(x)
      .dsvert_dp_irls_client_state_json(x$metadata$stages[[i]]), character(1L))
    if (length(unique(bytes)) != 1L)
      stop("Authorities returned inconsistent IRLS predecessors.", call. = FALSE)
  }
  first$metadata$stages <- first$metadata$stages[seq_len(common)]
  first$verified <- .dsvert_dp_irls_verify_chain(first$metadata$signed_contract,
    pins, first$metadata$stages, expected_root = first$metadata$root)
  first
}

.dsvert_dp_irls_client_publication <- function(analysis_id, datasources, metadata,
    pins, .aggregate = DSI::datashield.aggregate) {
  values <- .dsvert_dp_irls_client_call("dpIRLSPublicationDS", analysis_id,
                                      datasources, .aggregate)
  checked <- lapply(values, function(value) {
    if (!.dsvert_dp_has_exact_names(value,
        c("version", "analysis_id", "root", "stages")) ||
        !identical(value$version, "dsvert-dp-irls-publication-v1") ||
        !identical(value$analysis_id, analysis_id) ||
        !identical(value$root, metadata$root))
      stop("Invalid IRLS publication response.", call. = FALSE)
    .dsvert_dp_irls_verify_chain(metadata$signed_contract, pins, value$stages,
                                 expected_root = metadata$root)
    value
  })
  bytes <- vapply(checked, .dsvert_dp_irls_client_state_json, character(1L))
  if (length(unique(bytes)) != 1L)
    stop("IRLS publication is not durable at both authorities.", call. = FALSE)
  checked[[1L]]
}

# The internal stage driver is also the cost-gate boundary. It never accepts
# beta, predecessor, budget, family, seeds, bounds or any source values.
.dsvert_dp_irls_client_stage <- function(analysis_id, datasources, metadata,
    pins, .aggregate = DSI::datashield.aggregate) {
  old <- length(metadata$stages)
  contract <- metadata$signed_contract$contract
  if (old >= length(contract$accounting$stages))
    return(.dsvert_dp_irls_client_publication(analysis_id, datasources,
                                            metadata, pins, .aggregate))
  status <- .dsvert_dp_irls_client_call("dpIRLSStartDS", analysis_id,
                                       datasources, .aggregate)
  deadline <- proc.time()[["elapsed"]] + as.numeric(contract$operations$max_stage_seconds)
  frames <- stats::setNames(rep(list(""), length(datasources)), names(datasources))
  repeat {
    for (value in status) {
      if (!.dsvert_dp_has_exact_names(value, c("version", "analysis_id", "root",
                                               "status", "stage_index")) ||
          !identical(value$version, "dsvert-dp-irls-status-v1") ||
          !identical(value$analysis_id, analysis_id) ||
          !identical(value$root, metadata$root) ||
          !value$status %in% c("running", "stage-published", "complete") ||
          !identical(value$stage_index, as.character(old)))
        stop("Invalid IRLS stage status.", call. = FALSE)
    }
    if (all(vapply(status, function(x) x$status != "running", logical(1L)))) break
    if (proc.time()[["elapsed"]] > deadline)
      stop("IRLS stage exceeded its signed time envelope.", call. = FALSE)
    if (length(datasources) != 2L)
      stop("The admitted IRLS relay requires two authorities.", call. = FALSE)
    expressions <- stats::setNames(lapply(seq_along(datasources), function(i)
      call("dpIRLSRelayDS", analysis_id = analysis_id, frame = frames[[3L - i]])),
      names(datasources))
    relayed <- .dsvert_aggregate_strict(conns = datasources, expr = expressions,
      operation = "DP IRLS encrypted transport", .aggregate = .aggregate)
    for (peer in names(datasources)) {
      value <- relayed[[peer]]
      if (!.dsvert_dp_has_exact_names(value, c("version", "analysis_id", "root",
          "status", "stage_index", "frame")) ||
          !identical(value$version, "dsvert-dp-irls-relay-v1") ||
          !identical(value$analysis_id, analysis_id) ||
          !identical(value$root, metadata$root) ||
          !is.character(value$frame) || length(value$frame) != 1L || is.na(value$frame))
        stop("Invalid encrypted IRLS transport response.", call. = FALSE)
      frames[[peer]] <- value$frame
      value$frame <- NULL
      value$version <- "dsvert-dp-irls-status-v1"
      status[[peer]] <- value
    }
    if (!any(nzchar(unlist(frames)))) Sys.sleep(0.02)
  }
  publication <- .dsvert_dp_irls_client_publication(analysis_id, datasources,
                                                   metadata, pins, .aggregate)
  if (length(publication$stages) != old + 1L)
    stop("IRLS stage publication did not advance exactly once.", call. = FALSE)
  for (i in seq_len(old)) if (!identical(publication$stages[[i]], metadata$stages[[i]]))
    stop("IRLS publication changed a retained predecessor.", call. = FALSE)
  publication
}

.dsvert_dp_irls_client_decode <- function(stage, contract) {
  decoded <- list(); endpoint <- FALSE
  for (lane in c("g", "H")) {
    values <- stage$body$payload[[lane]]
    if (is.null(values)) next
    translation <- contract$numeric$translation[[lane]]
    limits <- contract$numeric$upper_limits[[lane]]
    unit <- .dsvert_dp_irls_client_bigq(contract$numeric$lattice[[lane]])
    decoded[[lane]] <- lapply(seq_along(values), function(j) {
      v <- gmp::as.bigz(values[[j]])
      endpoint <<- endpoint || v == 0 || v == gmp::as.bigz(limits[[j]])
      .dsvert_dp_irls_client_rational((v - gmp::as.bigz(translation[[j]])) * unit)
    })
  }
  list(released = decoded, endpoint_hit = endpoint)
}

#' Execute a signed differential-private IRLS path
#'
#' Selects a custodian-provisioned analysis. Set the session option
#' `dsvert.dp_irls.trusted_pins` to the independently verified named peer pinset
#' before use. All signed public updates are independently verified in R.
#' @param analysis_id Custodian analysis selector.
#' @param datasources Named DataSHIELD connections.
#' @param verbose Emit fixed public progress messages.
#' @return Verified released chain and coefficients. Sampling inference is
#'   unavailable in this B3 adapter.
#' @export
ds.vertDPIRLS <- function(analysis_id, datasources = NULL, verbose = TRUE) {
  analysis_id <- .dsvert_dp_irls_client_id(analysis_id)
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose))
    stop("verbose must be one non-missing logical value.", call. = FALSE)
  datasources <- .dsvert_dp_datasources(datasources)
  pins <- .dsvert_dp_irls_client_pinned(datasources)
  resolved <- .dsvert_dp_irls_client_metadata(analysis_id, datasources, pins)
  metadata <- resolved$metadata
  if (verbose) message("Executing the signed DP IRLS path.")
  while (length(metadata$stages) < length(metadata$signed_contract$contract$accounting$stages)) {
    publication <- .dsvert_dp_irls_client_stage(analysis_id, datasources, metadata, pins)
    metadata$stages <- publication$stages
  }
  verified <- .dsvert_dp_irls_verify_chain(metadata$signed_contract, pins,
    metadata$stages, expected_root = metadata$root)
  decoded <- lapply(metadata$stages, .dsvert_dp_irls_client_decode,
                      contract = metadata$signed_contract$contract)
  beta <- vapply(verified$beta, function(q)
    as.numeric(.dsvert_dp_irls_client_bigq(q)), numeric(1L))
  names(beta) <- unlist(metadata$signed_contract$contract$estimand$coefficient_order)
  structure(list(analysis_id = analysis_id, root = metadata$root,
    coefficients = beta, signed_contract = metadata$signed_contract,
    stages = metadata$stages, released = decoded, accounting = verified$accounting,
    inference = "unavailable_pending_B4", std_errors = NULL,
    endpoint_hit = any(vapply(decoded, `[[`, logical(1L), "endpoint_hit"))),
    class = c("dsvert_dp_irls", "list"))
}
