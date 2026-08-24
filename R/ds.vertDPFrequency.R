.dsvert_dp_frequency_accuracy_v1 <- function(config, confidence) {
  if (!is.numeric(confidence) || length(confidence) != 1L ||
      is.na(confidence) || !is.finite(confidence) ||
      confidence <= 0 || confidence >= 1) {
    stop("confidence must be one finite number in (0, 1)", call. = FALSE)
  }
  plan <- .dsvert_dp_frequency_client_plan_v1(config)
  selection <- config$backend_selection
  full <- selection$selected_plan
  certificate <- selection$selected_accuracy_certificate
  maximum <- config$coordinate_upper_bound
  maximum_big <- openssl::bignum(format(maximum, scientific = FALSE, trim = TRUE))
  capped <- function(value, positive = FALSE) {
    value <- .dsvert_dp_analysis_frequency_uint_v1(value, positive)
    if (value >= maximum_big) maximum else as.numeric(value)
  }
  signed_95 <- capped(certificate$simultaneous_95_abs)
  support <- capped(certificate$absolute_support, TRUE)
  tv <- min(1, .dsvert_dp_vector_fraction_upper(
    certificate$release_tv_upper_numerator,
    certificate$release_tv_upper_denominator))
  dimension <- as.integer(config$factor_domain$dimension)
  gaussian <- isTRUE(.dsvert_dp_analysis_frequency_profile_v1(
    selection$summary$selected_primitive)$gaussian)
  finite_support <- identical(certificate$method, "finite_support_v1")
  if (identical(as.numeric(confidence), 0.95)) {
    steps <- signed_95
    method <- paste("signed selected 95%", certificate$method,
                    "simultaneous certificate")
  } else if (gaussian) {
    accuracy_plan <- full
    accuracy_plan$maximum_noise_magnitude_two_peers <- format(
      support, scientific = FALSE, trim = TRUE)
    derived <- .dsvert_dp_vector_gaussian_accuracy_steps(
      accuracy_plan, dimension, confidence)
    steps <- min(derived$steps, maximum)
    finite_support <- isTRUE(derived$finite_support)
    if (confidence > 0.95) steps <- max(steps, signed_95)
    method <- if (finite_support) {
      "signed Gaussian finite-support simultaneous bound"
    } else {
      "confidence-specific signed Gaussian subgaussian/TV simultaneous bound"
    }
  } else {
    alpha <- .dsvert_dp_vector_next_down(.dsvert_dp_vector_next_down(
      1 - confidence - tv))
    steps <- support
    if (is.finite(alpha) && alpha > 0 && support > 0) {
      target <- .dsvert_dp_vector_next_down(log(.dsvert_dp_vector_next_down(
        alpha / dimension)))
      tail <- .dsvert_dp_vector_dyadic_tail_context(full)
      bound <- function(radius) .dsvert_dp_vector_plan_log_tail_upper(
        radius, tail, TRUE)
      if (bound(support - 1) <= target) {
        low <- -1
        high <- support - 1
        while (high - low > 1) {
          middle <- low + floor((high - low) / 2)
          if (bound(middle) <= target) high <- middle else low <- middle
        }
        steps <- high
        finite_support <- FALSE
      } else finite_support <- TRUE
    } else finite_support <- TRUE
    if (confidence > 0.95) steps <- max(steps, signed_95)
    method <- if (finite_support) {
      "signed convolution finite-support simultaneous bound"
    } else {
      "confidence-specific signed dyadic convolution-tail/TV simultaneous bound"
    }
  }
  list(radius = unname(min(steps, maximum)), confidence = confidence,
    method = method, implementation_tv_upper_bound = tv,
    finite_support_fallback = finite_support,
    selected_plan_sha256 = plan$full_plan_sha256)
}
.dsvert_dp_frequency_radius <- function(x, confidence) {
  x <- .dsvert_dp_frequency_contract(x)
  if (identical(x$proof$version,
                "dsvert-dp-frequency-synopsis-proof-v1")) {
    return(.dsvert_dp_frequency_synopsis_accuracy_v1(x, confidence)$radius)
  }
  .dsvert_dp_frequency_accuracy_v1(x$proof$config, confidence)$radius
}
.dsvert_dp_frequency_regions <- function(counts, radius, capacity) {
  if (!is.numeric(counts) || !length(counts) || anyNA(counts) ||
      any(!is.finite(counts)) || any(counts < 0) ||
      !.dsvert_dp_is_number(radius, 0) ||
      !.dsvert_dp_is_number(capacity, 1)) {
    stop("Invalid frequency-region inputs", call. = FALSE)
  }
  lower <- pmax(0, counts - radius)
  upper <- pmin(capacity, counts + radius)
  sum_lower <- sum(lower)
  sum_upper <- sum(upper)
  infeasible <- sum_lower > capacity || sum_upper <= 0
  if (infeasible) {
    proportions <- cbind(lower = rep(0, length(counts)),
                         upper = rep(1, length(counts)))
  } else {
    other_lower <- sum_lower - lower
    other_upper <- sum_upper - upper
    maximum_other <- pmin(other_upper, capacity - lower)
    lower_denominator <- lower + pmax(0, maximum_other)
    feasible_coordinate_upper <- pmax(lower, pmin(
      upper, capacity - other_lower))
    upper_denominator <- feasible_coordinate_upper + other_lower
    proportions <- cbind(
      lower = ifelse(lower_denominator > 0,
                     lower / lower_denominator, 0),
      upper = ifelse(upper_denominator > 0,
                     feasible_coordinate_upper / upper_denominator, 1))
  }
  rownames(proportions) <- names(counts)
  list(count = cbind(lower = lower, upper = upper),
    proportion = proportions, feasible = !infeasible,
    includes_zero_effective_count = sum_lower == 0,
    has_positive_effective_count = sum_upper > 0)
}
.dsvert_dp_frequency_public_execution_v1 <- function(execution) {
  fields <- c("version", "operation_id", "variable_name", "levels", "values",
              "source_owner", "finalizer_peer", "proof")
  proof_fields <- c("session_id", "claim", "config", "receipts", "contract",
                    "worker_static", "authorities", "authorizations", "release")
  valid <- .dsvert_dp_frequency_client_object_v1(execution, fields) &&
    identical(execution$version,
              .DSVERT_CLIENT_DP_FREQUENCY_RESULT_VERSION) &&
    is.character(execution$operation_id) && length(execution$operation_id) == 1L &&
    !is.na(execution$operation_id) &&
    grepl("^op_[0-9a-f]{32}$", execution$operation_id) &&
    .dsvert_dp_frequency_client_object_v1(execution$proof, proof_fields)
  if (!isTRUE(valid)) stop(
    "Invalid closed Frequency execution result.", call. = FALSE)
  proof <- execution$proof
  peers <- if (is.list(proof$receipts)) names(proof$receipts) else NULL
  if (is.null(peers) || length(peers) < 2L || length(peers) > 4096L ||
      anyNA(peers) || anyDuplicated(peers)) {
    stop("Invalid public Frequency proof.", call. = FALSE)
  }
  envelopes <- stats::setNames(lapply(peers, function(peer) list(
    config = proof$config, receipt = proof$receipts[[peer]])), peers)
  compiled <- .dsvert_dp_frequency_client_compile_v1(envelopes, peers, proof$claim)
  same <- identical(proof$claim, compiled$claim) &&
    identical(proof$config, compiled$config) &&
    identical(proof$receipts, compiled$receipts) &&
    identical(proof$contract, compiled$contract) &&
    identical(proof$worker_static, compiled$worker_static) &&
    identical(proof$authorities, compiled$authorities)
  session <- proof$session_id
  same <- same && is.character(session) && length(session) == 1L &&
    !is.na(session) && grepl(paste0(
      "^[0-9a-f]{8}-[0-9a-f]{4}-4[0-9a-f]{3}-",
      "[89ab][0-9a-f]{3}-[0-9a-f]{12}$"), session)
  if (!isTRUE(same)) stop(
    "Public Frequency proof failed recompilation.", call. = FALSE)
  authorizations <- .dsvert_dp_frequency_client_authorizations_v1(
    proof$authorizations, session, compiled)
  geometry <- .dsvert_dp_frequency_execution_geometry_v1(compiled$worker_static)
  context <- list(
    compiled = compiled, geometry = geometry,
    authorization_set_sha256 = .dsvert_dp_frequency_client_hash_v1(
      .DSVERT_CLIENT_DP_FREQUENCY_AUTH_SET_DOMAIN,
      unname(authorizations[compiled$authorities])))
  release <- .dsvert_dp_frequency_execution_release_v1(proof$release, context)
  levels <- unlist(compiled$config$factor_domain$levels, use.names = FALSE)
  upper <- compiled$config$coordinate_upper_bound
  values <- execution$values
  valid <- identical(release, proof$release) &&
    identical(execution$variable_name,
              compiled$config$factor_domain$variable_name) &&
    identical(execution$levels, levels) &&
    identical(execution$source_owner, compiled$authorities[[1L]]) &&
    identical(execution$finalizer_peer, compiled$authorities[[2L]]) &&
    is.numeric(values) && length(values) == geometry$d && !anyNA(values) &&
    all(is.finite(values)) && all(values == floor(values)) &&
    all(values >= 0 & values <= upper)
  if (!isTRUE(valid)) stop("Invalid public Frequency release values.", call. = FALSE)
  chunk_hashes <- vapply(seq_len(geometry$chunk_count) - 1L, function(index) {
    offset <- index * .DSVERT_CLIENT_DP_FREQUENCY_CHUNK_COORDINATES
    count <- min(.DSVERT_CLIENT_DP_FREQUENCY_CHUNK_COORDINATES,
                 geometry$d - offset)
    .dsvert_dp_frequency_execution_chunk_hash_v1(
      values[seq.int(offset + 1L, offset + count)], index, offset)
  }, character(1L))
  window_hashes <- vapply(seq_len(geometry$window_count) - 1L, function(index) {
    window <- .dsvert_dp_frequency_execution_window_v1(geometry, index)
    range <- seq.int(window$first_chunk + 1L,
                     window$first_chunk + window$chunks)
    .dsvert_dp_frequency_client_hash_v1(
      .DSVERT_CLIENT_DP_FREQUENCY_WINDOW_DOMAIN, list(
        version = "dsvert-dp-frequency-final-window-v1",
        window_index = index, coordinate_offset = window$offset,
        coordinate_count = window$count,
        chunk_hashes = as.list(chunk_hashes[range])))
  }, character(1L))
  committed <- identical(unname(chunk_hashes), unlist(
    release$final_chunk_hashes, use.names = FALSE)) &&
    identical(unname(window_hashes), unlist(
      release$final_window_hashes, use.names = FALSE)) &&
    identical(.dsvert_dp_frequency_execution_merkle_v1(chunk_hashes),
              release$final_vector_root)
  if (!isTRUE(committed)) stop(
    "Frequency values disagree with the signed release commitment.",
    call. = FALSE)
  list(compiled = compiled, release = release, values = unname(values),
    levels = levels, proof = proof, operation_id = execution$operation_id,
    source_owner = execution$source_owner, finalizer_peer = execution$finalizer_peer)
}
#' Differentially private one-way frequency distribution
#'
#' Reads one signed categorical marginal from the durable canonical Synopsis.
#' It performs no extra sampling or public opening; replay is post-processing
#' of that same per-artifact release. `server` is the required source owner.
#' @param data_name Name of the registered protected data frame.
#' @param variable Fixed-domain categorical variable.
#' @param server Required source-owner peer name.
#' @param datasources Named DataSHIELD connections.
#' @return A `ds.vertDPFrequency` object.
#' @export
ds.vertDPFrequency <- function(data_name, variable, server = NULL,
                               datasources = NULL) {
  resolved <- .dsvert_federation_argument(data_name, datasources)
  .dsvert_dp_frequency_impl(
    resolved$value, variable, server, resolved$datasources,
    DSI::datashield.aggregate)
}

.dsvert_dp_frequency_intercept_formula <- function(formula) {
  inherits(formula, "formula") && length(formula) == 3L &&
    is.symbol(formula[[2L]]) && {
      terms <- stats::terms(formula)
      identical(as.integer(attr(terms, "intercept")), 1L) &&
        length(attr(terms, "term.labels")) == 0L
    }
}

.dsvert_dp_frequency_intercept_artifact <- function(
    formula, data, server, datasources, method) {
  if (!inherits(formula, "formula") || length(formula) != 3L ||
      !is.symbol(formula[[2L]])) {
    stop("Automatic Frequency-backed ", method,
         " requires a simple outcome formula", call. = FALSE)
  }
  terms <- stats::terms(formula)
  if (!identical(as.integer(attr(terms, "intercept")), 1L) ||
      length(attr(terms, "term.labels")) != 0L) {
    stop("Automatic Frequency-backed ", method,
         " supports only an intercept-only y ~ 1 formula", call. = FALSE)
  }
  if (!is.character(data) || length(data) != 1L || is.na(data) ||
      !nzchar(data)) {
    stop("Automatic Frequency-backed ", method,
         " requires one non-empty data name", call. = FALSE)
  }
  if (is.null(server)) {
    stop("Automatic Frequency-backed ", method,
         " requires an explicit source owner", call. = FALSE)
  }
  ds.vertDPFrequency(data, as.character(formula[[2L]]), server = server,
                     datasources = datasources)
}

.dsvert_dp_frequency_synopsis_accuracy_v1 <- function(x, confidence) {
  table <- x
  table$table <- unname(x$counts)
  profile <- .dsvert_dp_table_vector_profile(table)
  if (!is.list(profile)) {
    stop("Invalid Synopsis Frequency profile", call. = FALSE)
  }
  total_coordinates <- table$synopsis_coordinate_count
  if (!.dsvert_dp_is_integer(total_coordinates, 1,
                             .DSVERT_DP_MAX_COORDINATES)) {
    stop("Invalid Synopsis Frequency coordinate count", call. = FALSE)
  }
  release <- list(
    epsilon = table$epsilon,
    implementation_delta = table$implementation_delta,
    mechanism = table$mechanism,
    backend = table$implementation,
    mechanism_plan = table$mechanism_plan,
    plan_sha256 = table$plan_sha256,
    backend_selection = table$backend_selection,
    backend_assessment = table$backend_assessment,
    manifest_sha256 = table$manifest_sha256)
  manifest <- list(workload = list(
    coordinate_count = total_coordinates,
    capsule_mechanism = table$synopsis_mechanism,
    mechanism_selection = table$mechanism_selection,
    release_lattice = list(
      output_lattice_bits = table$output_lattice_bits,
      output_lattice_scale = table$output_lattice_scale,
      natural_l1_sensitivity = table$l1_sensitivity,
      integer_l1_sensitivity_steps =
        table$l1_sensitivity * table$output_lattice_scale,
      natural_l2_sensitivity = table$l2_sensitivity,
      integer_l2_sensitivity_steps =
        table$l2_sensitivity * table$output_lattice_scale)))
  .dsvert_dp_vector_accuracy_radius(
    release, manifest, coordinate_count = total_coordinates,
    confidence = confidence, maximum_error = table$coordinate_maximum)
}

.dsvert_dp_frequency_synopsis_result_v1 <- function(
    data_name, variable, source, datasources, .aggregate,
    .run = .dsvert_dp_synopsis_vector_run) {
  run <- .run(datasources, .aggregate = .aggregate)
  context <- .dsvert_dp_vector_context(run, allow_synopsis = TRUE)
  block <- .dsvert_dp_capsule_single_block(
    context$layout, "categorical_marginals", dataset = data_name,
    owner_peer = source,
    predicate = function(candidate) identical(
      candidate$descriptor$column, variable),
    description = paste0("signed categorical marginal for '", variable, "'"))
  descriptor <- block$descriptor
  levels <- .dsvert_dp_capsule_manifest_strings(
    descriptor$levels, "categorical marginal levels", sorted = TRUE)
  capacity <- .dsvert_dp_vector_block_capacity(block)
  counts <- .dsvert_dp_capsule_vector_values(context$release, block)
  if (length(counts) != length(levels) || any(counts < 0) ||
      any(counts > capacity)) {
    stop("The released categorical marginal violates its signed domain",
         call. = FALSE)
  }
  metadata <- .dsvert_dp_vector_public_metadata(context)
  result <- c(metadata, list(
    status = if (sum(counts) > 0) "ok" else "dp_effective_count_zero",
    server = source, source_owner = source, variable = variable,
    levels = levels, counts = stats::setNames(counts, levels),
    effective_count_dp = sum(counts),
    proportions = if (sum(counts) > 0) stats::setNames(
      counts / sum(counts), levels) else
      stats::setNames(rep(NA_real_, length(counts)), levels),
    coordinate_family = "categorical_marginals",
    coordinate_descriptor = descriptor,
    coordinate_maximum = capacity,
    repeated_record_policy = descriptor$repeated_record_policy,
    missingness_policy = descriptor$missingness_policy,
    primitive = "signed_categorical_marginal_v1",
    public_openings = 1L,
    additional_privacy_cost = c(epsilon = 0, delta = 0),
    operation_id = metadata$execution_id,
    release_sha256 = metadata$final_vector_root,
    proof = list(
      version = "dsvert-dp-frequency-synopsis-proof-v1",
      coordinate_family = "categorical_marginals",
      coordinate_descriptor = descriptor,
      release_provenance = metadata$release_provenance),
    uncertainty_scope = paste(
      "Simultaneous regions cover DP mechanism noise only; population",
      "sampling uncertainty is excluded"),
    inferential_scope = paste(
      "Finite-snapshot DP frequencies; use ds.vertDPFrequencyInference()",
      "for conservative iid multinomial sampling regions")))
  accuracy <- .dsvert_dp_frequency_synopsis_accuracy_v1(result, 0.95)
  result$accuracy_simultaneous_95_abs <- accuracy$radius
  result$accuracy_simultaneous_confidence <- accuracy$confidence
  result$accuracy_simultaneous_method <- accuracy$method
  result$accuracy_implementation_tv_upper_bound <-
    accuracy$implementation_tv_upper_bound
  result$mechanism_regions <- .dsvert_dp_frequency_regions(
    result$counts, accuracy$radius, capacity)
  class(result) <- c("ds.vertDPFrequency", "list")
  .dsvert_dp_frequency_contract(result)
}

.dsvert_dp_frequency_legacy_impl <- function(
    data_name, variable, server = NULL, datasources = NULL, .aggregate,
    .prepare = .dsvert_dp_frequency_prepare_v1,
    .execute = .dsvert_dp_frequency_execute_v1) {
  values <- list(data_name = data_name, variable = variable)
  if (any(!vapply(values, function(value) is.character(value) &&
      length(value) == 1L && !is.na(value) && nzchar(value), logical(1L)))) {
    stop("data_name and variable must be non-empty strings", call. = FALSE)
  }
  if (is.null(server)) stop(
    "Frequency requires an explicit source owner before DSI.", call. = FALSE)
  if (!is.function(.prepare) || !is.function(.execute)) stop(
    "Invalid Frequency client dependency.", call. = FALSE)
  datasources <- .dsvert_dp_datasources(datasources)
  source <- .dsvert_dp_server(server, datasources)
  prepared <- .prepare(data_name, variable_name = variable, source_owner = source,
    datasources = datasources, .aggregate = .aggregate)
  closed <- .dsvert_dp_frequency_public_execution_v1(.execute(
    prepared, data_name, datasources, .aggregate = .aggregate))
  config <- closed$compiled$config
  release <- closed$release
  counts <- stats::setNames(closed$values, closed$levels)
  total <- sum(counts)
  proportions <- if (total > 0) counts / total else
    stats::setNames(rep(NA_real_, length(counts)), closed$levels)
  accuracy <- .dsvert_dp_frequency_accuracy_v1(config, 0.95)
  regions <- .dsvert_dp_frequency_regions(counts, accuracy$radius,
    config$coordinate_upper_bound)
  result <- list(
    released = TRUE, status = if (total > 0) "ok" else
      "dp_effective_count_zero", server = closed$source_owner,
    source_owner = closed$source_owner, finalizer_peer = closed$finalizer_peer,
    variable = variable, levels = closed$levels, counts = counts,
    effective_count_dp = total, proportions = proportions,
    coordinate_maximum = config$coordinate_upper_bound,
    repeated_record_policy = config$repeated_record_policy,
    missingness_policy = config$missingness_policy,
    primitive = release$primitive, mechanism = release$mechanism,
    sampler = release$sampler, implementation = "exact_mpc_ring128_v1",
    epsilon = config$privacy$epsilon, delta = config$privacy$delta,
    implementation_delta_upper = config$calibration$implementation_delta,
    adjacency = config$privacy$adjacency,
    accuracy_simultaneous_95_abs = accuracy$radius,
    accuracy_simultaneous_confidence = 0.95,
    accuracy_simultaneous_method = accuracy$method,
    accuracy_implementation_tv_upper_bound =
      accuracy$implementation_tv_upper_bound,
    mechanism_regions = regions, sticky_noise = TRUE,
    sticky_replay = TRUE, intermediate_values_exposed = FALSE,
    public_openings = 1L,
    composition_rule = "one_sticky_release_per_canonical_signed_artifact",
    additional_privacy_cost = c(epsilon = 0, delta = 0),
    artifact_key = release$artifact_key, release_sha256 = release$release_sha256,
    operation_id = closed$operation_id, proof = closed$proof,
    uncertainty_scope = paste(
      "Simultaneous regions cover DP mechanism noise only; population",
      "sampling uncertainty is excluded"),
    inferential_scope = paste(
      "Finite-snapshot DP frequencies; use ds.vertDPFrequencyInference()",
      "for conservative iid multinomial sampling regions"))
  class(result) <- c("ds.vertDPFrequency", "list")
  .dsvert_dp_frequency_contract(result)
}

.dsvert_dp_frequency_impl <- function(
    data_name, variable, server = NULL, datasources = NULL, .aggregate,
    .prepare = NULL, .execute = NULL,
    .run = .dsvert_dp_synopsis_vector_run) {
  values <- list(data_name = data_name, variable = variable)
  if (any(!vapply(values, function(value) is.character(value) &&
      length(value) == 1L && !is.na(value) && nzchar(value), logical(1L)))) {
    stop("data_name and variable must be non-empty strings", call. = FALSE)
  }
  if (is.null(server)) stop(
    "Frequency requires an explicit source owner before DSI.", call. = FALSE)
  datasources <- .dsvert_dp_datasources(datasources)
  source <- .dsvert_dp_server(server, datasources)
  if (is.null(.prepare) && is.null(.execute)) {
    return(.dsvert_dp_frequency_synopsis_result_v1(
      data_name, variable, source, datasources, .aggregate, .run = .run))
  }
  if (!is.function(.prepare) || !is.function(.execute)) stop(
    "Invalid Frequency client dependency.", call. = FALSE)
  .dsvert_dp_frequency_legacy_impl(
    data_name, variable, source, datasources, .aggregate,
    .prepare = .prepare, .execute = .execute)
}

.dsvert_dp_frequency_synopsis_contract_v1 <- function(x, fail) {
  proof <- x$proof
  proof_fields <- c(
    "version", "coordinate_family", "coordinate_descriptor",
    "release_provenance")
  roots <- c(
    "artifact_key", "execution_id", "contract_sha256", "attempt_sha256",
    "source_contract_sha256", "result_set_sha256", "final_vector_root")
  recursive_names <- function(value) {
    if (!is.list(value)) return(character())
    c(names(value), unlist(lapply(value, recursive_names), use.names = FALSE))
  }
  forbidden <- c(
    "capsule_id", "privacy_epoch", "noise_key_id", "history_gate",
    "request_limit", "operation_limit", "lifetime_budget",
    "lifetime_composition", "privacy_accountant", "release_instance",
    "allocation_certificate", "ledger", "reservation", "rate_limit",
    "catalog_limit", "quota")
  descriptor <- proof$coordinate_descriptor
  levels <- tryCatch(.dsvert_dp_capsule_manifest_strings(
    descriptor$levels, "Frequency synopsis levels", sorted = TRUE),
    error = function(error) NULL)
  capacity <- tryCatch(.dsvert_dp_vector_block_capacity(
    list(descriptor = descriptor)), error = function(error) NA_real_)
  accuracy <- tryCatch(
    .dsvert_dp_frequency_synopsis_accuracy_v1(x, 0.95),
    error = function(error) NULL)
  expected_proportions <- if (sum(x$counts) > 0) {
    x$counts / sum(x$counts)
  } else stats::setNames(rep(NA_real_, length(x$counts)), x$levels)
  regions <- tryCatch(.dsvert_dp_frequency_regions(
    x$counts, accuracy$radius, capacity), error = function(error) NULL)
  expected_status <- if (sum(x$counts) > 0) "ok" else "dp_effective_count_zero"
  profile <- .dsvert_dp_table_vector_profile(x)
  valid_roots <- all(vapply(roots, function(field) {
    .dsvert_vector_hex(x[[field]]) &&
      identical(x[[field]], proof$release_provenance[[field]])
  }, logical(1L)))
  provenance <- x$release_provenance
  peers <- tryCatch(names(provenance$ordered_peer_pinset),
                    error = function(error) NULL)
  pinset <- tryCatch(.dsvert_dp_synopsis_client_pinset_v1(
    provenance$ordered_peer_pinset, peers), error = function(error) NULL)
  designated <- tryCatch(.dsvert_dp_synopsis_client_string_list_v1(
    provenance$designated_noise_peers, "Frequency noise authority set", 2L),
    error = function(error) NULL)
  valid_authorities <- !is.null(pinset) && length(pinset) >= 2L &&
    !is.null(designated) && length(designated) == 2L &&
    !anyDuplicated(designated) && all(designated %in% names(pinset)) &&
    identical(provenance$peer_pinset_sha256,
              .dsvert_vector_hash(as.list(pinset)))
  privacy <- x$privacy
  valid_privacy <- is.list(privacy) && identical(names(privacy), c(
    "version", "per_artifact_epsilon", "per_artifact_delta",
    "sticky_noise", "unlimited_replay", "replay_is_postprocessing",
    "public_openings", "distinct_artifacts_compose",
    "finite_global_composition_claim")) &&
    identical(privacy$version, "dsvert-per-synopsis-dp-v1") &&
    .dsvert_dp_num_equal(privacy$per_artifact_epsilon, x$epsilon, 128) &&
    .dsvert_dp_num_equal(privacy$per_artifact_delta, x$delta, 128) &&
    identical(privacy$sticky_noise, TRUE) &&
    identical(privacy$unlimited_replay, TRUE) &&
    identical(privacy$replay_is_postprocessing, TRUE) &&
    identical(privacy$public_openings, 1L) &&
    identical(privacy$distinct_artifacts_compose, TRUE) &&
    identical(privacy$finite_global_composition_claim, FALSE)
  valid <- .dsvert_dp_has_exact_names(proof, proof_fields) &&
    identical(proof$version, "dsvert-dp-frequency-synopsis-proof-v1") &&
    identical(proof$coordinate_family, "categorical_marginals") &&
    identical(proof$coordinate_descriptor, x$coordinate_descriptor) &&
    identical(proof$release_provenance, x$release_provenance) &&
    valid_roots && is.list(provenance) && valid_authorities &&
    identical(provenance$version,
              "dsvert-stateless-synopsis-public-provenance-v1") &&
    identical(provenance$durable_replay, TRUE) &&
    is.character(x$variable) && length(x$variable) == 1L && !is.na(x$variable) &&
    identical(x$variable, descriptor$column) &&
    identical(x$source_owner, descriptor$owner_peer) &&
    identical(x$server, descriptor$owner_peer) &&
    identical(x$levels, levels) && identical(names(x$counts), levels) &&
    is.numeric(x$counts) && length(x$counts) == length(levels) &&
    !anyNA(x$counts) && all(is.finite(x$counts)) &&
    all(x$counts >= 0 & x$counts <= capacity) &&
    identical(as.numeric(x$coordinate_maximum), as.numeric(capacity)) &&
    identical(x$repeated_record_policy, descriptor$repeated_record_policy) &&
    identical(x$missingness_policy, descriptor$missingness_policy) &&
    identical(x$coordinate_family, "categorical_marginals") &&
    identical(x$primitive, "signed_categorical_marginal_v1") &&
    identical(x$status, expected_status) &&
    identical(as.numeric(x$effective_count_dp), sum(x$counts)) &&
    isTRUE(all.equal(x$proportions, expected_proportions,
                     tolerance = 64 * .Machine$double.eps)) &&
    is.list(profile) && identical(x$implementation, profile$backend) &&
    identical(x$sampler, profile$sampler) &&
    identical(x$mechanism, profile$release_mechanism) &&
    identical(x$sticky_noise, "one_immutable_canonical_synopsis_artifact") &&
    identical(x$sticky_replay, TRUE) && identical(x$unlimited_replay, TRUE) &&
    identical(x$source_values_exposed, FALSE) &&
    identical(x$intermediate_values_exposed, FALSE) &&
    identical(x$clipped_coordinates, NA_integer_) &&
    identical(x$clamp_activation_disclosed, FALSE) &&
    identical(as.numeric(x$public_openings), 1) &&
    identical(x$additional_privacy_cost, c(epsilon = 0, delta = 0)) &&
    identical(x$composition_rule,
              "one_sticky_release_per_canonical_signed_artifact") &&
    identical(x$operation_id, x$execution_id) &&
    identical(x$release_sha256, x$final_vector_root) &&
    identical(x$uncertainty_scope, paste(
      "Simultaneous regions cover DP mechanism noise only; population",
      "sampling uncertainty is excluded")) &&
    identical(x$inferential_scope, paste(
      "Finite-snapshot DP frequencies; use ds.vertDPFrequencyInference()",
      "for conservative iid multinomial sampling regions")) &&
    isTRUE(valid_privacy) && !any(recursive_names(x) %in% forbidden) &&
    !is.null(accuracy) && !is.null(regions) &&
    identical(as.numeric(x$accuracy_simultaneous_95_abs),
              as.numeric(accuracy$radius)) &&
    identical(as.numeric(x$accuracy_simultaneous_confidence), 0.95) &&
    identical(x$accuracy_simultaneous_method, accuracy$method) &&
    identical(as.numeric(x$accuracy_implementation_tv_upper_bound),
              as.numeric(accuracy$implementation_tv_upper_bound)) &&
    isTRUE(all.equal(x$mechanism_regions, regions,
                     tolerance = 64 * .Machine$double.eps))
  if (!isTRUE(valid)) fail()
  x
}

.dsvert_dp_frequency_contract <- function(x) {
  fail <- function() stop(
    "x must be a released, validated ds.vertDPFrequency object",
    call. = FALSE)
  if (!inherits(x, "ds.vertDPFrequency") || !is.list(x) ||
      !isTRUE(x$released) || !is.list(x$proof)) fail()
  if (identical(x$proof$version,
                "dsvert-dp-frequency-synopsis-proof-v1")) {
    return(.dsvert_dp_frequency_synopsis_contract_v1(x, fail))
  }
  execution <- list(
    version = .DSVERT_CLIENT_DP_FREQUENCY_RESULT_VERSION,
    operation_id = x$operation_id, variable_name = x$variable,
    levels = x$levels, values = unname(x$counts),
    source_owner = x$source_owner, finalizer_peer = x$finalizer_peer,
    proof = x$proof)
  closed <- tryCatch(.dsvert_dp_frequency_public_execution_v1(execution),
                     error = function(error) NULL)
  if (is.null(closed)) fail()
  config <- closed$compiled$config
  accuracy <- tryCatch(.dsvert_dp_frequency_accuracy_v1(config, 0.95),
                       error = function(error) NULL)
  expected_proportions <- if (sum(x$counts) > 0) x$counts / sum(x$counts) else
    stats::setNames(rep(NA_real_, length(x$counts)), x$levels)
  regions <- tryCatch(.dsvert_dp_frequency_regions(
    x$counts, accuracy$radius, config$coordinate_upper_bound),
    error = function(error) NULL)
  expected_status <- if (sum(x$counts) > 0) "ok" else "dp_effective_count_zero"
  valid <- !is.null(accuracy) && !is.null(regions) &&
    identical(names(x$counts), x$levels) &&
    isTRUE(all.equal(x$proportions, expected_proportions,
                     tolerance = 64 * .Machine$double.eps)) &&
    identical(as.numeric(x$effective_count_dp), sum(x$counts)) &&
    identical(as.numeric(x$coordinate_maximum),
              as.numeric(config$coordinate_upper_bound)) &&
    identical(x$status, expected_status) &&
    identical(x$server, closed$source_owner) &&
    identical(x$variable, config$factor_domain$variable_name) &&
    identical(x$repeated_record_policy, config$repeated_record_policy) &&
    identical(x$missingness_policy, config$missingness_policy) &&
    identical(x$primitive, closed$release$primitive) &&
    identical(x$mechanism, closed$release$mechanism) &&
    identical(x$sampler, closed$release$sampler) &&
    identical(x$implementation, "exact_mpc_ring128_v1") &&
    identical(as.numeric(x$epsilon), as.numeric(config$privacy$epsilon)) &&
    identical(as.numeric(x$delta), as.numeric(config$privacy$delta)) &&
    identical(as.numeric(x$implementation_delta_upper),
              as.numeric(config$calibration$implementation_delta)) &&
    identical(x$adjacency, config$privacy$adjacency) &&
    identical(x$sticky_noise, TRUE) && identical(x$sticky_replay, TRUE) &&
    identical(x$intermediate_values_exposed, FALSE) &&
    identical(as.numeric(x$public_openings), 1) &&
    identical(as.numeric(x$accuracy_simultaneous_95_abs),
              as.numeric(accuracy$radius)) &&
    identical(as.numeric(x$accuracy_simultaneous_confidence), 0.95) &&
    identical(x$accuracy_simultaneous_method, accuracy$method) &&
    identical(as.numeric(x$accuracy_implementation_tv_upper_bound),
              as.numeric(accuracy$implementation_tv_upper_bound)) &&
    identical(x$additional_privacy_cost, c(epsilon = 0, delta = 0)) &&
    identical(x$composition_rule,
              "one_sticky_release_per_canonical_signed_artifact") &&
    identical(x$uncertainty_scope, paste(
      "Simultaneous regions cover DP mechanism noise only; population",
      "sampling uncertainty is excluded")) &&
    identical(x$inferential_scope, paste(
      "Finite-snapshot DP frequencies; use ds.vertDPFrequencyInference()",
      "for conservative iid multinomial sampling regions")) &&
    identical(x$artifact_key, closed$release$artifact_key) &&
    identical(x$release_sha256, closed$release$release_sha256) &&
    isTRUE(all.equal(x$mechanism_regions, regions,
                     tolerance = 64 * .Machine$double.eps))
  if (!isTRUE(valid)) fail()
  x
}
.dsvert_dp_frequency_cp_regions_v1 <- function(
    event_lower, event_upper, capacity, alpha) {
  valid <- is.numeric(event_lower) && is.numeric(event_upper) &&
    length(event_lower) > 0L && length(event_lower) == length(event_upper) &&
    !anyNA(event_lower) && !anyNA(event_upper) &&
    all(is.finite(event_lower)) && all(is.finite(event_upper)) &&
    all(event_lower >= 0 & event_lower <= event_upper) &&
    .dsvert_dp_is_number(capacity, 1) &&
    .dsvert_dp_is_number(alpha, 0, 1, lower_open = TRUE) && alpha < 1
  if (!isTRUE(valid)) stop("Invalid Frequency sampling box.", call. = FALSE)
  margin <- 128 * .Machine$double.eps
  lower <- ceiling(pmax(0, event_lower - margin * pmax(1, event_lower)))
  upper <- floor(pmin(2^53 - 1, event_upper +
    margin * pmax(1, event_upper)))
  broad <- function() {
    result <- cbind(lower = rep(0, length(lower)),
                    upper = rep(1, length(lower)))
    rownames(result) <- names(event_lower)
    result
  }
  if (any(lower > upper) || sum(lower) > capacity) return(broad())
  total_lower <- sum(lower)
  total_upper <- sum(upper)
  nonevent_lower <- total_lower - lower
  nonevent_upper <- pmax(
    nonevent_lower, pmin(total_upper - upper, capacity - lower))
  result_lower <- numeric(length(lower))
  selected <- lower > 0
  result_lower[selected] <- stats::qbeta(
    alpha / 2, lower[selected], nonevent_upper[selected] + 1)
  result_upper <- rep(1, length(lower))
  selected <- nonevent_lower > 0
  result_upper[selected] <- stats::qbeta(
    1 - alpha / 2, upper[selected] + 1, nonevent_lower[selected])
  result <- cbind(lower = result_lower, upper = result_upper)
  rownames(result) <- names(event_lower)
  if (anyNA(result) || any(!is.finite(result)) ||
      any(result < 0 | result > 1) || any(result[, 1L] > result[, 2L])) {
    stop("Frequency sampling regions are not representable.", call. = FALSE)
  }
  result
}
#' Conservative sampling regions for a DP frequency distribution
#'
#' Pure post-processing with a DP count box and exact binomial regions.
#' @param x A validated `ds.vertDPFrequency` object.
#' @param level Requested joint coverage.
#' @param dp_fraction Fraction of total error assigned to the DP event.
#' @return A `ds.vertDPFrequencyInference` object.
#' @export
ds.vertDPFrequencyInference <- function(x, level = 0.95,
                                        dp_fraction = 0.5) {
  x <- .dsvert_dp_frequency_contract(x)
  if (!.dsvert_dp_is_number(level, 0, 1, lower_open = TRUE) || level >= 1 ||
      !.dsvert_dp_is_number(dp_fraction, 0, 1, lower_open = TRUE) ||
      dp_fraction >= 1) stop(
    "level and dp_fraction must be finite numbers in (0, 1)", call. = FALSE)
  alpha <- 1 - level
  dp_alpha <- alpha * dp_fraction
  sampling_alpha <- alpha - dp_alpha
  per_level_alpha <- sampling_alpha / length(x$counts)
  accuracy <- if (identical(x$proof$version,
                             "dsvert-dp-frequency-synopsis-proof-v1")) {
    .dsvert_dp_frequency_synopsis_accuracy_v1(x, 1 - dp_alpha)
  } else {
    .dsvert_dp_frequency_accuracy_v1(x$proof$config, 1 - dp_alpha)
  }
  regions <- .dsvert_dp_frequency_regions(
    x$counts, accuracy$radius, x$coordinate_maximum)
  intervals <- .dsvert_dp_frequency_cp_regions_v1(
    regions$count[, "lower"], regions$count[, "upper"],
    x$coordinate_maximum, per_level_alpha)
  result <- list(
    status = "ok", levels = x$levels, counts_dp = x$counts,
    proportions_dp = x$proportions, intervals = intervals, level = level,
    coverage_lower_bound = level,
    coverage_method = paste(
      "simultaneous DP mechanism count box plus Bonferroni union of",
      "Clopper-Pearson binomial intervals"),
    dp_event_confidence = 1 - dp_alpha,
    sampling_familywise_confidence = 1 - sampling_alpha,
    base_sampling_interval_confidence = 1 - per_level_alpha,
    mechanism_radius = accuracy$radius,
    mechanism_accuracy_method = accuracy$method,
    sampling_model = "iid multinomial privacy units in the public domain",
    p_values = NULL, hypothesis_tests = NULL,
    additional_server_calls = 0L,
    additional_privacy_cost = c(epsilon = 0, delta = 0),
    artifact_key = x$artifact_key, release_sha256 = x$release_sha256)
  class(result) <- c("ds.vertDPFrequencyInference", "list")
  result
}
#' @export
print.ds.vertDPFrequency <- function(x, ...) {
  print(data.frame(
    level = x$levels, count_dp = unname(x$counts),
    proportion_dp = unname(x$proportions), row.names = NULL), ...)
  cat("simultaneous 95% mechanism radius:",
      x$accuracy_simultaneous_95_abs, "\n")
  invisible(x)
}
#' @export
print.ds.vertDPFrequencyInference <- function(x, ...) {
  print(x$intervals, ...)
  cat("Conservative joint coverage >= ", format(100 * x$level), "%\n",
      sep = "")
  invisible(x)
}
