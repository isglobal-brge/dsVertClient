# Internal client compiler for server-authoritative Count analyses.  The
# analyst supplies only a protected object name and the complete federation;
# every statistical and cryptographic parameter comes from signed server
# receipts.

.DSVERT_CLIENT_DP_COUNT_CONFIG_VERSION <- "dsvert-dp-count-config-v1"
.DSVERT_CLIENT_DP_COUNT_RECEIPT_VERSION <- "dsvert-dp-count-receipt-v1"
.DSVERT_CLIENT_DP_COUNT_RECEIPT_DOMAIN <- "dsVert/dp-count/receipt/v1|"
.DSVERT_CLIENT_DP_COUNT_CONFIG_DOMAIN <- "dsVert/dp-count/config/v1|"
.DSVERT_CLIENT_DP_COUNT_RECEIPT_SET_DOMAIN <-
  "dsVert/dp-count/receipt-set/v1|"
.DSVERT_CLIENT_DP_COUNT_CONTRACT_DOMAIN <-
  "dsVert/dp-count/compiled-contract/v1|"
.DSVERT_CLIENT_DP_COUNT_WORKER_DOMAIN <-
  "dsVert/dp-count/worker-static/v1|"
.DSVERT_CLIENT_DP_COUNT_PUBLIC_AUTH_VERSION <-
  "dsvert-dp-count-public-authorization-v1"
.DSVERT_CLIENT_DP_COUNT_PUBLIC_AUTH_DOMAIN <-
  "dsVert/dp-count/public-authorization/v1|"
.DSVERT_CLIENT_DP_COUNT_PUBLIC_AUTH_SIGNATURE_DOMAIN <-
  "dsVert/dp-count/public-authorization-signature/v1|"
.DSVERT_CLIENT_DP_COUNT_WORKER_VERSION <-
  "dsvert-dp-count-worker-static-v1"
.DSVERT_CLIENT_DP_COUNT_CONFIG_MAX_BYTES <- 1024L * 1024L
.DSVERT_CLIENT_DP_COUNT_RECEIPTS_MAX_BYTES <- 16L * 1024L^2

.dsvert_dp_count_client_hash_v1 <- function(domain, value) {
  digest::digest(
    charToRaw(paste0(domain, .dsvert_joint_dp_client_json(value))),
    algo = "sha256", serialize = FALSE)
}

.dsvert_dp_count_client_hex_v1 <- function(value, what) {
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !grepl("^[0-9a-f]{64}$", value)) {
    stop("Invalid Count ", what, ".", call. = FALSE)
  }
  value
}

.dsvert_dp_count_client_peer_name_v1 <- function(value) {
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !grepl("^[A-Za-z0-9][A-Za-z0-9_.-]{0,127}$", value)) {
    stop("Invalid Count peer name.", call. = FALSE)
  }
  value
}

.dsvert_dp_count_client_positive_integer_v1 <- function(
    value, what, maximum = Inf) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      !is.finite(value) || value < 1 || value != floor(value) ||
      value > maximum) {
    stop("Invalid Count ", what, ".", call. = FALSE)
  }
  as.numeric(value)
}

.dsvert_dp_count_client_decimal_text_v1 <- function(value) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      !is.finite(value) || value < 0) {
    stop("Invalid Count exact decimal value.", call. = FALSE)
  }
  sprintf("%.17g", as.numeric(value))
}

.dsvert_dp_count_client_decimal_rational_v1 <- function(value) {
  text <- if (is.character(value) && length(value) == 1L && !is.na(value)) {
    value
  } else {
    .dsvert_dp_count_client_decimal_text_v1(value)
  }
  if (nchar(text, type = "bytes") > 1024L) {
    stop("Invalid Count exact decimal value.", call. = FALSE)
  }
  match <- regexec(
    "^([0-9]+)(?:\\.([0-9]+))?(?:[eE]([+-]?[0-9]+))?$",
    text, perl = TRUE)
  parts <- regmatches(text, match)[[1L]]
  if (!length(parts)) {
    stop("Invalid Count exact decimal value.", call. = FALSE)
  }
  fraction <- if (length(parts) >= 3L && nzchar(parts[[3L]])) {
    parts[[3L]]
  } else {
    ""
  }
  exponent <- if (length(parts) >= 4L && nzchar(parts[[4L]])) {
    suppressWarnings(as.integer(parts[[4L]]))
  } else {
    0L
  }
  if (is.na(exponent) || abs(exponent) > 4096L) {
    stop("Invalid Count exact decimal value.", call. = FALSE)
  }
  coefficient <- sub(
    "^0+(?=[0-9])", "", paste0(parts[[2L]], fraction), perl = TRUE)
  power <- exponent - nchar(fraction, type = "bytes")
  numerator <- openssl::bignum(coefficient)
  denominator <- openssl::bignum(1)
  if (power >= 0L) {
    numerator <- numerator * openssl::bignum(10)^power
  } else {
    denominator <- denominator * openssl::bignum(10)^(-power)
  }
  list(numerator = numerator, denominator = denominator)
}

.dsvert_dp_count_client_integer_text_v1 <- function(
    value, allow_zero = FALSE) {
  is.character(value) && length(value) == 1L && !is.na(value) &&
    nchar(value, type = "bytes") <= 4096L &&
    grepl(if (isTRUE(allow_zero)) "^(0|[1-9][0-9]*)$" else
      "^[1-9][0-9]*$", value)
}

.dsvert_dp_count_client_rational_leq_v1 <- function(
    numerator, denominator, bound) {
  openssl::bignum(numerator) * bound$denominator <=
    bound$numerator * openssl::bignum(denominator)
}

.dsvert_dp_count_client_plan_v1 <- function(certificate, config = NULL) {
  fields <- c(
    "version", "sampler", "bernoulli_bits", "stop_numerator",
    "max_geometric_steps", "sensitivity_steps", "coordinate_count",
    "epsilon_effective_upper_numerator",
    "epsilon_effective_upper_denominator",
    "implementation_delta_numerator", "implementation_delta_denominator",
    "implementation_delta_bound", "accounting", "bernoulli_trials",
    "aes_blocks")
  scalar_integer <- function(value, expected = NULL) {
    valid <- is.numeric(value) && length(value) == 1L && !is.na(value) &&
      is.finite(value) && value == floor(value)
    isTRUE(valid) &&
      (is.null(expected) || identical(as.numeric(value), expected))
  }
  if (!is.list(certificate) || is.null(names(certificate)) ||
      anyNA(names(certificate)) || anyDuplicated(names(certificate)) ||
      !setequal(names(certificate), fields) ||
      !identical(certificate$version,
                 "dsvert-joint-dp-laplace-plan-v2") ||
      !identical(certificate$sampler,
                 .DSVERT_DP_ANALYSIS_COUNT_TV_SAMPLER) ||
      !scalar_integer(certificate$bernoulli_bits, 8) ||
      !.dsvert_dp_count_client_integer_text_v1(
        certificate$stop_numerator) ||
      openssl::bignum(certificate$stop_numerator) >= openssl::bignum(2)^8L ||
      !.dsvert_dp_count_client_integer_text_v1(
        certificate$sensitivity_steps) ||
      !identical(certificate$sensitivity_steps, "1") ||
      !scalar_integer(certificate$coordinate_count, 1) ||
      !.dsvert_dp_count_client_integer_text_v1(
        certificate$epsilon_effective_upper_numerator, allow_zero = TRUE) ||
      !.dsvert_dp_count_client_integer_text_v1(
        certificate$epsilon_effective_upper_denominator) ||
      !.dsvert_dp_count_client_integer_text_v1(
        certificate$implementation_delta_numerator) ||
      !.dsvert_dp_count_client_integer_text_v1(
        certificate$implementation_delta_denominator) ||
      !is.character(certificate$implementation_delta_bound) ||
      length(certificate$implementation_delta_bound) != 1L ||
      is.na(certificate$implementation_delta_bound) ||
      !identical(certificate$implementation_delta_bound, paste0(
        certificate$implementation_delta_numerator, "/",
        certificate$implementation_delta_denominator)) ||
      !is.character(certificate$accounting) ||
      length(certificate$accounting) != 1L || is.na(certificate$accounting) ||
      !nzchar(certificate$accounting)) {
    stop("The Count planner returned an invalid exact certificate.",
         call. = FALSE)
  }
  for (field in c(
      "max_geometric_steps", "bernoulli_trials", "aes_blocks")) {
    .dsvert_dp_count_client_positive_integer_v1(
      certificate[[field]], paste("planner", field), 4096^2)
  }
  if (certificate$max_geometric_steps > 4096) {
    stop("The Count planner returned an invalid exact certificate.",
         call. = FALSE)
  }
  if (!is.null(config)) {
    config <- .dsvert_dp_count_client_config_v1(config)
    epsilon <- .dsvert_dp_count_client_decimal_rational_v1(
      .dsvert_dp_count_client_decimal_text_v1(config$privacy$epsilon))
    implementation <- .dsvert_dp_count_client_decimal_rational_v1(
      .dsvert_dp_count_client_decimal_text_v1(
        config$calibration$implementation_delta))
    total <- .dsvert_dp_count_client_decimal_rational_v1(
      .dsvert_dp_count_client_decimal_text_v1(config$privacy$delta))
    if (!.dsvert_dp_count_client_rational_leq_v1(
        certificate$epsilon_effective_upper_numerator,
        certificate$epsilon_effective_upper_denominator, epsilon)) {
      stop("The Count effective epsilon certificate exceeds epsilon.",
           call. = FALSE)
    }
    if (!.dsvert_dp_count_client_rational_leq_v1(
        certificate$implementation_delta_numerator,
        certificate$implementation_delta_denominator, implementation)) {
      stop("The Count implementation delta certificate exceeds calibration.",
           call. = FALSE)
    }
    if (!.dsvert_dp_count_client_rational_leq_v1(
        certificate$implementation_delta_numerator,
        certificate$implementation_delta_denominator, total)) {
      stop("The Count total delta certificate exceeds privacy delta.",
           call. = FALSE)
    }
  }
  .dsvert_dp_analysis_client_canonical_value_v1(certificate)
}

.dsvert_dp_count_client_config_v1 <- function(config) {
  fields <- c(
    "version", "domain", "cohort_id", "dataset_id", "dataset_version",
    "privacy_unit_column", "alignment_purpose", "count_upper_bound",
    "max_records_per_unit", "overflow_policy", "privacy", "calibration",
    "peer_pins", "backend_build_sha256", "transport_chunk_coordinates")
  if (!is.list(config) || is.null(names(config)) || anyNA(names(config)) ||
      anyDuplicated(names(config)) || !setequal(names(config), fields) ||
      !identical(config$version, .DSVERT_CLIENT_DP_COUNT_CONFIG_VERSION)) {
    stop("Invalid Count configuration.", call. = FALSE)
  }
  identifiers <- c(
    "domain", "cohort_id", "dataset_id", "dataset_version",
    "alignment_purpose")
  if (any(!vapply(config[identifiers],
                  .dsvert_dp_analysis_client_scalar_id, logical(1L)))) {
    stop("Invalid Count configuration.", call. = FALSE)
  }
  if (!is.character(config$privacy_unit_column) ||
      length(config$privacy_unit_column) != 1L ||
      is.na(config$privacy_unit_column) ||
      !grepl("^[A-Za-z._][A-Za-z0-9._]{0,127}$",
             config$privacy_unit_column)) {
    stop("Invalid Count privacy-unit column.", call. = FALSE)
  }
  records <- .dsvert_dp_count_client_positive_integer_v1(
    config$max_records_per_unit, "maximum records per privacy unit",
    .Machine$integer.max)
  upper <- .dsvert_dp_count_client_positive_integer_v1(
    config$count_upper_bound, "upper bound", 1000000)
  if (!identical(records, 1) ||
      !identical(config$overflow_policy, "reject_operation")) {
    stop("Count requires one aligned record per privacy unit.",
         call. = FALSE)
  }
  scalar <- function(value) {
    is.numeric(value) && length(value) == 1L && !is.na(value) &&
      is.finite(value)
  }
  privacy <- config$privacy
  calibration <- config$calibration
  if (!is.list(privacy) ||
      !identical(sort(names(privacy)), c("delta", "epsilon")) ||
      !is.list(calibration) ||
      !identical(names(calibration), "implementation_delta") ||
      !scalar(privacy$epsilon) || privacy$epsilon <= 0 ||
      privacy$epsilon > 8 ||
      !scalar(privacy$delta) || privacy$delta <= 0 || privacy$delta >= 1 ||
      !scalar(calibration$implementation_delta) ||
      calibration$implementation_delta <= 0) {
    stop("Invalid Count privacy parameters.", call. = FALSE)
  }
  implementation <- .dsvert_dp_count_client_decimal_rational_v1(
    .dsvert_dp_count_client_decimal_text_v1(
      calibration$implementation_delta))
  total <- .dsvert_dp_count_client_decimal_rational_v1(
    .dsvert_dp_count_client_decimal_text_v1(privacy$delta))
  if (!.dsvert_dp_count_client_rational_leq_v1(
      as.character(implementation$numerator),
      as.character(implementation$denominator), total)) {
    stop("Count implementation delta exceeds the total privacy delta.",
         call. = FALSE)
  }
  pins <- config$peer_pins
  if (!is.character(pins) || length(pins) < 2L || length(pins) > 4096L ||
      is.null(names(pins)) || anyNA(names(pins)) ||
      anyDuplicated(names(pins)) || any(!nzchar(names(pins))) ||
      any(!vapply(names(pins), function(peer) {
        tryCatch({
          .dsvert_dp_count_client_peer_name_v1(peer)
          TRUE
        }, error = function(error) FALSE)
      }, logical(1L)))) {
    stop("Invalid Count peer pins.", call. = FALSE)
  }
  normalized <- tryCatch(vapply(
    pins, .dsvert_dp_analysis_client_identity_pk, character(1L)),
    error = function(error) character())
  if (length(normalized) != length(pins) || anyDuplicated(normalized)) {
    stop("Invalid Count peer pins.", call. = FALSE)
  }
  names(normalized) <- names(pins)
  normalized <- normalized[order(names(normalized), method = "radix")]
  .dsvert_dp_count_client_hex_v1(
    config$backend_build_sha256, "backend build digest")
  chunk <- .dsvert_dp_count_client_positive_integer_v1(
    config$transport_chunk_coordinates, "transport chunk size",
    .Machine$integer.max)

  result <- config
  result$count_upper_bound <- upper
  result$max_records_per_unit <- records
  result$privacy <- list(
    delta = as.numeric(privacy$delta), epsilon = as.numeric(privacy$epsilon))
  result$calibration <- list(
    implementation_delta = as.numeric(calibration$implementation_delta))
  result$peer_pins <- normalized
  result$transport_chunk_coordinates <- chunk
  result[sort(names(result), method = "radix")]
}

.dsvert_dp_count_client_config_hash_v1 <- function(config) {
  config <- .dsvert_dp_count_client_config_v1(config)
  config$peer_pins <- as.list(config$peer_pins)
  .dsvert_dp_count_client_hash_v1(
    .DSVERT_CLIENT_DP_COUNT_CONFIG_DOMAIN, config)
}

.dsvert_dp_count_client_receipt_message_v1 <- function(unsigned) {
  charToRaw(paste0(
    .DSVERT_CLIENT_DP_COUNT_RECEIPT_DOMAIN,
    .dsvert_joint_dp_client_json(unsigned)))
}

.dsvert_dp_count_client_receipt_v1 <- function(receipt, config) {
  fields <- c(
    "version", "peer_name", "peer_identity_pk", "config_sha256",
    "psi_run_sha256", "snapshot_commitment", "sampler_plan", "signature")
  if (!is.list(receipt) || is.null(names(receipt)) || anyNA(names(receipt)) ||
      anyDuplicated(names(receipt)) || !setequal(names(receipt), fields) ||
      !identical(receipt$version, .DSVERT_CLIENT_DP_COUNT_RECEIPT_VERSION)) {
    stop("Invalid signed Count receipt fields.", call. = FALSE)
  }
  config <- .dsvert_dp_count_client_config_v1(config)
  unsigned <- receipt[setdiff(names(receipt), "signature")]
  unsigned$peer_name <- .dsvert_dp_count_client_peer_name_v1(
    unsigned$peer_name)
  unsigned$peer_identity_pk <- tryCatch(
    .dsvert_dp_analysis_client_identity_pk(unsigned$peer_identity_pk),
    error = function(error) stop("Invalid Count receipt identity.",
                                 call. = FALSE))
  for (field in c(
      "config_sha256", "psi_run_sha256", "snapshot_commitment")) {
    unsigned[[field]] <- .dsvert_dp_count_client_hex_v1(
      unsigned[[field]], paste("receipt", field))
  }
  unsigned$sampler_plan <- .dsvert_dp_count_client_plan_v1(
    unsigned$sampler_plan, config)
  unsigned <- .dsvert_dp_analysis_client_canonical_value_v1(unsigned)
  if (!identical(unsigned$config_sha256,
                 .dsvert_dp_count_client_config_hash_v1(config))) {
    stop("The Count receipt targets a different configuration.",
         call. = FALSE)
  }
  if (!unsigned$peer_name %in% names(config$peer_pins) ||
      !identical(unsigned$peer_identity_pk,
                 unname(config$peer_pins[[unsigned$peer_name]]))) {
    stop("The Count receipt identity is not pinned.", call. = FALSE)
  }
  signature <- .dsvert_joint_dp_client_b64url(
    receipt$signature, 64L, "Count receipt signature")
  public <- .dsvert_joint_dp_client_b64url(
    unsigned$peer_identity_pk, 32L, "Count receipt identity public key")
  verified <- tryCatch(openssl::ed25519_verify(
    .dsvert_dp_count_client_receipt_message_v1(unsigned), signature,
    openssl::read_ed25519_pubkey(public)), error = function(error) FALSE)
  if (!isTRUE(verified)) {
    stop("Count receipt signature verification failed.", call. = FALSE)
  }
  .dsvert_dp_analysis_client_canonical_value_v1(
    c(unsigned, list(signature = receipt$signature)))
}

.dsvert_dp_count_client_contract_v1 <- function(config, receipts) {
  config <- .dsvert_dp_count_client_config_v1(config)
  plan <- receipts[[1L]]$sampler_plan
  owner_snapshots <- stats::setNames(lapply(receipts, function(receipt) {
    list(
      version = .DSVERT_DP_ANALYSIS_SNAPSHOT_VERSION,
      dataset_id = config$dataset_id,
      dataset_version = config$dataset_version,
      snapshot_commitment = receipt$snapshot_commitment)
  }), vapply(receipts, `[[`, character(1L), "peer_identity_pk"))
  owner_snapshots <- owner_snapshots[
    order(names(owner_snapshots), method = "radix")]
  identities <- sort(names(owner_snapshots), method = "radix")
  constraints <- list(
    version = "dsvert-contribution-constraints-v1",
    policy_sha256 = .dsvert_dp_count_client_hash_v1(
      "dsVert/dp-count/contribution/v1|", list(
        max_records_per_unit = config$max_records_per_unit,
        overflow_policy = config$overflow_policy)))
  semantic <- list(
    version = .DSVERT_DP_ANALYSIS_SEMANTIC_VERSION,
    domain = config$domain,
    cohort_id = config$cohort_id,
    owner_snapshots = owner_snapshots,
    noise_authorities = as.list(unname(identities[1:2])),
    analysis = list(
      primitive = "joint-dp-laplace-v2",
      formula = NULL,
      effective_arguments = list(
        statistic = "aligned_privacy_unit_count",
        owner_combination = "vertical_membership_once_v1",
        count_bounds = list(lower = 0, upper = config$count_upper_bound),
        sampler_plan = plan)),
    privacy = list(
      version = "dsvert-per-analysis-dp-v1",
      adjacency = "add_remove_patient",
      privacy_unit = "patient",
      contribution = list(
        version = "dsvert-contribution-policy-v1",
        max_records_per_unit = config$max_records_per_unit,
        overflow_policy = config$overflow_policy,
        constraints = constraints),
      mechanism = list(
        family = "discrete_laplace",
        version = .DSVERT_DP_ANALYSIS_COUNT_TV_MECHANISM,
        sensitivity = list(
          version = "dsvert-sensitivity-v1", norm = "l1", value = 1),
        calibration = list(
          version = "dsvert-calibration-v1",
          noise_scale = 1 / config$privacy$epsilon,
          sampler = .DSVERT_DP_ANALYSIS_COUNT_TV_SAMPLER,
          implementation_delta =
            config$calibration$implementation_delta),
        randomness = list(
          version = "dsvert-randomness-plan-v1",
          lanes = list(final_noise = list(
            version = "dsvert-randomness-lane-v1",
            purpose = "privatize_final_vector",
            primitive = .DSVERT_DP_ANALYSIS_COUNT_TV_SAMPLER,
            coordinates = 1)))),
      epsilon = config$privacy$epsilon,
      delta = config$privacy$delta),
    numeric = list(
      version = "dsvert-numeric-semantics-v1",
      value_bits = 127,
      fractional_bits = 0,
      rounding = "toward_zero",
      overflow = "reject",
      sampler_encoding = "aes128ctr_integer_coordinate_v2",
      output_encoding = "twos_complement_integer_v1"),
    public_shape = list(count = 1))
  semantic <- .dsvert_dp_analysis_semantic_validate_v1(semantic)
  contract <- list(
    version = .DSVERT_DP_ANALYSIS_CONTRACT_VERSION,
    artifact_key = .dsvert_dp_analysis_artifact_key_v1(semantic),
    semantic = semantic,
    execution = list(
      version = .DSVERT_DP_ANALYSIS_EXECUTION_VERSION,
      peer_pins = as.list(config$peer_pins),
      backend = list(
        kernel = "joint-dp-laplace-v2",
        ring = "ring127",
        build_sha256 = config$backend_build_sha256),
      transport = list(
        chunk_coordinates = config$transport_chunk_coordinates)))
  .dsvert_dp_analysis_contract_validate_v1(contract)
}

.dsvert_dp_count_client_hex_raw_v1 <- function(value, what) {
  value <- .dsvert_dp_count_client_hex_v1(value, what)
  as.raw(strtoi(
    substring(value, seq.int(1L, 63L, 2L), seq.int(2L, 64L, 2L)),
    base = 16L))
}

.dsvert_dp_count_client_commitment_context_v1 <- function(
    transcript, role, peer_id) {
  if (!role %in% c("garbler", "evaluator") ||
      !is.character(peer_id) || length(peer_id) != 1L || is.na(peer_id) ||
      !grepl(.DSVERT_CLIENT_EXACT_GC_PEER_RE, peer_id)) {
    stop("Invalid Count worker commitment role.", call. = FALSE)
  }
  digest::digest(c(
    charToRaw("dsvert-joint-dp-private-seed-commitment-v2"), as.raw(0L),
    .dsvert_dp_count_client_hex_raw_v1(transcript, "transcript hash"),
    as.raw(0L), charToRaw(role), as.raw(0L), charToRaw(peer_id)),
    algo = "sha256", serialize = FALSE)
}

.dsvert_dp_count_client_worker_static_v1 <- function(contract) {
  contract <- .dsvert_dp_analysis_contract_validate_v1(contract)
  binding <- .dsvert_exact_gc_analysis_binding(contract)
  semantic <- contract$semantic
  arguments <- semantic$analysis$effective_arguments
  bounds <- arguments$count_bounds
  if (!is.list(arguments) || is.null(names(arguments)) ||
      !setequal(names(arguments), c(
        "statistic", "owner_combination", "count_bounds", "sampler_plan")) ||
      !identical(arguments$statistic, "aligned_privacy_unit_count") ||
      !identical(arguments$owner_combination,
                 "vertical_membership_once_v1") ||
      !is.list(bounds) ||
      !identical(sort(names(bounds)), c("lower", "upper")) ||
      !is.numeric(bounds$lower) || length(bounds$lower) != 1L ||
      is.na(bounds$lower) || !is.finite(bounds$lower) ||
      !identical(as.numeric(bounds$lower), 0)) {
    stop("Invalid Count worker semantic contract.", call. = FALSE)
  }
  upper <- .dsvert_dp_count_client_positive_integer_v1(
    bounds$upper, "worker upper bound", 1000000)
  certificate <- .dsvert_dp_count_client_plan_v1(arguments$sampler_plan)
  mechanism <- semantic$privacy$mechanism
  roles <- binding$binding$authority_roles
  peer_ids <- stats::setNames(vapply(
    roles, .dsvert_exact_gc_identity_peer_id, character(1L)), names(roles))
  .dsvert_dp_analysis_client_canonical_value_v1(list(
    version = .DSVERT_CLIENT_DP_COUNT_WORKER_VERSION,
    ring_bits = 127L,
    frac_bits = 0L,
    coordinate_count = 1L,
    sampler = mechanism$calibration$sampler,
    epsilon = .dsvert_dp_count_client_decimal_text_v1(
      semantic$privacy$epsilon),
    allocated_delta = .dsvert_dp_count_client_decimal_text_v1(
      mechanism$calibration$implementation_delta),
    sensitivity_steps = "1",
    bernoulli_bits = as.integer(certificate$bernoulli_bits),
    stop_numerator = certificate$stop_numerator,
    max_geometric_steps = as.integer(certificate$max_geometric_steps),
    implementation_delta_numerator =
      certificate$implementation_delta_numerator,
    implementation_delta_denominator =
      certificate$implementation_delta_denominator,
    encoded_lower = "0",
    encoded_upper = as.character(as.integer(upper)),
    transcript_hash = binding$sha256,
    garbler_commitment_context =
      .dsvert_dp_count_client_commitment_context_v1(
        binding$sha256, "garbler", peer_ids[["garbler"]]),
    evaluator_commitment_context =
      .dsvert_dp_count_client_commitment_context_v1(
        binding$sha256, "evaluator", peer_ids[["evaluator"]])))
}

.dsvert_dp_count_client_compile_v1 <- function(envelopes, server_names) {
  if (!is.character(server_names) || length(server_names) < 2L ||
      length(server_names) > 4096L || anyNA(server_names) ||
      any(!nzchar(server_names)) || anyDuplicated(server_names) ||
      !is.list(envelopes) || is.null(names(envelopes)) ||
      anyNA(names(envelopes)) || anyDuplicated(names(envelopes)) ||
      !setequal(names(envelopes), server_names)) {
    stop("Count compile results do not cover the complete federation.",
         call. = FALSE)
  }
  envelopes <- envelopes[server_names]
  configs <- lapply(envelopes, function(envelope) {
    if (!is.list(envelope) || is.null(names(envelope)) ||
        anyNA(names(envelope)) || anyDuplicated(names(envelope)) ||
        !setequal(names(envelope), c("config", "receipt"))) {
      stop("A peer returned an invalid Count compile envelope.",
           call. = FALSE)
    }
    .dsvert_dp_count_client_config_v1(envelope$config)
  })
  reference <- configs[[1L]]
  if (!all(vapply(configs, identical, logical(1L), reference)) ||
      !setequal(names(reference$peer_pins), server_names)) {
    stop("Count peers disagree on one full-K configuration.",
         call. = FALSE)
  }
  receipts <- Map(function(envelope, site) {
    receipt <- .dsvert_dp_count_client_receipt_v1(
      envelope$receipt, reference)
    if (!identical(receipt$peer_name, site)) {
      stop("Count compile receipts do not match their connected peers.",
           call. = FALSE)
    }
    receipt
  }, envelopes, server_names)
  peers <- vapply(receipts, `[[`, character(1L), "peer_name")
  if (anyDuplicated(peers) || !setequal(peers, names(reference$peer_pins))) {
    stop("Invalid Count receipt coverage.", call. = FALSE)
  }
  names(receipts) <- peers
  receipts <- receipts[names(reference$peer_pins)]
  if (length(unique(vapply(
      receipts, `[[`, character(1L), "psi_run_sha256"))) != 1L) {
    stop("Count receipts disagree on the PSI run.", call. = FALSE)
  }
  plan <- receipts[[1L]]$sampler_plan
  if (!all(vapply(receipts, function(receipt) {
    identical(receipt$sampler_plan, plan)
  }, logical(1L)))) {
    stop("Count receipts disagree on the sampler certificate.",
         call. = FALSE)
  }
  contract <- .dsvert_dp_count_client_contract_v1(reference, receipts)
  binding <- .dsvert_exact_gc_analysis_binding(contract)
  worker <- .dsvert_dp_count_client_worker_static_v1(contract)
  roles <- unlist(
    binding$binding$authority_roles[c("garbler", "evaluator")],
    use.names = TRUE)
  authorities <- names(reference$peer_pins)[match(
    unname(roles), unname(reference$peer_pins))]
  if (length(authorities) != 2L || anyNA(authorities) ||
      anyDuplicated(authorities)) {
    stop("Count noise authorities do not match the full-K pinset.",
         call. = FALSE)
  }
  list(
    config = reference,
    receipts = receipts,
    contract = contract,
    binding = binding,
    worker_static = worker,
    authorities = unname(authorities),
    config_sha256 = .dsvert_dp_count_client_config_hash_v1(reference),
    receipt_set_sha256 = .dsvert_dp_count_client_hash_v1(
      .DSVERT_CLIENT_DP_COUNT_RECEIPT_SET_DOMAIN, receipts),
    psi_run_sha256 = receipts[[1L]]$psi_run_sha256,
    contract_sha256 = .dsvert_dp_count_client_hash_v1(
      .DSVERT_CLIENT_DP_COUNT_CONTRACT_DOMAIN, contract),
    worker_static_sha256 = .dsvert_dp_count_client_hash_v1(
      .DSVERT_CLIENT_DP_COUNT_WORKER_DOMAIN, worker))
}

.dsvert_dp_count_client_public_authorization_v1 <- function(
    value, peer, session_id, compiled) {
  fields <- c(
    "version", "session_id", "artifact_key", "config_sha256",
    "receipt_set_sha256", "psi_run_sha256", "contract_sha256",
    "analysis_binding_sha256", "worker_static_sha256", "local_authority",
    "commitment_context", "seed_commitment", "authorization_sha256",
    "signature")
  if (!is.list(value) || is.null(names(value)) || anyNA(names(value)) ||
      anyDuplicated(names(value)) || !setequal(names(value), fields)) {
    stop("A Count authority returned an invalid public authorization.",
         call. = FALSE)
  }
  peer <- .dsvert_dp_count_client_peer_name_v1(peer)
  expected_role <- names(compiled$binding$binding$authority_roles)[match(
    unname(compiled$config$peer_pins[[peer]]),
    unlist(compiled$binding$binding$authority_roles, use.names = FALSE))]
  authority <- value$local_authority
  expected_authority <- list(
    peer_name = peer,
    identity_pk = unname(compiled$config$peer_pins[[peer]]),
    role = unname(expected_role))
  expected_context <- compiled$worker_static[[paste0(
    expected_role, "_commitment_context")]]
  valid <- identical(value$version,
                     .DSVERT_CLIENT_DP_COUNT_PUBLIC_AUTH_VERSION) &&
    identical(value$session_id, session_id) &&
    identical(value$artifact_key, compiled$contract$artifact_key) &&
    identical(value$config_sha256, compiled$config_sha256) &&
    identical(value$receipt_set_sha256, compiled$receipt_set_sha256) &&
    identical(value$psi_run_sha256, compiled$psi_run_sha256) &&
    identical(value$contract_sha256, compiled$contract_sha256) &&
    identical(value$analysis_binding_sha256, compiled$binding$sha256) &&
    identical(value$worker_static_sha256,
              compiled$worker_static_sha256) &&
    is.list(authority) && !is.null(names(authority)) &&
    setequal(names(authority), names(expected_authority)) &&
    identical(authority[names(expected_authority)], expected_authority) &&
    identical(value$commitment_context, expected_context)
  if (!isTRUE(valid)) {
    stop("A Count authority returned a misbound public authorization.",
         call. = FALSE)
  }
  for (field in c(
      "artifact_key", "config_sha256", "receipt_set_sha256",
      "psi_run_sha256", "contract_sha256", "analysis_binding_sha256",
      "worker_static_sha256", "commitment_context", "seed_commitment",
      "authorization_sha256")) {
    .dsvert_dp_count_client_hex_v1(value[[field]],
                                   paste("authorization", field))
  }
  core <- value[setdiff(names(value), c("authorization_sha256", "signature"))]
  expected_hash <- .dsvert_dp_count_client_hash_v1(
    .DSVERT_CLIENT_DP_COUNT_PUBLIC_AUTH_DOMAIN, core)
  if (!identical(value$authorization_sha256, expected_hash)) {
    stop("A Count authority returned a corrupt public authorization.",
         call. = FALSE)
  }
  signed <- c(core, list(authorization_sha256 = expected_hash))
  signature <- .dsvert_joint_dp_client_b64url(
    value$signature, 64L, "Count authorization signature")
  public <- .dsvert_joint_dp_client_b64url(
    expected_authority$identity_pk, 32L,
    "Count authorization identity public key")
  verified <- tryCatch(openssl::ed25519_verify(
    charToRaw(paste0(
      .DSVERT_CLIENT_DP_COUNT_PUBLIC_AUTH_SIGNATURE_DOMAIN,
      .dsvert_joint_dp_client_json(signed))),
    signature, openssl::read_ed25519_pubkey(public)),
    error = function(error) FALSE)
  if (!isTRUE(verified)) {
    stop("Count public authorization signature verification failed.",
         call. = FALSE)
  }
  .dsvert_dp_analysis_client_canonical_value_v1(value)
}

.dsvert_dp_count_client_authorization_set_v1 <- function(
    values, session_id, compiled) {
  authorities <- compiled$authorities
  if (!is.list(values) || is.null(names(values)) || anyNA(names(values)) ||
      anyDuplicated(names(values)) || !setequal(names(values), authorities)) {
    stop("Count public authorizations do not cover both noise authorities.",
         call. = FALSE)
  }
  values <- values[authorities]
  verified <- Map(function(value, peer) {
    .dsvert_dp_count_client_public_authorization_v1(
      value, peer, session_id, compiled)
  }, values, authorities)
  names(verified) <- authorities
  roles <- vapply(verified, function(value) {
    value$local_authority$role
  }, character(1L))
  commitments <- vapply(verified, `[[`, character(1L), "seed_commitment")
  if (!identical(unname(roles), c("garbler", "evaluator")) ||
      anyDuplicated(commitments)) {
    stop("Count public authorizations have invalid authority coverage.",
         call. = FALSE)
  }
  verified
}

.dsvert_dp_count_compile_authorize_bind_v1 <- function(
    data_name, datasources,
    .aggregate = DSI::datashield.aggregate,
    .session_id = .dsvert_uuid4,
    .setup_exact = .dsvert_setup_exact_gc_transport) {
  if (!is.character(data_name) || length(data_name) != 1L ||
      is.na(data_name) ||
      !grepl("^[a-zA-Z._][a-zA-Z0-9._]*$", data_name)) {
    stop("Invalid Count protected data name.", call. = FALSE)
  }
  if (!is.list(datasources) || length(datasources) < 2L ||
      is.null(names(datasources)) || anyNA(names(datasources)) ||
      any(!nzchar(names(datasources))) || anyDuplicated(names(datasources))) {
    stop("Count requires the complete named K>=2 datasource set.",
         call. = FALSE)
  }
  if (!is.function(.aggregate) || !is.function(.session_id) ||
      !is.function(.setup_exact)) {
    stop("Invalid Count client orchestration dependency.", call. = FALSE)
  }
  server_names <- names(datasources)
  envelopes <- .dsvert_aggregate_strict(
    datasources,
    call(name = "dsvertDPCountCompileDS", data_name = data_name),
    operation = "Count signed analysis compilation",
    .aggregate = .aggregate)
  compiled <- .dsvert_dp_count_client_compile_v1(envelopes, server_names)
  session_id <- .session_id()
  if (!is.character(session_id) || length(session_id) != 1L ||
      is.na(session_id) || !grepl(
        "^[0-9a-f]{8}-[0-9a-f]{4}-4[0-9a-f]{3}-[89ab][0-9a-f]{3}-[0-9a-f]{12}$",
        session_id)) {
    stop("Invalid Count client session identifier.", call. = FALSE)
  }
  authorities <- compiled$authorities
  authority_conns <- datasources[authorities]
  wire_config <- compiled$config
  wire_config$peer_pins <- as.list(wire_config$peer_pins)
  config_json <- .dsvert_joint_dp_client_json(wire_config)
  receipts_json <- .dsvert_joint_dp_client_json(
    unname(compiled$receipts))
  if (nchar(config_json, type = "bytes") >
        .DSVERT_CLIENT_DP_COUNT_CONFIG_MAX_BYTES ||
      nchar(receipts_json, type = "bytes") >
        .DSVERT_CLIENT_DP_COUNT_RECEIPTS_MAX_BYTES) {
    stop("Count authorization input exceeds its bounded wire contract.",
         call. = FALSE)
  }
  authorize_call <- call(
    name = "dsvertDPCountAuthorizeDS",
    config_json = .dsvert_dsi_text_encode(
      config_json, "Count authorization configuration"),
    receipts_json = .dsvert_dsi_text_encode(
      receipts_json, "Count authorization receipts"),
    session_id = session_id)
  authorized <- .dsvert_aggregate_strict(
    authority_conns, authorize_call,
    operation = "Count two-authority authorization",
    .aggregate = .aggregate)
  authorized <- .dsvert_dp_count_client_authorization_set_v1(
    authorized, session_id, compiled)
  selected <- match(authorities, server_names)
  transport <- .setup_exact(
    datasources = datasources,
    server_names = server_names,
    servers = selected,
    session_id = session_id,
    cleanup_purpose = .DSVERT_CLIENT_EXACT_GC_CROSS_CLEANUP_PURPOSE,
    analysis_contract = compiled$contract,
    .aggregate = .aggregate)
  list(
    session_id = session_id,
    contract = compiled$contract,
    authorities = authorities,
    authorizations = authorized,
    transport = transport)
}
