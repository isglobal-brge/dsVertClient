.client_analysis_identity_pk <- function(index) {
  encoded <- jsonlite::base64_enc(as.raw(rep(as.integer(index), 32L)))
  chartr("+/", "-_", sub("=+$", "", encoded))
}

.client_analysis_contract_fixture <- function(k = 2L) {
  owners <- paste0("site_", seq_len(k))
  pins <- setNames(vapply(
    seq_along(owners), .client_analysis_identity_pk, character(1L)),
    owners)
  list(
    version = "dsvert-analysis-contract-v1",
    artifact_key = strrep("0", 64L),
    semantic = list(
      version = "dsvert-analysis-semantic-v1",
      domain = "study-domain",
      cohort_id = "cohort-v1",
      owner_snapshots = setNames(lapply(seq_along(owners), function(index) {
        list(
          version = "dsvert-analysis-snapshot-v1",
          dataset_id = "cohort_table",
          dataset_version = "v1",
          snapshot_commitment = strrep(sprintf("%x", index), 64L))
      }), unname(pins)),
      noise_authorities = unname(pins[seq_len(2L)]),
      analysis = list(
        primitive = "glm-binomial-logit-v1",
        formula = list(
          response = "outcome", intercept = TRUE,
          terms = list("age", "treatment")),
        effective_arguments = list(
          link = "logit", missing = "complete_case")),
      privacy = list(
        version = "dsvert-per-analysis-dp-v1",
        adjacency = "add_remove_patient",
        privacy_unit = "patient",
        contribution = list(
          version = "dsvert-contribution-policy-v1",
          max_records_per_unit = 1,
          overflow_policy = "reject_operation",
          constraints = list(
            version = "dsvert-contribution-constraints-v1",
            policy_sha256 = strrep("c", 64L))),
        mechanism = list(
          family = "gaussian",
          version = "gaussian-output-perturbation-v1",
          sensitivity = list(
            version = "dsvert-sensitivity-v1",
            norm = "l2",
            value = 1),
          calibration = list(
            version = "dsvert-calibration-v1",
            noise_scale = 5,
            sampler = "gaussian-one-draw-v1",
            implementation_delta = 1e-9),
          randomness = list(
            version = "dsvert-randomness-plan-v1",
            lanes = list(
              final_noise = list(
                version = "dsvert-randomness-lane-v1",
                purpose = "privatize_final_vector",
                primitive = "gaussian-one-draw-v1",
                coordinates = 3)))),
        epsilon = 1,
        delta = 1e-6),
      numeric = list(
        version = "dsvert-numeric-semantics-v1",
        value_bits = 120,
        fractional_bits = 32,
        rounding = "nearest_even",
        overflow = "reject",
        sampler_encoding = "chacha20_absolute_coordinate_v1",
        output_encoding = "fixed_point_v1"),
      public_shape = list(coefficients = 3, covariance = c(3, 3))),
    execution = list(
      version = "dsvert-analysis-execution-v1",
      peer_pins = as.list(pins),
      backend = list(
        kernel = "glm-binomial-logit-v1",
        ring = "ring127",
        build_sha256 = strrep("a", 64L)),
      transport = list(chunk_coordinates = 4096)))
}

.client_analysis_count_contract_fixture <- function(k = 2L) {
  contract <- .client_analysis_contract_fixture(k)
  contract$semantic$analysis$primitive <- "joint-dp-laplace-v2"
  contract$semantic$analysis["formula"] <- list(NULL)
  contract$semantic$analysis$effective_arguments <- list(
    statistic = "admitted_privacy_unit_count")
  contract$semantic$privacy$mechanism <- list(
    family = "discrete_laplace",
    version = "discrete-laplace-output-perturbation-tv-v2",
    sensitivity = list(
      version = "dsvert-sensitivity-v1", norm = "l1", value = 1),
    calibration = list(
      version = "dsvert-calibration-v1",
      noise_scale = 1,
      sampler = "hkdf-sha256-aes128ctr-two-geometric-tv-v2",
      implementation_delta = 1e-9),
    randomness = list(
      version = "dsvert-randomness-plan-v1",
      lanes = list(
        final_noise = list(
          version = "dsvert-randomness-lane-v1",
          purpose = "privatize_final_vector",
          primitive = "hkdf-sha256-aes128ctr-two-geometric-tv-v2",
          coordinates = 1))))
  contract$semantic$privacy$epsilon <- 1
  contract$semantic$privacy$delta <- 1e-6
  contract$semantic$numeric <- list(
    version = "dsvert-numeric-semantics-v1",
    value_bits = 127,
    fractional_bits = 0,
    rounding = "toward_zero",
    overflow = "reject",
    sampler_encoding = "aes128ctr_integer_coordinate_v2",
    output_encoding = "twos_complement_integer_v1")
  contract$semantic$public_shape <- list(count = 1)
  contract$execution$backend$kernel <- "joint-dp-laplace-v2"
  contract$execution$backend$ring <- "ring127"
  contract$artifact_key <-
    .dsvert_dp_analysis_artifact_key_v1(contract$semantic)
  contract
}

.client_analysis_frequency_contract_fixture <- function(
    k = 2L, profile = c("convolution", "gaussian"),
    levels = c("control", "caf\u00e9", "case"),
    chunk_coordinates = min(8192L, length(levels))) {
  profile <- match.arg(profile)
  contract <- .client_analysis_contract_fixture(k)
  owner_ids <- sort(names(contract$semantic$owner_snapshots), method = "radix")
  source_owner <- owner_ids[[min(2L, length(owner_ids))]]
  secondary <- sort(setdiff(owner_ids, source_owner), method = "radix")[[1L]]
  dimension <- length(levels)
  convolution <- identical(profile, "convolution")
  primitive <- if (convolution) {
    "independent_full_global_draw_convolution_ring128_v3"
  } else {
    paste0("independent_full_global_dyadic_discrete_gaussian_",
           "tv_bounded_ring128_v2")
  }
  registry <- .dsvert_dp_analysis_frequency_profile_v1(primitive)
  plan_version <- registry$plan
  sampler <- registry$sampler
  mechanism_version <- if (convolution) {
    "two-independent-complete-vector-discrete-laplace-draws-v3"
  } else {
    paste0("two-independent-complete-vector-dyadic-discrete-gaussian-",
           "tv-bounded-draws-v2")
  }
  full_plan_sha256 <- strrep(if (convolution) "4" else "5", 64L)
  planner_request <- .dsvert_dp_analysis_frequency_planner_request_v1(
    registry, list(epsilon = 1), list(implementation_delta = 0.001),
    list(value = 1), dimension)
  planner_request_sha256 <- .dsvert_dp_analysis_frequency_hash_v1(
    registry$request_domain, planner_request)
  allocated <- .dsvert_dp_analysis_frequency_decimal_fraction_v1(
    planner_request$delta)
  allocated_delta <- list(
    numerator = as.character(allocated$numerator),
    denominator = as.character(allocated$denominator))
  no_wrap_sha256 <- .dsvert_dp_analysis_frequency_hash_v1(
    "dsVert/frequency/ring128-no-wrap/v1|", list(
      version = "dsvert-frequency-ring128-no-wrap-v1",
      coordinate_upper_bound = "1000",
      maximum_noise_per_peer = "100",
      maximum_noise_release = "200"))
  contract$semantic$version <-
    "dsvert-analysis-semantic-fixed-categorical-vector-v1"
  contract$semantic$noise_authorities <- NULL
  contract$semantic$noise_authority_roles <- list(
    version = "dsvert-frequency-noise-authority-roles-v1",
    role_order = list("source_owner", "secondary_noise_authority"),
    authority_ids = list(source_owner, secondary))
  contract$semantic$analysis <- list(
    primitive = primitive,
    formula = NULL,
    effective_arguments = list(
      version = "dsvert-fixed-domain-categorical-frequency-v1",
      statistic = "aligned_fixed_domain_categorical_frequency",
      source_owner = source_owner,
      dataset_id = "cohort_table",
      dataset_version = "v1",
      variable_id = "status",
      levels = as.list(levels),
      dimension = dimension,
      repeated_record_policy = "consistent_level_else_exclude_v1",
      missingness_policy = "missing_or_out_of_domain_rows_are_ignored",
      coordinate_bounds = list(lower = 0, upper = 1000),
      sampler_plan = list(
        version = "dsvert-frequency-plan-summary-v1",
        physical_plan_version = plan_version,
        full_plan_sha256 = full_plan_sha256,
        planner_request_sha256 = planner_request_sha256,
        coordinate_order_sha256 =
          .dsvert_dp_analysis_frequency_coordinate_order_sha256_v1(
            as.list(levels)),
        d = dimension,
        chunk_coordinates = chunk_coordinates,
        implementation_delta = list(
          numerator = "1", denominator = "1000000"),
        allocated_delta = allocated_delta,
        core_delta = list(
          numerator = if (convolution) "0" else "1",
          denominator = if (convolution) "1" else "2000"),
        maximum_noise_per_peer = "100",
        no_wrap_sha256 = no_wrap_sha256,
        profile_sha256 = .dsvert_dp_analysis_frequency_hash_v1(
          "dsVert/frequency/physical-profile/v1|", registry),
        backend_selection = list(
          version = "dsvert-frequency-backend-selection-v1",
          policy_sha256 =
            .dsvert_dp_analysis_frequency_policy_sha256_v1(registry)))))
  contract$semantic$privacy$adjacency <- "add_remove_patient"
  contract$semantic$privacy$mechanism <- list(
    family = if (convolution) "discrete_laplace" else "gaussian",
    version = mechanism_version,
    sensitivity = list(
      version = "dsvert-sensitivity-v1",
      norm = if (convolution) "l1" else "l2",
      value = 1),
    calibration = list(
      version = "dsvert-calibration-v1",
      sampler = sampler,
      implementation_delta = 0.001),
    randomness = list(
      version = "dsvert-randomness-plan-v1",
      lanes = list(final_noise = list(
        version = "dsvert-randomness-lane-v1",
        purpose = "privatize_final_vector",
        primitive = sampler,
        coordinates = dimension))))
  contract$semantic$privacy$epsilon <- 1
  contract$semantic$privacy$delta <- 0.01
  contract$semantic$numeric <- list(
    version = "dsvert-numeric-semantics-v1",
    value_bits = 128,
    fractional_bits = 0,
    rounding = "toward_zero",
    overflow = "reject",
    output_encoding = "twos_complement_integer_v1")
  contract$semantic$public_shape <- list(counts = dimension)
  contract$execution$backend$kernel <- primitive
  contract$execution$backend$ring <- "ring128"
  contract$execution$transport$chunk_coordinates <- chunk_coordinates
  contract$artifact_key <-
    .dsvert_dp_analysis_artifact_key_v1(contract$semantic)
  contract
}

test_that("client validates fixed categorical Frequency contracts for K=2,3,5", {
  contracts <- lapply(c(2L, 3L, 5L), function(k) {
    .dsvert_dp_analysis_contract_validate_v1(
      .client_analysis_frequency_contract_fixture(k))
  })
  expect_true(all(vapply(contracts, function(contract) {
    arguments <- contract$semantic$analysis$effective_arguments
    identical(arguments$dimension, 3) &&
      identical(arguments$source_owner,
                contract$semantic$noise_authority_roles$authority_ids[[1L]]) &&
      identical(contract$semantic$privacy$mechanism$randomness$lanes$
                  final_noise$coordinates, 3)
  }, logical(1L))))
  expect_identical(
    contracts[[2L]]$artifact_key,
    "c3d3b47bbd9f8b9f9d4b0ebfd967332496df3c564ec618c54f3007d464e83e4c")

  singleton <- .dsvert_dp_analysis_contract_validate_v1(
    .client_analysis_frequency_contract_fixture(
      3L, levels = "only", chunk_coordinates = 1L))
  expect_identical(singleton$semantic$public_shape, list(counts = 1))
  expect_identical(singleton$semantic$analysis$primitive,
                   "independent_full_global_draw_convolution_ring128_v3")
  expect_identical(
    singleton$semantic$analysis$effective_arguments$sampler_plan$allocated_delta,
    list(denominator = "1000", numerator = "1"))
  expect_silent(.dsvert_dp_analysis_contract_validate_v1(
    .client_analysis_frequency_contract_fixture(
      3L, "gaussian", levels = "only", chunk_coordinates = 1L)))

  convolution_profile <- .dsvert_dp_analysis_frequency_profile_v1(
    singleton$semantic$analysis$primitive)
  expect_identical(convolution_profile$chunk_partition_version,
                   "contiguous_full_chunks_except_last_v1")
  expect_identical(
    .dsvert_dp_analysis_frequency_planner_request_v1(
      convolution_profile, list(epsilon = 1),
      list(implementation_delta = 0.001), list(value = 1), 1),
    list(epsilon = "1e+00", delta = "1e-03", sensitivity_steps = "1",
         total_coordinate_count = 1L))

  gaussian <- .dsvert_dp_analysis_contract_validate_v1(
    .client_analysis_frequency_contract_fixture(5L, "gaussian"))
  expect_identical(
    gaussian$semantic$privacy$mechanism$calibration$sampler,
    paste0("cks-target-outward-rational-dyadic-cdf-hkdf-sha256-",
           "chacha20-coordinate-domain-v2"))
  expect_identical(
    gaussian$semantic$analysis$effective_arguments$sampler_plan$allocated_delta,
    list(denominator = "100000000000000000000",
         numerator = "99999999999997161"))
  expect_identical(
    gaussian$artifact_key,
    "2864afbf9954701c31b074352cc68159457d0ca5e47752712964ed03217bb17b")
  gaussian_profile <- .dsvert_dp_analysis_frequency_profile_v1(
    gaussian$semantic$analysis$primitive)
  expect_identical(
    .dsvert_dp_analysis_frequency_planner_request_v1(
      gaussian_profile, list(epsilon = 1),
      list(implementation_delta = 0.001), list(value = 1), 3),
    list(
      epsilon = "9.9999999999997158e-01",
      delta = "9.9999999999997161e-04",
      l2_sensitivity_steps = "1.0000000000000284e+00",
      total_coordinate_count = 3L))

  original <- contracts[[2L]]
  source <- original$semantic$analysis$effective_arguments$source_owner
  expect_false(identical(
    source, sort(names(original$semantic$owner_snapshots),
                 method = "radix")[[1L]]))
  expect_identical(original$semantic$noise_authority_roles$role_order,
                   list("source_owner", "secondary_noise_authority"))

  for (kind in c("convolution", "gaussian")) {
    replace_one <- .client_analysis_frequency_contract_fixture(3L, kind)
    replace_one$semantic$privacy$adjacency <- "replace_one_fixed_cohort"
    sensitivity <- replace_one$semantic$privacy$mechanism$sensitivity
    sensitivity$value <- if (identical(kind, "gaussian")) sqrt(2) else 2
    replace_one$semantic$privacy$mechanism$sensitivity <- sensitivity
    registry <- .dsvert_dp_analysis_frequency_profile_v1(
      replace_one$semantic$analysis$primitive)
    request <- .dsvert_dp_analysis_frequency_planner_request_v1(
      registry, replace_one$semantic$privacy,
      replace_one$semantic$privacy$mechanism$calibration, sensitivity, 3)
    plan <- replace_one$semantic$analysis$effective_arguments$sampler_plan
    plan$planner_request_sha256 <- .dsvert_dp_analysis_frequency_hash_v1(
      registry$request_domain, request)
    plan$full_plan_sha256 <- strrep("d", 64L)
    replace_one$semantic$analysis$effective_arguments$sampler_plan <- plan
    replace_one$artifact_key <-
      .dsvert_dp_analysis_artifact_key_v1(replace_one$semantic)
    expect_silent(.dsvert_dp_analysis_contract_validate_v1(replace_one))
  }

  for (kind in c("convolution", "gaussian")) {
    under <- .client_analysis_frequency_contract_fixture(3L, kind)
    sensitivity <- under$semantic$privacy$mechanism$sensitivity
    sensitivity$value <- 1 - 512 * .Machine$double.eps
    under$semantic$privacy$mechanism$sensitivity <- sensitivity
    registry <- .dsvert_dp_analysis_frequency_profile_v1(
      under$semantic$analysis$primitive)
    request <- .dsvert_dp_analysis_frequency_planner_request_v1(
      registry, under$semantic$privacy,
      under$semantic$privacy$mechanism$calibration, sensitivity, 3)
    under$semantic$analysis$effective_arguments$sampler_plan$
      planner_request_sha256 <- .dsvert_dp_analysis_frequency_hash_v1(
        registry$request_domain, request)
    expect_error(.dsvert_dp_analysis_contract_validate_v1(under),
                 "Frequency", info = kind)
  }
})

test_that("client Frequency identity is canonical and fully semantic", {
  original <- .client_analysis_frequency_contract_fixture(3L)
  validated <- .dsvert_dp_analysis_contract_validate_v1(original)

  aliased <- original
  source <- aliased$semantic$analysis$effective_arguments$source_owner
  alias <- paste0(" \n", chartr("-_", "+/", source), "=\t")
  source_index <- match(source, names(aliased$semantic$owner_snapshots))
  names(aliased$semantic$owner_snapshots)[[source_index]] <- alias
  aliased$semantic$analysis$effective_arguments$source_owner <- alias
  aliased$semantic$noise_authority_roles$authority_ids[[1L]] <- alias
  pin_index <- which(vapply(
    aliased$execution$peer_pins, identical, logical(1L), source))
  aliased$execution$peer_pins[[pin_index]] <- alias
  aliased$artifact_key <-
    .dsvert_dp_analysis_artifact_key_v1(aliased$semantic)
  expect_identical(
    .dsvert_dp_analysis_contract_validate_v1(aliased)$artifact_key,
    validated$artifact_key)

  operational <- original
  operational$execution$backend$build_sha256 <- strrep("b", 64L)
  names(operational$execution$peer_pins) <- paste0(
    "connection_", seq_along(operational$execution$peer_pins))
  expect_identical(
    .dsvert_dp_analysis_contract_validate_v1(operational)$artifact_key,
    validated$artifact_key)

  expect_error(.client_analysis_frequency_contract_fixture(
    3L, chunk_coordinates = 1L), "Frequency")
  expect_error(.client_analysis_frequency_contract_fixture(
    3L, "gaussian", chunk_coordinates = 2L), "Frequency")

  replanned <- original
  replanned$semantic$analysis$effective_arguments$sampler_plan$
    full_plan_sha256 <- strrep("e", 64L)
  replanned$artifact_key <-
    .dsvert_dp_analysis_artifact_key_v1(replanned$semantic)
  expect_false(identical(replanned$artifact_key, validated$artifact_key))

  looser_total_delta <- original
  looser_total_delta$semantic$privacy$delta <- 0.02
  looser_total_delta$artifact_key <- .dsvert_dp_analysis_artifact_key_v1(
    looser_total_delta$semantic)
  expect_silent(
    .dsvert_dp_analysis_contract_validate_v1(looser_total_delta))
  expect_false(identical(
    looser_total_delta$artifact_key, validated$artifact_key))

  nfc <- .client_analysis_frequency_contract_fixture(
    3L, levels = c("\u00e9", "x"), chunk_coordinates = 2L)
  nfd <- .client_analysis_frequency_contract_fixture(
    3L, levels = c("e\u0301", "x"), chunk_coordinates = 2L)
  expect_false(identical(nfc$artifact_key, nfd$artifact_key))

  latin1_label <- iconv("caf\u00e9", from = "UTF-8", to = "latin1")
  Encoding(latin1_label) <- "latin1"
  expect_error(.client_analysis_frequency_contract_fixture(
    3L, levels = c(latin1_label, "x"), chunk_coordinates = 2L), "Frequency")

  large <- original
  large_plan <- large$semantic$analysis$effective_arguments$sampler_plan
  large_plan$implementation_delta <- list(
    numerator = "9007199254740993",
    denominator = "9007199254740993000000001")
  large$semantic$analysis$effective_arguments$sampler_plan <- large_plan
  large$artifact_key <- .dsvert_dp_analysis_artifact_key_v1(large$semantic)
  expect_silent(.dsvert_dp_analysis_contract_validate_v1(large))

  boundary <- original
  boundary_plan <- boundary$semantic$analysis$effective_arguments$sampler_plan
  boundary_plan$maximum_noise_per_peer <-
    "85070591730234615865843651857942052363"
  boundary_plan$no_wrap_sha256 <- .dsvert_dp_analysis_frequency_hash_v1(
    "dsVert/frequency/ring128-no-wrap/v1|", list(
      version = "dsvert-frequency-ring128-no-wrap-v1",
      coordinate_upper_bound = "1000",
      maximum_noise_per_peer =
        "85070591730234615865843651857942052363",
      maximum_noise_release =
        "170141183460469231731687303715884104726"))
  boundary$semantic$analysis$effective_arguments$sampler_plan <- boundary_plan
  boundary$artifact_key <-
    .dsvert_dp_analysis_artifact_key_v1(boundary$semantic)
  expect_silent(.dsvert_dp_analysis_contract_validate_v1(boundary))

  million_bound <- original
  million_plan <-
    million_bound$semantic$analysis$effective_arguments$sampler_plan
  million_bound$semantic$analysis$effective_arguments$coordinate_bounds$upper <-
    1000000
  million_plan$no_wrap_sha256 <- .dsvert_dp_analysis_frequency_hash_v1(
    "dsVert/frequency/ring128-no-wrap/v1|", list(
      version = "dsvert-frequency-ring128-no-wrap-v1",
      coordinate_upper_bound = "1000000",
      maximum_noise_per_peer = "100", maximum_noise_release = "200"))
  million_bound$semantic$analysis$effective_arguments$sampler_plan <-
    million_plan
  million_bound$artifact_key <-
    .dsvert_dp_analysis_artifact_key_v1(million_bound$semantic)
  expect_silent(.dsvert_dp_analysis_contract_validate_v1(million_bound))

  overflow <- boundary
  overflow_plan <- overflow$semantic$analysis$effective_arguments$sampler_plan
  overflow_plan$maximum_noise_per_peer <-
    "85070591730234615865843651857942052364"
  overflow_plan$no_wrap_sha256 <- .dsvert_dp_analysis_frequency_hash_v1(
    "dsVert/frequency/ring128-no-wrap/v1|", list(
      version = "dsvert-frequency-ring128-no-wrap-v1",
      coordinate_upper_bound = "1000",
      maximum_noise_per_peer =
        "85070591730234615865843651857942052364",
      maximum_noise_release =
        "170141183460469231731687303715884104728"))
  overflow$semantic$analysis$effective_arguments$sampler_plan <- overflow_plan
  expect_error(.dsvert_dp_analysis_artifact_key_v1(overflow$semantic))
})

test_that("client Frequency contracts fail closed on every bound dimension", {
  original <- .client_analysis_frequency_contract_fixture(3L)
  invalid_semantic <- list(
    roles = function(x) {
      x$noise_authority_roles$role_order <-
        list("secondary_noise_authority", "source_owner"); x
    },
    role_identity = function(x) {
      x$noise_authority_roles$authority_ids <-
        rev(x$noise_authority_roles$authority_ids); x
    },
    source = function(x) {
      x$analysis$effective_arguments$source_owner <-
        .client_analysis_identity_pk(99L); x
    },
    source_dataset = function(x) {
      x$analysis$effective_arguments$dataset_version <- "v2"; x
    },
    dimension = function(x) {
      x$analysis$effective_arguments$dimension <- 2; x
    },
    duplicate_level = function(x) {
      x$analysis$effective_arguments$levels[[2L]] <-
        x$analysis$effective_arguments$levels[[1L]]; x
    },
    oversized_level = function(x) {
      x$analysis$effective_arguments$levels[[1L]] <- strrep("x", 1025L); x
    },
    lane = function(x) {
      x$privacy$mechanism$randomness$lanes$final_noise$coordinates <- 2; x
    },
    shape = function(x) { x$public_shape$counts <- 2; x },
    primitive = function(x) {
      x$analysis$primitive <- "exact_gc_one_joint_discrete_laplace_draw_ring128_v3"; x
    },
    numeric_primitive = function(x) { x$analysis$primitive <- 1; x },
    mechanism = function(x) {
      x$privacy$mechanism$version <-
        "discrete-laplace-output-perturbation-v1"; x
    },
    sampler = function(x) {
      x$privacy$mechanism$calibration$sampler <-
        "discrete-laplace-one-draw-v1"; x
    },
    legacy_noise_scale = function(x) {
      x$privacy$mechanism$calibration$noise_scale <- 1; x
    },
    plan_dimension = function(x) {
      x$analysis$effective_arguments$sampler_plan$d <- 2; x
    },
    physical_plan = function(x) {
      x$analysis$effective_arguments$sampler_plan$physical_plan_version <-
        "dsvert-joint-dp-vector-laplace-plan-v3"; x
    },
    planner_request = function(x) {
      x$analysis$effective_arguments$sampler_plan$planner_request_sha256 <-
        strrep("f", 64L); x
    },
    profile = function(x) {
      x$analysis$effective_arguments$sampler_plan$profile_sha256 <-
        strrep("f", 64L); x
    },
    order = function(x) {
      x$analysis$effective_arguments$sampler_plan$coordinate_order_sha256 <-
        strrep("f", 64L); x
    },
    derived_lattice = function(x) {
      x$analysis$effective_arguments$sampler_plan$output_lattice_bits <- 2; x
    },
    raw_bound = function(x) {
      x$analysis$effective_arguments$coordinate_bounds$upper <- 999; x
    },
    partition = function(x) {
      x$analysis$effective_arguments$sampler_plan$chunk_coordinates <- 1; x
    },
    stream = function(x) {
      x$analysis$effective_arguments$sampler_plan$stream_binding <-
        list(mode = "absolute_coordinate_per_peer"); x
    },
    stream_role = function(x) {
      x$analysis$effective_arguments$sampler_plan$stream_binding$
        peer_role_binding <- "connection_alias"; x
    },
    runtime_claim = function(x) {
      x$analysis$effective_arguments$sampler_plan$stream_binding$
        runtime_sticky_claimed <- TRUE; x
    },
    fraction = function(x) {
      x$analysis$effective_arguments$sampler_plan$implementation_delta <- list(
          numerator = "2", denominator = "2000000"); x
    },
    allocation = function(x) {
      x$analysis$effective_arguments$sampler_plan$allocated_delta <-
        list(numerator = "1", denominator = "10"); x
    },
    selection = function(x) {
      x$analysis$effective_arguments$sampler_plan$backend_selection$
        runtime_failure_consulted <- TRUE; x
    },
    selection_backend = function(x) {
      x$analysis$effective_arguments$sampler_plan$backend_selection$
        policy_sha256 <- strrep("f", 64L); x
    },
    output = function(x) {
      x$analysis$effective_arguments$sampler_plan$output_transform <-
        "unclamped"; x
    },
    preimage = function(x) {
      x$analysis$effective_arguments$sampler_plan$preimage_validation$
        server_replan_required <- FALSE; x
    },
    no_wrap = function(x) {
      x$analysis$effective_arguments$sampler_plan$maximum_noise_per_peer <-
          "170141183460469231731687303715884105728"; x
    },
    no_wrap_hash = function(x) {
      x$analysis$effective_arguments$sampler_plan$no_wrap_sha256 <-
        strrep("f", 64L); x
    },
    zero_noise = function(x) {
      plan <- x$analysis$effective_arguments$sampler_plan
      plan$maximum_noise_per_peer <- "0"
      plan$no_wrap_sha256 <- .dsvert_dp_analysis_frequency_hash_v1(
        "dsVert/frequency/ring128-no-wrap/v1|", list(
          version = "dsvert-frequency-ring128-no-wrap-v1",
          coordinate_upper_bound = "1000",
          maximum_noise_per_peer = "0", maximum_noise_release = "0"))
      x$analysis$effective_arguments$sampler_plan <- plan
      x
    },
    excessive_epsilon = function(x) { x$privacy$epsilon <- 8.1; x },
    excessive_bound = function(x) {
      x$analysis$effective_arguments$coordinate_bounds$upper <- 1000001; x
    },
    vector_adjacency = function(x) {
      x$privacy$adjacency <- c("add_remove_patient", "add_remove_patient"); x
    },
    vector_snapshot_hash = function(x) {
      x$owner_snapshots[[1L]]$snapshot_commitment <-
        c(strrep("a", 64L), strrep("b", 64L)); x
    },
    vector_policy_hash = function(x) {
      x$privacy$contribution$constraints$policy_sha256 <-
        c(strrep("a", 64L), strrep("b", 64L)); x
    },
    numeric = function(x) { x$numeric$value_bits <- 127; x })
  for (mutate in invalid_semantic) {
    expect_error(.dsvert_dp_analysis_artifact_key_v1(
      mutate(original$semantic)))
  }

  wrong_kernel <- original
  wrong_kernel$execution$backend$kernel <- "other-backend-v1"
  expect_error(.dsvert_dp_analysis_contract_validate_v1(wrong_kernel),
               "execution backend|execution contract")
  wrong_chunk <- original
  wrong_chunk$execution$transport$chunk_coordinates <- 1
  expect_error(.dsvert_dp_analysis_contract_validate_v1(wrong_chunk),
               "execution transport|execution contract")
})

test_that("client accepts only fail-closed canonical Count TV contracts", {
  contracts <- lapply(c(2L, 3L, 5L), function(k) {
    .dsvert_dp_analysis_contract_validate_v1(
      .client_analysis_count_contract_fixture(k))
  })
  expect_true(all(vapply(contracts, function(contract) {
    identical(contract$semantic$numeric$fractional_bits, 0) &&
      identical(contract$semantic$numeric$value_bits, 127) &&
      identical(contract$semantic$numeric$overflow, "reject")
  }, logical(1L))))
  calibration <- contracts[[1L]]$semantic$privacy$mechanism$calibration
  expect_lt(calibration$implementation_delta,
            contracts[[1L]]$semantic$privacy$delta)

  original <- .client_analysis_count_contract_fixture(3L)
  expect_identical(
    original$artifact_key,
    "f930efe30f04bd6c2d118c4c69ae82a43cd15f62703b91b43ae99d464133e6b0")
  reordered <- original
  reordered$semantic <- reordered$semantic[rev(names(reordered$semantic))]
  reordered$execution <- reordered$execution[rev(names(reordered$execution))]
  expect_identical(
    .dsvert_dp_analysis_contract_validate_v1(reordered),
    .dsvert_dp_analysis_contract_validate_v1(original))

  invalid <- list(
    sampler = function(x) {
      x$semantic$privacy$mechanism$calibration$sampler <-
        "discrete-laplace-one-draw-v1"
      x$semantic$privacy$mechanism$randomness$lanes$final_noise$primitive <-
        "discrete-laplace-one-draw-v1"
      x
    },
    mechanism = function(x) {
      x$semantic$privacy$mechanism$version <-
        "discrete-laplace-output-perturbation-v1"
      x
    },
    legacy_pair = function(x) {
      x$semantic$privacy$mechanism$version <-
        "discrete-laplace-output-perturbation-v1"
      x$semantic$privacy$mechanism$calibration$sampler <-
        "discrete-laplace-one-draw-v1"
      x$semantic$privacy$mechanism$calibration$implementation_delta <- 0
      x$semantic$privacy$mechanism$randomness$lanes$final_noise$primitive <-
        "discrete-laplace-one-draw-v1"
      x
    },
    zero_implementation_delta = function(x) {
      x$semantic$privacy$mechanism$calibration$implementation_delta <- 0
      x
    },
    zero_delta = function(x) {
      x$semantic$privacy$delta <- 0
      x$semantic$privacy$mechanism$calibration$implementation_delta <- 0
      x
    },
    excessive_implementation_delta = function(x) {
      x$semantic$privacy$mechanism$calibration$implementation_delta <- 2e-6
      x
    },
    insufficient_scale = function(x) {
      x$semantic$privacy$mechanism$calibration$noise_scale <- 1 - 1e-12
      x
    },
    overflowing_scale = function(x) {
      x$semantic$privacy$mechanism$sensitivity$value <- 1e308
      x$semantic$privacy$mechanism$calibration$noise_scale <- 1e308
      x$semantic$privacy$epsilon <- 1e-308
      x
    },
    nonnumeric_certificate = function(x) {
      x$semantic$privacy$mechanism$calibration$implementation_delta <- "1e-9"
      x
    },
    saturating_overflow = function(x) {
      x$semantic$numeric$overflow <- "saturate"
      x
    },
    vector_noise = function(x) {
      x$semantic$privacy$mechanism$randomness$lanes$final_noise$coordinates <-
        2
      x
    },
    nonunit_sensitivity = function(x) {
      x$semantic$privacy$mechanism$sensitivity$value <- 2
      x$semantic$privacy$mechanism$calibration$noise_scale <- 2
      x
    },
    extra_lane = function(x) {
      x$semantic$privacy$mechanism$randomness$lanes$internal <- list(
        version = "dsvert-randomness-lane-v1",
        purpose = "confidential_internal_randomness",
        primitive = "aes128ctr-internal-v1",
        coordinates = 1)
      x
    },
    non_count_primitive = function(x) {
      x$semantic$analysis$primitive <- "other-analysis-v1"
      x$execution$backend$kernel <- "other-analysis-v1"
      x
    },
    unsupported_ring = function(x) {
      x$execution$backend$ring <- "ring128"
      x
    },
    narrow_value_bits = function(x) {
      x$semantic$numeric$value_bits <- 126
      x
    },
    wrong_rounding = function(x) {
      x$semantic$numeric$rounding <- "nearest_even"
      x
    },
    wrong_sampler_encoding = function(x) {
      x$semantic$numeric$sampler_encoding <- "other_encoding_v1"
      x
    },
    wrong_output_encoding = function(x) {
      x$semantic$numeric$output_encoding <- "other_encoding_v1"
      x
    },
    wrong_shape = function(x) {
      x$semantic$public_shape <- list(value = 1)
      x
    },
    negative_fractional_bits = function(x) {
      x$semantic$numeric$fractional_bits <- -1
      x
    },
    full_width_fractional_bits = function(x) {
      x$semantic$numeric$fractional_bits <- 127
      x
    })
  for (mutate in invalid) {
    expect_error(.dsvert_dp_analysis_contract_validate_v1(
      mutate(original)))
  }

  generic_zero <- .client_analysis_contract_fixture()
  generic_zero$semantic$numeric$fractional_bits <- 0
  generic_zero$artifact_key <-
    .dsvert_dp_analysis_artifact_key_v1(generic_zero$semantic)
  validated_zero <- .dsvert_dp_analysis_contract_validate_v1(generic_zero)
  expect_identical(validated_zero$semantic$numeric$fractional_bits, 0)
  bad_lane <- generic_zero
  bad_lane$semantic$privacy$mechanism$randomness$lanes$final_noise$coordinates <-
    0
  expect_error(.dsvert_dp_analysis_contract_validate_v1(bad_lane))
})

test_that("client validates the server semantic artifact key for K=2,3,5", {
  contracts <- lapply(c(2L, 3L, 5L), function(k) {
    contract <- .client_analysis_contract_fixture(k)
    contract$artifact_key <-
      .dsvert_dp_analysis_artifact_key_v1(contract$semantic)
    .dsvert_dp_analysis_contract_validate_v1(contract)
  })
  expect_true(all(vapply(contracts, function(contract) {
    identical(contract$artifact_key,
              .dsvert_dp_analysis_artifact_key_v1(contract$semantic))
  }, logical(1L))))

  original <- contracts[[2L]]
  expect_identical(
    original$artifact_key,
    "051e83176ab341f6d9461d71c97b9c14bb765fdb4e7f0220fd8a1d3579de4709")
  reordered <- original[rev(names(original))]
  reordered$semantic <- reordered$semantic[rev(names(reordered$semantic))]
  reordered$semantic$owner_snapshots <-
    reordered$semantic$owner_snapshots[
      rev(names(reordered$semantic$owner_snapshots))]
  reordered$semantic$noise_authorities <-
    rev(reordered$semantic$noise_authorities)
  expect_identical(
    .dsvert_dp_analysis_contract_validate_v1(reordered), original)

  vector_arguments <- original
  vector_arguments$semantic$analysis$effective_arguments$opaque <- c(1, 2)
  vector_arguments$artifact_key <- .dsvert_dp_analysis_artifact_key_v1(
    vector_arguments$semantic)
  list_arguments <- original
  list_arguments$semantic$analysis$effective_arguments$opaque <- list(1, 2)
  list_arguments$artifact_key <- .dsvert_dp_analysis_artifact_key_v1(
    list_arguments$semantic)
  expect_identical(
    .dsvert_dp_analysis_contract_validate_v1(vector_arguments),
    .dsvert_dp_analysis_contract_validate_v1(list_arguments))

  no_formula <- original
  no_formula$semantic$analysis["formula"] <- list(NULL)
  no_formula$artifact_key <- .dsvert_dp_analysis_artifact_key_v1(
    no_formula$semantic)
  expect_null(.dsvert_dp_analysis_contract_validate_v1(
    no_formula)$semantic$analysis$formula)

  for (ambiguous in list(
      stats::setNames(c(1, 2), c("age", "treatment")),
      matrix(1:4, nrow = 2L),
      structure(list(1, 2), dim = c(1L, 2L)))) {
    bad <- original
    bad$semantic$analysis$effective_arguments$ambiguous <- ambiguous
    expect_error(
      .dsvert_dp_analysis_contract_validate_v1(bad),
      "canonical analysis contract")
  }

  changed_execution <- original
  changed_execution$execution$backend$build_sha256 <- strrep("b", 64L)
  changed_execution$execution$backend$ring <- "ring128"
  changed_execution$execution$transport$chunk_coordinates <- 8192
  names(changed_execution$execution$peer_pins) <- paste0(
    "connection_", seq_along(changed_execution$execution$peer_pins))
  expect_identical(
    .dsvert_dp_analysis_contract_validate_v1(changed_execution)$artifact_key,
    original$artifact_key)

  tampered <- original
  tampered$semantic$analysis$effective_arguments$link <- "probit"
  expect_error(
    .dsvert_dp_analysis_contract_validate_v1(tampered),
    "artifact key")
  downgraded <- original
  downgraded$version <- "dsvert-analysis-contract-v0"
  expect_error(
    .dsvert_dp_analysis_contract_validate_v1(downgraded),
    "analysis contract")

  numeric <- original
  numeric$semantic$numeric$fractional_bits <- 40
  numeric$artifact_key <-
    .dsvert_dp_analysis_artifact_key_v1(numeric$semantic)
  expect_false(identical(numeric$artifact_key, original$artifact_key))

  authorities <- original
  authorities$semantic$noise_authorities <- unname(
    names(authorities$semantic$owner_snapshots)[2:3])
  expect_error(
    .dsvert_dp_analysis_artifact_key_v1(authorities$semantic),
    "semantic contract")

  edge <- .client_analysis_contract_fixture(2L)
  edge_pk <- .client_analysis_identity_pk(255L)
  names(edge$semantic$owner_snapshots)[1L] <- edge_pk
  edge$execution$peer_pins[[1L]] <- edge_pk
  edge$semantic$noise_authorities <- unname(sort(
    names(edge$semantic$owner_snapshots), method = "radix"))
  edge$artifact_key <- .dsvert_dp_analysis_artifact_key_v1(edge$semantic)
  expect_identical(
    sort(names(.dsvert_dp_analysis_contract_validate_v1(
      edge)$semantic$owner_snapshots), method = "radix"),
    sort(names(edge$semantic$owner_snapshots), method = "radix"))

  malformed <- .client_analysis_contract_fixture(2L)
  original_pk <- names(malformed$semantic$owner_snapshots)[1L]
  malformed_pk <- paste0(
    substr(original_pk, 1L, 10L), "=",
    substr(original_pk, 11L, nchar(original_pk)))
  names(malformed$semantic$owner_snapshots)[1L] <- malformed_pk
  malformed$execution$peer_pins[[1L]] <- malformed_pk
  malformed$semantic$noise_authorities[[1L]] <- malformed_pk
  expect_error(
    .dsvert_dp_analysis_contract_validate_v1(malformed),
    "semantic contract")

  aliased <- .client_analysis_contract_fixture(2L)
  canonical_pk <- .client_analysis_identity_pk(255L)
  alias_pk <- paste0(" \n", chartr("-_", "+/", canonical_pk), "=\t")
  names(aliased$semantic$owner_snapshots)[1L] <- alias_pk
  aliased$execution$peer_pins[[1L]] <- alias_pk
  aliased$semantic$noise_authorities[[1L]] <- alias_pk
  aliased$artifact_key <- .dsvert_dp_analysis_artifact_key_v1(
    aliased$semantic)
  expect_true(canonical_pk %in% names(
    .dsvert_dp_analysis_contract_validate_v1(
      aliased)$semantic$owner_snapshots))

  overpadded <- aliased
  overpadded_pk <- paste0(chartr("-_", "+/", canonical_pk), "==")
  names(overpadded$semantic$owner_snapshots)[1L] <- overpadded_pk
  overpadded$execution$peer_pins[[1L]] <- overpadded_pk
  overpadded$semantic$noise_authorities[[1L]] <- overpadded_pk
  expect_error(
    .dsvert_dp_analysis_contract_validate_v1(overpadded),
    "semantic contract")

  bad <- original
  bad$semantic$analysis$effective_arguments$session_id <- "session-1"
  expect_error(
    .dsvert_dp_analysis_contract_validate_v1(bad),
    "Operational request fields")
  for (field in c(
      "data_name", "peer_name", "frontdoor", "route", "ledger_path",
      "lifetime_limit", "privacy_epoch", "noise_epoch", "noise_key_id",
      "connection_order", "format", "postprocessing")) {
    bad <- original
    bad$semantic$analysis$effective_arguments[[field]] <- "operational"
    expect_error(
      .dsvert_dp_analysis_contract_validate_v1(bad),
      "Operational request fields")
  }

  bad <- original
  bad$execution$peer_pins[[2L]] <- bad$execution$peer_pins[[1L]]
  expect_error(
    .dsvert_dp_analysis_contract_validate_v1(bad),
    "execution contract")

  bad <- original
  bad$semantic$privacy$mechanism$calibration$sampler <-
    "gaussian-evil-v1"
  bad$semantic$privacy$mechanism$randomness$lanes$final_noise$primitive <-
    "gaussian-evil-v1"
  expect_error(
    .dsvert_dp_analysis_contract_validate_v1(bad),
    "semantic contract")

  bad <- original
  bad$semantic$privacy$mechanism$calibration$noise_scale <-
    .Machine$double.xmin
  expect_error(
    .dsvert_dp_analysis_contract_validate_v1(bad),
    "semantic contract")

  bad <- original
  bad$semantic$privacy$mechanism$family <- "laplace"
  bad$semantic$privacy$mechanism$version <-
    "laplace-output-perturbation-v1"
  bad$semantic$privacy$mechanism$sensitivity$norm <- "l1"
  bad$semantic$privacy$mechanism$sensitivity$value <- 1e308
  bad$semantic$privacy$mechanism$calibration$noise_scale <- 1
  bad$semantic$privacy$mechanism$calibration$sampler <-
    "laplace-one-draw-v1"
  bad$semantic$privacy$mechanism$calibration$implementation_delta <- 0
  bad$semantic$privacy$mechanism$randomness$lanes$final_noise$primitive <-
    "laplace-one-draw-v1"
  bad$semantic$privacy$epsilon <- 1e-308
  expect_error(
    .dsvert_dp_analysis_contract_validate_v1(bad),
    "semantic contract")

  bad <- original
  bad$semantic$privacy$mechanism$sensitivity$value <- 1e308
  bad$semantic$privacy$mechanism$calibration$noise_scale <- 1e308
  bad$semantic$privacy$mechanism$calibration$implementation_delta <- 1e-9
  bad$semantic$privacy$epsilon <- 1e-308
  bad$semantic$privacy$delta <- 0.1
  expect_error(
    .dsvert_dp_analysis_contract_validate_v1(bad),
    "semantic contract")

  bad <- original
  bad$semantic$numeric$irrelevant <- "accepted"
  expect_error(
    .dsvert_dp_analysis_contract_validate_v1(bad),
    "semantic contract")

  bad <- original
  bad$semantic$numeric$value_bits <- 128
  bad$artifact_key <- .dsvert_dp_analysis_artifact_key_v1(bad$semantic)
  expect_error(
    .dsvert_dp_analysis_contract_validate_v1(bad),
    "execution backend")

  bad <- original
  bad$semantic$privacy$mechanism$randomness$lanes$final_noise$primitive <-
    "laplace-one-draw-v1"
  expect_error(
    .dsvert_dp_analysis_contract_validate_v1(bad),
    "semantic contract")
})
