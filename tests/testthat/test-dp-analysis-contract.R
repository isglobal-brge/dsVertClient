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
