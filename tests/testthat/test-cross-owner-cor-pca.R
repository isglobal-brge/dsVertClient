# Focused current-certificate tests. The shared helper keeps real signed
# RELEASE/replay and two-authority evidence validation, with explicit seams for
# authenticated bundle/compile reconstruction; this is not federation evidence.
.vertcor_mixed_fixture <- function(k) {
  seed <- .gaussian_current_fixture(k)
  artifact <- seed$raw_artifact
  peers <- names(seed$conns)
  variables <- c(paste0(peers, "$x"), paste0(peers[[1L]], "$z"))
  outcome <- variables[[1L]]
  predictors <- sort(variables[-1L], method = "radix")
  descriptor <- function(reference) {
    ref <- .dsvert_dp_gaussian_reference(reference)
    list(column = ref$column, dataset = paste0("data_", ref$server),
         owner_peer = ref$server, lower = 0, upper = 1)
  }
  artifact$outcome <- descriptor(outcome)
  artifact$dataset <- artifact$outcome$dataset
  artifact$owner_peer <- artifact$outcome$owner_peer
  artifact$predictors <- stats::setNames(lapply(predictors, descriptor), predictors)
  artifact$predictor_order <- as.list(predictors)
  artifact$input_variable_order <- as.list(c(predictors, outcome))
  artifact$participating_peers <- as.list(peers)
  artifact$design_terms <- as.list(c("(Intercept)", predictors))
  m <- length(variables)
  dimension <- m * (m + 1) / 2 + m + 2
  artifact$coordinate_count <- dimension
  artifact$statistic_maximum <- as.list(rep(seed$capacity * seed$scale, dimension))
  artifact$source_raw_l1_sensitivity <- dimension * seed$scale
  artifact$source_raw_l2_sensitivity <- sqrt(dimension) * seed$scale
  artifact$natural_l1_sensitivity <- dimension
  artifact$natural_l2_sensitivity <- sqrt(dimension)
  artifact$transcript$variable_count <- m
  artifact$transcript$validity_product_rounds <- m - 1L
  # Binary factorial columns make every sufficient statistic lattice-exact.
  # The last variable is the fixed outcome; at least one predictor shares its
  # owner, and each remaining peer owns one predictor (including K=5).
  values <- vapply(seq_len(m), function(bit) {
    as.numeric(bitwAnd(0:63, bitwShiftL(1L, bit - 1L)) != 0L)
  }, numeric(64L))
  colnames(values) <- c(predictors, outcome)
  design <- cbind(1, values[, predictors, drop = FALSE])
  gram <- crossprod(design)
  coordinates <- c(nrow(values), gram[upper.tri(gram, diag = TRUE)],
    crossprod(design, values[, outcome]), sum(values[, outcome]^2))
  fixture <- .gaussian_current_fixture(k, artifact, coordinates)
  fixture$oracle <- stats::cor(values)
  fixture
}

.vertcor_no_transport <- function() {
  function(...) stop("offline path called transport or discovery", call. = FALSE)
}

.vertcor_reseal <- function(certificate) {
  certificate$certificate_sha256 <- .dsvert_dp_gaussian_certificate_hash(
    certificate[setdiff(names(certificate), "certificate_sha256")])
  certificate
}

test_that("current cross-owner certificates give all mixed pairs and PCA at K2/3/5", {
  for (k in c(2L, 3L, 5L)) {
    fixture <- .vertcor_mixed_fixture(k)
    .with_gaussian_current_fixture(fixture, {
      source <- .gaussian_current_source(fixture)
      testthat::local_mocked_bindings(
        .dsvert_dp_datasources = .vertcor_no_transport(),
        .dsvert_federation_argument = .vertcor_no_transport(),
        .dsvert_dp_synopsis_vector_run = .vertcor_no_transport(),
        .package = "dsVertClient")
      testthat::local_mocked_bindings(
        datashield.aggregate = .vertcor_no_transport(),
        datashield.connections_find = .vertcor_no_transport(), .package = "DSI")
      certificate_bytes <- serialize(source$certificate, NULL, version = 3L)
      artifact_bytes <- serialize(fixture$raw_artifact, NULL, version = 3L)
      release_bytes <- serialize(fixture$vector_release, NULL, version = 3L)
      result <- ds.vertCor(gaussian_result = source$certificate,
        trusted_pinset = fixture$pins, verbose = FALSE)
      expect_equal(result$correlation, fixture$oracle, tolerance = 1e-12)
      expect_identical(result$var_names, colnames(fixture$oracle))
      expect_identical(result$estimand_missingness, "complete_case_joint")
      expect_identical(result$cross_owner_state, "exact_gc_to_joint_dp_vector_v1")
      expect_identical(result$additional_server_calls_after_synopsis, 0L)
      expect_identical(result$additional_privacy_cost, c(epsilon = 0, delta = 0))
      explicit <- ds.vertDPCor(artifact_family = "gaussian_models",
        gaussian_result = source$certificate, trusted_pinset = fixture$pins)
      expect_identical(explicit, result)
      alias <- ds.vert.cor(gaussian_result = source$certificate,
        trusted_pinset = fixture$pins, verbose = FALSE)
      expect_identical(alias$correlation, result$correlation)
      expect_identical(alias$artifact_key, result$artifact_key)
      replay <- ds.vertCor(gaussian_result = source$certificate,
        trusted_pinset = fixture$pins, verbose = FALSE)
      expect_identical(serialize(replay, NULL, version = 3L),
                       serialize(result, NULL, version = 3L))
      rebuilt <- .gaussian_current_source(fixture)
      expect_identical(serialize(rebuilt$certificate, NULL, version = 3L), certificate_bytes)
      expect_identical(serialize(rebuilt$context$release, NULL, version = 3L),
                       serialize(source$context$release, NULL, version = 3L))
      expect_identical(serialize(fixture$raw_artifact, NULL, version = 3L), artifact_bytes)
      expect_identical(serialize(fixture$vector_release, NULL, version = 3L), release_bytes)
      pca <- ds.vertPCA(cor_result = result, trusted_pinset = fixture$pins,
                       verbose = FALSE)
      direct <- ds.vertPCA(gaussian_result = source$certificate,
                          trusted_pinset = fixture$pins, verbose = FALSE)
      expect_identical(direct, pca)
      alias_pca <- ds.vert.pca(gaussian_result = source$certificate,
        trusted_pinset = fixture$pins, verbose = FALSE)
      expect_identical(alias_pca$loadings, pca$loadings)
      expect_identical(alias_pca$additional_server_calls, 0L)
      expect_equal(sum(pca$eigenvalues), ncol(fixture$oracle), tolerance = 1e-12)
      expect_equal(unname(pca$loadings %*% diag(pca$eigenvalues) %*%
        t(pca$loadings)), unname(fixture$oracle), tolerance = 1e-12)
      expect_null(pca$scores)
      expect_false(pca$scores_available)
      expect_false(pca$synopsis_read_performed)
      expect_identical(pca$additional_server_calls, 0L)
      # Repeated eigenvalues: assert the invariant eigenspace reconstruction,
      # never a particular eigenvector sign or basis.
      chosen <- rev(result$var_names[c(1L, length(result$var_names))])
      subset <- ds.vertCor(variables = chosen, gaussian_result = source$certificate,
                           trusted_pinset = fixture$pins, verbose = FALSE)
      expect_identical(subset$var_names, chosen)
      expect_identical(subset$artifact_key, result$artifact_key)
      expect_identical(subset$final_vector_root, result$final_vector_root)
      expect_equal(subset$n_obs, result$n_obs)
      by_owner <- stats::setNames(rep(list("x"), k), names(fixture$conns))
      named <- ds.vertCor(variables = by_owner, gaussian_result = source$certificate,
                          trusted_pinset = fixture$pins, verbose = FALSE)
      expect_identical(named$var_names, paste0(names(fixture$conns), "$x"))
      expect_error(ds.vertCor(variables = c("x", "z"),
        gaussian_result = source$certificate, trusted_pinset = fixture$pins,
        verbose = FALSE), "resolve uniquely")
    })
  }
})

test_that("cross quantization is revalidated with a nonvacuous interval", {
  # Five custodians, two variable owners: site_a computes without owning a
  # model column, and site_d/site_e are membership-only participants.
  fixture <- .gaussian_current_fixture(5L,
    coordinates = c(80, 80, 40, 32.5, 40, 29, 26.75))
  .with_gaussian_current_fixture(fixture, {
    source <- .gaussian_current_source(fixture)
    result <- ds.vertCor(gaussian_result = source$certificate,
      trusted_pinset = fixture$pins, verbose = FALSE)
    expect_gt(result$correlation_95_interval_complete_case$lower[1L, 2L], -1)
    expect_lt(result$correlation_95_interval_complete_case$lower[1L, 2L], 0.9)
    expect_no_error(ds.vertPCA(cor_result = result,
      trusted_pinset = fixture$pins, verbose = FALSE))
    tampered <- result
    tampered$correlation_95_interval_complete_case$lower[1L, 2L] <- 0.99
    expect_error(ds.vertPCA(cor_result = tampered,
      trusted_pinset = fixture$pins, verbose = FALSE), "do not match")
  })
})

test_that("offline correlation refuses unanchored, stale, or detached evidence", {
  fixture <- .gaussian_current_fixture(2L)
  .with_gaussian_current_fixture(fixture, {
    source <- .gaussian_current_source(fixture)
    certificate <- source$certificate
    pin_hash <- certificate$peer_context$peer_pinset_sha256
    if (exists(pin_hash, envir = .dsvert_dp_gaussian_trust_cache, inherits = FALSE)) {
      rm(list = pin_hash, envir = .dsvert_dp_gaussian_trust_cache)
    }
    expect_error(ds.vertCor(gaussian_result = certificate, verbose = FALSE), "externally anchored")
    pins <- fixture$pins; pins[[1L]] <- pins[[2L]]
    expect_error(ds.vertCor(gaussian_result = certificate, trusted_pinset = pins,
      verbose = FALSE), "trusted pinset")
    expect_error(ds.vertCor(data_name = "wrong", gaussian_result = certificate,
      trusted_pinset = fixture$pins, verbose = FALSE), "data_name disagrees")
    expect_error(ds.vertCor(analysis_id = "wrong", gaussian_result = certificate,
      trusted_pinset = fixture$pins, verbose = FALSE), "analysis_id disagrees")
    expect_error(ds.vertDPCor(gaussian_result = certificate), "artifact_family")
    expect_error(ds.vertCor(gaussian_result = list(version = "legacy")), "current signed")
    expect_false(exists(pin_hash, envir = .dsvert_dp_gaussian_trust_cache, inherits = FALSE))
    for (field in c("block_values_sha256", "source_contract_sha256", "final_vector_root")) {
      bad <- certificate; bad[[field]] <- strrep("0", 64)
      expect_error(ds.vertCor(gaussian_result = .vertcor_reseal(bad),
        trusted_pinset = fixture$pins, verbose = FALSE))
    }
    for (field in c("signature", "peer_identity_pk", "source_contract_sha256",
      "private_layout_sha256", "transcript_sha256", "numeric_certificate_sha256",
      "exact_intermediates_exposed")) {
      bad <- certificate
      evidence <- .dsvert_dp_gaussian_synopsis_evidence_decode(
        bad$signed_evidence$cross_gaussian_evidence_json, "test evidence")
      evidence[[1L]][[field]] <- if (field == "exact_intermediates_exposed") TRUE else strrep("0", 64)
      bad$signed_evidence$cross_gaussian_evidence_json <- .dsvert_joint_dp_client_json(evidence)
      expect_error(ds.vertCor(gaussian_result = .vertcor_reseal(bad),
        trusted_pinset = fixture$pins, verbose = FALSE))
    }
    for (shape in c("missing", "nested", "wrong_peer")) {
      bad <- certificate
      evidence <- .dsvert_dp_gaussian_synopsis_evidence_decode(
        bad$signed_evidence$cross_gaussian_evidence_json, "test evidence")
      evidence <- switch(shape, missing = evidence[1L],
        nested = stats::setNames(list(evidence), certificate$analysis_id),
        wrong_peer = stats::setNames(evidence, c("site_a", "impostor")))
      bad$signed_evidence$cross_gaussian_evidence_json <- .dsvert_joint_dp_client_json(evidence)
      expect_error(ds.vertCor(gaussian_result = .vertcor_reseal(bad),
        trusted_pinset = fixture$pins, verbose = FALSE))
    }
    fit <- .gaussian_current_fit(fixture)
    accepted <- ds.vertCor(gaussian_result = fit, trusted_pinset = fixture$pins,
                           verbose = FALSE)
    fit$sufficient_statistics_dp$gram_projected[1L, 1L] <- 100000
    fit$provenance_integrity <- TRUE
    fit$provenance_authenticity <- "caller_anchored"
    expect_error(ds.vertCor(gaussian_result = fit, trusted_pinset = fixture$pins,
      verbose = FALSE), "does not match")
    bad <- accepted; bad$correlation[1L, 2L] <- 0
    expect_error(ds.vertPCA(cor_result = bad, trusted_pinset = fixture$pins,
      verbose = FALSE), "do not match")
  })
})

test_that("singular Gaussian regression still permits certificate-only PCA", {
  original <- .vertcor_mixed_fixture(2L)
  x <- rep(c(0, 1), 32L)
  y <- rep(c(0, 0, 1, 1), 16L)
  design <- cbind(1, x, x)
  gram <- crossprod(design)
  coordinates <- c(64, gram[upper.tri(gram, diag = TRUE)],
    crossprod(design, y), sum(y^2))
  fixture <- .gaussian_current_fixture(2L, original$raw_artifact, coordinates)
  .with_gaussian_current_fixture(fixture, {
    expect_error(.gaussian_current_fit(fixture), "singular")
    source <- .gaussian_current_source(fixture)
    pca <- ds.vertPCA(gaussian_result = source$certificate,
      trusted_pinset = fixture$pins, verbose = FALSE)
    expect_equal(unname(pca$eigenvalues), c(2, 1, 0), tolerance = 1e-12)
    expect_null(pca$scores)
    expect_false(pca$scores_available)
  })
})

test_that("qualified duplicate names work through both public Gaussian front doors", {
  fixture <- .vertcor_mixed_fixture(3L)
  .with_gaussian_current_fixture(fixture, {
    source <- .gaussian_current_source(fixture)
    testthat::local_mocked_bindings(
      .dsvert_dp_gaussian_synopsis_release = function(...) source,
      .package = "dsVertClient")
    fit <- ds.vertDPGaussian(fixture$artifact$dataset,
      fixture$artifact$analysis_id, datasources = fixture$conns)
    formula <- stats::as.formula(paste(
      tail(.dsvert_dp_gaussian_variable_order(fixture$artifact), 1L), "~",
      paste(fixture$artifact$predictor_order, collapse = " + ")))
    adapter <- ds.vertGLM(formula, data = fixture$artifact$dataset,
      family = "gaussian", dp_analysis_id = fixture$artifact$analysis_id,
      verbose = FALSE, datasources = fixture$conns)
    expect_identical(adapter$coefficients, fit$coefficients)
    expect_identical(adapter$provenance_certificate, fit$provenance_certificate)
    expect_named(fit$coefficients, fixture$artifact$design_terms)
    bad <- fixture$raw_artifact
    bad$predictors[[1L]]$owner_peer <- "site_z"
    expect_error(.dsvert_dp_gaussian_cross_artifact(bad,
      bad$dataset, bad$analysis_id, NULL, bad$adjacency,
      fixture$scale, fixture$capacity), "predictor order changed")
  })
})

test_that("cross complete-case certificates refuse absent intercepts and zero variance", {
  original <- .gaussian_current_fixture(2L)
  artifact <- original$raw_artifact
  artifact$intercept <- FALSE
  artifact$design_terms <- artifact$predictor_order
  artifact$coordinate_count <- 4
  artifact$statistic_maximum <- as.list(rep(original$capacity * original$scale, 4L))
  artifact$source_raw_l1_sensitivity <- 4 * original$scale
  artifact$source_raw_l2_sensitivity <- 2 * original$scale
  artifact$natural_l1_sensitivity <- 4
  artifact$natural_l2_sensitivity <- 2
  fixture <- .gaussian_current_fixture(2L, artifact, c(80, 30, 24, 30))
  .with_gaussian_current_fixture(fixture, {
    source <- .gaussian_current_source(fixture)
    condition <- tryCatch(ds.vertCor(gaussian_result = source$certificate,
      trusted_pinset = fixture$pins, verbose = FALSE), error = identity)
    expect_s3_class(condition, "non_identifiable")
    expect_identical(condition$reason, "pearson_requires_intercept_marginal_sums")
  })
  fixture <- .gaussian_current_fixture(2L, coordinates = c(80, 80, 0, 0, 40, 0, 30))
  .with_gaussian_current_fixture(fixture, {
    source <- .gaussian_current_source(fixture)
    condition <- tryCatch(ds.vertCor(gaussian_result = source$certificate,
      trusted_pinset = fixture$pins, verbose = FALSE), error = identity)
    expect_s3_class(condition, "non_identifiable")
    expect_identical(condition$reason, "zero_dp_complete_case_variance")
  })
})

test_that("signed bare identifiers remain distinct from physical-column aliases", {
  original <- .vertcor_mixed_fixture(3L)
  artifact <- original$raw_artifact
  # A signed catalog can retain one bare key while qualifying its homonyms.
  names(artifact$predictors)[names(artifact$predictors) == "site_b$x"] <- "x"
  order <- sort(names(artifact$predictors), method = "radix")
  artifact$predictors <- artifact$predictors[order]
  artifact$predictor_order <- as.list(order)
  artifact$design_terms <- as.list(c("(Intercept)", order))
  artifact$input_variable_order <- as.list(c(order, "site_a$x"))
  fixture <- .gaussian_current_fixture(3L, artifact, original$coordinates)
  .with_gaussian_current_fixture(fixture, {
    source <- .gaussian_current_source(fixture)
    correlation <- ds.vertCor(gaussian_result = source$certificate,
      trusted_pinset = fixture$pins, verbose = FALSE)
    expect_identical(correlation$var_names, c(order, "site_a$x"))
    expect_equal(unname(correlation$correlation), diag(4), tolerance = 1e-12)
    testthat::local_mocked_bindings(
      .dsvert_dp_gaussian_synopsis_release = function(...) source,
      .package = "dsVertClient")
    formula <- stats::as.formula(paste("site_a$x ~", paste(order, collapse = " + ")))
    fit <- ds.vertGLM(formula, data = artifact$dataset, family = "gaussian",
      dp_analysis_id = artifact$analysis_id, verbose = FALSE,
      datasources = fixture$conns)
    expect_named(fit$coefficients, c("(Intercept)", order))
  })
})
