formal_glm_authority <- function() {
  list(
    consortium_id = "consortium-test",
    capsule_id = "capsule-test",
    logical_snapshot = "snapshot-test",
    policy_contract_sha256 = strrep("1", 64L),
    custodian_peers = c("site_c", "site_a", "site_b"),
    designated_peers = c("site_b", "site_a"),
    pinset_sha256 = strrep("2", 64L))
}

formal_glm_model <- function(family = "binomial", formula = "y ~ x + group",
                             capacity = 4L, adjacency = "add_remove",
                             offset_mode = "none") {
  columns <- list(
    y = if (identical(family, "binomial")) {
      list(kind = "binary", owner = "site_a")
    } else {
      list(kind = "count", owner = "site_a", upper = "8")
    },
    x = list(kind = "numeric", owner = "site_b", lower = "-1", upper = "1"),
    z = list(kind = "numeric", owner = "site_c", lower = "-1", upper = "1"),
    group = list(kind = "categorical", owner = "site_b",
                 levels = c("control", "treated"), reference = "control",
                 contrast = "treatment"),
    w = list(kind = "weight", owner = "site_a", upper = "1"),
    off = list(kind = "offset", owner = "site_a", lower = "-1", upper = "1"),
    exposure = list(kind = "exposure", owner = "site_a",
                    lower = "0.5", upper = "2"))
  if (identical(offset_mode, "none")) {
    offset <- list(mode = "none", column = NULL)
  } else if (identical(offset_mode, "bounded_offset")) {
    offset <- list(mode = "bounded_offset", column = "off")
  } else {
    offset <- list(mode = "log_exposure", column = "exposure")
  }
  list(
    formula = formula,
    family = family,
    capacity = capacity,
    columns = columns,
    adjacency = adjacency,
    patient_collapse = list(
      unit = "aligned_patient", repeated_records = "reject_duplicates",
      row_order_invariant = TRUE, max_records_per_unit = 1L,
      conflict_policy = "zero_weight"),
    missingness = "complete_tuple_zero_weight",
    clipping = list(
      numeric = "clamp_then_quantize",
      categorical = "registered_level_or_zero_weight",
      binomial = "binary_or_zero_weight",
      poisson_count = "integer_then_patient_cap",
      offset = "clamp_then_quantize",
      weight = "clamp_then_quantize"),
    weights = list(mode = "unit", column = NULL),
    offset = offset,
    coefficient_box = "1",
    ridge = "1",
    optimizer = list(alpha = "0.1", iterations = 8L, start = "zero"),
    numeric = list(
      x_fraction_bits = 12L, offset_fraction_bits = 12L,
      beta_fraction_bits = 16L, link_fraction_bits = 20L,
      working_fraction_bits = 32L, reference_precision_bits = 64L),
    link = list(segments = 4L,
                construction = "uniform_monotone_dyadic_linear"),
    privacy = list(
      mechanism = "joint_discrete_gaussian_one_global_draw",
      epsilon = "1", delta = "1e-6",
      allocation = "one_stacked_capsule_vector"))
}

test_that("formal GLM rational primitives are exact and sign-safe", {
  expect_identical(
    .dsvert_glm_rat_json(.dsvert_glm_rat("-0.125")),
    list(numerator = "-1", denominator = "8"))
  huge <- paste0("1", strrep("0", 200L))
  exact <- .dsvert_glm_rat_add(huge, "1")
  expect_identical(
    .dsvert_glm_rat_json(exact)$numerator,
    paste0("1", strrep("0", 199L), "1"))
  expect_identical(.dsvert_glm_rat_cmp(
    .dsvert_glm_rat_round_dyadic("0.5", 0L), "0"), 0L)
  expect_identical(.dsvert_glm_rat_cmp(
    .dsvert_glm_rat_round_dyadic("1.5", 0L), "2"), 0L)
  expect_identical(.dsvert_glm_rat_cmp(
    .dsvert_glm_rat_round_dyadic("-1.5", 0L), "-2"), 0L)
  expect_equal(.dsvert_glm_rat_double(.dsvert_glm_rat_div(huge, huge)), 1)
  expect_identical(
    .dsvert_glm_rat_cmp(.dsvert_glm_rat_exp_interval("0", 64L)$lower,
                        "1"), 0L)
  log_one <- .dsvert_glm_rat_log_interval("1", 64L)
  expect_identical(.dsvert_glm_rat_cmp(log_one$lower, "0"), 0L)
  expect_identical(.dsvert_glm_rat_cmp(log_one$upper, "0"), 0L)
})

test_that("formal GLM compilation is canonical and signature-ready", {
  first <- formal_glm_model()
  second <- formal_glm_model(formula = "y ~ group + x")
  second$columns <- second$columns[rev(names(second$columns))]
  second$columns$z$upper <- "0.75"
  second$clipping <- second$clipping[rev(names(second$clipping))]

  compiled_first <- .dsvert_formal_glm_compile(
    first, formal_glm_authority())
  compiled_second <- .dsvert_formal_glm_compile(
    second, formal_glm_authority())

  expect_s3_class(compiled_first, "dsvert_formal_glm_compilation")
  expect_identical(compiled_first$canonical_json,
                   compiled_second$canonical_json)
  expect_identical(compiled_first$sha256, compiled_second$sha256)
  expect_true(startsWith(
    rawToChar(compiled_first$signature_payload),
    .DSVERT_FORMAL_GLM_SIGNATURE_DOMAIN))
  expect_false(compiled_first$production_release_ready)
  expect_true(.dsvert_formal_glm_validate_compilation(compiled_first))
  expect_identical(
    compiled_first$unsigned_schema$estimand$coefficient_order,
    as.list(c("(Intercept)", "group[treated]", "x")))
  expect_identical(
    compiled_first$unsigned_schema$authority$designated_peers,
    as.list(c("site_a", "site_b")))
})

test_that("formal GLM authority supports K=2 and K>=3 custodians", {
  model <- formal_glm_model(formula = "y ~ x")
  model$columns$z <- NULL
  authority <- formal_glm_authority()
  authority$custodian_peers <- c("site_b", "site_a")
  compiled <- .dsvert_formal_glm_compile(model, authority)
  expect_identical(
    compiled$unsigned_schema$authority$custodian_peers,
    as.list(c("site_a", "site_b")))

  bad <- authority
  bad$custodian_peers <- "site_a"
  expect_error(
    .dsvert_formal_glm_compile(model, bad),
    class = "dsvert_formal_glm_specification_error")
})

test_that("intercept and treatment contrasts compile without ambiguity", {
  with_intercept <- .dsvert_formal_glm_compile(
    formal_glm_model(formula = "y ~ group"), formal_glm_authority())
  without_intercept <- .dsvert_formal_glm_compile(
    formal_glm_model(formula = "y ~ 0 + group"), formal_glm_authority())
  expect_identical(
    with_intercept$unsigned_schema$estimand$coefficient_order,
    as.list(c("(Intercept)", "group[treated]")))
  expect_identical(
    without_intercept$unsigned_schema$estimand$coefficient_order,
    as.list(c("group[control]", "group[treated]")))
  expect_true(with_intercept$unsigned_schema$estimand$intercept)
  expect_false(without_intercept$unsigned_schema$estimand$intercept)
})

test_that("formal GLM compilation detects commitment tampering", {
  compiled <- .dsvert_formal_glm_compile(
    formal_glm_model(), formal_glm_authority())
  tampered <- compiled
  tampered$unsigned_schema$optimizer$iterations <- 9
  expect_error(
    .dsvert_formal_glm_validate_compilation(tampered),
    class = "dsvert_formal_glm_specification_error")
  tampered <- compiled
  tampered$sha256 <- strrep("0", 64L)
  expect_error(
    .dsvert_formal_glm_validate_compilation(tampered),
    class = "dsvert_formal_glm_specification_error")
})

test_that("formal GLM schema requires two valid pinned signatures", {
  compiled <- .dsvert_formal_glm_compile(
    formal_glm_model(), formal_glm_authority())
  peers <- c("site_a", "site_b")
  custodians <- c(peers, "site_c")
  keys <- stats::setNames(lapply(custodians, function(peer) {
    openssl::ed25519_keygen()
  }), custodians)
  pinset <- vapply(keys, function(key) {
    .dsvert_exact_gc_b64url_encode(tail(as.raw(as.list(key)$pubkey), 32L))
  }, character(1L))
  signatures <- lapply(keys[peers], function(key) {
    .dsvert_exact_gc_b64url_encode(openssl::ed25519_sign(
      compiled$signature_payload, key))
  })
  signed <- .dsvert_formal_glm_verify_cross_signatures(
    compiled, signatures, pinset)
  expect_s3_class(signed, "dsvert_formal_glm_signed_schema")
  expect_identical(signed$schema_sha256, compiled$sha256)

  bad <- signatures
  bad$site_b <- bad$site_a
  expect_error(
    .dsvert_formal_glm_verify_cross_signatures(compiled, bad, pinset),
    class = "dsvert_formal_glm_specification_error")
  expect_error(
    .dsvert_formal_glm_verify_cross_signatures(
      compiled, signatures["site_a"], pinset),
    class = "dsvert_formal_glm_specification_error")
})

test_that("unsupported formulas and public domains fail before execution", {
  formulas <- c("y ~ x:z", "y ~ log(x)", "y ~ .", "cbind(y, y) ~ x")
  for (formula in formulas) {
    model <- formal_glm_model(formula = formula)
    expect_error(
      .dsvert_formal_glm_compile(model, formal_glm_authority()),
      class = "dsvert_formal_glm_specification_error")
  }

  bad <- formal_glm_model()
  bad$columns$group$levels <- c("control", "control")
  expect_error(.dsvert_formal_glm_compile(bad, formal_glm_authority()),
               class = "dsvert_formal_glm_specification_error")
  bad <- formal_glm_model()
  bad$columns$group$contrast <- "sum"
  expect_error(.dsvert_formal_glm_compile(bad, formal_glm_authority()),
               class = "dsvert_formal_glm_specification_error")
  bad <- formal_glm_model()
  bad$ridge <- "0"
  expect_error(.dsvert_formal_glm_compile(bad, formal_glm_authority()),
               class = "dsvert_formal_glm_specification_error")
  bad <- formal_glm_model()
  bad$optimizer$alpha <- "2"
  expect_error(.dsvert_formal_glm_compile(bad, formal_glm_authority()),
               class = "dsvert_formal_glm_specification_error")
  bad <- formal_glm_model()
  bad$missingness <- "available_case"
  expect_error(.dsvert_formal_glm_compile(bad, formal_glm_authority()),
               class = "dsvert_formal_glm_specification_error")
  bad <- formal_glm_model()
  bad$patient_collapse$repeated_records <-
    "deterministic_patient_aggregate_then_clip"
  expect_error(.dsvert_formal_glm_compile(bad, formal_glm_authority()),
               class = "dsvert_formal_glm_specification_error")
})

test_that("link tables are exact-rational, continuous and monotone", {
  for (family in c("binomial", "poisson")) {
    model <- formal_glm_model(family = family, formula = "y ~ x")
    compiled <- .dsvert_formal_glm_compile(model, formal_glm_authority())
    table <- compiled$unsigned_schema$link_surrogate
    values <- lapply(table$values, .dsvert_glm_rat)
    slopes <- lapply(table$slopes, .dsvert_glm_rat)
    expect_true(all(vapply(seq_len(length(values) - 1L), function(index) {
      .dsvert_glm_rat_cmp(values[[index]], values[[index + 1L]]) <= 0L
    }, logical(1L))))
    expect_true(all(vapply(slopes, function(value) {
      .dsvert_glm_rat_cmp(value, "0") >= 0L
    }, logical(1L))))
    zero_index <- which(vapply(table$knots, function(value) {
      .dsvert_glm_rat_cmp(value, "0") == 0L
    }, logical(1L)))
    expect_length(zero_index, 1L)
    expected <- if (family == "binomial") "0.5" else "1"
    expect_identical(
      .dsvert_glm_rat_cmp(table$values[[zero_index]], expected), 0L)
    expect_true(table$error_proof$continuity)
    expect_true(table$error_proof$monotone)
    expect_true(.dsvert_glm_rat_cmp(
      table$uniform_mean_error_upper, "0") > 0L)
    domain <- c(.dsvert_glm_rat_double(table$domain_lower),
                .dsvert_glm_rat_double(table$domain_upper))
    eta <- seq(domain[[1L]], domain[[2L]], length.out = 101L)
    approximation <- vapply(eta, function(value) {
      .dsvert_glm_rat_double(.dsvert_formal_glm_pwl_mean(value, table))
    }, numeric(1L))
    truth <- if (family == "binomial") stats::plogis(eta) else exp(eta)
    error_bound <- .dsvert_glm_rat_double(table$uniform_mean_error_upper)
    expect_true(max(abs(approximation - truth)) <=
                  error_bound * (1 + 1e-12))
  }
})

test_that("regularized rational oracle handles separation and fixed transcript", {
  model <- formal_glm_model(formula = "y ~ x", capacity = 4L)
  compiled <- .dsvert_formal_glm_compile(model, formal_glm_authority())
  separated <- data.frame(y = c(0, 0, 1, 1),
                          x = c(-1, -0.5, 0.5, 1))
  missing <- data.frame(y = rep(NA_real_, 4), x = rep(NA_real_, 4))
  fit <- .dsvert_formal_glm_rational_oracle(compiled, separated)
  empty_fit <- .dsvert_formal_glm_rational_oracle(compiled, missing)

  expect_s3_class(fit, "dsvert_formal_glm_reference")
  expect_identical(fit$status, "estimable_regularized_reference")
  expect_true(all(is.finite(fit$coefficients)))
  expect_true(all(abs(fit$coefficients) <= 1))
  expect_identical(fit$logical_work, empty_fit$logical_work)
  expect_identical(fit$logical_work$data_dependent_stopping, FALSE)
  expect_error(
    .dsvert_formal_glm_ordinary_reference(compiled, separated),
    class = "non_identifiable")
})

test_that("full ridge keeps collinear and all-zero Poisson targets unique", {
  binomial <- formal_glm_model(formula = "y ~ x + z", capacity = 4L)
  compiled_binomial <- .dsvert_formal_glm_compile(
    binomial, formal_glm_authority())
  collinear <- data.frame(y = c(0, 1, 0, 1), x = c(-1, 0, 0.5, 1),
                          z = c(-1, 0, 0.5, 1))
  expect_true(all(is.finite(
    .dsvert_formal_glm_rational_oracle(
      compiled_binomial, collinear)$coefficients)))
  expect_error(
    .dsvert_formal_glm_ordinary_reference(compiled_binomial, collinear),
    class = "non_identifiable")

  poisson <- formal_glm_model(
    family = "poisson", formula = "y ~ x", capacity = 4L)
  compiled_poisson <- .dsvert_formal_glm_compile(
    poisson, formal_glm_authority())
  zeros <- data.frame(y = rep(0, 4), x = c(-1, 0, 0.5, 1))
  fit <- .dsvert_formal_glm_rational_oracle(compiled_poisson, zeros)
  expect_true(all(is.finite(fit$coefficients)))
  expect_error(
    .dsvert_formal_glm_ordinary_reference(compiled_poisson, zeros),
    class = "non_identifiable")
})

test_that("ordinary central comparator is explicitly separate when estimable", {
  binomial_model <- formal_glm_model(
    formula = "y ~ x", capacity = 8L)
  binomial <- .dsvert_formal_glm_compile(
    binomial_model, formal_glm_authority())
  binomial_data <- data.frame(
    y = c(0, 1, 0, 1, 1, 0, 1, 0),
    x = c(-1, -0.75, -0.5, -0.25, 0.25, 0.5, 0.75, 1))
  ordinary_binomial <- .dsvert_formal_glm_ordinary_reference(
    binomial, binomial_data)
  surrogate_binomial <- .dsvert_formal_glm_rational_oracle(
    binomial, binomial_data)
  expect_identical(ordinary_binomial$status, "ordinary_reference_estimable")
  expect_identical(surrogate_binomial$status,
                   "estimable_regularized_reference")
  expect_true(ordinary_binomial$certificate$not_the_private_release_estimand)

  poisson_model <- formal_glm_model(
    family = "poisson", formula = "y ~ x", capacity = 8L)
  poisson <- .dsvert_formal_glm_compile(
    poisson_model, formal_glm_authority())
  poisson_data <- data.frame(
    y = c(1, 2, 1, 3, 2, 4, 3, 5),
    x = c(-1, -0.75, -0.5, -0.25, 0.25, 0.5, 0.75, 1))
  ordinary_poisson <- .dsvert_formal_glm_ordinary_reference(
    poisson, poisson_data)
  expect_identical(ordinary_poisson$status, "ordinary_reference_estimable")
  expect_true(all(is.finite(ordinary_poisson$coefficients)))
})

test_that("offset, exposure, weights and clipping have fixed semantics", {
  model <- formal_glm_model(
    family = "poisson", formula = "y ~ x", offset_mode = "log_exposure")
  model$weights <- list(mode = "bounded_column", column = "w")
  compiled <- .dsvert_formal_glm_compile(model, formal_glm_authority())
  data <- data.frame(
    y = c(0, 20, 1, 2), x = c(-20, 20, 0, NA),
    exposure = c(0.25, 4, 1, 1), w = c(-1, 2, 1, 1))
  fit <- .dsvert_formal_glm_rational_oracle(compiled, data)
  expect_true(all(is.finite(fit$coefficients)))
  expect_identical(
    compiled$unsigned_schema$estimand$offset$mode, "log_exposure")
  expect_true(.dsvert_glm_rat_cmp(
    compiled$unsigned_schema$estimand$offset$log_quantization_error_upper,
    "0") >= 0L)
})

test_that("sensitivity uses quantized rather than raw decimal bounds", {
  model <- formal_glm_model(formula = "y ~ x", offset_mode = "bounded_offset")
  model$columns$x$lower <- "-0.2"
  model$columns$x$upper <- "0.2"
  model$columns$w$upper <- "0.2"
  model$columns$off$lower <- "-0.2"
  model$columns$off$upper <- "0.2"
  model$weights <- list(mode = "bounded_column", column = "w")
  model$numeric$x_fraction_bits <- 2L
  model$numeric$offset_fraction_bits <- 2L
  compiled <- .dsvert_formal_glm_compile(model, formal_glm_authority())
  x_term <- compiled$unsigned_schema$estimand$coefficients[[2L]]$term
  expect_identical(.dsvert_glm_rat_cmp(x_term$abs_bound, "0.25"), 0L)
  expect_identical(.dsvert_glm_rat_cmp(
    compiled$unsigned_schema$estimand$weights$maximum_patient_weight,
    "0.25"), 0L)
  expect_identical(.dsvert_glm_rat_cmp(
    compiled$unsigned_schema$estimand$offset$upper, "0.25"), 0L)
  rows <- .dsvert_formal_glm_fixture_rows(
    compiled,
    data.frame(y = c(0, 1, 0, 1), x = c(-100, 100, 0, 0),
               w = rep(100, 4), off = c(-100, 100, 0, 0)))
  expect_identical(.dsvert_glm_rat_cmp(
    rows[[2L]]$design[[2L]], "0.25"), 0L)
  expect_identical(.dsvert_glm_rat_cmp(rows[[2L]]$weight, "0.25"), 0L)
  expect_identical(.dsvert_glm_rat_cmp(rows[[2L]]$offset, "0.25"), 0L)
})

test_that("certified sensitivity bounds enumerated adjacent reference fits", {
  base <- formal_glm_model(formula = "y ~ 1", capacity = 2L)
  for (adjacency in c("add_remove", "replace_one")) {
    base$adjacency <- adjacency
    compiled <- .dsvert_formal_glm_compile(base, formal_glm_authority())
    left <- if (adjacency == "add_remove") {
      data.frame(y = c(0, NA_real_))
    } else {
      data.frame(y = c(0, 0))
    }
    right <- if (adjacency == "add_remove") {
      data.frame(y = c(0, 1))
    } else {
      data.frame(y = c(0, 1))
    }
    fit_left <- .dsvert_formal_glm_rational_oracle(compiled, left)
    fit_right <- .dsvert_formal_glm_rational_oracle(compiled, right)
    observed <- .dsvert_glm_rat_abs(.dsvert_glm_rat_sub(
      fit_left$coefficients_rational[[1L]],
      fit_right$coefficients_rational[[1L]]))
    bound <- .dsvert_glm_rat(
      compiled$unsigned_schema$theorem_certificate$
        global_l2_sensitivity_reference)
    expect_lte(.dsvert_glm_rat_cmp(observed, bound), 0L)
  }
})

test_that("joint L2 sensitivity covers simultaneous outcome and feature change", {
  for (adjacency in c("add_remove", "replace_one")) {
    model <- formal_glm_model(
      formula = "y ~ x", capacity = 2L, adjacency = adjacency)
    compiled <- .dsvert_formal_glm_compile(model, formal_glm_authority())
    left <- if (adjacency == "add_remove") {
      data.frame(y = c(0, NA_real_), x = c(-1, NA_real_))
    } else {
      data.frame(y = c(0, 0), x = c(-1, -1))
    }
    right <- data.frame(y = c(0, 1), x = c(-1, 1))
    beta_left <- .dsvert_formal_glm_rational_oracle(
      compiled, left)$coefficients_rational
    beta_right <- .dsvert_formal_glm_rational_oracle(
      compiled, right)$coefficients_rational
    squared_distance <- .dsvert_formal_glm_sum(Map(function(a, b) {
      .dsvert_glm_rat_pow(.dsvert_glm_rat_sub(a, b), 2L)
    }, beta_left, beta_right))
    bound <- .dsvert_glm_rat(
      compiled$unsigned_schema$theorem_certificate$
        global_l2_sensitivity_reference)
    expect_lte(.dsvert_glm_rat_cmp(
      squared_distance, .dsvert_glm_rat_pow(bound, 2L)), 0L)
  }
})

test_that("Phase-0 artifact contract does not pretend moments are sufficient", {
  compiled <- .dsvert_formal_glm_compile(
    formal_glm_model(), formal_glm_authority())
  requirements <- .dsvert_formal_glm_artifact_requirements(compiled)
  expect_false(requirements$protected_input_artifact$finite_sufficient_statistic)
  expect_identical(
    requirements$public_dp_output_artifact$intermediate_release_count, 0)
  expect_false(
    requirements$public_dp_output_artifact$score_hessian_deviance_residual_release)
  expect_true("signed_segment_comparison_and_mux" %in%
                unlist(requirements$exact_gc_required, use.names = FALSE))
  expect_identical(requirements$production_phase, "phase1_not_connected")
})

test_that("unrepresentable public link domains fail typed and never fallback", {
  model <- formal_glm_model(formula = "y ~ x")
  model$coefficient_box <- "100"
  expect_error(
    .dsvert_formal_glm_compile(model, formal_glm_authority()),
    class = "dsvert_numeric_backend_unrepresentable")
})
