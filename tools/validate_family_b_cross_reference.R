#!/usr/bin/env Rscript
# SOURCE-TREE TEST HARNESS ONLY. All rows and RNG seeds below are public
# synthetic fixtures. The custom DSLite method and mocked release boundary
# are never registered by either package. This is not production MPC/PSI
# validation; that gate belongs to the parallel fused-producer integration.

.family_b_reference_run <- function(root, output = NULL) {
  root <- normalizePath(root, mustWork = TRUE)
  server_dir <- file.path(root, "dsVert")
  client_dir <- file.path(root, "dsVertClient")
  pkgload::load_all(server_dir, quiet = TRUE)
  pkgload::load_all(client_dir, quiet = TRUE)
  fixture_env <- new.env(parent = asNamespace("dsVert"))
  sys.source(file.path(server_dir, "tests/testthat/helper-crossgrid-family-b.R"), fixture_env)
  ns <- asNamespace("dsVertClient")
  mpc_dir <- file.path(server_dir, "inst/dsvert-mpc")
  scratch <- tempfile("family-b-reference-")
  dir.create(scratch)
  on.exit(unlink(scratch, recursive = TRUE), add = TRUE)
  binary <- file.path(scratch, "family-b-reference.test")
  oldwd <- setwd(mpc_dir)
  family_sources <- c("k2_exact_gc_multinomial_loss.go",
    "k2_exact_gc_multinomial_loss_profile.go", "k2_exact_gc_ordinal_loss.go",
    "k2_exact_gc_multinomial_loss_test.go", "k2_exact_gc_ordinal_loss_test.go",
    "k2_exact_gc_multinomial_loss_bridge_test.go")
  status <- system2("go", c("test", "-c", "-tags=dsvert_family_reference_test",
    "-o", shQuote(binary), family_sources),
    stdout = file.path(scratch, "build.log"), stderr = file.path(scratch, "build.log"))
  setwd(oldwd)
  if (status != 0L) stop(paste(readLines(file.path(scratch, "build.log")), collapse = "\n"))

  # Each peer has an independent DSLite session and only its vertical columns.
  # The test method intentionally exports PUBLIC synthetic rows to the tagged
  # reference harness. It must never be added to a production allowlist.
  config <- list(AggregateMethods = data.frame(name = "familyBReference",
    value = paste0("function(D, seed, scale, candidates) { set.seed(seed); ",
      "list(source=D, noise=stats::rgamma(candidates,0.5,scale=scale)-",
      "stats::rgamma(candidates,0.5,scale=scale)) }"),
    class = "script", stringsAsFactors = FALSE), AssignMethods = data.frame(), Options = list())
  seed <- 9182026L
  set.seed(seed)
  n <- 12000L
  x <- sample(c(0, 1), n, replace = TRUE)
  z <- sample(c(0, 1), n, replace = TRUE)
  ids <- sprintf("public-%05d", seq_len(n))
  results <- list()
  for (family in c("multinomial", "ordinal")) {
    if (family == "multinomial") {
      truth <- c(0, 2, -2, 0, -2, 2)
      eta <- cbind(0, 2*x-2*z, -2*x+2*z)
      probability <- exp(eta) / rowSums(exp(eta))
      grid <- list(truth, rep(0, 6), -truth)
      fixture <- fixture_env$.family_b_contract_fixture(family, beta_grid = grid, capacity = n)
    } else {
      truth <- c(0, 2, -2)
      eta <- 2*x-2*z
      cuts <- c(-0.75, 0.75)
      cumulative <- cbind(0, plogis(cuts[1]-eta), plogis(cuts[2]-eta), 1)
      probability <- cumulative[, -1] - cumulative[, -4]
      grid <- lapply(c(1, 0, -1), function(multiplier) {
        list(beta = multiplier*truth, thresholds = cuts)
      })
      fixture <- fixture_env$.family_b_contract_fixture(family, candidate_grid = grid, capacity = n)
    }
    y <- vapply(seq_len(n), function(i) sample.int(3, 1, prob = probability[i, ]), integer(1)) - 1L
    pooled <- data.frame(x = x, z = z, y = factor(letters[y+1L], levels = letters[1:3]))
    server_validate <- get(paste0(".dsvert_dp_", family, "_grid_cross_contract_validate"), asNamespace("dsVert"))
    client_validate <- get(paste0(".dsvert_dp_", family, "_grid_cross_contract_validate"), ns)
    server_contract <- server_validate(fixture$contract, fixture$policy, fixture$schema)
    contract <- client_validate(fixture$contract, fixture$policy, fixture$schema)
    stopifnot(identical(server_contract, contract))
    spec <- contract$spec
    candidates <- if (family == "multinomial") spec$beta_grid else spec$candidate_grid
    m <- length(candidates)
    owner_a <- data.frame(id = ids, x = x, y = y)
    owner_b <- data.frame(id = ids, z = z)
    owner_a <- owner_a[sample.int(n), ]
    owner_b <- owner_b[sample.int(n), ]
    peers <- lapply(list(owner_a, owner_b), function(rows) {
      s <- DSLite::newDSLiteServer(tables = list(D = rows), config = config,
        strict = TRUE, home = tempfile("dslite-", tmpdir = scratch))
      sid <- s$newSession(profile = "default")
      s$assignTable(sid, "D", "D")
      list(server = s, sid = sid)
    })
    materialize <- function(epsilon) lapply(seq_along(peers), function(i) {
      peers[[i]]$server$aggregate(peers[[i]]$sid, call("familyBReference", quote(D),
        seed + 100L*i + as.integer(epsilon),
        spec$sensitivity$raw_l1_sensitivity / epsilon, m))
    })
    parts <- materialize(1)
    # Public fixture alignment only. Production must consume private PSI shares.
    aligned <- merge(parts[[1]]$source, parts[[2]]$source, by = "id", sort = TRUE)
    stopifnot(nrow(aligned) == n, identical(aligned$id, ids))
    unique_rows <- unique(aligned[, c("x", "z", "y")])
    key <- function(d) paste(d$x, d$z, d$y, sep = ":")
    counts <- tabulate(match(key(aligned), key(unique_rows)), nbins = nrow(unique_rows))
    encoded_beta <- if (family == "multinomial") spec$beta_encoded else
      lapply(spec$candidate_encoded, `[[`, "beta")
    encoded_thresholds <- if (family == "ordinal") lapply(spec$candidate_encoded, `[[`, "thresholds") else NULL
    caps <- vapply(spec$sensitivity$candidate_bounds, `[[`, numeric(1), "per_patient_cap")
    request <- list(family = family, classes = 3, g = spec$numeric_grid_bits,
      features_encoded = lapply(seq_len(nrow(unique_rows)), function(i) {
        as.list(unname(as.numeric(unique_rows[i, c("x", "z")])) * 2^50)
      }), beta_encoded = encoded_beta, thresholds_encoded = encoded_thresholds,
      outcomes = as.list(unique_rows$y), valid = as.list(rep(TRUE, nrow(unique_rows))),
      caps = as.list(caps))
    input <- file.path(scratch, paste0(family, ".json"))
    out <- file.path(scratch, paste0(family, "-losses.json"))
    # Every numeric request field is an exact integer. Default double JSON
    # formatting emits scientific notation and loses f50 integer digits.
    jsonlite::write_json(request, input, auto_unbox = TRUE, digits = 0, null = "null")
    roundtrip <- jsonlite::read_json(input)
    feature_numbers <- function(rows) lapply(rows, function(row) {
      as.numeric(unlist(row, use.names = FALSE))
    })
    stopifnot(identical(feature_numbers(roundtrip$features_encoded),
                        feature_numbers(request$features_encoded)))
    status <- system2(binary, c("-test.run=^TestFamilyReferenceBridge$", "-test.v"),
      env = c(paste0("DSVERT_FAMILY_REFERENCE_INPUT=", shQuote(input)),
              paste0("DSVERT_FAMILY_REFERENCE_OUTPUT=", shQuote(out))),
      stdout = file.path(scratch, "reference.log"),
      stderr = file.path(scratch, "reference.log"))
    if (status != 0L || !file.exists(out)) stop(paste(readLines(file.path(scratch, "reference.log")), collapse = "\n"))
    response <- jsonlite::read_json(out, simplifyVector = TRUE)
    quantized <- response$rows
    stopifnot(is.matrix(quantized), identical(dim(quantized), c(nrow(unique_rows), m)))
    # JSON may simplify small row coordinates to R's 32-bit integer type.
    # Signed accumulated coordinates fit the exact binary64 integer domain.
    storage.mode(quantized) <- "double"
    integer_totals <- colSums(quantized * counts)
    stopifnot(all(is.finite(integer_totals)),
      all(integer_totals == floor(integer_totals)), all(integer_totals >= 0),
      all(integer_totals <= unlist(spec$sensitivity$maximum_coordinates)))
    real_rows <- sapply(candidates, function(candidate) {
      design <- cbind(1, unique_rows$x, unique_rows$z)
      if (family == "multinomial") {
        eta <- cbind(0, design %*% matrix(unlist(candidate), nrow = 3))
        log(rowSums(exp(eta))) - eta[cbind(seq_len(nrow(eta)), unique_rows$y+1L)]
      } else {
        eta <- drop(design %*% unlist(candidate$beta))
        thresholds <- unlist(candidate$thresholds)
        cdf <- cbind(0, plogis(thresholds[1]-eta), plogis(thresholds[2]-eta), 1)
        -log((cdf[, -1]-cdf[, -4])[cbind(seq_along(eta), unique_rows$y+1L)])
      }
    })
    real_totals <- colSums(real_rows * counts)
    maximum_error <- max(abs(quantized / 2^spec$numeric_grid_bits - real_rows))
    profile_error <- as.numeric(spec$numeric_contract$certified_uniform_error)
    output_error <- if (spec$numeric_grid_bits < 16) {
      2^(-spec$numeric_grid_bits - 1)
    } else 0
    row_error_bound <- profile_error + output_error
    stopifnot(maximum_error <= row_error_bound)
    exact_best <- which.min(real_totals)
    profile_best <- which.min(integer_totals)
    if (family == "multinomial") {
      fit <- nnet::multinom(y ~ x + z, data = pooled, trace = FALSE, maxit = 1000)
      central_parameters <- as.numeric(t(stats::coef(fit)))
      grid_parameters <- lapply(candidates, unlist)
    } else {
      fit <- MASS::polr(y ~ x + z, data = pooled, method = "logistic", Hess = FALSE)
      central_parameters <- c(0, stats::coef(fit), fit$zeta)
      grid_parameters <- lapply(candidates, function(candidate) unlist(c(candidate$beta, candidate$thresholds)))
    }
    central_nearest <- which.min(vapply(grid_parameters, function(value) sum((value-central_parameters)^2), numeric(1)))
    stopifnot(central_nearest == exact_best)
    cache <- new.env(parent = emptyenv())
    for (epsilon in c(1, 4, 8)) {
      release <- function(contract, datasources) {
        key <- paste(contract$artifact$spec_sha256, epsilon, sep = ":")
        if (exists(key, cache, inherits = FALSE)) return(get(key, cache))
        noise <- materialize(epsilon)
        noisy <- integer_totals + round(noise[[1]]$noise + noise[[2]]$noise)
        value <- get(paste0(".dsvert_dp_", family, "_grid_cross_postprocess"), ns)(contract, noisy)
        assign(key, value, cache)
        value
      }
      # This explicit testthat namespace mock is SOURCE-ONLY. Outside this
      # lexical test block the public frontdoor must still fail closed.
      frontdoor <- get(if (family == "multinomial") "dp_multinomial_grid" else "dp_ordinal_grid", ns)
      invoke <- function() frontdoor(peer_a$y ~ peer_a$x + peer_b$z, "cohort", "family-b-grid",
        fixture$contract, fixture$policy, fixture$schema)
      bindings <- list(code = quote({ first <- invoke(); replay <- invoke(); list(first, replay) }),
        .package = "dsVertClient")
      bindings[[paste0(".dsvert_dp_", family, "_grid_cross_release")]] <- release
      selected <- do.call(testthat::with_mocked_bindings, bindings, envir = environment())
      stopifnot(identical(selected[[1]], selected[[2]]),
        selected[[1]]$selected_candidate == exact_best,
        is.null(selected[[1]]$standard_errors))
      results[[length(results)+1L]] <- list(family = family, peers = 2, epsilon = epsilon,
        seed = seed, n = n, candidates = m, exact_best = exact_best,
        profile_best = profile_best,
        dp_selected = selected[[1]]$selected_candidate, central_nearest = central_nearest,
        central_fit = if (family == "multinomial") "nnet::multinom" else "MASS::polr",
        central_negative_log_likelihood = -as.numeric(stats::logLik(fit)),
        exact_grid_negative_log_likelihood = real_totals[exact_best],
        gap = sort(real_totals)[2]-min(real_totals), maximum_row_error = maximum_error,
        numeric_profile = spec$numeric_contract$version,
        reference_noise = "two_independent_gamma_half_differences_summed_then_rounded",
        certified_row_error = row_error_bound,
        certified_total_error = n * row_error_bound,
        observed_maximum_total_error = max(abs(integer_totals / 2^spec$numeric_grid_bits - real_totals)),
        natural_noise_scale = spec$sensitivity$natural_l1_sensitivity / epsilon,
        certified_error_to_noise_ratio = n * row_error_bound * epsilon /
          spec$sensitivity$natural_l1_sensitivity,
        sticky_replay_identical = TRUE, plaintext_reference_test_tag = "dsvert_family_reference_test",
        production_mpc_psi_release_validated = FALSE)
    }
    unavailable <- tryCatch(invoke(), error = identity)
    stopifnot(inherits(unavailable, "dsvert_dp_public_failure"))
  }
  report <- list(scope = "two-peer DSLite public-synthetic reference harness; test-only release mock",
    alignment = "public_synthetic_id_join_not_private_PSI",
    noise = "reference_continuous_Laplace_postprocessing_not_production_discrete_sampler",
    production_release_gate = "pending fused-producer integration", cases = results)
  if (!is.null(output)) jsonlite::write_json(report, output, pretty = TRUE, auto_unbox = TRUE, digits = NA)
  report
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if ("--help" %in% args) {
    cat("Usage: Rscript tools/validate_family_b_cross_reference.R ROOT [OUTPUT.json]\n",
      "Public synthetic two-peer DSLite reference test; not production release validation.\n")
  } else {
    if (!length(args)) stop("A sibling-checkout ROOT is required")
    print(.family_b_reference_run(args[1], if (length(args)>1) args[2] else NULL))
  }
}
