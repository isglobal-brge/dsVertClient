#!/usr/bin/env Rscript
# Source-tree synthetic integration test. This deliberately uses a plaintext
# Go TEST executable, never the production MPC command surface. It validates
# family arithmetic, signed R contracts, DSLite transport and DP selection.
# It does not certify the pending fused producer/PSI/result-evidence wiring.
args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) {
  cat("Source-tree synthetic NB/LASSO DSLite reference harness.\n",
      "Usage: Rscript tools/validate_fama_reference_dslite.R [--output=PATH]\n",
      "Builds Go tests with dsvert_family_reference; no production fallback.\n")
  quit(save = "no")
}
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L])
root <- normalizePath(file.path(dirname(script), "..", ".."), mustWork = TRUE)
server_root <- file.path(root, "dsVert")
client_root <- file.path(root, "dsVertClient")
go_root <- file.path(server_root, "inst", "dsvert-mpc")
for (package in c("DSLite", "DSI", "pkgload", "jsonlite", "MASS", "glmnet", "processx")) {
  if (!requireNamespace(package, quietly = TRUE)) stop("Missing reference-test dependency: ", package)
}
pkgload::load_all(server_root, quiet = TRUE)
pkgload::load_all(client_root, quiet = TRUE)
message("Family A synthetic reference: package namespaces loaded")
source(file.path(server_root, "tests", "testthat", "helper-nb-grid-cross.R"))
lasso_helpers <- new.env(parent = asNamespace("dsVertClient"))
sys.source(file.path(client_root, "tests", "testthat", "helper-lasso-cross.R"), envir = lasso_helpers)
ns <- function(name, package = "dsVertClient") getFromNamespace(name, package)
work <- tempfile("fama-reference-")
dir.create(work)
run <- function() {
  binary <- file.path(work, "family-reference.test")
  # Build a minimal test package so Go applies the build tag normally. Passing
  # explicitly named .go files would bypass build constraints for those files.
  reference_source <- file.path(work, "source")
  dir.create(reference_source)
  files <- c("go.mod", "go.sum", "k2_exact_gc_nb_loss.go", "k2_exact_gc_nb_loss_test.go",
    "cross_grid_contract_v1.go", "cross_grid_reference_v1_test.go",
    "k2_exact_gc_lasso_loss_reference_cli_test.go")
  stopifnot(all(file.copy(file.path(go_root, files), reference_source)))
  build <- processx::run("go", c("test", "-c", "-tags", "dsvert_family_reference", "-o", binary, "."),
    wd = reference_source, error_on_status = FALSE, timeout = 600)
  if (build$status != 0L) stop(build$stderr)
  # _test.go remains excluded even when the tag is enabled in a normal build.
  production_files <- processx::run("go", c("list", "-tags", "dsvert_family_reference",
    "-f", "{{join .GoFiles \"\\n\"}}"), wd = go_root)$stdout
  stopifnot(!grepl("loss_reference_cli_test", production_files, fixed = TRUE))
  untagged_tests <- processx::run("go", c("list", "-f",
    "{{join .TestGoFiles \"\\n\"}}"), wd = go_root)$stdout
  stopifnot(!grepl("loss_reference_cli_test", untagged_tests, fixed = TRUE))
  fails_closed <- function(expression) {
    failure <- tryCatch({ force(expression); NULL }, error = identity)
    stopifnot(inherits(failure, "dsvert_dp_public_failure"))
    TRUE
  }
  reference <- function(request) {
    input <- file.path(work, "input.json"); output <- file.path(work, "output.json")
    unlink(output)
    jsonlite::write_json(request, input, auto_unbox = TRUE, digits = NA)
    env <- c(DSVERT_FAMA_SYNTHETIC_INPUT = input, DSVERT_FAMA_SYNTHETIC_OUTPUT = output)
    result <- processx::run(binary, "-test.run=^TestFamilyAReferenceCommand$", env = env,
      wd = go_root, error_on_status = FALSE, timeout = 900)
    if (result$status != 0L || !file.exists(output)) stop("Tagged reference failed: ", result$stdout)
    value <- jsonlite::read_json(output, simplifyVector = TRUE)
    stopifnot(isTRUE(value$reference_only), identical(value$family, request$family))
    value$coordinates
  }
  n <- 4096L; bits <- 18; scale <- 2^bits
  set.seed(20260918)
  x <- matrix(sample(c(0, 0.25, 0.5, 0.75, 1), n * 2, replace = TRUE), ncol = 2)
  colnames(x) <- c("x", "z")
  design <- cbind(1, x)
  beta <- list(c(-0.5, 1, 0.5), c(0, 0, 0), c(1, -1, -0.5))
  records <- list()
  for (family in c("nb", "binomial", "poisson", "gaussian")) {
    message("Family A synthetic reference: ", family)
    eta <- as.vector(design %*% beta[[1L]])
    maximum <- if (family %in% c("binomial", "gaussian")) 1 else 32
    y <- switch(family,
      nb = pmin(rnbinom(n, mu = exp(eta), size = 2), maximum),
      binomial = rbinom(n, 1, plogis(eta)),
      poisson = pmin(rpois(n, exp(eta)), maximum),
      gaussian = pmin(1, pmax(0, 0.2 + 0.4*x[,1] + 0.2*x[,2] + rnorm(n, sd = 0.025))))
    grid <- if (family == "gaussian") list(c(0,0,0),c(0.2,0.4,0.2),c(0.8,-0.4,-0.2)) else beta
    fixture <- .nb_grid_cross_fixture(capacity = n, bits = bits, beta_grid = grid,
      theta_grid = rep(2, length(grid)), max_outcome = maximum)
    if (family %in% c("binomial", "poisson")) {
      raw <- fixture$raw; raw$theta_grid <- NULL
      raw$version <- paste0(family, "_grid_cross_v1")
      spec <- ns(".dsvert_dp_glm_grid_cross_spec", "dsVert")(raw, fixture$policy, fixture$authenticated)
      artifact <- ns(".dsvert_dp_glm_grid_cross_artifact", "dsVert")(spec)
      fixture$contract <- fixture$sign(list(version = fixture$contract$version, spec = spec,
        artifact = artifact, source_contract = ns(".dsvert_dp_glm_grid_cross_source_contract", "dsVert")(spec, artifact)))
    }
    spec <- fixture$contract$spec
    grid <- lapply(spec$beta_grid, unlist, use.names = FALSE)
    # Both DSLite custodians validate NB/binomial/Poisson signed contracts, then
    # persist only their synthetic local columns for the tagged test oracle.
    # Gaussian uses the NB fixture for this transport stage only; its existing
    # Gaussian artifact validator and moment geometry are checked below. Its
    # signed Gaussian server adapter is still an explicit integration gate.
    # No such aggregate method is added to either package's runtime allowlist.
    references <- character()
    conns <- NULL
    peers <- c("peer_a", "peer_b")
    for (index in seq_along(peers)) {
      server <- DSLite::newDSLiteServer(tables = list(synthetic = data.frame(id = seq_len(n))))
      # DSLite rebinds environments and serializes method source. Preserve exact
      # binary64 metadata with an R serialization inside a constant string;
      # deparsing a numeric list can otherwise change a signed floating value.
      local_values <- if (index == 1L) list(id = seq_along(y),
        feature = sprintf("%.0f", round(x[,1] * 2^50)),
        outcome = sprintf("%.0f", if (family == "gaussian") round(y*2^50) else y)) else
        list(id = seq_along(y), feature = sprintf("%.0f", round(x[,2] * 2^50)))
      path <- file.path(work, paste0("synthetic-peer-", index, ".json"))
      payload <- jsonlite::base64_enc(serialize(list(
        validator = if (family %in% c("nb", "gaussian"))
          ".dsvert_dp_nb_grid_cross_contract_validate" else
          ".dsvert_dp_glm_grid_cross_contract_validate",
        contract = fixture$contract, policy = fixture$policy, schema = fixture$schema,
        local_values = local_values, path = path, owner = index), NULL, xdr = TRUE))
      method <- eval(substitute(
        function() {
          config <- unserialize(jsonlite::base64_dec(PAYLOAD))
          validator <- getFromNamespace(config$validator, "dsVert")
          validator(config$contract, config$policy, config$schema)
          jsonlite::write_json(config$local_values, config$path, auto_unbox = TRUE, digits = NA)
          list(reference_only = TRUE, owner = config$owner,
               source_sha256 = digest::digest(file = config$path, algo = "sha256"))
        }, list(PAYLOAD = payload)))
      server$aggregateMethod("familySyntheticReferenceStageDS", method)
      object <- paste0("fama_reference_", Sys.getpid(), "_", index)
      assign(object, server, envir = .GlobalEnv)
      references[[peers[[index]]]] <- object
    }
    stage <- tryCatch({
      builder <- DSI::newDSLoginBuilder()
      for (peer in peers) builder$append(server = peer, url = references[[peer]],
        table = "synthetic", driver = "DSLiteDriver")
      conns <- DSI::datashield.login(builder$build(), assign = FALSE)
      DSI::datashield.aggregate(conns, quote(familySyntheticReferenceStageDS()))
    }, error = function(error) {
      # Every object in this source-only harness is generated synthetic data.
      # Preserve DSLite diagnostics before logout clears the failed session.
      diagnostic <- paste(capture.output(print(DSI::datashield.errors())), collapse = "\n")
      stop("Synthetic DSLite reference staging failed:\n", diagnostic, call. = FALSE)
    }, finally = {
      if (!is.null(conns)) DSI::datashield.logout(conns)
      rm(list = unname(references), envir = .GlobalEnv)
    })
    stopifnot(length(stage) == 2L, all(vapply(stage, function(value) isTRUE(value$reference_only), logical(1L))))
    first <- jsonlite::read_json(file.path(work, "synthetic-peer-1.json"), simplifyVector = TRUE)
    second <- jsonlite::read_json(file.path(work, "synthetic-peer-2.json"), simplifyVector = TRUE)
    stopifnot(identical(first$id, second$id))
    features <- cbind(first$feature, second$feature)
    caps <- vapply(spec$sensitivity$candidate_bounds, `[[`, numeric(1L), "per_patient_cap")
    request <- list(family = family, g = bits,
      features = lapply(seq_len(n), function(i) as.list(features[i,])),
      outcomes = as.list(first$outcome),
      validity = rep(list(as.list(rep(1L, 4))), n), beta = spec$beta_encoded,
      caps = as.list(caps), theta_exponent = as.list(rep(1L, length(grid))))
    coordinates <- reference(request)
    real_loss <- function(b) {
      linear <- as.vector(design %*% b)
      switch(family, nb = sum(-dnbinom(y, mu = exp(linear), size = 2, log = TRUE)),
        binomial = sum(pmax(linear,0)+log1p(exp(-abs(linear)))-y*linear),
        poisson = sum(exp(linear)-y*linear+lgamma(y+1)),
        gaussian = 0.5*sum((y-linear)^2))
    }
    exact <- vapply(grid, real_loss, numeric(1L))
    lambda <- if (family == "nb") 0 else 0.01
    penalty <- vapply(grid, function(b) ns(".dsvert_dp_lasso_cross_penalty")(b, lambda, n, bits), numeric(1L))
    if (family == "gaussian") {
      dx <- round(design*scale); dy <- round(y*scale)
      expected <- c(n*scale, unlist(lapply(seq_len(ncol(dx)), function(right) {
        vapply(seq_len(right), function(left) sum(floor(dx[,left]*dx[,right]/scale)), numeric(1L))
      })), vapply(seq_len(ncol(dx)), function(k) sum(floor(dx[,k]*dy/scale)), numeric(1L)),
        sum(floor(dy^2/scale)))
      stopifnot(identical(as.numeric(coordinates), as.numeric(expected)))
      artifact <- lasso_helpers$.lasso_cross_gaussian_artifact(n, bits)
      ns(".dsvert_dp_gaussian_cross_artifact")(artifact, "cohort", "loss",
        "peer_a", "add_remove_patient", scale, n)
      descriptor <- ns(".dsvert_dp_lasso_cross_base")(artifact, fixture$policy, "gaussian")
      delta1 <- length(coordinates)*scale
    } else {
      row_error <- if (family == "nb") (maximum + 2) * 0.00003936 +
        1024 * 1.4655e-14 + 2^-63 + 0.5/scale else 0.75/scale
      stopifnot(max(abs(coordinates/scale-exact)) <= n*row_error)
      descriptor <- if (family == "nb") NULL else
        ns(".dsvert_dp_lasso_cross_base")(fixture$contract, fixture$policy)
      delta1 <- spec$sensitivity$raw_l1_sensitivity
    }
    lasso_spec <- if (family == "nb") NULL else ns(".dsvert_dp_lasso_cross_spec")(
      list(version = "lasso_grid_cross_v1", analysis_id = "synthetic_lasso",
        candidate_grid = lapply(grid, function(b) list(lambda = lambda, beta = as.list(b)))), descriptor)
    if (family == "nb") {
      ns(".dsvert_dp_nb_grid_cross_contract_validate")(fixture$contract,
        fixture$policy, fixture$schema)
      production_failed_closed <- fails_closed(ns("dp_nb_grid")(
        peer_a$y ~ peer_a$x + peer_b$z, "cohort", spec$analysis_id,
        fixture$contract, fixture$policy, fixture$schema))
      missing_signature <- fixture$contract
      missing_signature$signatures$peer_b <- NULL
      missing_signature_failed_closed <- fails_closed(
        ns(".dsvert_dp_nb_grid_cross_contract_validate")(
          missing_signature, fixture$policy, fixture$schema))
    } else {
      unsigned <- list(version = "dsvert-cross-owner-lasso-signed-contract-v1",
                       spec = lasso_spec)
      message <- ns(".dsvert_dp_lasso_cross_message")(unsigned)
      keys <- get("keys", envir = environment(fixture$sign), inherits = FALSE)
      unsigned$signatures <- lapply(keys, function(key) sub("=+$", "", chartr("+/", "-_",
        gsub("[\r\n]", "", jsonlite::base64_enc(openssl::ed25519_sign(message, key))))))
      lasso_contract <- unsigned
      production_failed_closed <- fails_closed(ns("dp_lasso_grid")(
        c("peer_a$x", "peer_b$z"), lasso_contract, fixture$contract,
        fixture$policy, fixture$schema))
      # Gaussian production wiring is deliberately unavailable. The synthetic
      # descriptor above has already passed its existing artifact validator;
      # only this source test injects that public descriptor for signed L1 checks.
      base_validator <- if (family == "gaussian") function(...) descriptor else
        ns(".dsvert_dp_lasso_cross_validate_base")
      ns(".dsvert_dp_lasso_cross_contract_validate")(lasso_contract,
        fixture$policy, fixture$schema, fixture$contract,
        .base_validator = base_validator)
      missing_signature <- lasso_contract
      missing_signature$signatures$peer_b <- NULL
      missing_signature_failed_closed <- fails_closed(
        ns(".dsvert_dp_lasso_cross_contract_validate")(
          missing_signature, fixture$policy, fixture$schema, fixture$contract,
          .base_validator = base_validator))
    }
    central <- if (family == "nb") MASS::glm.nb(y ~ x[,1]+x[,2]) else
      glmnet::glmnet(x, y, family = family, alpha = 1, lambda = lambda,
        standardize = FALSE, intercept = TRUE, thresh = 1e-12)
    central_beta <- if (family == "nb") unname(coef(central)) else as.numeric(coef(central))
    common_objective <- real_loss(central_beta)/n+lambda*sum(abs(central_beta[-1L]))
    central_fit_objective <- if (family == "nb") -mean(dnbinom(y,
      mu = exp(as.vector(design %*% central_beta)), size = central$theta, log = TRUE)) else
      common_objective
    exact_best <- which.min(exact+penalty/scale)
    for (epsilon in c(1,4,8)) {
      # Each independent synthetic authority contributes a complete calibrated
      # discrete-Laplace draw. Their sum is conservative epsilon-DP; this is a
      # statistical reference, not a test of the production joint sampler.
      noise <- function(seed) {
        set.seed(seed)
        prob <- -expm1(-epsilon/delta1)
        rgeom(length(coordinates), prob)-rgeom(length(coordinates), prob)
      }
      noisy <- coordinates+noise(701L+epsilon)+noise(1709L+epsilon)
      maxima <- if (family == "gaussian") unlist(descriptor$maximum_coordinates) else
        unlist(spec$sensitivity$maximum_coordinates)
      released <- pmax(0, pmin(maxima, noisy))
      selected <- if (family == "nb")
        ns(".dsvert_dp_nb_grid_cross_moment")(released, spec) else
        ns(".dsvert_dp_lasso_cross_postprocess")(released, lasso_spec)
      chosen <- if (family == "nb") selected$selected_candidate else
        selected$selected_candidates[[1L]]
      stopifnot(length(chosen) == 1L, chosen %in% seq_along(grid))
      stopifnot(is.null(selected$standard_errors))
      records[[length(records)+1L]] <- list(family = family, epsilon = epsilon,
        delta_cap = 2^-100, reference_noise_delta = 0, peers = 2L, capacity = n,
        selected = chosen, exact_best = exact_best, agrees_with_exact = chosen == exact_best,
        production_failed_closed = production_failed_closed,
        missing_signature_failed_closed = missing_signature_failed_closed,
        profile_route = if (family == "nb") "certified_q16_piecewise_quadratic64" else
          if (family == "gaussian") "existing_cross_owner_sufficient_statistics" else
            "frozen_step1_reference_pending_step2_piecewise",
        server_contract_route = if (family == "gaussian")
          "nb_synthetic_stage_only_gaussian_server_adapter_pending" else
            paste0(family, "_signed_server_contract_validator"),
        maximum_profile_loss_error = if (family == "gaussian") NULL else n*row_error,
        profile_error_to_noise_scale = if (family == "gaussian") NULL else
          n*row_error/(delta1/scale/epsilon),
        central_coefficients = as.list(central_beta),
        central_theta = if (family == "nb") central$theta else NULL,
        selected_coefficients = as.list(grid[[chosen]]),
        central_l1_distance = sum(abs(grid[[chosen]]-central_beta)),
        central_objective = central_fit_objective,
        common_signed_theta_objective = if (family == "nb") common_objective else NULL,
        selected_exact_objective = (exact[chosen]+penalty[chosen]/scale)/n,
        replay_identical = identical(noisy, coordinates+noise(701L+epsilon)+noise(1709L+epsilon)))
    }
  }
  list(reference_only = TRUE, production_release_validated = FALSE,
    validation_context = list(r_version = R.version.string, platform = R.version$platform,
      logical_cores = parallel::detectCores(),
      go_version = trimws(processx::run("go", "version")$stdout),
      test_build_tag = "dsvert_family_reference",
      f50_json_encoding = "canonical_decimal_strings",
      noise_validation = "fixed_seed_synthetic_reference_only",
      nb_profile_sha256 = ns(".dsvert_dp_nb_grid_cross_profile", "dsVert")()$profile_sha256),
    pending = c("fused producer", "PSI runtime", "authenticated joint DP evidence",
      "binomial/Poisson piecewise producer", "Gaussian signed server adapter"), cases = records)
}
result <- tryCatch(run(), finally = unlink(work, recursive = TRUE))
output <- sub("^--output=", "", grep("^--output=", args, value = TRUE))
if (length(output)) jsonlite::write_json(result, output[[1L]], auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null")
cat(jsonlite::toJSON(result, auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null"), "\n")
