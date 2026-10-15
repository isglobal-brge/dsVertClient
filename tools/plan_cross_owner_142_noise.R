# Public synthetic planning only: no DSI calls, custodians, protocol seams or
# release artifacts. Run after sourcing the documented isolated env.sh:
# timeout 900 Rscript tools/plan_cross_owner_142_noise.R [output.json]
# This estimates finite Cor/PCA availability, not statistical precision.

cross_owner_142_noise_parameters <- function() {
  scale <- 256
  sensitivity <- 29 * scale
  epsilon <- 8
  q <- exp(-epsilon / sensitivity)
  p <- -expm1(-epsilon / sensitivity)
  list(scale = scale, sensitivity = sensitivity, epsilon = epsilon,
       q = q, geometric_probability = p)
}

cross_owner_142_noise_scope_check <- function() {
  references <- c("a", "b", "site_b$x", "site_c$x", "unused_d", "unused_e")
  owners <- c("site_a", "site_a", "site_b", "site_c", "site_d", "site_e")
  columns <- stats::setNames(lapply(seq_along(references), function(i) list(
    kind = "numeric", column = sub(".*[$]", "", references[[i]]),
    dataset = "DA", owner_peer = owners[[i]], lower = 0, upper = 1)), references)
  policy <- list(datasets = list(DA = list()), numeric_bounds = stats::setNames(
    rep(list(c(0, 1)), length(references)), references))
  scope <- get(".dsvert_dp_capsule_workload_scope", asNamespace("dsVert"))(
    list(mode = "catalog_v1", numeric_moments = character(),
      categorical_marginals = character(), categorical_pairs = list(), correlations = list()),
    columns, policy, list(), list(), list(cor_primary = list(
      version = "v2", dataset = "DA", outcome = "a",
      predictors = references[2:4], intercept = TRUE)), list())
  stopifnot(identical(scope$numeric_moments, references[1:4]))
  list(kind = "offline production scope function on public metadata; no signature or release claim",
    included_numeric = scope$numeric_moments,
    excluded_unreferenced_numeric = setdiff(references, scope$numeric_moments),
    resulting_coordinate_count = 1L + 3L * length(scope$numeric_moments) + 16L)
}

cross_owner_142_noise_summary <- function(results, elapsed_seconds,
    origin = "fresh public Monte Carlo") {
  parameters <- cross_owner_142_noise_parameters()
  ns <- asNamespace("dsVertClient")
  tail <- get(".dsvert_dp_vector_convolution_log_tail", ns)
  for (index in seq_along(results)) {
    row <- results[[index]]
    row$finite_fraction <- row$finite / row$trials
    row$finite_fraction_binomial_95_interval <- unname(
      binom.test(row$finite, row$trials)$conf.int)
    # Symmetry converts the exact two-sided tail into P(noise <= -n*scale).
    row$count_nonpositive_probability <- exp(tail(
      row$n * parameters$scale - 1, parameters$epsilon,
      parameters$sensitivity)) / 2
    # A strict radius below n/10 gives positive 2x2 principal moment blocks.
    # The ten one-sided events are two count/intercept lower tails, four
    # second-moment lower tails and four first-moment upper tails.
    radius_steps <- floor((row$n * parameters$scale - 1) / 10)
    row$sufficient_identifiability_failure_union_bound <- min(1, 5 * exp(
      tail(radius_steps, parameters$epsilon, parameters$sensitivity)))
    results[[index]] <- row
  }
  radius <- function(dimension) {
    low <- 0; high <- 100000
    while (high - low > 1) {
      middle <- floor((low + high) / 2)
      if (tail(middle, parameters$epsilon, parameters$sensitivity) <=
          log(0.05 / dimension)) high <- middle else low <- middle
    }
    high / parameters$scale
  }
  list(
    version = "dsvert-1.4.2-public-noise-planning-v1",
    classification = "public synthetic planning; no federation evidence",
    evidence_origin = origin,
    package_version = as.character(utils::packageVersion("dsVertClient")),
    r_version = R.version.string, seed = 20261015L, rng_kind = RNGkind(),
    scope = list(variables = c("site_a$a", "site_a$b", "site_b$x", "site_c$x"),
      data_owners = 3L, variable_count = 4L, bounds = c(0, 1),
      adjacency = "add_remove_patient", epsilon = 8, delta = 0,
      numeric_grid_bits = 8L, explicit_catalog = "empty except Gaussian spec",
      automatic_numeric_moment_coordinates = 12L,
      gaussian_coordinates = 16L, admitted_count_coordinates = 1L,
      total_coordinates = 29L, natural_l1_sensitivity = 29,
      integer_l1_sensitivity_steps = parameters$sensitivity),
    noise = list(
      mechanism = "two-independent-complete-vector-discrete-laplace-draws-v4",
      per_authority_natural_scale = 29 / 8,
      per_coordinate_standard_deviation =
        sqrt(4 * parameters$q) / parameters$geometric_probability / parameters$scale,
      q = parameters$q, geometric_probability = parameters$geometric_probability,
      simulation = "G1-G2+G3-G4; independent zero-based geometric variables; divide by 256",
      sampler_implementation_delta = "0/1",
      note = "Ideal exact law, up to the negligible signed Ring128 representability bound; public R RNG is used only for planning"),
    observed_signed_reference = list(
      file = ".tmp/vertcor/smoke-k3-certificate.rds", variable_count = 3L,
      capacity = 4L, total_coordinates = 21L,
      natural_l1_sensitivity = 21, integer_l1_sensitivity_steps = 5376,
      explanation = "The signed reference includes nine automatic numeric-moment coordinates plus eleven Gaussian coordinates plus admitted count"),
    simulation = list(trials_per_n = 50000L, order = c(16L, 32L, 64L, 128L, 256L),
      exact_input = "Four balanced mutually orthogonal binary columns, bits 0 through 3 of row index minus one",
      implementation = "Unmodified .dsvert_dp_gaussian_unpack, then the production Cor projected-mass and centered-diagonal tolerance checks",
      unpack_body_sha256 = digest::digest(paste(deparse(body(
        get(".dsvert_dp_gaussian_unpack", ns))), collapse = "\n"),
        algo = "sha256", serialize = FALSE),
      elapsed_seconds = elapsed_seconds, results = results),
    uncertainty = list(
      gaussian_block_simultaneous_95_sum_radius = radius(16L),
      entire_release_simultaneous_95_sum_radius = radius(29L),
      conservative_quantization_sum_allowance_at_n64 = 64 * (2 + 1 / (4 * 256)) / 256,
      actual_binary_fixture_quantization_error = 0,
      sufficient_condition = "For strict r<n/10, count/intercept >= n-r, each sum <= n/2+r and each square sum >= n/2-r imply positive principal determinants: n^2/4-(5/2)*n*r>0. Eigenvalue clipping adds a PSD matrix, preserving those positive principal blocks.",
      interpretation = "Binomial intervals quantify Monte Carlo sampling uncertainty under this public model. Finite output does not establish accurate correlations, informative confidence regions, stable PCA directions or a guarantee for a future sticky release."),
    decision = list(n = 64L, unit_capacity = 64L,
      reason = "No finite-output failures in 50000 public simulations; stronger margin than n32 while avoiding the larger exact-MPC cost of n128/n256"),
    k5_membership_extension = list(data_owners = 3L, federation_members = 5L,
      local_schema_requirement = "site_d/e each need an unused bounded public synthetic column because bootstrap rejects empty local schemas",
      release_scope = "Do not list site_d/e columns in catalog, correlations or Gaussian spec; catalog_v1 includes only explicit or referenced numeric columns",
      expected_total_coordinates = 29L, expected_integer_l1_sensitivity_steps = 7424L,
      offline_scope_check = cross_owner_142_noise_scope_check(),
      requirement = "Check the actual signed manifest before release; this is not a five-data-owner fit"))
}

cross_owner_142_noise_simulate <- function() {
  parameters <- cross_owner_142_noise_parameters()
  unpack <- get(".dsvert_dp_gaussian_unpack", asNamespace("dsVertClient"))
  artifact <- list(design_terms = c("(Intercept)", "b", "site_b$x", "site_c$x"))
  trials <- 50000L
  set.seed(20261015)
  started <- proc.time()[["elapsed"]]
  results <- list()
  for (n in c(16L, 32L, 64L, 128L, 256L)) {
    index <- seq_len(n) - 1L
    variables <- sapply(0:3, function(bit)
      as.numeric(bitwAnd(index, bitwShiftL(1L, bit)) != 0L))
    y <- variables[, 1L]; design <- cbind(1, variables[, -1L])
    gram <- crossprod(design)
    truth <- c(n, gram[upper.tri(gram, diag = TRUE)],
      drop(crossprod(design, y)), sum(y * y))
    stopifnot(length(truth) == 16L)
    p <- parameters$geometric_probability
    noises <- matrix(rgeom(trials * 16L, p) - rgeom(trials * 16L, p) +
      rgeom(trials * 16L, p) - rgeom(trials * 16L, p), nrow = trials) / parameters$scale
    finite <- bad_count <- bad_mass <- bad_variance <- 0L
    minimum_variance <- numeric(trials)
    for (i in seq_len(trials)) {
      coordinates <- truth + noises[i, ]
      if (coordinates[[1L]] <= 0) {
        bad_count <- bad_count + 1L; minimum_variance[[i]] <- NA_real_; next
      }
      augmented <- unpack(coordinates, artifact, n)$augmented_projected
      mass <- augmented[1L, 1L]
      tolerance <- 256 * .Machine$double.eps * max(1, abs(mass)) * 4L
      if (!is.finite(mass) || mass <= tolerance) {
        bad_mass <- bad_mass + 1L; minimum_variance[[i]] <- NA_real_; next
      }
      selected <- 2:5; sums <- augmented[1L, selected]
      centered <- augmented[selected, selected, drop = FALSE] - outer(sums, sums) / mass
      centered <- (centered + t(centered)) / 2
      variance_tolerance <- 256 * .Machine$double.eps * max(1, max(abs(centered))) * 4L
      minimum_variance[[i]] <- min(diag(centered))
      if (any(!is.finite(diag(centered))) || any(diag(centered) <= variance_tolerance)) {
        bad_variance <- bad_variance + 1L; next
      }
      finite <- finite + 1L
    }
    results[[length(results) + 1L]] <- list(n = n, trials = trials, finite = finite,
      count_failures = bad_count, mass_failures = bad_mass, variance_failures = bad_variance,
      minimum_variance_quantiles = unname(quantile(minimum_variance,
        c(0, 0.01, 0.05, 0.5), na.rm = TRUE)),
      cumulative_elapsed_seconds = proc.time()[["elapsed"]] - started)
    message("Public noise planning n=", n, ": ", finite, "/", trials, " finite")
  }
  cross_owner_142_noise_summary(results, proc.time()[["elapsed"]] - started)
}

if (sys.nframe() == 0L) {
  suppressPackageStartupMessages(library(dsVertClient))
  args <- commandArgs(TRUE)
  output <- if (length(args)) args[[1L]] else
    "inst/validation/v1.4.2/vertcor04-noise-planning.json"
  result <- cross_owner_142_noise_simulate()
  writeLines(jsonlite::toJSON(result, auto_unbox = TRUE, pretty = TRUE,
    digits = 16, null = "null"), output)
}
