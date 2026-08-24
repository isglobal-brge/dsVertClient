.synopsis_describe_real_e2e_server <- function() {
  skip_if_not_installed("pkgload")
  configured_server <- Sys.getenv("DSVERT_SERVER_SOURCE", unset = "")
  server_path <- if (nzchar(configured_server)) {
    normalizePath(configured_server, mustWork = FALSE)
  } else {
    normalizePath(file.path(
      test_path(), "..", "..", "..", "dsVert"), mustWork = FALSE)
  }
  skip_if_not(dir.exists(server_path), "requires the sibling dsVert source")
  pkgload::load_all(server_path, quiet = TRUE)
  asNamespace("dsVert")
}

# The empty default runs the complete release battery.  A named family is a
# developer focal gate and never changes the production route.
.synopsis_real_e2e_family <- Sys.getenv(
  "DSVERT_TEST_SYNOPSIS_E2E_FAMILY", unset = "")
.synopsis_real_e2e_families <- c(
  "describe", "same_owner_contingency", "cross_owner_contingency",
  "stratified_epidemiology", "causal_standardization", "frequency",
  "mi", "mantel_haenszel", "roc", "survival", "correlation", "gaussian", "cross_owner_tamper",
  "gaussian_lasso_focal", "lmm")
if (nzchar(.synopsis_real_e2e_family) &&
    !.synopsis_real_e2e_family %in% .synopsis_real_e2e_families) {
  stop("unknown DSVERT_TEST_SYNOPSIS_E2E_FAMILY", call. = FALSE)
}

.synopsis_real_e2e_only <- function(family) {
  if (nzchar(.synopsis_real_e2e_family) &&
      !identical(.synopsis_real_e2e_family, family)) {
    skip(paste("focused on", .synopsis_real_e2e_family))
  }
}

.synopsis_real_e2e_focal_only <- function() {
  if (!identical(.synopsis_real_e2e_family, "gaussian_lasso_focal")) {
    skip("set DSVERT_TEST_SYNOPSIS_E2E_FAMILY=gaussian_lasso_focal")
  }
}

# The default remains the complete K=2/K=3/K=5 release gate.  This selector
# only lets a developer isolate one topology while diagnosing a slow failure.
.synopsis_real_e2e_peer_counts <- function(
    focused = Sys.getenv("DSVERT_TEST_SYNOPSIS_E2E_K", unset = "")) {
  if (!is.character(focused) || length(focused) != 1L || is.na(focused)) {
    stop("DSVERT_TEST_SYNOPSIS_E2E_K must be a string", call. = FALSE)
  }
  if (!nzchar(focused)) return(c(2L, 3L, 5L))
  values <- suppressWarnings(as.integer(
    strsplit(focused, ",", fixed = TRUE)[[1L]]))
  if (!length(values) || anyNA(values) || anyDuplicated(values) ||
      any(!values %in% c(2L, 3L, 5L))) {
    stop("DSVERT_TEST_SYNOPSIS_E2E_K must be a unique subset of 2,3,5",
         call. = FALSE)
  }
  values
}

test_that("the Synopsis real-E2E topology selector preserves the full gate", {
  expect_identical(.synopsis_real_e2e_peer_counts(""), c(2L, 3L, 5L))
  expect_identical(.synopsis_real_e2e_peer_counts("2,5"), c(2L, 5L))
  expect_error(.synopsis_real_e2e_peer_counts("4"), "unique subset")
})

.synopsis_describe_real_e2e_fixture <- function(k, server_ns) {
  get_server <- function(name) get(name, envir = server_ns, inherits = FALSE)
  b64url <- function(raw) chartr(
    "+/", "-_", sub("=+$", "", jsonlite::base64_enc(raw)))
  peers <- c("peer_a", "peer_b", if (k > 2L) {
    paste0("witness_", seq_len(k - 2L))
  } else character())
  identities <- stats::setNames(lapply(seq_along(peers), function(index) {
    get_server(".callMpcTool")("derive-identity", list(
      seed = jsonlite::base64_enc(as.raw((seq_len(32L) + 37L * index) %% 256L))))
  }), peers)
  pins <- vapply(identities, function(value) {
    b64url(jsonlite::base64_dec(value$identity_pk))
  }, character(1L))
  root <- tempfile("synopsis-describe-rock-e2e-")
  dir.create(root, mode = "0700", recursive = TRUE)
  common <- list(
    schema_version = 1L,
    mechanism_version = "dsvert-dp-v7-contingency-unit-aggregation-1",
    policy_contract = get_server(".DSVERT_DP_SYNOPSIS_POLICY_CONTRACT"),
    domain = "synopsis-describe-rock-e2e", cohort_id = "cohort-v1",
    logical_peers = peers, peer_pinset = pins,
    peer_pinset_sha256 = get_server(".dsvert_dp_synopsis_pinset_hash_v1")(pins),
    peer_count = as.integer(k), designated_noise_peers = peers[1:2],
    global_total_epsilon = 1, global_total_delta = 1e-6,
    adjacency = "add_remove_patient", patient_column = "patient_id",
    unit_capacity = 200L, fixed_cohort_size = NULL,
    max_records_per_unit = 1L, overflow_policy = "reject_snapshot",
    contingency_unit_aggregation_policy = "consistent_cell_else_exclude_v1",
    numeric_grid_bits = 8L,
    noise_selection = list(version = "synopsis-describe-real-e2e-v1"),
    transcript_privacy = "test-only", snapshot_binding = "internal_test_bypass",
    alignment_binding = "internal_test_bypass", require_snapshot_digest = FALSE,
    require_alignment_manifest = FALSE, categorical_levels = list(),
    capsule_dataset_mapping = NULL, lock_timeout_ms = 30000L,
    state_private = TRUE)
  policies <- stats::setNames(lapply(peers, function(peer) {
    variable <- paste0("x_", peer)
    data_name <- paste0("data_", peer)
    specs <- list(describe = list(), survival = list(), gaussian = list(),
                  vertical_cross = list())
    if (identical(peer, "peer_a")) specs$describe$primary <- list(
      version = "v1", dataset = data_name, variables = variable,
      histogram_grids = stats::setNames(list(c(5, 10)), variable),
      allocation = c(count = 0.25, sum = 0.25, sumsq = 0.25,
                     histogram = 0.25))
    value <- c(common, list(
      peer_name = peer, own_identity_pk = unname(pins[[peer]]),
      datasets = stats::setNames(list(list(
        id = paste0("cohort-", peer), version = "v1",
        snapshot_sha256 = NULL, alignment_manifest_hash = NULL,
        alignment_manifest_version = 1L)), data_name),
      numeric_bounds = stats::setNames(list(c(0, 10)), variable),
      capsule_workload_scope = list(
        mode = "catalog_v1", numeric_moments = "x_peer_a",
        categorical_marginals = character(), categorical_pairs = list(),
        correlations = list()),
      synopsis_state_path = file.path(root, peer, "synopsis")))
    value$capsule_dataset_mapping <- stats::setNames(list(variable), data_name)
    value$capsule_workload_specs <- specs
    value
  }), peers)
  for (peer in peers) dir.create(dirname(policies[[peer]]$synopsis_state_path),
                                  mode = "0700", recursive = TRUE)
  snapshots <- stats::setNames(lapply(peers, function(peer) {
    variable <- paste0("x_", peer)
    data_name <- paste0("data_", peer)
    data <- data.frame(patient_id = paste0("u", seq_len(100L)),
                       stringsAsFactors = FALSE)
    data[[variable]] <- rep(c(0, 10), 50L)
    stats::setNames(list(list(
      data = data,
      dataset = list(public = list(
        data_name = data_name, id = paste0("cohort-", peer), version = "v1",
        alignment_manifest_hash = NULL, alignment_manifest_version = 1L),
        fingerprint = strrep(if (identical(peer, "peer_a")) "a" else "b", 64L)))),
      data_name)
  }), peers)
  state <- new.env(parent = emptyenv())
  state$storage <- stats::setNames(lapply(peers, function(...) {
    new.env(parent = emptyenv())
  }), peers)
  state$source_prepare <- 0L
  state$start <- 0L
  list(root = root, peers = peers, identities = identities, pins = pins,
       policies = policies, snapshots = snapshots,
       secrets = stats::setNames(lapply(seq_along(peers), function(index) {
         as.raw(rep(16L + index, 32L))
       }), peers), state = state, server_ns = server_ns)
}

.synopsis_describe_real_e2e_session <- function(fixture, peer, session_id) {
  storage <- fixture$state$storage[[peer]]
  if (!exists(session_id, envir = storage, inherits = FALSE)) {
    value <- new.env(parent = emptyenv())
    assign(".session_id", session_id, envir = value)
    assign(session_id, value, envir = storage)
  }
  get(session_id, envir = storage, inherits = FALSE)
}

.synopsis_contingency_real_e2e_fixture <- function(k, server_ns) {
  fixture <- .synopsis_describe_real_e2e_fixture(k, server_ns)
  pair_scope <- list(
    mode = "catalog_v1", numeric_moments = "x_peer_a",
    categorical_marginals = character(),
    categorical_pairs = list(c("exposure", "outcome")),
    correlations = list())
  for (peer in fixture$peers) {
    policy <- fixture$policies[[peer]]
    policy$capsule_workload_scope <- pair_scope
    if (identical(peer, "peer_a")) {
      policy$categorical_levels <- list(
        exposure = c("no", "yes"), outcome = c("no", "yes"))
      policy$capsule_dataset_mapping[["data_peer_a"]] <- c(
        "x_peer_a", "exposure", "outcome")
    }
    fixture$policies[[peer]] <- policy
  }
  data <- fixture$snapshots$peer_a[["data_peer_a"]]$data
  data$exposure <- rep(c("no", "yes", "no", "yes"), each = 25L)
  data$outcome <- rep(c("no", "no", "yes", "yes"), each = 25L)
  fixture$snapshots$peer_a[["data_peer_a"]]$data <- data
  fixture
}

.synopsis_stratified_contingency_real_e2e_fixture <- function(k, server_ns) {
  fixture <- .synopsis_describe_real_e2e_fixture(k, server_ns)
  pair_scope <- list(
    mode = "catalog_v1", numeric_moments = "x_peer_a",
    categorical_marginals = character(),
    categorical_pairs = list(c("stratum", "outcome")),
    correlations = list())
  for (peer in fixture$peers) {
    policy <- fixture$policies[[peer]]
    policy$capsule_workload_scope <- pair_scope
    if (identical(peer, "peer_a")) {
      policy$categorical_levels <- list(
        stratum = c("young", "middle", "old"),
        outcome = c("no", "yes"))
      policy$capsule_dataset_mapping[["data_peer_a"]] <- c(
        "x_peer_a", "stratum", "outcome")
    }
    fixture$policies[[peer]] <- policy
  }
  data <- data.frame(
    patient_id = paste0("u", seq_len(120L)),
    x_peer_a = rep(c(0, 10), length.out = 120L),
    stratum = rep(c("young", "middle", "old"), each = 40L),
    outcome = c(rep("no", 30L), rep("yes", 10L),
                rep("no", 25L), rep("yes", 15L),
                rep("no", 20L), rep("yes", 20L)),
    stringsAsFactors = FALSE)
  fixture$snapshots$peer_a[["data_peer_a"]]$data <- data
  fixture
}

.synopsis_mantel_haenszel_real_e2e_fixture <- function(k, server_ns) {
  fixture <- .synopsis_describe_real_e2e_fixture(k, server_ns)
  cells <- c(
    "exposed_event", "exposed_nonevent",
    "unexposed_event", "unexposed_nonevent")
  pair_scope <- list(
    mode = "catalog_v1", numeric_moments = "x_peer_a",
    categorical_marginals = character(),
    categorical_pairs = list(c("stratum", "cell")),
    correlations = list())
  for (peer in fixture$peers) {
    policy <- fixture$policies[[peer]]
    policy$unit_capacity <- 200L
    policy$capsule_workload_scope <- pair_scope
    if (identical(peer, "peer_a")) {
      policy$categorical_levels <- list(
        stratum = c("young", "middle", "old"), cell = cells)
      policy$capsule_dataset_mapping[["data_peer_a"]] <- c(
        "x_peer_a", "stratum", "cell")
    }
    fixture$policies[[peer]] <- policy
  }
  counts <- matrix(
    c(16L, 8L, 6L, 10L,
      20L, 10L, 8L, 12L,
      12L, 16L, 4L, 18L),
    nrow = 3L, byrow = TRUE,
    dimnames = list(c("young", "middle", "old"), cells))
  rows <- lapply(rownames(counts), function(stratum) {
    data.frame(
      stratum = rep(stratum, sum(counts[stratum, ])),
      cell = rep(colnames(counts), counts[stratum, ]),
      stringsAsFactors = FALSE)
  })
  data <- do.call(rbind, rows)
  data$patient_id <- paste0("u", seq_len(nrow(data)))
  data$x_peer_a <- rep(c(0, 10), length.out = nrow(data))
  fixture$snapshots$peer_a[["data_peer_a"]]$data <-
    data[, c("patient_id", "x_peer_a", "stratum", "cell")]
  fixture
}

.synopsis_roc_real_e2e_fixture <- function(k, server_ns) {
  fixture <- .synopsis_describe_real_e2e_fixture(k, server_ns)
  bins <- c("low", "mid_low", "mid_high", "high")
  pair_scope <- list(
    mode = "catalog_v1", numeric_moments = "x_peer_a",
    categorical_marginals = character(),
    categorical_pairs = list(c("disease", "score_bin")),
    correlations = list())
  for (peer in fixture$peers) {
    policy <- fixture$policies[[peer]]
    policy$unit_capacity <- 200L
    policy$capsule_workload_scope <- pair_scope
    if (identical(peer, "peer_a")) {
      policy$categorical_levels <- list(
        disease = c("case", "control"), score_bin = bins)
      policy$capsule_dataset_mapping[["data_peer_a"]] <- c(
        "x_peer_a", "disease", "score_bin")
    }
    fixture$policies[[peer]] <- policy
  }
  counts <- rbind(case = c(5L, 10L, 20L, 35L),
                  control = c(35L, 20L, 10L, 5L))
  rows <- lapply(rownames(counts), function(disease) {
    data.frame(
      disease = rep(disease, sum(counts[disease, ])),
      score_bin = rep(bins, counts[disease, ]),
      stringsAsFactors = FALSE)
  })
  data <- do.call(rbind, rows)
  data$patient_id <- paste0("u", seq_len(nrow(data)))
  data$x_peer_a <- rep(c(0, 10), length.out = nrow(data))
  fixture$snapshots$peer_a[["data_peer_a"]]$data <-
    data[, c("patient_id", "x_peer_a", "disease", "score_bin")]
  fixture
}

.synopsis_causal_contingency_real_e2e_fixture <- function(k, server_ns) {
  fixture <- .synopsis_describe_real_e2e_fixture(k, server_ns)
  arms <- c("young_control", "young_treated", "middle_control",
            "middle_treated", "old_control", "old_treated")
  for (peer in fixture$peers) {
    policy <- fixture$policies[[peer]]
    policy$unit_capacity <- 256L
    policy$capsule_workload_scope <- list(
      mode = "catalog_v1", numeric_moments = "x_peer_a",
      categorical_marginals = character(),
      categorical_pairs = list(c("arm", "outcome")),
      correlations = list())
    if (identical(peer, "peer_a")) {
      policy$categorical_levels <- list(
        arm = arms, outcome = c("no", "yes"))
      policy$capsule_dataset_mapping[["data_peer_a"]] <- c(
        "x_peer_a", "arm", "outcome")
    }
    fixture$policies[[peer]] <- policy
  }
  events <- c(8L, 16L, 10L, 20L, 12L, 24L)
  arm <- rep(arms, each = 40L)
  outcome <- unlist(lapply(events, function(event_count) c(
    rep("yes", event_count), rep("no", 40L - event_count))),
    use.names = FALSE)
  fixture$snapshots$peer_a[["data_peer_a"]]$data <- data.frame(
    patient_id = paste0("u", seq_along(arm)),
    x_peer_a = rep(c(0, 10), length.out = length(arm)),
    arm = arm, outcome = outcome, stringsAsFactors = FALSE)
  fixture
}

.synopsis_cross_contingency_real_e2e_fixture <- function(k, server_ns) {
  fixture <- .synopsis_describe_real_e2e_fixture(k, server_ns)
  get_server <- function(name) get(name, envir = server_ns, inherits = FALSE)
  scope <- list(
    mode = "catalog_v1", numeric_moments = character(),
    categorical_marginals = character(), categorical_pairs = list(),
    correlations = list())
  cross_spec <- list(
    version = "v2", left_dataset = "data_peer_a",
    right_dataset = "data_peer_b", left = "disease", right = "exposure",
    family = "categorical_pair")
  for (peer in fixture$peers) {
    policy <- fixture$policies[[peer]]
    policy$unit_capacity <- 16L
    policy$capsule_workload_scope <- scope
    policy$capsule_workload_specs <- list(
      describe = list(), survival = list(), gaussian = list(),
      vertical_cross = list())
    if (identical(peer, "peer_a")) {
      policy$categorical_levels <- list(disease = c("no", "yes"))
      policy$capsule_dataset_mapping[["data_peer_a"]] <- c(
        "x_peer_a", "disease")
      policy$capsule_workload_specs$vertical_cross$cross_table <- cross_spec
    } else if (identical(peer, "peer_b")) {
      policy$categorical_levels <- list(exposure = c("high", "low"))
      policy$capsule_dataset_mapping[["data_peer_b"]] <- c(
        "x_peer_b", "exposure")
    }
    fixture$policies[[peer]] <- policy
  }
  for (peer in fixture$peers) {
    data_name <- paste0("data_", peer)
    fixture$snapshots[[peer]][[data_name]]$data <-
      fixture$snapshots[[peer]][[data_name]]$data[seq_len(16L), , drop = FALSE]
  }
  fixture$snapshots$peer_a[["data_peer_a"]]$data$disease <-
    rep(c("no", "yes"), each = 8L)
  fixture$snapshots$peer_b[["data_peer_b"]]$data$exposure <-
    rep(c("high", "low"), 8L)
  token <- chartr(
    "+/", "-_", sub("=+$", "",
                       jsonlite::base64_enc(as.raw(seq_len(32L) - 1L))))
  for (peer in fixture$peers) {
    data_name <- paste0("data_", peer)
    aligned <- get_server(".psi_attach_alignment_manifest")(
      fixture$snapshots[[peer]][[data_name]]$data, "patient_id", token)
    alignment <- get_server(".psi_validate_alignment_manifest")(aligned)
    fixture$snapshots[[peer]][[data_name]]$data <- aligned
    fixture$policies[[peer]]$datasets[[data_name]]$alignment_manifest_hash <-
      alignment$hash
    fixture$policies[[peer]]$datasets[[data_name]]$alignment_manifest_version <-
      alignment$version
    fixture$snapshots[[peer]][[data_name]]$dataset$public$
      alignment_manifest_hash <- alignment$hash
    fixture$snapshots[[peer]][[data_name]]$dataset$public$
      alignment_manifest_version <- alignment$version
  }
  fixture
}

.synopsis_frequency_real_e2e_fixture <- function(k, server_ns) {
  fixture <- .synopsis_describe_real_e2e_fixture(k, server_ns)
  scope <- list(
    mode = "catalog_v1", numeric_moments = "x_peer_a",
    categorical_marginals = "status", categorical_pairs = list(),
    correlations = list(), strict_missing_categorical = "status")
  for (peer in fixture$peers) {
    policy <- fixture$policies[[peer]]
    policy$capsule_workload_scope <- scope
    if (identical(peer, "peer_a")) {
      policy$categorical_levels <- list(status = c("case", "control"))
      policy$capsule_dataset_mapping[["data_peer_a"]] <- c(
        "x_peer_a", "status")
    }
    fixture$policies[[peer]] <- policy
  }
  fixture$snapshots$peer_a[["data_peer_a"]]$data$status <-
    rep(c("case", "control"), each = 50L)
  fixture
}

.synopsis_mi_real_e2e_fixture <- function(k, server_ns) {
  fixture <- .synopsis_frequency_real_e2e_fixture(k, server_ns)
  fixture$snapshots$peer_a[["data_peer_a"]]$data$status <- c(
    rep("case", 45L), rep("control", 45L), rep(NA_character_, 10L))
  fixture
}

.synopsis_survival_real_e2e_fixture <- function(k, server_ns) {
  fixture <- .synopsis_describe_real_e2e_fixture(k, server_ns)
  policy <- fixture$policies$peer_a
  policy$numeric_bounds$time_peer_a <- c(0, 10)
  policy$categorical_levels <- list(status_peer_a = c("censor", "event"))
  policy$capsule_dataset_mapping[["data_peer_a"]] <- c(
    "x_peer_a", "time_peer_a", "status_peer_a")
  policy$capsule_workload_specs$survival$primary <- list(
    version = "v1", dataset = "data_peer_a", time = "time_peer_a",
    event = "status_peer_a", censor = "censor", time_grid = c(5, 10),
    entry = NULL)
  fixture$policies$peer_a <- policy
  data <- fixture$snapshots$peer_a[["data_peer_a"]]$data
  data$time_peer_a <- rep(c(2, 8), each = 50L)
  data$status_peer_a <- rep(c("censor", "event"), each = 50L)
  fixture$snapshots$peer_a[["data_peer_a"]]$data <- data
  fixture
}

.synopsis_correlation_real_e2e_fixture <- function(k, server_ns) {
  fixture <- .synopsis_describe_real_e2e_fixture(k, server_ns)
  scope <- list(
    mode = "catalog_v1", numeric_moments = c("x_peer_a", "y_peer_a"),
    categorical_marginals = character(), categorical_pairs = list(),
    correlations = list(c("x_peer_a", "y_peer_a")))
  for (peer in fixture$peers) {
    policy <- fixture$policies[[peer]]
    policy$capsule_workload_scope <- scope
    if (identical(peer, "peer_a")) {
      policy$numeric_bounds$y_peer_a <- c(0, 10)
      policy$capsule_dataset_mapping[["data_peer_a"]] <- c(
        "x_peer_a", "y_peer_a")
    }
    fixture$policies[[peer]] <- policy
  }
  data <- fixture$snapshots$peer_a[["data_peer_a"]]$data
  data$y_peer_a <- rep(c(10, 0), 50L)
  fixture$snapshots$peer_a[["data_peer_a"]]$data <- data
  fixture
}

.synopsis_gaussian_real_e2e_fixture <- function(k, server_ns, n = 10000L) {
  fixture <- .synopsis_correlation_real_e2e_fixture(k, server_ns)
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 64L ||
      n != as.integer(n)) {
    stop("n must be one integer of at least 64", call. = FALSE)
  }
  n <- as.integer(n)
  for (peer in fixture$peers) {
    policy <- fixture$policies[[peer]]
    policy$unit_capacity <- n
    fixture$policies[[peer]] <- policy
  }
  policy <- fixture$policies$peer_a
  policy$capsule_workload_specs$gaussian$gaussian_primary <- list(
    version = "v1", dataset = "data_peer_a", outcome = "y_peer_a",
    predictors = "x_peer_a", intercept = TRUE)
  fixture$policies$peer_a <- policy
  fixture$snapshots$peer_a[["data_peer_a"]]$data <- data.frame(
    patient_id = paste0("u", seq_len(n)),
    x_peer_a = rep(c(0, 10), length.out = n),
    y_peer_a = rep(c(10, 0), length.out = n),
    stringsAsFactors = FALSE)
  fixture
}

.synopsis_lmm_real_e2e_fixture <- function(k, server_ns, n = 10000L) {
  fixture <- .synopsis_describe_real_e2e_fixture(k, server_ns)
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 400L ||
      n != as.integer(n) || n %% 100L != 0L) {
    stop("n must be an integer multiple of 100 of at least 400", call. = FALSE)
  }
  n <- as.integer(n)
  sites <- sprintf("site_%03d", seq_len(n %/% 100L))
  for (peer in fixture$peers) {
    policy <- fixture$policies[[peer]]
    policy$unit_capacity <- n
    policy$capsule_workload_scope <- list(
      mode = "catalog_v1", numeric_moments = character(),
      categorical_marginals = character(), categorical_pairs = list(),
      correlations = list())
    policy$capsule_workload_specs$describe <- list()
    policy$capsule_workload_specs$survival <- list()
    if (identical(peer, "peer_a")) {
      policy$numeric_bounds$y_peer_a <- c(0, 10)
      policy$categorical_levels$site_peer_a <- sites
      policy$capsule_dataset_mapping[["data_peer_a"]] <- c(
        "x_peer_a", "y_peer_a", "site_peer_a")
      policy$capsule_workload_specs$gaussian$lmm_primary <- list(
        version = "random_intercept_v1", dataset = "data_peer_a",
        outcome = "y_peer_a", cluster = "site_peer_a",
        max_patients_per_cluster = 100L)
    }
    fixture$policies[[peer]] <- policy
  }
  site <- rep(sites, each = 100L)
  effect <- rep(seq(1.5, 8.5, length.out = length(sites)), each = 100L)
  within <- rep(c(-0.30, -0.10, 0.10, 0.30), length.out = n)
  fixture$snapshots$peer_a[["data_peer_a"]]$data <- data.frame(
    patient_id = paste0("u", seq_len(n)),
    x_peer_a = rep(c(0, 10), length.out = n),
    y_peer_a = pmin(10, pmax(0, effect + within)),
    site_peer_a = site, stringsAsFactors = FALSE)
  fixture
}

.synopsis_glmm_real_e2e_fixture <- function(k, server_ns, n = 10000L) {
  fixture <- .synopsis_lmm_real_e2e_fixture(k, server_ns, n = n)
  sites <- fixture$policies$peer_a$categorical_levels$site_peer_a
  within_site <- rep(seq_len(100L), length(sites))
  probability <- rep(seq(0.20, 0.80, length.out = length(sites)), each = 100L)
  fixture$policies$peer_a$numeric_bounds$y_peer_a <- c(0, 1)
  fixture$snapshots$peer_a[["data_peer_a"]]$data$y_peer_a <- as.numeric(
    within_site <= round(100 * probability))
  fixture
}

.synopsis_describe_real_e2e_dispatch <- function(fixture) {
  get_server <- function(name) get(name, envir = fixture$server_ns,
                                   inherits = FALSE)
  call_mpc <- get_server(".callMpcTool")
  function(conns, expr, error = NULL, ...) {
    calls <- if (is.list(expr) && !is.call(expr)) expr else {
      stats::setNames(rep(list(expr), length(conns)), names(conns))
    }
    stats::setNames(lapply(names(conns), function(peer) {
      call <- calls[[peer]]
      method <- as.character(call[[1L]])
      args <- as.list(call)[-1L]
      if (identical(method, "dsvertIdentityPkDS")) {
        return(list(identity_pk = unname(fixture$pins[[peer]])))
      }
      if (identical(method, "dsvertDPSynopsisSourcePrepareDS")) {
        fixture$state$source_prepare <- fixture$state$source_prepare + 1L
      }
      if (identical(method, "dsvertDPSynopsisStartDS")) {
        fixture$state$start <- fixture$state$start + 1L
      }
      value <- tryCatch(testthat::with_mocked_bindings(
        do.call(get_server(method), args),
        .dsvert_dp_synopsis_policy_v1 = function() fixture$policies[[peer]],
        .dsvert_dp_policy = function() fixture$policies[[peer]],
        .dsvert_dp_synopsis_state_path_v1 = function() {
          fixture$policies[[peer]]$synopsis_state_path
        },
        .dsvert_session_storage_root = function() {
          root <- file.path(fixture$root, peer, "rock")
          if (!dir.exists(root)) {
            dir.create(root, recursive = TRUE, mode = "0700")
          }
          Sys.chmod(root, mode = "0700")
          normalizePath(root, winslash = "/", mustWork = TRUE)
        },
        .dsvert_require_configured_local_peer_name = function() peer,
        .get_trusted_peers = function() {
          designated <- fixture$peers[1:2]
          fixture$pins[setdiff(designated, peer)]
        },
        .exact_gc_designated_policy_context = function() {
          designated <- sort(fixture$peers[1:2], method = "radix")
          list(
            peer_name = peer, designated = designated,
            pins = fixture$pins[designated],
            full_pinset_sha256 =
              fixture$policies[[peer]]$peer_pinset_sha256,
            consortium_id = "synopsis-real-e2e-exact-gc")
        },
        .dsvert_dp_secret = function() fixture$secrets[[peer]],
        .get_identity_keypair = function() fixture$identities[[peer]],
        .session_storage = function() fixture$state$storage[[peer]],
        .S = function(id) .synopsis_describe_real_e2e_session(
          fixture, peer, id),
        .dsvert_dp_resolve_snapshot = function(policy, data_name, ...) {
          fixture$snapshots[[peer]][[data_name]]
        },
        .callMpcTool = function(command, input) call_mpc(command, input),
        .package = "dsVert"), error = function(condition) {
          if (is.function(error)) error(peer, conditionMessage(condition))
          NULL
        })
      value
    }), names(conns))
  }
}

test_that("real Synopsis Describe is plausible and Rock-replayable at K=2/3/5", {
  .synopsis_real_e2e_only("describe")
  server_ns <- .synopsis_describe_real_e2e_server()
  describe <- get(".dsvert_dp_describe_impl", asNamespace("dsVertClient"),
                  inherits = FALSE)
  meanvar <- get(".dsvert_dp_meanvar_impl", asNamespace("dsVertClient"),
                 inherits = FALSE)
  status_impl <- get(".dsvert_dp_status_impl", asNamespace("dsVertClient"),
                     inherits = FALSE)
  plan_impl <- get(".dsvert_dp_synopsis_plan_impl",
                   asNamespace("dsVertClient"), inherits = FALSE)
  quantile <- get("ds.vertDPQuantile", asNamespace("dsVertClient"),
                  inherits = FALSE)
  median <- get("ds.vertDPMedian", asNamespace("dsVertClient"),
                inherits = FALSE)
  for (k in .synopsis_real_e2e_peer_counts()) {
    fixture <- .synopsis_describe_real_e2e_fixture(k, server_ns)
    on.exit(unlink(fixture$root, recursive = TRUE, force = TRUE), add = TRUE)
    conns <- stats::setNames(lapply(fixture$peers, function(peer) {
      structure(list(peer = peer), class = "dsvert_synopsis_real_e2e_connection")
    }), fixture$peers)
    dispatch <- .synopsis_describe_real_e2e_dispatch(fixture)
    status <- status_impl(conns, dispatch)
    plan <- plan_impl(conns, status, dispatch)
    expect_s3_class(status, "ds.vertDPStatus")
    expect_s3_class(plan, "ds.vertDPCapsulePlan")
    expect_identical(plan$consortium$K, as.integer(k))
    expect_identical(plan$guarantees$data_access, FALSE)
    expect_identical(plan$guarantees$release_created, FALSE)
    expect_identical(plan$guarantees$operation_limit, FALSE)
    expect_identical(plan$guarantees$request_limit, FALSE)
    expect_identical(plan$guarantees$rate_limit, FALSE)
    expect_identical(plan$guarantees$catalog_limit, FALSE)
    expect_identical(fixture$state$source_prepare, 0L)
    expect_identical(fixture$state$start, 0L)
    first <- describe("data_peer_a", "primary", 0.5, "peer_a", conns,
                      dispatch)
    expect_s3_class(first, "ds.vertDPDescribe")
    expect_true(isTRUE(first$released))
    expect_identical(fixture$state$source_prepare, 1L)
    expect_identical(fixture$state$start, 2L)
    expect_identical(first$release_provenance$designated_noise_peers,
                     as.list(fixture$peers[1:2]))
    expect_length(first$release_provenance$ordered_peer_pinset, k)
    expect_true(is.finite(first$descriptives$n_dp[[1L]]))
    expect_gt(first$descriptives$n_dp[[1L]], 0)
    expect_true(is.finite(first$descriptives$mean[[1L]]))
    expect_gte(first$descriptives$mean[[1L]], 0)
    expect_lte(first$descriptives$mean[[1L]], 10)
    expect_true(is.finite(first$descriptives$variance[[1L]]))
    expect_gte(first$descriptives$variance[[1L]], 0)
    expect_lte(first$descriptives$variance[[1L]], 25)
    source_values <- fixture$snapshots$peer_a$data_peer_a$data$x_peer_a
    source_mean <- mean(source_values)
    source_variance <- mean((source_values - source_mean)^2)
    describe_row <- first$descriptives[1L, , drop = FALSE]
    expect_lte(describe_row$n_simultaneous_95_lower, length(source_values))
    expect_gte(
      describe_row$n_dp + describe_row$count_noise_radius_simultaneous_95,
      length(source_values))
    expect_lte(describe_row$mean_mechanism_grid_lower_95, source_mean)
    expect_gte(describe_row$mean_mechanism_grid_upper_95, source_mean)
    expect_lte(describe_row$variance_mechanism_grid_lower_95, source_variance)
    expect_gte(describe_row$variance_mechanism_grid_upper_95, source_variance)

    route_describe <- function(data_name, analysis_id, probs,
                               server = NULL, datasources = NULL,
                               .aggregate, resume = NULL) {
      describe(data_name, analysis_id, probs, server, datasources, dispatch,
               resume = resume)
    }
    public_describe <- testthat::with_mocked_bindings(
      .dsvert_dp_describe_impl = route_describe,
      ds.vertDPDescribe(
        "data_peer_a", "primary", probs = 0.5, server = "peer_a",
        datasources = conns),
      .package = "dsVertClient")
    expect_s3_class(public_describe, "ds.vertDPDescribe")
    expect_identical(public_describe$statistics, first$statistics)
    expect_identical(public_describe$final_vector_root, first$final_vector_root)
    expect_identical(c(fixture$state$source_prepare, fixture$state$start),
                     c(1L, 2L))

    route_meanvar <- function(data_name, variable, server = NULL,
                              datasources = NULL, .aggregate) {
      meanvar(data_name, variable, server, datasources, dispatch)
    }
    public_meanvar <- testthat::with_mocked_bindings(
      .dsvert_dp_meanvar_impl = route_meanvar,
      ds.vertDPMeanVar("data_peer_a", "x_peer_a", "peer_a", conns),
      .package = "dsVertClient")
    expect_s3_class(public_meanvar, "ds.vertDPMeanVar")
    expect_true(is.finite(public_meanvar$mean))
    expect_true(is.finite(public_meanvar$variance))
    expect_identical(c(fixture$state$source_prepare, fixture$state$start),
                     c(1L, 2L))

    postprocess_before <- c(fixture$state$source_prepare, fixture$state$start)
    quantiles <- quantile(first, probs = c(0.25, 0.5, 0.75))
    medians <- median(first)
    expect_s3_class(quantiles, "ds.vertDPQuantile")
    expect_s3_class(medians, "ds.vertDPMedian")
    expect_true(all(quantiles$status == "ok_binned_postprocessed_estimate"))
    expect_true(all(is.finite(quantiles$bin_lower)))
    expect_true(all(is.finite(quantiles$bin_upper)))
    expect_true(all(quantiles$bin_lower >= 0 & quantiles$bin_upper <= 10))
    expect_identical(attr(quantiles, "additional_server_calls"), 0L)
    expect_identical(attr(medians, "additional_server_calls"), 0L)
    expect_identical(c(fixture$state$source_prepare, fixture$state$start),
                     postprocess_before)

    legacy_describe <- function(data_name, analysis_id, probs, server = NULL,
                                datasources = NULL, resume = NULL) {
      describe(data_name, analysis_id, probs, server, datasources, dispatch,
               resume = resume)
    }
    legacy <- testthat::with_mocked_bindings(
      ds.vertDPDescribe = legacy_describe,
      ds.vertDesc("data_peer_a", analysis_id = "primary",
                  probs = c(0.25, 0.5, 0.75), verbose = FALSE,
                  datasources = conns),
      .package = "dsVertClient")
    legacy_alias <- testthat::with_mocked_bindings(
      ds.vertDPDescribe = legacy_describe,
      ds.vert.desc("data_peer_a", analysis_id = "primary",
                   probs = c(0.25, 0.5, 0.75), verbose = FALSE,
                   datasources = conns),
      .package = "dsVertClient")
    expect_s3_class(legacy, "ds.vertDesc")
    expect_s3_class(legacy_alias, "ds.vertDesc")
    expect_identical(legacy, legacy_alias)
    expect_true(all(is.finite(legacy$n)))
    expect_true(all(legacy$mean >= legacy$range_low &
                    legacy$mean <= legacy$range_high))
    expect_true(all(c("q25", "q50", "q75") %in% names(legacy)))
    expect_identical(c(fixture$state$source_prepare, fixture$state$start),
                     postprocess_before)

    before <- c(fixture$state$source_prepare, fixture$state$start)
    fixture$state$storage <- stats::setNames(lapply(fixture$peers, function(...) {
      new.env(parent = emptyenv())
    }), fixture$peers)
    moments <- meanvar("data_peer_a", "x_peer_a", "peer_a", conns, dispatch)
    expect_s3_class(moments, "ds.vertDPMeanVar")
    expect_true(is.finite(moments$n))
    expect_gt(moments$n, 0)
    expect_true(is.finite(moments$mean))
    expect_gte(moments$mean, 0)
    expect_lte(moments$mean, 10)
    expect_true(is.finite(moments$variance))
    expect_gte(moments$variance, 0)
    expect_lte(moments$variance, 25)
    expect_identical(c(fixture$state$source_prepare, fixture$state$start), before)

    fixture$state$storage <- stats::setNames(lapply(fixture$peers, function(...) {
      new.env(parent = emptyenv())
    }), fixture$peers)
    replay <- describe("data_peer_a", "primary", 0.5, "peer_a", conns,
                       dispatch, resume = first$resume)
    expect_identical(replay$statistics, first$statistics)
    expect_identical(replay$final_vector_root, first$final_vector_root)
    expect_identical(c(fixture$state$source_prepare, fixture$state$start), before)
  }
})

test_that("real same-owner Synopsis contingency is plausible and Rock-replayable at K=2/3/5", {
  .synopsis_real_e2e_only("same_owner_contingency")
  server_ns <- .synopsis_describe_real_e2e_server()
  contingency <- get(".dsvert_dp_contingency_impl",
                     asNamespace("dsVertClient"), inherits = FALSE)
  for (k in .synopsis_real_e2e_peer_counts()) {
    fixture <- .synopsis_contingency_real_e2e_fixture(k, server_ns)
    on.exit(unlink(fixture$root, recursive = TRUE, force = TRUE), add = TRUE)
    conns <- stats::setNames(lapply(fixture$peers, function(peer) {
      structure(list(peer = peer), class = "dsvert_synopsis_real_e2e_connection")
    }), fixture$peers)
    dispatch <- .synopsis_describe_real_e2e_dispatch(fixture)
    first <- contingency("data_peer_a", "exposure", "outcome", "peer_a",
                         conns, dispatch)
    expect_s3_class(first, "ds.vertDPContingency")
    expect_true(isTRUE(first$released))
    expect_false(isTRUE(first$cross_owner))
    expect_identical(fixture$state$source_prepare, 1L)
    expect_identical(fixture$state$start, 2L)
    expect_identical(first$release_provenance$designated_noise_peers,
                     as.list(fixture$peers[1:2]))
    expect_length(first$release_provenance$ordered_peer_pinset, k)
    expect_identical(dim(first$table), c(2L, 2L))
    expect_true(all(is.finite(first$table)))
    expect_true(all(first$table >= 0 & first$table <= 200))
    expect_gt(sum(first$table), 0)

    route_contingency <- function(data_name, row_var, col_var, server = NULL,
                                  datasources = NULL, .aggregate) {
      contingency(data_name, row_var, col_var, server, datasources, dispatch)
    }
    public <- testthat::with_mocked_bindings(
      .dsvert_dp_contingency_impl = route_contingency,
      ds.vertDPContingency(
        "data_peer_a", "exposure", "outcome", server = "peer_a",
        datasources = conns),
      .package = "dsVertClient")
    expect_s3_class(public, "ds.vertDPContingency")
    expect_identical(public$table, first$table)
    expect_identical(public$final_vector_root, first$final_vector_root)
    expect_identical(c(fixture$state$source_prepare, fixture$state$start),
                     c(1L, 2L))

    before <- c(fixture$state$source_prepare, fixture$state$start)
    chisq <- ds.vertChisq(first, simulations = 128L)
    fisher <- ds.vertFisher(first, simulations = 128L)
    chisq_alias <- ds.vert.chisq(first, simulations = 128L)
    fisher_alias <- ds.vert.fisher(first, simulations = 128L)
    expect_s3_class(chisq, "ds.vertChisq")
    expect_s3_class(fisher, "ds.vertFisher")
    expect_identical(chisq_alias$frontdoor, "ds.vert.chisq")
    expect_identical(chisq_alias$route, "ds.vertChisq")
    expect_identical(fisher_alias$frontdoor, "ds.vert.fisher")
    expect_identical(fisher_alias$route, "ds.vertFisher")
    expect_identical(chisq_alias$statistic, chisq$statistic)
    expect_identical(chisq_alias$p_value, chisq$p_value)
    expect_identical(fisher_alias$statistic, fisher$statistic)
    expect_identical(fisher_alias$p_value, fisher$p_value)
    expect_identical(chisq$status, "ok")
    expect_identical(fisher$status, "ok")
    expect_true(is.finite(chisq$statistic) && chisq$statistic >= 0)
    expect_true(is.finite(fisher$statistic) && fisher$statistic >= 0)
    expect_true(chisq$p_value >= 0 && chisq$p_value <= 1)
    expect_true(fisher$p_value >= 0 && fisher$p_value <= 1)
    for (inference in list(chisq, fisher)) {
      expect_identical(inference$additional_server_calls, 0L)
      expect_identical(inference$additional_privacy_cost,
                       c(epsilon = 0, delta = 0))
    }

    # These historical epidemiology views must all be deterministic
    # post-processing of this one signed 2-by-2 release: no second query, no
    # fresh noise, and no hidden privacy cost. The balanced fixture makes the
    # resulting finite-snapshot quantities plausibly estimable at K=2, K=3,
    # and K=5 without treating the DP result as a sampling claim.
    epi <- ds.vertDPEpi2x2(first, exposed = "yes", event = "yes")
    epi_inference <- ds.vertDPEpi2x2Inference(
      first, exposed = "yes", event = "yes")
    # The historical frontdoor must recognise a validated DP contingency
    # result and delegate to the same zero-call, mechanism-aware view rather
    # than applying an ordinary Wald calculation to noisy cells.
    epi_historical <- ds.vertEpi2x2(
      chisq, exposed = "yes", event = "yes")
    prevalence <- ds.vertDPPrevalenceRatio(
      first, exposed = "yes", prevalent = "yes")
    prevalence_inference <- ds.vertDPPrevalenceRatioInference(
      first, exposed = "yes", prevalent = "yes")
    diagnostic <- ds.vertDPDiagnostic2x2(
      first, disease_positive = "yes", test_positive = "yes")
    diagnostic_inference <- ds.vertDPDiagnostic2x2Inference(
      first, disease_positive = "yes", test_positive = "yes")
    direct <- ds.vertDPDirectStandardization(
      first, standard_weights = c(no = 0.4, yes = 0.6), event = "yes")
    direct_inference <- ds.vertDPDirectStandardizationInference(
      first, standard_weights = c(no = 0.4, yes = 0.6), event = "yes")
    indirect <- ds.vertDPIndirectStandardization(
      first, expected_rates = c(no = 0.2, yes = 0.3), event = "yes")
    indirect_inference <- ds.vertDPIndirectStandardizationInference(
      first, expected_rates = c(no = 0.2, yes = 0.3), event = "yes")

    for (value in list(
      epi, epi_inference, epi_historical, prevalence, prevalence_inference,
      diagnostic, diagnostic_inference, direct, direct_inference,
      indirect, indirect_inference)) {
      expect_identical(value$additional_server_calls, 0L)
      expect_identical(value$additional_privacy_cost,
                       c(epsilon = 0, delta = 0))
    }
    expect_s3_class(epi_historical, "ds.vertDPEpi2x2")
    expect_identical(epi_historical$point_estimates, epi$point_estimates)
    expect_identical(epi_historical$mechanism_regions,
                     epi$mechanism_regions)
    epi_risks <- unlist(epi$point_estimates[c(
      "risk_exposed", "risk_unexposed", "population_risk")],
      use.names = FALSE)
    expect_true(all(is.finite(epi_risks)))
    expect_true(all(epi_risks >= 0 & epi_risks <= 1))
    expect_identical(prevalence$prevalence_point_estimates[[
      "prevalence_ratio"]], epi$point_estimates[["risk_ratio"]])
    expect_true(all(is.finite(diagnostic$estimates[c(
      "sensitivity", "specificity", "accuracy", "balanced_accuracy")])))
    expect_true(all(diagnostic$estimates[c(
      "sensitivity", "specificity", "accuracy", "balanced_accuracy")] >= 0 &
      diagnostic$estimates[c(
        "sensitivity", "specificity", "accuracy", "balanced_accuracy")] <= 1))
    for (result in list(direct, direct_inference)) {
      region <- if (inherits(result,
                             "ds.vertDPDirectStandardizationInference")) {
        result$combined_region
      } else {
        result$mechanism_region
      }
      expect_true(is.finite(result$estimate))
      expect_gte(result$estimate, region[["lower"]] - 1e-12)
      expect_lte(result$estimate, region[["upper"]] + 1e-12)
      expect_gte(region[["lower"]], 0)
      expect_lte(region[["upper"]], 1)
    }
    for (result in list(indirect, indirect_inference)) {
      region <- if (inherits(result,
                             "ds.vertDPIndirectStandardizationInference")) {
        result$combined_region
      } else {
        result$mechanism_region
      }
      expect_true(is.finite(result$estimate))
      expect_gte(result$estimate, region[["lower"]] - 1e-12)
      expect_lte(result$estimate, region[["upper"]] + 1e-12)
      expect_gte(region[["lower"]], 0)
      expect_true(is.finite(region[["upper"]]))
    }
    expect_identical(c(fixture$state$source_prepare, fixture$state$start), before)

    fixture$state$storage <- stats::setNames(lapply(fixture$peers, function(...) {
      new.env(parent = emptyenv())
    }), fixture$peers)
    replay <- contingency("data_peer_a", "exposure", "outcome", "peer_a",
                          conns, dispatch)
    expect_identical(replay$table, first$table)
    expect_identical(replay$final_vector_root, first$final_vector_root)
    expect_identical(c(fixture$state$source_prepare, fixture$state$start), before)

    # The historical IPW name has a deliberately narrow, exact special case:
    # treatment ~ 1 reduces normalized IPW to the arm risk difference.  It
    # consumes this already-released signed table; it does not create a
    # propensity fit, individual weights, or another private draw.
    ipw <- testthat::with_mocked_bindings(
      ds.vertDPContingency = function(data, row_var, col_var, server = NULL,
                                      datasources = NULL) {
        expect_identical(c(data, row_var, col_var, server), c(
          "data_peer_a", "exposure", "outcome", "peer_a"))
        replay
      },
      ds.vertIPW(
        outcome ~ exposure, exposure ~ 1, data = "data_peer_a",
        treated = "yes", event = "yes", server = "peer_a",
        datasources = conns, verbose = FALSE),
      .package = "dsVertClient")
    expect_s3_class(ipw, "ds.vertIPW")
    expect_true(is.finite(ipw$estimate))
    expect_gte(ipw$estimate, -1)
    expect_lte(ipw$estimate, 1)
    expect_true(all(is.finite(ipw$confidence_region)))
    expect_true(all(ipw$confidence_region >= -1 & ipw$confidence_region <= 1))
    expect_identical(ipw$source_artifact_key, replay$artifact_key)
    expect_identical(ipw$final_vector_root, replay$final_vector_root)
    expect_identical(ipw$weights_released, FALSE)
    expect_identical(ipw$additional_server_calls_after_artifact, 0L)
    expect_identical(c(fixture$state$source_prepare, fixture$state$start), before)
  }
})

test_that("real stratified Synopsis supports sticky standardisation at K=2/3/5", {
  .synopsis_real_e2e_only("stratified_epidemiology")
  server_ns <- .synopsis_describe_real_e2e_server()
  contingency <- get(".dsvert_dp_contingency_impl",
                     asNamespace("dsVertClient"), inherits = FALSE)
  for (k in .synopsis_real_e2e_peer_counts()) {
    fixture <- .synopsis_stratified_contingency_real_e2e_fixture(k, server_ns)
    on.exit(unlink(fixture$root, recursive = TRUE, force = TRUE), add = TRUE)
    conns <- stats::setNames(lapply(fixture$peers, function(peer) {
      structure(list(peer = peer), class = "dsvert_synopsis_real_e2e_connection")
    }), fixture$peers)
    dispatch <- .synopsis_describe_real_e2e_dispatch(fixture)
    first <- contingency("data_peer_a", "stratum", "outcome", "peer_a",
                         conns, dispatch)
    expect_s3_class(first, "ds.vertDPContingency")
    expect_true(isTRUE(first$released))
    expect_identical(dim(first$table), c(3L, 2L))
    expect_setequal(rownames(first$table), c("young", "middle", "old"))
    expect_identical(colnames(first$table), c("no", "yes"))
    expect_true(all(is.finite(first$table) & first$table >= 0))
    expect_identical(fixture$state$source_prepare, 1L)
    expect_identical(fixture$state$start, 2L)

    direct <- ds.vertDPDirectStandardization(
      first, c(young = 0.2, middle = 0.3, old = 0.5), event = "yes")
    direct_inference <- ds.vertDPDirectStandardizationInference(
      first, c(young = 0.2, middle = 0.3, old = 0.5), event = "yes")
    indirect <- ds.vertDPIndirectStandardization(
      first, c(young = 0.1, middle = 0.2, old = 0.4), event = "yes")
    indirect_inference <- ds.vertDPIndirectStandardizationInference(
      first, c(young = 0.1, middle = 0.2, old = 0.4), event = "yes")
    for (result in list(direct, direct_inference, indirect,
                        indirect_inference)) {
      expect_identical(result$additional_server_calls, 0L)
      expect_identical(result$additional_privacy_cost,
                       c(epsilon = 0, delta = 0))
    }
    expect_equal(
      direct$weights,
      unname(c(young = 0.2, middle = 0.3, old = 0.5)[rownames(first$table)]),
      tolerance = 0)
    expect_equal(
      indirect$expected_rates,
      unname(c(young = 0.1, middle = 0.2, old = 0.4)[rownames(first$table)]),
      tolerance = 0)
    for (result in list(direct, direct_inference)) {
      region <- if (inherits(result,
                             "ds.vertDPDirectStandardizationInference")) {
        result$combined_region
      } else {
        result$mechanism_region
      }
      expect_true(is.finite(result$estimate))
      expect_gte(result$estimate, region[["lower"]] - 1e-12)
      expect_lte(result$estimate, region[["upper"]] + 1e-12)
      expect_gte(region[["lower"]], 0)
      expect_lte(region[["upper"]], 1)
      expect_gt(result$estimate, 0.1)
      expect_lt(result$estimate, 0.75)
    }
    for (result in list(indirect, indirect_inference)) {
      region <- if (inherits(result,
                             "ds.vertDPIndirectStandardizationInference")) {
        result$combined_region
      } else {
        result$mechanism_region
      }
      expect_true(is.finite(result$estimate))
      expect_gte(result$estimate, region[["lower"]] - 1e-12)
      expect_lte(result$estimate, region[["upper"]] + 1e-12)
      expect_gte(region[["lower"]], 0)
      expect_true(is.finite(region[["upper"]]))
      expect_gt(result$estimate, 0.5)
      expect_lt(result$estimate, 3)
    }

    before <- c(fixture$state$source_prepare, fixture$state$start)
    fixture$state$storage <- stats::setNames(lapply(fixture$peers, function(...) {
      new.env(parent = emptyenv())
    }), fixture$peers)
    replay <- contingency("data_peer_a", "stratum", "outcome", "peer_a",
                          conns, dispatch)
    expect_identical(replay$table, first$table)
    expect_identical(replay$final_vector_root, first$final_vector_root)
    expect_identical(c(fixture$state$source_prepare, fixture$state$start), before)
  }
})

test_that("real Synopsis Mantel-Haenszel is plausible and Rock-replayable at K=2/3/5", {
  .synopsis_real_e2e_only("mantel_haenszel")
  server_ns <- .synopsis_describe_real_e2e_server()
  contingency <- get(".dsvert_dp_contingency_impl",
                     asNamespace("dsVertClient"), inherits = FALSE)
  cells <- c(
    "exposed_event", "exposed_nonevent",
    "unexposed_event", "unexposed_nonevent")
  for (k in .synopsis_real_e2e_peer_counts()) {
    fixture <- .synopsis_mantel_haenszel_real_e2e_fixture(k, server_ns)
    on.exit(unlink(fixture$root, recursive = TRUE, force = TRUE), add = TRUE)
    conns <- stats::setNames(lapply(fixture$peers, function(peer) {
      structure(list(peer = peer), class = "dsvert_synopsis_real_e2e_connection")
    }), fixture$peers)
    dispatch <- .synopsis_describe_real_e2e_dispatch(fixture)
    first <- contingency("data_peer_a", "stratum", "cell", "peer_a",
                         conns, dispatch)
    expect_s3_class(first, "ds.vertDPContingency")
    expect_true(isTRUE(first$released))
    expect_identical(dim(first$table), c(3L, 4L))
    expect_identical(rownames(first$table),
                     sort(c("young", "middle", "old")))
    expect_identical(colnames(first$table), cells)
    expect_identical(first$unit_aggregation_policy,
                     "consistent_joint_cell_else_exclude_v1")
    expect_equal(first$artifact_l1_sensitivity, 1, tolerance = 0)
    expect_true(all(is.finite(first$table) & first$table >= 0))

    before <- c(fixture$state$source_prepare, fixture$state$start)
    fit <- ds.vertDPMantelHaenszel(first)
    expect_s3_class(fit, "ds.vertDPMantelHaenszel")
    expect_identical(fit$status, "ok")
    expect_identical(fit$estimate_type, "finite")
    expect_true(is.finite(fit$estimate) && fit$estimate > 0)
    expect_true(is.finite(fit$mechanism_region[["lower"]]))
    expect_false(is.na(fit$mechanism_region[["upper"]]))
    expect_lte(fit$mechanism_region[["lower"]], fit$estimate)
    expect_gte(fit$mechanism_region[["upper"]], fit$estimate)
    expect_identical(fit$additional_server_calls, 0L)
    expect_identical(fit$additional_privacy_cost,
                     c(epsilon = 0, delta = 0))
    expect_identical(c(fixture$state$source_prepare, fixture$state$start), before)

    fixture$state$storage <- stats::setNames(lapply(fixture$peers, function(...) {
      new.env(parent = emptyenv())
    }), fixture$peers)
    replay <- contingency("data_peer_a", "stratum", "cell", "peer_a",
                          conns, dispatch)
    expect_identical(replay$table, first$table)
    expect_identical(replay$final_vector_root, first$final_vector_root)
    expect_identical(c(fixture$state$source_prepare, fixture$state$start), before)
  }
})

test_that("real Synopsis ROC is plausible and Rock-replayable at K=2/3/5", {
  .synopsis_real_e2e_only("roc")
  server_ns <- .synopsis_describe_real_e2e_server()
  contingency <- get(".dsvert_dp_contingency_impl",
                     asNamespace("dsVertClient"), inherits = FALSE)
  bins <- c("low", "mid_low", "mid_high", "high")
  for (k in .synopsis_real_e2e_peer_counts()) {
    fixture <- .synopsis_roc_real_e2e_fixture(k, server_ns)
    on.exit(unlink(fixture$root, recursive = TRUE, force = TRUE), add = TRUE)
    conns <- stats::setNames(lapply(fixture$peers, function(peer) {
      structure(list(peer = peer), class = "dsvert_synopsis_real_e2e_connection")
    }), fixture$peers)
    dispatch <- .synopsis_describe_real_e2e_dispatch(fixture)
    first <- contingency("data_peer_a", "disease", "score_bin", "peer_a",
                         conns, dispatch)
    expect_s3_class(first, "ds.vertDPContingency")
    expect_true(isTRUE(first$released))
    expect_identical(dim(first$table), c(2L, 4L))
    expect_setequal(rownames(first$table), c("case", "control"))
    expect_setequal(colnames(first$table), bins)
    expect_true(all(is.finite(first$table) & first$table >= 0))

    before <- c(fixture$state$source_prepare, fixture$state$start)
    fit <- ds.vertDPROC(first, disease_positive = "case", score_order = bins,
                        direction = "higher")
    expect_s3_class(fit, "ds.vertDPROC")
    expect_identical(fit$status, "ok")
    expect_true(is.finite(fit$auc) && fit$auc >= 0 && fit$auc <= 1)
    expect_equal(nrow(fit$curve), length(bins) + 1L)
    expect_true(all(is.finite(as.matrix(fit$curve[, c(
      "sensitivity", "specificity", "false_positive_rate")]))))
    expect_true(all(fit$curve$sensitivity >= 0 & fit$curve$sensitivity <= 1))
    expect_true(all(fit$curve$specificity >= 0 & fit$curve$specificity <= 1))
    expect_lte(fit$auc_mechanism_region[["lower"]], fit$auc)
    expect_gte(fit$auc_mechanism_region[["upper"]], fit$auc)
    expect_identical(fit$additional_server_calls, 0L)
    expect_identical(fit$additional_privacy_cost,
                     c(epsilon = 0, delta = 0))
    expect_identical(c(fixture$state$source_prepare, fixture$state$start), before)

    fixture$state$storage <- stats::setNames(lapply(fixture$peers, function(...) {
      new.env(parent = emptyenv())
    }), fixture$peers)
    replay <- contingency("data_peer_a", "disease", "score_bin", "peer_a",
                          conns, dispatch)
    expect_identical(replay$table, first$table)
    expect_identical(replay$final_vector_root, first$final_vector_root)
    expect_identical(c(fixture$state$source_prepare, fixture$state$start), before)
  }
})

test_that("real Synopsis causal standardisation is plausible and Rock-replayable at K=2/3/5", {
  .synopsis_real_e2e_only("causal_standardization")
  server_ns <- .synopsis_describe_real_e2e_server()
  contingency <- get(".dsvert_dp_contingency_impl",
                     asNamespace("dsVertClient"), inherits = FALSE)
  strata <- c("young", "young", "middle", "middle", "old", "old")
  treatment <- rep(c("control", "treated"), 3L)
  weights <- c(young = 0.2, middle = 0.3, old = 0.5)
  for (k in .synopsis_real_e2e_peer_counts()) {
    fixture <- .synopsis_causal_contingency_real_e2e_fixture(k, server_ns)
    on.exit(unlink(fixture$root, recursive = TRUE, force = TRUE), add = TRUE)
    conns <- stats::setNames(lapply(fixture$peers, function(peer) {
      structure(list(peer = peer), class = "dsvert_synopsis_real_e2e_connection")
    }), fixture$peers)
    dispatch <- .synopsis_describe_real_e2e_dispatch(fixture)
    first <- contingency("data_peer_a", "arm", "outcome", "peer_a",
                         conns, dispatch)
    expect_s3_class(first, "ds.vertDPContingency")
    expect_true(isTRUE(first$released))
    expect_identical(dim(first$table), c(6L, 2L))
    expect_identical(rownames(first$table), c(
      "middle_control", "middle_treated", "old_control",
      "old_treated", "young_control", "young_treated"))
    expect_identical(colnames(first$table), c("no", "yes"))
    expect_true(all(is.finite(first$table) & first$table >= 0 &
                    first$table <= 256))
    expect_identical(fixture$state$source_prepare, 1L)
    expect_identical(fixture$state$start, 2L)

    causal <- ds.vertDPCausalStandardization(
      first, strata = strata, treatment = treatment, treated = "treated",
      standard_weights = weights, event = "yes")
    causal_inference <- ds.vertDPCausalStandardizationInference(
      first, strata = strata, treatment = treatment, treated = "treated",
      standard_weights = weights, event = "yes")
    expect_s3_class(causal, "ds.vertDPCausalStandardization")
    expect_s3_class(causal_inference,
                    "ds.vertDPCausalStandardizationInference")
    for (result in list(causal, causal_inference)) {
      expect_identical(result$additional_server_calls, 0L)
      expect_identical(result$additional_privacy_cost,
                       c(epsilon = 0, delta = 0))
      expect_equal(result$standard_weights, weights, tolerance = 0)
    }
    for (name in c("risk_treated", "risk_control")) {
      expect_true(all(is.finite(causal$mechanism_regions[[name]])))
      expect_true(all(causal$mechanism_regions[[name]] >= 0 &
                      causal$mechanism_regions[[name]] <= 1))
      expect_true(all(is.finite(causal_inference$combined_regions[[name]])))
      expect_true(all(causal_inference$combined_regions[[name]] >= 0 &
                      causal_inference$combined_regions[[name]] <= 1))
    }
    expect_identical(causal_inference$coverage_lower_bound, 0.95)

    before <- c(fixture$state$source_prepare, fixture$state$start)
    fixture$state$storage <- stats::setNames(lapply(fixture$peers, function(...) {
      new.env(parent = emptyenv())
    }), fixture$peers)
    replay <- contingency("data_peer_a", "arm", "outcome", "peer_a",
                          conns, dispatch)
    expect_identical(replay$table, first$table)
    expect_identical(replay$final_vector_root, first$final_vector_root)
    expect_identical(c(fixture$state$source_prepare, fixture$state$start), before)
  }
})

test_that("real cross-owner Synopsis contingency is plausible and Rock-replayable at K=2/3/5", {
  .synopsis_real_e2e_only("cross_owner_contingency")
  server_ns <- .synopsis_describe_real_e2e_server()
  contingency <- get(".dsvert_dp_contingency_impl",
                     asNamespace("dsVertClient"), inherits = FALSE)
  for (k in .synopsis_real_e2e_peer_counts()) {
    fixture <- .synopsis_cross_contingency_real_e2e_fixture(k, server_ns)
    on.exit(unlink(fixture$root, recursive = TRUE, force = TRUE), add = TRUE)
    conns <- stats::setNames(lapply(fixture$peers, function(peer) {
      structure(list(peer = peer), class = "dsvert_synopsis_real_e2e_connection")
    }), fixture$peers)
    dispatch <- .synopsis_describe_real_e2e_dispatch(fixture)

    first <- contingency(
      "data_peer_a", "peer_a$disease", "peer_b$exposure", NULL,
      conns, dispatch)
    expect_s3_class(first, "ds.vertDPContingency")
    expect_true(isTRUE(first$released))
    expect_true(isTRUE(first$cross_owner))
    expect_identical(first$servers, c("peer_a", "peer_b"))
    expect_identical(fixture$state$source_prepare, 2L)
    expect_identical(fixture$state$start, 2L)
    expect_identical(first$release_provenance$designated_noise_peers,
                     as.list(fixture$peers[1:2]))
    expect_length(first$release_provenance$ordered_peer_pinset, k)
    expect_identical(dim(first$table), c(2L, 2L))
    expect_true(all(is.finite(first$table)))
    expect_true(all(first$table >= 0 & first$table <= 200))
    expect_gt(sum(first$table), 0)

    route_contingency <- function(data_name, row_var, col_var, server = NULL,
                                  datasources = NULL, .aggregate) {
      contingency(data_name, row_var, col_var, server, datasources, dispatch)
    }
    public <- testthat::with_mocked_bindings(
      .dsvert_dp_contingency_impl = route_contingency,
      ds.vertDPContingency(
        "data_peer_a", "peer_a$disease", "peer_b$exposure",
        datasources = conns),
      .package = "dsVertClient")
    expect_s3_class(public, "ds.vertDPContingency")
    expect_identical(public$table, first$table)
    expect_identical(public$final_vector_root, first$final_vector_root)
    expect_identical(c(fixture$state$source_prepare, fixture$state$start),
                     c(2L, 2L))

    inference <- ds.vertChisqCross(
      first, correct = TRUE, fisher = TRUE, simulations = 128L,
      verbose = FALSE)
    expect_s3_class(inference, "ds.vertChisq")
    expect_true(isTRUE(inference$cross_owner))
    expect_identical(inference$source_dp_release, first)
    if (isTRUE(inference$decision_available)) {
      expect_true(is.finite(inference$p_value) &&
                  inference$p_value >= 0 && inference$p_value <= 1)
    } else {
      expect_true(inference$status %in% c(
        "not_tested_degenerate_dp_table",
        "not_tested_minimum_expected_count"))
      expect_true(is.na(inference$p_value))
    }
    if (is.finite(inference$fisher_p)) {
      expect_true(inference$fisher_p >= 0 && inference$fisher_p <= 1)
    } else {
      expect_false(inference$fisher$decision_available)
      expect_true(inference$fisher$status %in% c(
        "not_tested_degenerate_dp_table",
        "not_tested_minimum_expected_count",
        "not_tested_degenerate_conditional_support"))
      expect_true(is.na(inference$fisher$p_value))
    }

    historical <- ds.vert.chisq_cross(
      first, correct = TRUE, fisher = TRUE, simulations = 128L,
      verbose = FALSE)
    expect_s3_class(historical, "ds.vertChisq")
    expect_identical(historical$frontdoor, "ds.vert.chisq_cross")
    expect_identical(historical$route, "ds.vertChisqCross")
    expect_identical(historical$statistic, inference$statistic)
    expect_identical(historical$p_value, inference$p_value)
    expect_identical(historical$fisher_p, inference$fisher_p)
    expect_identical(historical$source_dp_release, first)
    expect_identical(historical$additional_server_calls, 0L)
    expect_identical(historical$additional_privacy_cost,
                     c(epsilon = 0, delta = 0))
    expect_identical(c(fixture$state$source_prepare, fixture$state$start),
                     c(2L, 2L))

    before <- c(fixture$state$source_prepare, fixture$state$start)
    fixture$state$storage <- stats::setNames(lapply(fixture$peers, function(...) {
      new.env(parent = emptyenv())
    }), fixture$peers)
    replay <- contingency(
      "data_peer_a", "peer_a$disease", "peer_b$exposure", NULL,
      conns, dispatch)
    expect_identical(replay$table, first$table)
    expect_identical(replay$final_vector_root, first$final_vector_root)
    expect_identical(c(fixture$state$source_prepare, fixture$state$start), before)
  }
})

test_that("real Synopsis Frequency is plausible and Rock-replayable at K=2/3/5", {
  .synopsis_real_e2e_only("frequency")
  server_ns <- .synopsis_describe_real_e2e_server()
  frequency <- get(".dsvert_dp_frequency_impl",
                   asNamespace("dsVertClient"), inherits = FALSE)
  for (k in .synopsis_real_e2e_peer_counts()) {
    fixture <- .synopsis_frequency_real_e2e_fixture(k, server_ns)
    on.exit(unlink(fixture$root, recursive = TRUE, force = TRUE), add = TRUE)
    conns <- stats::setNames(lapply(fixture$peers, function(peer) {
      structure(list(peer = peer), class = "dsvert_synopsis_real_e2e_connection")
    }), fixture$peers)
    dispatch <- .synopsis_describe_real_e2e_dispatch(fixture)
    first <- frequency("data_peer_a", "status", "peer_a", conns, dispatch)
    expect_s3_class(first, "ds.vertDPFrequency")
    expect_true(isTRUE(first$released))
    expect_identical(first$proof$version,
                     "dsvert-dp-frequency-synopsis-proof-v1")
    expect_identical(first$levels, c("case", "control"))
    expect_identical(fixture$state$source_prepare, 1L)
    expect_identical(fixture$state$start, 2L)
    expect_identical(first$release_provenance$designated_noise_peers,
                     as.list(fixture$peers[1:2]))
    expect_length(first$release_provenance$ordered_peer_pinset, k)
    expect_true(all(is.finite(first$counts)))
    expect_true(all(first$counts >= 0 & first$counts <= 200))
    expect_gt(sum(first$counts), 0)
    expect_true(all(first$mechanism_regions$lower >= 0))
    expect_true(all(first$mechanism_regions$upper <= 200))
    source_levels <- fixture$snapshots$peer_a$data_peer_a$data$status
    source_counts <- table(factor(source_levels, levels = first$levels))
    expect_lte(first$mechanism_regions$lower, as.numeric(source_counts))
    expect_gte(first$mechanism_regions$upper, as.numeric(source_counts))

    route_frequency <- function(data_name, variable, server = NULL,
                                datasources = NULL, .aggregate) {
      frequency(data_name, variable, server, datasources, dispatch)
    }
    public <- testthat::with_mocked_bindings(
      .dsvert_dp_frequency_impl = route_frequency,
      ds.vertDPFrequency("data_peer_a", "status", "peer_a", conns),
      .package = "dsVertClient")
    expect_s3_class(public, "ds.vertDPFrequency")
    expect_identical(public$counts, first$counts)
    expect_identical(public$final_vector_root, first$final_vector_root)
    expect_identical(c(fixture$state$source_prepare, fixture$state$start),
                     c(1L, 2L))

    inference <- ds.vertDPFrequencyInference(first, level = 0.95)
    expect_s3_class(inference, "ds.vertDPFrequencyInference")
    expect_true(all(is.finite(inference$intervals)))
    expect_true(all(inference$intervals >= 0 & inference$intervals <= 1))

    before <- c(fixture$state$source_prepare, fixture$state$start)
    fixture$state$storage <- stats::setNames(lapply(fixture$peers, function(...) {
      new.env(parent = emptyenv())
    }), fixture$peers)
    replay <- frequency("data_peer_a", "status", "peer_a", conns, dispatch)
    expect_identical(replay$counts, first$counts)
    expect_identical(replay$final_vector_root, first$final_vector_root)
    expect_identical(c(fixture$state$source_prepare, fixture$state$start), before)

    tampered <- first
    tampered$counts[[1L]] <- tampered$counts[[1L]] + 1
    expect_error(ds.vertDPFrequencyInference(tampered), "validated")
  }
})

test_that("real Synopsis categorical MI is plausible and Rock-replayable at K=2/3/5", {
  .synopsis_real_e2e_only("mi")
  server_ns <- .synopsis_describe_real_e2e_server()
  mi <- get(".dsvert_mi_synopsis_result_v1",
            asNamespace("dsVertClient"), inherits = FALSE)
  for (k in .synopsis_real_e2e_peer_counts()) {
    fixture <- .synopsis_mi_real_e2e_fixture(k, server_ns)
    on.exit(unlink(fixture$root, recursive = TRUE, force = TRUE), add = TRUE)
    conns <- stats::setNames(lapply(fixture$peers, function(peer) {
      structure(list(peer = peer), class = "dsvert_synopsis_real_e2e_connection")
    }), fixture$peers)
    dispatch <- .synopsis_describe_real_e2e_dispatch(fixture)
    run <- function(datasources, .aggregate) {
      dsVertClient:::.dsvert_dp_synopsis_vector_run(
        datasources, .aggregate = dispatch)
    }
    first <- mi(
      status ~ 1, "data_peer_a", NULL, 8L, "binomial", conns, dispatch,
      .run = run)
    expect_s3_class(first, "ds.vertMI")
    expect_identical(first$status, "ok")
    expect_identical(first$family, "binomial")
    expect_identical(first$outcome_levels, c("case", "control"))
    expect_true(is.finite(first$coefficients[["(Intercept)"]]))
    expect_true(all(is.finite(first$probabilities)))
    expect_equal(sum(first$probabilities), 1, tolerance = 1e-12)
    expect_gte(first$completed_count_dp, first$admitted_count_dp)
    expect_gte(first$missing_count_dp, 0)
    expect_false("completed_counts" %in% names(first))
    expect_identical(first$additional_privacy_cost,
                     c(epsilon = 0, delta = 0))
    expect_identical(c(fixture$state$source_prepare, fixture$state$start),
                     c(1L, 2L))
    before <- c(fixture$state$source_prepare, fixture$state$start)
    fixture$state$storage <- stats::setNames(lapply(fixture$peers, function(...) {
      new.env(parent = emptyenv())
    }), fixture$peers)
    replay <- mi(
      status ~ 1, "data_peer_a", NULL, 8L, "binomial", conns, dispatch,
      .run = run)
    expect_identical(replay$completed_draws_sha256,
                     first$completed_draws_sha256)
    expect_identical(replay$final_vector_root, first$final_vector_root)
    expect_identical(c(fixture$state$source_prepare, fixture$state$start), before)
  }
})

test_that("real Synopsis survival is plausible and Rock-replayable at K=2/3/5", {
  .synopsis_real_e2e_only("survival")
  server_ns <- .synopsis_describe_real_e2e_server()
  survival <- get(".dsvert_dp_survival_impl", asNamespace("dsVertClient"),
                  inherits = FALSE)
  for (k in .synopsis_real_e2e_peer_counts()) {
    fixture <- .synopsis_survival_real_e2e_fixture(k, server_ns)
    on.exit(unlink(fixture$root, recursive = TRUE, force = TRUE), add = TRUE)
    conns <- stats::setNames(lapply(fixture$peers, function(peer) {
      structure(list(peer = peer), class = "dsvert_synopsis_real_e2e_connection")
    }), fixture$peers)
    dispatch <- .synopsis_describe_real_e2e_dispatch(fixture)
    first <- survival("data_peer_a", "primary", "peer_a", conns, dispatch)
    expect_s3_class(first, "ds.vertDPSurvival")
    expect_true(isTRUE(first$released))
    expect_identical(fixture$state$source_prepare, 1L)
    expect_identical(fixture$state$start, 2L)
    expect_identical(first$release_provenance$designated_noise_peers,
                     as.list(fixture$peers[1:2]))
    expect_length(first$release_provenance$ordered_peer_pinset, k)
    expect_equal(first$curve$time, c(5, 10), tolerance = 0)
    expect_true(all(is.finite(first$curve$kaplan_meier)))
    expect_true(all(first$curve$kaplan_meier >= 0 &
                    first$curve$kaplan_meier <= 1))
    expect_true(all(is.finite(first$curve$nelson_aalen)))
    expect_true(all(first$curve$nelson_aalen >= 0))

    route_survival <- function(data_name, analysis_id, server = NULL,
                               datasources = NULL, .aggregate) {
      survival(data_name, analysis_id, server, datasources, dispatch)
    }
    public <- testthat::with_mocked_bindings(
      .dsvert_dp_survival_impl = route_survival,
      ds.vertDPSurvival("data_peer_a", "primary", "peer_a", conns),
      .package = "dsVertClient")
    expect_s3_class(public, "ds.vertDPSurvival")
    expect_identical(public$curve, first$curve)
    expect_identical(public$final_vector_root, first$final_vector_root)
    expect_identical(c(fixture$state$source_prepare, fixture$state$start),
                     c(1L, 2L))

    before <- c(fixture$state$source_prepare, fixture$state$start)
    kaplan_meier <- ds.vertDPKaplanMeier(first)
    nelson_aalen <- ds.vertDPNelsonAalen(first)
    cumulative_incidence <- ds.vertDPCumulativeIncidence(first, "event")
    rmst <- ds.vertDPRMST(first)
    rmtl <- ds.vertDPRMTL(first)
    quantiles <- ds.vertDPSurvivalQuantile(first, c(0.25, 0.5, 0.75))
    median <- ds.vertDPMedianSurvival(first)
    # A self-contrast retains the one signed artifact's simultaneous
    # mechanism event; it is not a second release or an independence claim.
    survival_contrast <- ds.vertDPSurvivalContrast(
      first, first, comparison_label = "same", reference_label = "baseline")
    rmst_contrast <- ds.vertDPRMSTContrast(
      first, first, tau = c(5, 10), comparison_label = "same",
      reference_label = "baseline")
    expect_true(all(is.finite(kaplan_meier$kaplan_meier)))
    expect_true(all(is.finite(nelson_aalen$nelson_aalen)))
    expect_true(all(is.finite(cumulative_incidence$cumulative_incidence)))
    expect_true(all(cumulative_incidence$cumulative_incidence >= 0 &
                    cumulative_incidence$cumulative_incidence <= 1))
    expect_true(all(is.finite(rmst$rmst) & rmst$rmst >= 0 &
                    rmst$rmst <= 10))
    expect_true(all(is.finite(rmtl$rmtl) & rmtl$rmtl >= 0 &
                    rmtl$rmtl <= 10))
    expect_identical(quantiles$probability, c(0.25, 0.5, 0.75))
    expect_identical(median$probability, 0.5)
    for (view in list(kaplan_meier, nelson_aalen, cumulative_incidence,
                      rmst, rmtl, quantiles, median)) {
      expect_identical(attr(view, "additional_server_calls"), 0L)
      expect_identical(attr(view, "additional_privacy_cost"),
                       c(epsilon = 0, delta = 0))
    }
    for (contrast in list(survival_contrast, rmst_contrast)) {
      expect_identical(attr(contrast, "additional_server_calls"), 0L)
      expect_identical(attr(contrast, "additional_privacy_cost"),
                       c(epsilon = 0, delta = 0))
      expect_identical(attr(contrast, "joint_event"),
                       "same_signed_survival_artifact")
    }
    expect_equal(survival_contrast$survival_difference,
                 rep(0, nrow(survival_contrast)), tolerance = 0)
    finite_survival_ratios <- is.finite(survival_contrast$survival_ratio)
    expect_true(all(survival_contrast$survival_ratio[
      finite_survival_ratios] == 1))
    expect_equal(rmst_contrast$rmst_difference,
                 rep(0, nrow(rmst_contrast)), tolerance = 0)
    finite_rmst_ratios <- is.finite(rmst_contrast$rmst_ratio)
    expect_true(all(rmst_contrast$rmst_ratio[finite_rmst_ratios] == 1))
    expect_identical(c(fixture$state$source_prepare, fixture$state$start), before)

    fixture$state$storage <- stats::setNames(lapply(fixture$peers, function(...) {
      new.env(parent = emptyenv())
    }), fixture$peers)
    replay <- survival("data_peer_a", "primary", "peer_a", conns, dispatch)
    expect_identical(replay$curve, first$curve)
    expect_identical(replay$final_vector_root, first$final_vector_root)
    expect_identical(c(fixture$state$source_prepare, fixture$state$start), before)
  }
})

test_that("real same-owner Synopsis correlation is plausible and Rock-replayable at K=2/3/5", {
  .synopsis_real_e2e_only("correlation")
  server_ns <- .synopsis_describe_real_e2e_server()
  correlation <- get(".dsvert_dp_cor_impl", asNamespace("dsVertClient"),
                     inherits = FALSE)
  for (k in .synopsis_real_e2e_peer_counts()) {
    fixture <- .synopsis_correlation_real_e2e_fixture(k, server_ns)
    on.exit(unlink(fixture$root, recursive = TRUE, force = TRUE), add = TRUE)
    conns <- stats::setNames(lapply(fixture$peers, function(peer) {
      structure(list(peer = peer), class = "dsvert_synopsis_real_e2e_connection")
    }), fixture$peers)
    dispatch <- .synopsis_describe_real_e2e_dispatch(fixture)
    first <- correlation(
      "data_peer_a", "data_peer_a::peer_a", c("x_peer_a", "y_peer_a"),
      "peer_a", conns, dispatch)
    expect_s3_class(first, "ds.vertDPCor")
    expect_identical(fixture$state$source_prepare, 1L)
    expect_identical(fixture$state$start, 2L)
    expect_identical(first$release_provenance$designated_noise_peers,
                     as.list(fixture$peers[1:2]))
    expect_length(first$release_provenance$ordered_peer_pinset, k)
    expect_true(all(is.finite(first$correlation_raw_pairwise)))
    expect_true(all(first$correlation_raw_pairwise >= -1 &
                    first$correlation_raw_pairwise <= 1))
    expect_true(all(is.finite(first$correlation)))
    expect_equal(first$correlation, t(first$correlation), tolerance = 1e-12)
    expect_equal(unname(diag(first$correlation)), c(1, 1), tolerance = 1e-12)
    expect_gte(min(eigen(first$correlation, symmetric = TRUE,
                         only.values = TRUE)$values), -1e-12)
    expect_identical(first$cross_owner_state, "reserved_not_materialized")
    expect_identical(first$additional_server_calls_after_synopsis, 0L)
    expect_identical(first$additional_privacy_cost, c(epsilon = 0, delta = 0))

    route_correlation <- function(data_name, analysis_id, variables = NULL,
                                  server = NULL, datasources = NULL,
                                  .aggregate) {
      correlation(data_name, analysis_id, variables, server, datasources,
                  dispatch)
    }
    public <- testthat::with_mocked_bindings(
      .dsvert_dp_cor_impl = route_correlation,
      ds.vertDPCor(
        "data_peer_a", "data_peer_a::peer_a", c("x_peer_a", "y_peer_a"),
        "peer_a", conns),
      .package = "dsVertClient")
    expect_s3_class(public, "ds.vertDPCor")
    expect_identical(public$correlation, first$correlation)
    expect_identical(c(fixture$state$source_prepare, fixture$state$start),
                     c(1L, 2L))

    before <- c(fixture$state$source_prepare, fixture$state$start)
    fixture$state$storage <- stats::setNames(lapply(fixture$peers, function(...) {
      new.env(parent = emptyenv())
    }), fixture$peers)
    replay <- correlation(
      "data_peer_a", "data_peer_a::peer_a", c("x_peer_a", "y_peer_a"),
      "peer_a", conns, dispatch)
    expect_identical(replay$correlation_raw_pairwise,
                     first$correlation_raw_pairwise)
    expect_identical(replay$correlation, first$correlation)
    expect_identical(replay$final_vector_root, first$final_vector_root)
    expect_identical(c(fixture$state$source_prepare, fixture$state$start), before)
  }
})

test_that("real same-owner Gaussian Synopsis and correlation are plausible and Rock-replayable at K=2/3/5", {
  .synopsis_real_e2e_only("gaussian")
  server_ns <- .synopsis_describe_real_e2e_server()
  gaussian <- get(".dsvert_dp_gaussian_impl", asNamespace("dsVertClient"),
                  inherits = FALSE)
  correlation <- get(".dsvert_dp_cor_gaussian_impl", asNamespace("dsVertClient"),
                     inherits = FALSE)
  for (k in .synopsis_real_e2e_peer_counts()) {
    fixture <- .synopsis_gaussian_real_e2e_fixture(k, server_ns)
    on.exit(unlink(fixture$root, recursive = TRUE, force = TRUE), add = TRUE)
    conns <- stats::setNames(lapply(fixture$peers, function(peer) {
      structure(list(peer = peer), class = "dsvert_synopsis_real_e2e_connection")
    }), fixture$peers)
    dispatch <- .synopsis_describe_real_e2e_dispatch(fixture)
    fit <- gaussian(
      "data_peer_a", "gaussian_primary", 0, "peer_a", conns, dispatch)

    expect_s3_class(fit, "ds.vertDPGaussian")
    expect_identical(fixture$state$source_prepare, 1L)
    expect_identical(fixture$state$start, 2L)
    expect_identical(fit$release_provenance$designated_noise_peers,
                     as.list(fixture$peers[1:2]))
    expect_length(fit$release_provenance$ordered_peer_pinset, k)
    expect_identical(fit$cross_owner_state, "reserved_not_materialized")
    expect_true(fit$provenance_integrity)
    expect_identical(fit$provenance_authenticity,
                     "session_transport_anchored")
    expect_true(is.finite(fit$n_obs) && fit$n_obs > 9000 && fit$n_obs < 11000)
    expect_true(all(is.finite(fit$coefficients_original_scale)))
    expect_equal(fit$coefficients_original_scale,
                 c(`(Intercept)` = 10, x_peer_a = -1), tolerance = 0.1)
    mechanism_region <- ds.vertConfint(fit, type = "mechanism")
    expect_identical(attr(mechanism_region, "interval_scope"),
                     "simultaneous_dp_mechanism")
    expect_identical(attr(mechanism_region, "sampling_inference"), FALSE)
    expect_true(all(is.finite(as.matrix(mechanism_region))))
    expect_true(all(mechanism_region$lower <= mechanism_region$estimate))
    expect_true(all(mechanism_region$estimate <= mechanism_region$upper))
    expect_true(all(mechanism_region$mechanism_radius >= 0))
    expect_identical(
      ds.vertConfint(fit, type = "mechanism"), mechanism_region)
    mechanism_wald <- ds.vertWald(
      fit, parm = "x_peer_a", null = 0, type = "mechanism")
    expect_identical(mechanism_wald$distribution,
                     "simultaneous_dp_mechanism_region")
    expect_false(mechanism_wald$sampling_inference)
    expect_null(mechanism_wald$p_value)
    expect_equal(mechanism_wald$lower,
                 mechanism_region["x_peer_a", "lower"])
    expect_equal(mechanism_wald$upper,
                 mechanism_region["x_peer_a", "upper"])
    mechanism_contrast <- ds.vertContrast(
      fit, K = c(0, 1), type = "mechanism")
    expect_identical(mechanism_contrast$distribution,
                     "simultaneous_dp_mechanism_region")
    expect_false(mechanism_contrast$sampling_inference)
    expect_null(mechanism_contrast$p_value)
    expect_equal(mechanism_contrast$estimate,
                 mechanism_region["x_peer_a", "estimate"])
    expect_equal(mechanism_contrast$lower,
                 mechanism_region["x_peer_a", "lower"])
    expect_equal(mechanism_contrast$upper,
                 mechanism_region["x_peer_a", "upper"])
    mechanism_lr <- ds.vertLR(
      reduced = "(Intercept)", full = fit, type = "mechanism")
    expect_identical(mechanism_lr$distribution,
                     "simultaneous_dp_mechanism_rss_reduction")
    expect_false(mechanism_lr$sampling_inference)
    expect_null(mechanism_lr$p_value)
    expect_true(is.finite(mechanism_lr$lower))
    expect_true(is.finite(mechanism_lr$upper))
    expect_gte(mechanism_lr$lower, 0)
    expect_lte(mechanism_lr$lower, mechanism_lr$statistic)
    expect_gte(mechanism_lr$upper, mechanism_lr$statistic)

    route_gaussian <- function(data_name, analysis_id, ridge = 0,
                               server = NULL, datasources = NULL,
                               .aggregate) {
      gaussian(data_name, analysis_id, ridge, server, datasources, dispatch)
    }
    public <- testthat::with_mocked_bindings(
      .dsvert_dp_gaussian_impl = route_gaussian,
      ds.vertDPGaussian(
        "data_peer_a", "gaussian_primary", ridge = 0, server = "peer_a",
        datasources = conns),
      .package = "dsVertClient")
    expect_s3_class(public, "ds.vertDPGaussian")
    expect_identical(public$coefficients_original_scale,
                     fit$coefficients_original_scale)
    expect_identical(public$final_vector_root, fit$final_vector_root)
    expect_identical(c(fixture$state$source_prepare, fixture$state$start),
                     c(1L, 2L))

    adapter <- testthat::with_mocked_bindings(
      ds.vertDPGaussian = function(
          data_name, analysis_id, ridge = 0, server = NULL,
          datasources = NULL) {
        gaussian(data_name, analysis_id, ridge, server, datasources, dispatch)
      },
      ds.vertGLM(
        y_peer_a ~ x_peer_a, data = "data_peer_a", family = "gaussian",
        dp_analysis_id = "gaussian_primary",
        missing = "complete_case_capsule", verbose = FALSE,
        datasources = conns),
      .package = "dsVertClient")
    expect_s3_class(adapter, "ds.vertDPGaussian")
    expect_identical(adapter$called_via,
                     "ds.vertGLM_explicit_dp_analysis_id")
    expect_identical(adapter$legacy_glm_estimand, FALSE)
    expect_identical(adapter$provenance_certificate$certificate_sha256,
                     fit$provenance_certificate$certificate_sha256)
    expect_equal(adapter$coefficients_original_scale,
                 fit$coefficients_original_scale, tolerance = 0)
    expect_identical(c(fixture$state$source_prepare, fixture$state$start),
                     c(1L, 2L))

    legacy_glm <- testthat::with_mocked_bindings(
      ds.vertDPGaussian = function(
          data_name, analysis_id, ridge = 0, server = NULL,
          datasources = NULL) {
        gaussian(data_name, analysis_id, ridge, server, datasources, dispatch)
      },
      ds.vert.glm(
        y_peer_a ~ x_peer_a, data = "data_peer_a", family = "gaussian",
        dp_analysis_id = "gaussian_primary",
        missing = "complete_case_capsule", verbose = FALSE,
        datasources = conns),
      .package = "dsVertClient")
    expect_s3_class(legacy_glm, "ds.vertDPGaussian")
    expect_identical(legacy_glm$frontdoor, "ds.vert.glm")
    expect_identical(legacy_glm$route, "ds.vertGLM")
    expect_equal(legacy_glm$coefficients_original_scale,
                 fit$coefficients_original_scale, tolerance = 0)
    expect_identical(c(fixture$state$source_prepare, fixture$state$start),
                     c(1L, 2L))

    # The protected Gaussian release, its public GLM route, correlation, and
    # replay are verified for every supported topology below.  The remaining
    # LASSO calls are deterministic post-processing of that already signed
    # artifact: running their full 2,000-iteration alias matrix again at K=3
    # and K=5 adds no MPC, privacy, or topology coverage.
    if (identical(k, 2L)) {
    iterative <- testthat::with_mocked_bindings(
      ds.vertDPGaussian = function(
          data_name, analysis_id, ridge = 0, server = NULL,
          datasources = NULL) {
        gaussian(data_name, analysis_id, ridge, server, datasources, dispatch)
      },
      ds.vertLASSOIter(
        y_peer_a ~ x_peer_a, data = "data_peer_a", family = "gaussian",
        lambda = c(0.1, 0.05, 0.01), max_outer = 20L, tol = 1e-10,
        dp_analysis_id = "gaussian_primary", verbose = FALSE,
        datasources = conns),
      .package = "dsVertClient")
    expect_s3_class(iterative, "ds.vertLASSOIter")
    expect_s3_class(iterative, "ds.vertDPLASSOPath")
    expect_identical(iterative$input_provenance,
                     "signed_dp_gaussian_synopsis")
    expect_identical(iterative$additional_privacy_cost,
                     c(epsilon = 0, delta = 0))
    expect_true(all(vapply(iterative$path_certificates, function(value) {
      isTRUE(value$kkt$satisfied)
    }, logical(1L))))
    expect_identical(c(fixture$state$source_prepare, fixture$state$start),
                     c(1L, 2L))

    iterative_alias <- testthat::with_mocked_bindings(
      ds.vertDPGaussian = function(
          data_name, analysis_id, ridge = 0, server = NULL,
          datasources = NULL) {
        gaussian(data_name, analysis_id, ridge, server, datasources, dispatch)
      },
      ds.vert.lasso_iter(
        y_peer_a ~ x_peer_a, data = "data_peer_a", family = "gaussian",
        lambda = c(0.1, 0.05, 0.01), max_outer = 20L, tol = 1e-10,
        dp_analysis_id = "gaussian_primary", verbose = FALSE,
        datasources = conns),
      .package = "dsVertClient")
    expect_s3_class(iterative_alias, "ds.vertDPLASSOPath")
    expect_identical(iterative_alias$frontdoor, "ds.vert.lasso_iter")
    expect_identical(iterative_alias$route,
                     "ds.vertLASSOIter(signed-gaussian-synopsis)")
    expect_equal(iterative_alias$paths, iterative$paths, tolerance = 1e-12)
    expect_identical(c(fixture$state$source_prepare, fixture$state$start),
                     c(1L, 2L))

    lasso <- ds.vertLASSOProximal(
      fit, lambda = 0.05, max_iter = 500L, tol = 1e-10)
    expect_s3_class(lasso, "ds.vertLASSOProximal")
    expect_true(lasso$converged)
    expect_true(lasso$kkt$satisfied)
    expect_identical(lasso$input_provenance,
                     "signed_dp_gaussian_synopsis")
    expect_identical(lasso$additional_server_calls_after_synopsis, 0L)
    expect_identical(lasso$additional_privacy_cost,
                     c(epsilon = 0, delta = 0))
    expect_true(all(is.finite(lasso$coefficients)))
    expect_lt(lasso$coefficients[["x_peer_a"]], 0)
    expect_gte(abs(lasso$coefficients[["x_peer_a"]]), 0.25)
    expect_lte(abs(lasso$coefficients[["x_peer_a"]]),
               abs(fit$coefficients_original_scale[["x_peer_a"]]) + 1e-8)
    expect_lte(abs(lasso$coefficients[["(Intercept)"]] -
                   fit$coefficients_original_scale[["(Intercept)"]]), 1.5)
    expect_identical(c(fixture$state$source_prepare, fixture$state$start),
                     c(1L, 2L))

    selection <- ds.vertLASSOCV(
      fit, lambda_grid = c(0.1, 0.05, 0.01), criterion = "BIC")
    expect_s3_class(selection, "ds.vertLASSOCV")
    expect_false(selection$cross_validation)
    expect_false(selection$one_standard_error_rule)
    expect_identical(selection$input_provenance,
                     "signed_dp_gaussian_synopsis")
    expect_identical(selection$additional_server_calls_after_synopsis, 0L)
    expect_identical(selection$additional_privacy_cost,
                     c(epsilon = 0, delta = 0))
    expect_identical(c(fixture$state$source_prepare, fixture$state$start),
                     c(1L, 2L))

    lasso_path <- ds.vertLASSO(
      fit, lambda_1 = 0.05, alpha_grid = c(1, 0.5, 0.1),
      max_iter = 500L, tol = 1e-10)
    one_step <- ds.vertLASSO1Step(
      fit, lambda = c(0.1, 0.05, 0.01), max_iter = 500L, tol = 1e-10)
    lasso_alias <- ds.vert.lasso(
      fit, lambda_1 = 0.05, alpha_grid = c(1, 0.5, 0.1),
      max_iter = 500L, tol = 1e-10)
    proximal_alias <- ds.vert.lasso_proximal(
      fit, lambda = 0.05, max_iter = 500L, tol = 1e-10)
    one_step_alias <- ds.vert.lasso_1step(
      fit, lambda = c(0.1, 0.05, 0.01), max_iter = 500L, tol = 1e-10)
    selection_alias <- ds.vert.lasso_cv(
      fit, lambda_grid = c(0.1, 0.05, 0.01), criterion = "BIC")
    expect_identical(lasso_alias$frontdoor, "ds.vert.lasso")
    expect_identical(proximal_alias$frontdoor, "ds.vert.lasso_proximal")
    expect_identical(one_step_alias$frontdoor, "ds.vert.lasso_1step")
    expect_identical(selection_alias$frontdoor, "ds.vert.lasso_cv")
    expect_equal(lasso_alias$paths, lasso_path$paths, tolerance = 1e-12)
    expect_equal(proximal_alias$coefficients, lasso$coefficients,
                 tolerance = 1e-12)
    expect_equal(one_step_alias$paths, one_step$paths, tolerance = 1e-12)
    expect_equal(selection_alias$lambda.min, selection$lambda.min,
                 tolerance = 1e-12)
    expect_identical(c(fixture$state$source_prepare, fixture$state$start),
                     c(1L, 2L))
    }

    cor <- correlation(
      "data_peer_a", "gaussian_primary", c("x_peer_a", "y_peer_a"),
      "peer_a", conns, dispatch)
    expect_equal(unname(diag(cor$correlation)), c(1, 1), tolerance = 1e-12)
    expect_equal(cor$correlation["x_peer_a", "y_peer_a"], -1,
                 tolerance = 0.05)
    expect_identical(c(fixture$state$source_prepare, fixture$state$start),
                     c(1L, 2L))

    pca <- ds.vertPCA(cor_result = cor, n_components = 2L, verbose = FALSE)

    # Exercise the user-facing correlation and PCA calls over the same real
    # Synopsis artifact.  The harness supplies its authenticated test
    # dispatch at the DSI boundary; the public functions must not create a
    # second release or route through a legacy correlation implementation.
    route_correlation <- function(data_name, analysis_id, variables = NULL,
                                  server = NULL, datasources = NULL,
                                  .aggregate) {
      correlation(data_name, analysis_id, variables, server, datasources,
                  dispatch)
    }
    public <- testthat::with_mocked_bindings(
      .dsvert_dp_cor_gaussian_impl = route_correlation,
      list(
        cor = ds.vertCor(
          "data_peer_a", c("x_peer_a", "y_peer_a"),
          analysis_id = "gaussian_primary", verbose = FALSE,
          datasources = conns),
        cor_alias = ds.vert.cor(
          "data_peer_a", c("x_peer_a", "y_peer_a"),
          analysis_id = "gaussian_primary", verbose = FALSE,
          datasources = conns),
        pca = ds.vertPCA(
          "data_peer_a", c("x_peer_a", "y_peer_a"), n_components = 2L,
          analysis_id = "gaussian_primary", verbose = FALSE,
          datasources = conns),
        pca_alias = ds.vert.pca(
          "data_peer_a", c("x_peer_a", "y_peer_a"), n_components = 2L,
          analysis_id = "gaussian_primary", verbose = FALSE,
          datasources = conns)),
      .package = "dsVertClient")
    expect_s3_class(public$cor, "ds.vertDPCor")
    expect_identical(public$cor$correlation, cor$correlation)
    expect_identical(public$cor_alias$frontdoor, "ds.vert.cor")
    expect_identical(public$cor_alias$route, "ds.vertCor")
    expect_identical(public$cor_alias$correlation, cor$correlation)
    expect_s3_class(public$pca, "ds.pca")
    expect_identical(public$pca$eigenvalues, pca$eigenvalues)
    expect_identical(public$pca_alias$frontdoor, "ds.vert.pca")
    expect_identical(public$pca_alias$route, "ds.vertPCA")
    expect_identical(public$pca_alias$eigenvalues, public$pca$eigenvalues)
    expect_identical(c(fixture$state$source_prepare, fixture$state$start),
                     c(1L, 2L))

    expect_s3_class(pca, "ds.pca")
    expect_true(all(is.finite(pca$eigenvalues) & pca$eigenvalues >= 0))
    expect_gt(pca$variance_pct[[1L]], 95)
    expect_identical(pca$additional_server_calls, 0L)
    expect_identical(pca$additional_server_calls_after_synopsis, 0L)
    expect_identical(pca$additional_privacy_cost,
                     c(epsilon = 0, delta = 0))
    pca_alias <- ds.vert.pca(
      cor_result = cor, n_components = 2L, verbose = FALSE)
    expect_identical(pca_alias$frontdoor, "ds.vert.pca")
    expect_identical(pca_alias$route, "ds.vertPCA")
    expect_identical(pca_alias$eigenvalues, pca$eigenvalues)
    expect_identical(pca_alias$loadings, pca$loadings)
    expect_identical(pca_alias$additional_server_calls_after_synopsis, 0L)
    expect_identical(pca_alias$additional_privacy_cost,
                     c(epsilon = 0, delta = 0))
    expect_identical(c(fixture$state$source_prepare, fixture$state$start),
                     c(1L, 2L))

    before <- c(fixture$state$source_prepare, fixture$state$start)
    fixture$state$storage <- stats::setNames(lapply(fixture$peers, function(...) {
      new.env(parent = emptyenv())
    }), fixture$peers)
    replay <- gaussian(
      "data_peer_a", "gaussian_primary", 0, "peer_a", conns, dispatch)
    expect_identical(serialize(replay, NULL, version = 3L),
                     serialize(fit, NULL, version = 3L))
    expect_identical(c(fixture$state$source_prepare, fixture$state$start), before)
  }
})

test_that("real random-intercept LMM Synopsis is plausible and Rock-replayable at K=2/3/5", {
  .synopsis_real_e2e_only("lmm")
  server_ns <- .synopsis_describe_real_e2e_server()
  lmm <- get(".dsvert_dp_lmm_impl", asNamespace("dsVertClient"),
             inherits = FALSE)
  for (k in .synopsis_real_e2e_peer_counts()) {
    fixture <- .synopsis_lmm_real_e2e_fixture(k, server_ns)
    on.exit(unlink(fixture$root, recursive = TRUE, force = TRUE), add = TRUE)
    conns <- stats::setNames(lapply(fixture$peers, function(peer) {
      structure(list(peer = peer), class = "dsvert_synopsis_real_e2e_connection")
    }), fixture$peers)
    dispatch <- .synopsis_describe_real_e2e_dispatch(fixture)
    fit <- lmm("data_peer_a", "lmm_primary", "peer_a", conns, dispatch)

    expect_s3_class(fit, "ds.vertDPLMM")
    expect_identical(fixture$state$source_prepare, 1L)
    expect_identical(fixture$state$start, 2L)
    expect_identical(fit$cross_owner_state, "reserved_not_materialized")
    expect_identical(fit$additional_server_calls_after_synopsis, 0L)
    expect_identical(fit$additional_privacy_cost, c(epsilon = 0, delta = 0))
    expect_identical(fit$status, "ok")
    expect_true(all(is.finite(c(fit$coefficients, fit$sigma2,
                                fit$sigma_b2, fit$icc))))
    expect_true(fit$coefficients[["(Intercept)"]] > 3 &&
                fit$coefficients[["(Intercept)"]] < 7)
    expect_true(fit$sigma2 >= 0 && fit$sigma_b2 >= 0 &&
                fit$icc >= 0 && fit$icc <= 1)
    expect_true(isTRUE(fit$provenance_certificate$
      public_dp_coordinates_only))

    route_lmm <- function(data_name, analysis_id, server = NULL,
                          datasources = NULL, .aggregate) {
      lmm(data_name, analysis_id, server, datasources, dispatch)
    }
    public <- testthat::with_mocked_bindings(
      .dsvert_dp_lmm_impl = route_lmm,
      list(
        typed = ds.vertDPLMM(
          "data_peer_a", "lmm_primary", server = "peer_a",
          datasources = conns),
        legacy = ds.vertLMM(
          y_peer_a ~ 1, data = "data_peer_a", cluster_col = "site_peer_a",
          analysis_id = "lmm_primary", reml = FALSE, datasources = conns),
        alias = ds.vert.lmm(
          y_peer_a ~ 1, data = "data_peer_a", cluster_col = "site_peer_a",
          analysis_id = "lmm_primary", datasources = conns),
        k3 = if (length(conns) >= 3L) ds.vertLMM.k3(
          y_peer_a ~ 1, data = "data_peer_a", cluster_col = "site_peer_a",
          analysis_id = "lmm_primary", datasources = conns) else NULL),
      .package = "dsVertClient")
    expect_s3_class(public$typed, "ds.vertDPLMM")
    expect_s3_class(public$legacy, "ds.vertLMM")
    expect_s3_class(public$alias, "ds.vertLMM")
    expect_identical(public$typed$coefficients, fit$coefficients)
    expect_identical(public$legacy$coefficients, fit$coefficients)
    expect_identical(public$alias$coefficients, fit$coefficients)
    expect_identical(public$typed$final_vector_root, fit$final_vector_root)
    expect_identical(public$legacy$reml, FALSE)
    expect_identical(public$legacy$cluster_sizes, NULL)
    if (k >= 3L) {
      expect_s3_class(public$k3, "ds.vertLMM")
      expect_identical(public$k3$coefficients, fit$coefficients)
      expect_identical(public$k3$frontdoor, "ds.vertLMM.k3")
    } else {
      expect_null(public$k3)
    }
    expect_identical(c(fixture$state$source_prepare, fixture$state$start),
                     c(1L, 2L))

    before <- c(fixture$state$source_prepare, fixture$state$start)
    fixture$state$storage <- stats::setNames(lapply(fixture$peers, function(...) {
      new.env(parent = emptyenv())
    }), fixture$peers)
    replay <- lmm("data_peer_a", "lmm_primary", "peer_a", conns, dispatch)
    expect_identical(replay$coefficients, fit$coefficients)
    expect_identical(replay$final_vector_root, fit$final_vector_root)
    expect_identical(c(fixture$state$source_prepare, fixture$state$start), before)
  }
})

test_that("real binary random-intercept GLMM moment route is plausible and Rock-replayable at K=2/3/5", {
  .synopsis_real_e2e_only("lmm")
  server_ns <- .synopsis_describe_real_e2e_server()
  lmm <- get(".dsvert_dp_lmm_impl", asNamespace("dsVertClient"),
             inherits = FALSE)
  for (k in .synopsis_real_e2e_peer_counts()) {
    fixture <- .synopsis_glmm_real_e2e_fixture(k, server_ns)
    on.exit(unlink(fixture$root, recursive = TRUE, force = TRUE), add = TRUE)
    conns <- stats::setNames(lapply(fixture$peers, function(peer) {
      structure(list(peer = peer), class = "dsvert_synopsis_real_e2e_connection")
    }), fixture$peers)
    dispatch <- .synopsis_describe_real_e2e_dispatch(fixture)
    route_lmm <- function(data_name, analysis_id, server = NULL,
                          datasources = NULL, .aggregate) {
      lmm(data_name, analysis_id, server, datasources, dispatch)
    }
    public <- testthat::with_mocked_bindings(
      .dsvert_dp_lmm_impl = route_lmm,
      list(
        direct = ds.vertGLMM(
          y_peer_a ~ 1, data = "data_peer_a", cluster_col = "site_peer_a",
          analysis_id = "lmm_primary", datasources = conns),
        alias = ds.vert.glmm(
          y_peer_a ~ 1, data = "data_peer_a", cluster_col = "site_peer_a",
          analysis_id = "lmm_primary", datasources = conns)),
      .package = "dsVertClient")
    fit <- public$direct
    expect_s3_class(fit, "ds.vertGLMM")
    expect_identical(fit$family, "binomial")
    expect_identical(
      fit$estimand,
      "binary_random_intercept_population_average_moment_approximation")
    expect_true(is.finite(fit$coefficients[["(Intercept)"]]))
    expect_true(fit$marginal_probability > 0.40 &&
                fit$marginal_probability < 0.60)
    expect_true(fit$icc_observed > 0 && fit$icc_observed < 1)
    expect_true(is.finite(fit$sigma_b2) && fit$sigma_b2 > 0)
    expect_null(fit$standard_errors)
    expect_null(fit$cluster_sizes)
    expect_identical(public$alias$coefficients, fit$coefficients)
    expect_identical(public$alias$frontdoor, "ds.vert.glmm")
    expect_identical(c(fixture$state$source_prepare, fixture$state$start),
                     c(1L, 2L))

    before <- c(fixture$state$source_prepare, fixture$state$start)
    fixture$state$storage <- stats::setNames(lapply(fixture$peers, function(...) {
      new.env(parent = emptyenv())
    }), fixture$peers)
    replay <- testthat::with_mocked_bindings(
      .dsvert_dp_lmm_impl = route_lmm,
      ds.vertGLMM(
        y_peer_a ~ 1, data = "data_peer_a", cluster_col = "site_peer_a",
        analysis_id = "lmm_primary", datasources = conns),
      .package = "dsVertClient")
    expect_identical(replay$coefficients, fit$coefficients)
    expect_identical(replay$provenance_certificate$certificate_sha256,
                     fit$provenance_certificate$certificate_sha256)
    expect_identical(c(fixture$state$source_prepare, fixture$state$start), before)
  }
})

test_that("focused Gaussian LASSO pseudo-IC is plausible and replayable at K=2/3/5", {
  .synopsis_real_e2e_focal_only()
  server_ns <- .synopsis_describe_real_e2e_server()
  gaussian <- get(".dsvert_dp_gaussian_impl", asNamespace("dsVertClient"),
                  inherits = FALSE)
  for (k in .synopsis_real_e2e_peer_counts()) {
    fixture <- .synopsis_gaussian_real_e2e_fixture(k, server_ns, n = 512L)
    on.exit(unlink(fixture$root, recursive = TRUE, force = TRUE), add = TRUE)
    conns <- stats::setNames(lapply(fixture$peers, function(peer) {
      structure(list(peer = peer), class = "dsvert_synopsis_real_e2e_connection")
    }), fixture$peers)
    fit <- gaussian(
      "data_peer_a", "gaussian_primary", 0, "peer_a", conns,
      .synopsis_describe_real_e2e_dispatch(fixture))
    selection <- ds.vertLASSOCV(
      fit, lambda_grid = c(0.1, 0.05, 0.01), criterion = "BIC")
    path <- ds.vertLASSO(
      fit, lambda_1 = 0.1, alpha_grid = c(1, 0.5, 0.1),
      max_iter = 2000L, tol = 1e-10)
    one_step <- ds.vertLASSO1Step(
      fit, lambda = c(0.1, 0.05, 0.01),
      max_iter = 2000L, tol = 1e-10)

    expect_s3_class(selection, "ds.vertLASSOCV")
    expect_s3_class(path, "ds.vertDPLASSOPath")
    expect_s3_class(one_step, "ds.vertLASSO1Step")
    expect_s3_class(one_step, "ds.vertDPLASSOPath")
    expect_identical(fit$release_provenance$designated_noise_peers,
                     as.list(fixture$peers[1:2]))
    expect_length(fit$release_provenance$ordered_peer_pinset, k)
    expect_true(selection$selection_available, info = paste("K =", k))
    expect_false(selection$cross_validation, info = paste("K =", k))
    expect_false(selection$one_standard_error_rule, info = paste("K =", k))
    expect_equal(selection$lambda.min, 0.01, tolerance = 1e-12,
                 info = paste("K =", k))
    expect_equal(selection$lambda.parsimonious, 0.01, tolerance = 1e-12,
                 info = paste("K =", k))
    expect_true(all(vapply(selection$path_certificates, function(value) {
      isTRUE(value$kkt$satisfied)
    }, logical(1L))), info = paste("K =", k))
    expect_gt(selection$beta.min[["(Intercept)"]], 5)
    expect_lt(selection$beta.min[["(Intercept)"]], 15)
    expect_lt(selection$beta.min[["x_peer_a"]], -0.25)
    expect_identical(selection$additional_server_calls_after_synopsis, 0L)
    expect_identical(selection$additional_privacy_cost,
                     c(epsilon = 0, delta = 0))
    expect_identical(path$input_provenance,
                     "signed_dp_gaussian_synopsis")
    expect_true(all(vapply(path$path_certificates, function(value) {
      isTRUE(value$kkt$satisfied)
    }, logical(1L))), info = paste("K =", k))
    expect_lt(path$paths[[1L]][["x_peer_a"]], -0.25)
    expect_identical(path$additional_server_calls_after_synopsis, 0L)
    expect_identical(path$additional_privacy_cost,
                     c(epsilon = 0, delta = 0))
    expect_identical(one_step$input_provenance,
                     "signed_dp_gaussian_synopsis")
    expect_equal(one_step$lambda, c(0.1, 0.05, 0.01), tolerance = 1e-12)
    expect_equal(unname(one_step$paths), unname(path$paths),
                 tolerance = 1e-12)
    expect_true(all(vapply(one_step$path_certificates, function(value) {
      isTRUE(value$kkt$satisfied)
    }, logical(1L))), info = paste("K =", k))
    expect_identical(one_step$additional_server_calls_after_synopsis, 0L)

    before <- c(fixture$state$source_prepare, fixture$state$start)
    fixture$state$storage <- stats::setNames(lapply(fixture$peers, function(...) {
      new.env(parent = emptyenv())
    }), fixture$peers)
    replay <- gaussian(
      "data_peer_a", "gaussian_primary", 0, "peer_a", conns,
      .synopsis_describe_real_e2e_dispatch(fixture))
    replay_selection <- ds.vertLASSOCV(
      replay, lambda_grid = c(0.1, 0.05, 0.01), criterion = "BIC")
    replay_path <- ds.vertLASSO(
      replay, lambda_1 = 0.1, alpha_grid = c(1, 0.5, 0.1),
      max_iter = 2000L, tol = 1e-10)
    replay_one_step <- ds.vertLASSO1Step(
      replay, lambda = c(0.1, 0.05, 0.01),
      max_iter = 2000L, tol = 1e-10)
    expect_identical(serialize(replay_selection, NULL, version = 3L),
                     serialize(selection, NULL, version = 3L))
    expect_identical(serialize(replay_path, NULL, version = 3L),
                     serialize(path, NULL, version = 3L))
    expect_identical(serialize(replay_one_step, NULL, version = 3L),
                     serialize(one_step, NULL, version = 3L))
    expect_identical(c(fixture$state$source_prepare, fixture$state$start),
                     before)
  }
})

test_that("cross-owner Synopsis rejects a tampered witness before mutation", {
  .synopsis_real_e2e_only("cross_owner_tamper")
  server_ns <- .synopsis_describe_real_e2e_server()
  contingency <- get(".dsvert_dp_contingency_impl",
                     asNamespace("dsVertClient"), inherits = FALSE)
  peer_counts <- intersect(.synopsis_real_e2e_peer_counts(), c(3L, 5L))
  if (!length(peer_counts)) skip("cross-owner witness tamper requires K=3 or K=5")
  for (k in peer_counts) {
    fixture <- .synopsis_cross_contingency_real_e2e_fixture(k, server_ns)
    on.exit(unlink(fixture$root, recursive = TRUE, force = TRUE), add = TRUE)
    witness <- fixture$peers[[3L]]
    fixture$pins[[witness]] <- fixture$pins[["peer_a"]]
    conns <- stats::setNames(lapply(fixture$peers, function(peer) {
      structure(list(peer = peer), class = "dsvert_synopsis_real_e2e_connection")
    }), fixture$peers)

    expect_error(.dsvert_dp_contingency_impl(
      "data_peer_a", "peer_a$disease", "peer_b$exposure", NULL,
      conns, .synopsis_describe_real_e2e_dispatch(fixture)))
    expect_identical(fixture$state$source_prepare, 0L)
    expect_identical(fixture$state$start, 0L)
  }
})
