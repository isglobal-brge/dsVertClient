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
  "frequency", "survival", "correlation", "gaussian",
  "cross_owner_tamper", "gaussian_lasso_focal")
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
    correlations = list())
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
  for (k in c(2L, 3L, 5L)) {
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
  for (k in c(2L, 3L, 5L)) {
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

    before <- c(fixture$state$source_prepare, fixture$state$start)
    chisq <- ds.vertChisq(first, simulations = 128L)
    fisher <- ds.vertFisher(first, simulations = 128L)
    expect_s3_class(chisq, "ds.vertChisq")
    expect_s3_class(fisher, "ds.vertFisher")
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
    expect_identical(c(fixture$state$source_prepare, fixture$state$start), before)

    fixture$state$storage <- stats::setNames(lapply(fixture$peers, function(...) {
      new.env(parent = emptyenv())
    }), fixture$peers)
    replay <- contingency("data_peer_a", "exposure", "outcome", "peer_a",
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
  for (k in c(2L, 3L, 5L)) {
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
  for (k in c(2L, 3L, 5L)) {
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

test_that("real Synopsis survival is plausible and Rock-replayable at K=2/3/5", {
  .synopsis_real_e2e_only("survival")
  server_ns <- .synopsis_describe_real_e2e_server()
  survival <- get(".dsvert_dp_survival_impl", asNamespace("dsVertClient"),
                  inherits = FALSE)
  for (k in c(2L, 3L, 5L)) {
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

    before <- c(fixture$state$source_prepare, fixture$state$start)
    kaplan_meier <- ds.vertDPKaplanMeier(first)
    nelson_aalen <- ds.vertDPNelsonAalen(first)
    cumulative_incidence <- ds.vertDPCumulativeIncidence(first, "event")
    rmst <- ds.vertDPRMST(first)
    rmtl <- ds.vertDPRMTL(first)
    quantiles <- ds.vertDPSurvivalQuantile(first, c(0.25, 0.5, 0.75))
    median <- ds.vertDPMedianSurvival(first)
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
  for (k in c(2L, 3L, 5L)) {
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
  for (k in c(2L, 3L, 5L)) {
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

    lasso <- ds.vertLASSOProximal(
      fit, lambda = 0.05, max_iter = 2000L, tol = 1e-10)
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

    cor <- correlation(
      "data_peer_a", "gaussian_primary", c("x_peer_a", "y_peer_a"),
      "peer_a", conns, dispatch)
    expect_equal(unname(diag(cor$correlation)), c(1, 1), tolerance = 1e-12)
    expect_equal(cor$correlation["x_peer_a", "y_peer_a"], -1,
                 tolerance = 0.05)
    expect_identical(c(fixture$state$source_prepare, fixture$state$start),
                     c(1L, 2L))

    pca <- ds.vertPCA(cor_result = cor, n_components = 2L, verbose = FALSE)
    expect_s3_class(pca, "ds.pca")
    expect_true(all(is.finite(pca$eigenvalues) & pca$eigenvalues >= 0))
    expect_gt(pca$variance_pct[[1L]], 95)
    expect_identical(pca$additional_server_calls, 0L)
    expect_identical(pca$additional_server_calls_after_synopsis, 0L)
    expect_identical(pca$additional_privacy_cost,
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

test_that("focused Gaussian LASSO pseudo-IC is plausible and replayable at K=2/3/5", {
  .synopsis_real_e2e_focal_only()
  server_ns <- .synopsis_describe_real_e2e_server()
  gaussian <- get(".dsvert_dp_gaussian_impl", asNamespace("dsVertClient"),
                  inherits = FALSE)
  for (k in c(2L, 3L, 5L)) {
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
  for (k in c(3L, 5L)) {
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
