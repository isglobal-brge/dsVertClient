# Process-isolated DSLite smoke. Run under an external timeout of 900 s.
# Synthetic data only. The vertcor03 aggregate wall-time cap is 45 minutes.
args <- commandArgs(TRUE)
root <- if (length(args)) normalizePath(args[[1L]]) else normalizePath("..")
client <- file.path(root, "dsVertClient")
server <- file.path(root, "dsVert")
out <- file.path(root, ".tmp", "vertcor")
source(file.path(client, "inst/validation/v1.3.0/worked_example_custodian.R"))
pkgload::load_all(client, quiet = TRUE)

# Preserve the already-public signed DP evidence even if a downstream
# certificate check rejects it. This observer never substitutes a return value.
options(vertcor142.evidence_path = file.path(out, "smoke-latest-certificate.rds"))
trace(".dsvert_dp_gaussian_synopsis_trusted", where = asNamespace("dsVertClient"),
  tracer = quote(saveRDS(certificate, getOption("vertcor142.evidence_path"))), print = FALSE)

# Record the actual endpoint error inside its custodian process. The previous
# transcript-stop probe could retain a benign publication polling condition.
# This observes the original allowlisted operation; it never replaces a check.
.we_worker_submit <- function(worker, operation, args, async) {
  if (!identical(worker$get_state(), "idle"))
    stop("The isolated DSLite peer already has a pending command", call. = FALSE)
  job <- switch(operation,
    aggregate = function(expr) withCallingHandlers(
      get(".dsvert_e2e_server", .GlobalEnv)$aggregate(
        get(".dsvert_e2e_sid", .GlobalEnv), expr),
      error = function(e) assign(".synthetic_endpoint_failure",
        list(operation = deparse(expr[[1L]]), message = conditionMessage(e)), .GlobalEnv)),
    assign = function(symbol, expr) {
      get(".dsvert_e2e_server", .GlobalEnv)$assignExpr(
        get(".dsvert_e2e_sid", .GlobalEnv), symbol, expr)
      NULL
    })
  if (isTRUE(async)) { worker$call(job, args = args); return(NULL) }
  worker$run(job, args = args)
}
started <- Sys.time()
prior_seconds <- as.numeric(Sys.getenv("DSVERT_SMOKE_PRIOR_SECONDS", "0"))
prior_analyses <- as.integer(Sys.getenv("DSVERT_SMOKE_PRIOR_ANALYSES", "0"))
stopifnot(is.finite(prior_seconds), prior_seconds >= 0, prior_seconds < 2700,
          !is.na(prior_analyses), prior_analyses >= 0L, prior_analyses < 100L)
invocations <- prior_analyses
event_path <- file.path(out, "vertcor03-analysis-events.jsonl")
record_event <- function(label, status, seconds, detail = NULL) {
  event <- list(label = label, status = status, seconds = seconds, detail = detail)
  cat(jsonlite::toJSON(event, auto_unbox = TRUE, null = "null"), "\n",
      file = event_path, append = TRUE)
}
run_analysis <- function(label, code) {
  stopifnot(prior_seconds + as.numeric(difftime(Sys.time(), started, units = "secs")) < 2680)
  invocations <<- invocations + 1L
  cat("ANALYSIS", invocations, label, "\n"); flush.console()
  began <- Sys.time()
  value <- tryCatch(force(code), error = identity)
  seconds <- as.numeric(difftime(Sys.time(), began, units = "secs"))
  if (inherits(value, "error")) {
    record_event(label, "failed", seconds, conditionMessage(value))
    stop(value)
  }
  record_event(label, "passed", seconds)
  cat("PASS", label, "seconds", seconds, "\n"); flush.console()
  value
}
observe_peers <- function(peers, capacity = 64L) {
  for (peer in peers) peer$worker$run(function(capacity) {
    options(dsvert.psi.max_input_ids = max(64L, capacity))
    invisible(NULL)
  }, args = list(capacity = capacity))
}
report_peer_errors <- function(peers) {
  for (owner in names(peers)) {
    diagnostic <- try(peers[[owner]]$worker$run(function() {
      get0(".synthetic_endpoint_failure", .GlobalEnv, ifnotfound = "none")
    }), silent = TRUE)
    cat("SYNTHETIC LAST ENDPOINT CONDITION (may be expected polling)",
        owner, as.character(diagnostic), "\n")
  }
}
run_k <- function(k) {
  case_started <- Sys.time()
  state <- tempfile(paste0("smoke-k", k, "-"), tmpdir = out)
  dir.create(state, mode = "0700", recursive = TRUE)
  peers <- list()
  on.exit({ we_close_peers(peers); unlink(state, recursive = TRUE) }, add = TRUE)
  owners <- paste0("site_", letters[seq_len(k)])
  refs <- c("site_a$a", paste0(owners[-1L], "$x"))
  physical <- sub("^[^$]+[$]", "", refs)
  duplicate <- duplicated(physical) | duplicated(physical, fromLast = TRUE)
  registered <- ifelse(duplicate, refs, physical)
  # One balanced variable per owner keeps this a bounded lifecycle smoke.
  # Mixed same/cross-owner pairs and larger workloads have separate unit coverage.
  row_setting <- Sys.getenv("DSVERT_SMOKE_ROWS", "")
  n <- if (nzchar(row_setting)) as.integer(row_setting) else if (k <= 3L) 8L else 16L
  stopifnot(!is.na(n), n >= (if (k <= 3L) 4L else 16L),
            n <= 128L, n %% 4L == 0L)
  specs <- stats::setNames(lapply(owners, function(owner) {
    path <- file.path(state, owner)
    dir.create(path, mode = "0700")
    list(name = owner, state_dir = path, dslite_home = file.path(path, "dslite"),
         mpc_binary = Sys.getenv("DS_VERT_MPC"), dataset_id = paste0("cor142-k", k),
         dataset_version = "v1")
  }), owners)
  raw <- stats::setNames(lapply(seq_along(owners), function(i) {
    data <- data.frame(patient_id = sprintf("u%03d", seq_len(n)))
    mask <- (if (k == 3L) c(1L, 2L, 3L) else c(1L, 2L, 4L, 8L, 3L))[[i]]
    value <- vapply(seq_len(n) - 1L, function(row)
      sum(as.integer(intToBits(bitwAnd(row, mask)))) %% 2, numeric(1L))
    if (i == 1L) data$a <- value else data$x <- value
    data
  }), owners)
  policies <- stats::setNames(lapply(owners, function(owner) {
    columns <- setdiff(names(raw[[owner]]), "patient_id")
    list(total_epsilon = 8, total_delta = 0,
         domain = paste0("cor142-k", k), cohort_id = "synthetic-v1",
         adjacency = "add_remove_patient", patient_column = "patient_id",
         unit_capacity = n, numeric_grid_bits = 8L, max_records_per_unit = 1L,
         overflow_policy = "reject_snapshot", designated_noise_peers = owners[1:2],
         datasets = list(DA = list(id = specs[[owner]]$dataset_id, version = "v1")),
         workload_scope = list(mode = "catalog_v1", numeric_moments = character(),
           categorical_marginals = character(), categorical_pairs = list(),
           correlations = list()),
         numeric_bounds = stats::setNames(rep(list(c(0, 1)), length(columns)), columns),
         capsule_dataset_mapping = list(DA = columns),
         gaussian_specs = if (owner == owners[[1L]]) list(cor_primary = list(
           version = "v2", dataset = "DA", outcome = registered[[1L]],
           predictors = registered[-1L], intercept = TRUE)) else list(),
         synopsis_state_path = file.path(specs[[owner]]$state_dir, "synopsis"))
  }), owners)
  for (owner in owners) peers[[owner]] <- we_boot_peer(
    specs[[owner]], server, raw = raw[[owner]], policy = policies[[owner]])
  identities <- ds.getIdentityPks(we_connections(peers))
  pins <- stats::setNames(unname(unlist(identities)), owners)
  we_close_peers(peers); peers <- list()
  for (owner in owners) peers[[owner]] <- we_boot_peer(
    specs[[owner]], server, pins = pins[setdiff(owners, owner)],
    raw = raw[[owner]], policy = policies[[owner]])
  observe_peers(peers, n)
  on.exit(report_peer_errors(peers), add = TRUE, after = FALSE)
  conns <- we_connections(peers)
  cat("PEERS", k, "PIDS", vapply(peers, function(p) p$worker$get_pid(), integer(1L)), "\n")
  ds.psiAlign("D", "patient_id", "DA", datasources = conns, verbose = FALSE)
  record_event(paste0("K", k, " setup and PSI"), "passed",
    as.numeric(difftime(Sys.time(), case_started, units = "secs")))
  cor <- run_analysis(paste0("K", k, " Cor"), ds.vertCor(
    "DA", refs, analysis_id = "cor_primary", datasources = conns, verbose = FALSE))
  saveRDS(cor, file.path(out, paste0("smoke-k", k, "-cor.rds")))
  pca <- run_analysis(paste0("K", k, " PCA"), {
    value <- ds.vertPCA(cor_result = cor, verbose = FALSE)
    stopifnot(is.null(value$scores), nrow(value$loadings) == length(refs))
    value
  })
  if (k == 2L) {
    replay <- run_analysis("K2 sticky Cor replay", {
      value <- ds.vertCor("DA", refs, analysis_id = "cor_primary",
        datasources = conns, verbose = FALSE)
      saveRDS(value, file.path(out, "smoke-k2-cor-replay.rds"))
      stopifnot(identical(serialize(cor$correlation, NULL), serialize(value$correlation, NULL)),
        identical(serialize(cor$signed_artifact, NULL), serialize(value$signed_artifact, NULL)),
        identical(serialize(cor$provenance_certificate, NULL),
          serialize(value$provenance_certificate, NULL)))
      for (field in c("release_set_json", "replay_set_json")) {
        original <- cor$provenance_certificate$signed_evidence[[field]]
        stopifnot(!is.null(original), identical(original,
          value$provenance_certificate$signed_evidence[[field]]))
      }
      value
    })
    cat("REPLAY_IDENTICAL artifact certificate release correlation\n")
    offline <- run_analysis("K2 offline Cor to PCA zero DSI", {
      audit <- new.env(parent = emptyenv()); audit$calls <- 0L
      assign(".vertcor142_offline_audit", audit, .GlobalEnv)
      functions <- c("datashield.aggregate", "datashield.assign", "datashield.connections_find")
      for (fn in functions) trace(fn, where = asNamespace("DSI"),
        tracer = quote({ .GlobalEnv$.vertcor142_offline_audit$calls <-
          .GlobalEnv$.vertcor142_offline_audit$calls + 1L }), print = FALSE)
      converted <- tryCatch({
        offline_cor <- ds.vertCor(variables = refs,
          gaussian_result = cor$provenance_certificate, trusted_pinset = pins, verbose = FALSE)
        offline_pca <- ds.vertPCA(cor_result = offline_cor, trusted_pinset = pins, verbose = FALSE)
        stopifnot(identical(offline_cor$correlation, cor$correlation),
          identical(offline_pca$loadings, pca$loadings), is.null(offline_pca$scores),
          identical(audit$calls, 0L))
        list(correlation = offline_cor, pca = offline_pca, dsi_calls = audit$calls)
      }, finally = {
        for (fn in functions) untrace(fn, where = asNamespace("DSI"))
        rm(".vertcor142_offline_audit", envir = .GlobalEnv)
      })
      converted
    })
    saveRDS(offline, file.path(out, "smoke-k2-offline.rds"))
  }
  cat("COMPLETE K", k, "elapsed", as.numeric(difftime(Sys.time(), started, units = "secs")), "\n")
}
cases <- strsplit(Sys.getenv("DSVERT_SMOKE_CASES", "cor2,cor3"), ",", fixed = TRUE)[[1L]]
stopifnot(length(cases) > 0L, all(cases %in% c("cor2", "cor3", "cor5", "binomial", "poisson", "lasso")))
source(file.path(client, "tools/validate_cross_owner_142_lasso.R"))
status <- tryCatch({
  for (case in cases) {
    if (startsWith(case, "cor")) run_k(as.integer(sub("cor", "", case)))
    else if (identical(case, "lasso"))
      vertcor142_lasso_smoke(c("binomial", "poisson"), root, run_analysis, observe_peers)
    else vertcor142_lasso_smoke(case, root, run_analysis, observe_peers)
  }
  "complete"
}, error = function(e) paste("incomplete:", conditionMessage(e)))
cat("SMOKE", status, "invocations", invocations, "elapsed_seconds",
    as.numeric(difftime(Sys.time(), started, units = "secs")), "\n")
if (!identical(status, "complete")) quit(status = 1L)
