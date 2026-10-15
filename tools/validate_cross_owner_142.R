# Short process-isolated DSLite smoke. Run under an external timeout of 900 s.
# Synthetic data only. At most ten analysis/post-processing invocations.
args <- commandArgs(TRUE)
root <- if (length(args)) normalizePath(args[[1L]]) else normalizePath("..")
client <- file.path(root, "dsVertClient")
server <- file.path(root, "dsVert")
out <- file.path(root, ".tmp", "vertcor")
source(file.path(client, "inst/validation/v1.3.0/worked_example_custodian.R"))
pkgload::load_all(client, quiet = TRUE)
started <- Sys.time()
prior_seconds <- as.numeric(Sys.getenv("DSVERT_SMOKE_PRIOR_SECONDS", "0"))
prior_analyses <- as.integer(Sys.getenv("DSVERT_SMOKE_PRIOR_ANALYSES", "0"))
stopifnot(is.finite(prior_seconds), prior_seconds >= 0, prior_seconds < 900,
          !is.na(prior_analyses), prior_analyses >= 0L, prior_analyses < 10L)
invocations <- prior_analyses
run_analysis <- function(label, code) {
  stopifnot(invocations < 10L,
            prior_seconds + as.numeric(difftime(Sys.time(), started, units = "secs")) < 880)
  invocations <<- invocations + 1L
  cat("ANALYSIS", invocations, label, "\n"); flush.console()
  value <- withCallingHandlers(force(code), error = function(e) {
    if (grepl("invalid synopsis schema binding", conditionMessage(e), fixed = TRUE)) {
      candidates <- Filter(function(frame) exists("schema_raw", frame, inherits = FALSE), sys.frames())
      if (length(candidates)) {
        frame <- tail(candidates, 1L)[[1L]]
        saveRDS(list(schema = get("schema", frame),
          schema_raw = get("schema_raw", frame),
          bootstrap_set_sha256 = get("bootstrap_set_sha256", frame),
          context_policy = get("context", frame)$policy),
          file.path(out, "smoke-schema-binding-diagnostic.rds"))
      }
    }
  })
  cat("PASS", label, "\n"); flush.console()
  value
}
run_k <- function(k) {
  state <- tempfile(paste0("smoke-k", k, "-"), tmpdir = out)
  dir.create(state, mode = "0700", recursive = TRUE)
  peers <- list()
  on.exit({ we_close_peers(peers); unlink(state, recursive = TRUE) }, add = TRUE)
  owners <- paste0("site_", letters[seq_len(k)])
  refs <- c("site_a$a", "site_a$b", paste0(owners[-1L], "$x"))
  physical <- sub("^[^$]+[$]", "", refs)
  duplicate <- duplicated(physical) | duplicated(physical, fromLast = TRUE)
  registered <- ifelse(duplicate, refs, physical)
  n <- 128L
  specs <- stats::setNames(lapply(owners, function(owner) {
    path <- file.path(state, owner)
    dir.create(path, mode = "0700")
    list(name = owner, state_dir = path, dslite_home = file.path(path, "dslite"),
         mpc_binary = Sys.getenv("DS_VERT_MPC"), dataset_id = paste0("cor142-k", k),
         dataset_version = "v1")
  }), owners)
  raw <- stats::setNames(lapply(seq_along(owners), function(i) {
    data <- data.frame(patient_id = sprintf("u%03d", seq_len(n)))
    if (i == 1L) {
      data$a <- rep(c(0, 1), length.out = n)
      data$b <- rep(c(0, 0, 1, 1), length.out = n)
    } else data$x <- ((seq_len(n) * (i + 1L)) %% 11L) / 10
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
  for (owner in owners) peers[[owner]]$worker$run(function() {
    trace(".dsvert_dp_transcript_stop", where = asNamespace("dsVert"),
      tracer = quote(assign(".synthetic_failure", conditionMessage(error), .GlobalEnv)),
      print = FALSE)
    invisible(NULL)
  })
  on.exit({
    for (owner in names(peers)) {
      diagnostic <- try(peers[[owner]]$worker$run(function() {
        get0(".synthetic_failure", .GlobalEnv, ifnotfound = "none")
      }), silent = TRUE)
      cat("SYNTHETIC CUSTODIAN DIAGNOSTIC", owner, as.character(diagnostic), "\n")
    }
  }, add = TRUE, after = FALSE)
  conns <- we_connections(peers)
  cat("PEERS", k, "PIDS", vapply(peers, function(p) p$worker$get_pid(), integer(1L)), "\n")
  ds.psiAlign("D", "patient_id", "DA", datasources = conns, verbose = FALSE)
  cor <- run_analysis(paste0("K", k, " Cor"), ds.vertCor(
    "DA", refs, analysis_id = "cor_primary", datasources = conns, verbose = FALSE))
  saveRDS(cor, file.path(out, paste0("smoke-k", k, "-cor.rds")))
  pca <- run_analysis(paste0("K", k, " PCA"), ds.vertPCA(cor_result = cor, verbose = FALSE))
  stopifnot(is.null(pca$scores), nrow(pca$loadings) == length(refs))
  if (k == 2L && prior_analyses == 0L) {
    replay <- run_analysis("K2 sticky Cor replay", ds.vertCor(
      "DA", refs, analysis_id = "cor_primary", datasources = conns, verbose = FALSE))
    stopifnot(identical(cor$correlation, replay$correlation))
  }
  cat("COMPLETE K", k, "elapsed", as.numeric(difftime(Sys.time(), started, units = "secs")), "\n")
}
topologies <- if (prior_analyses >= 5L) c(2L, 3L) else c(2L, 3L, 5L)
status <- tryCatch({ for (k in topologies) run_k(k); "complete" },
                   error = function(e) paste("incomplete:", conditionMessage(e)))
cat("SMOKE", status, "cumulative_invocations", invocations, "elapsed_seconds",
    as.numeric(difftime(Sys.time(), started, units = "secs")), "\n")
if (!identical(status, "complete")) quit(status = 1L)
