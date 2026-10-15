# One bounded real B3 run. Launch only through v3-dsVert/tools/run_cross_owner_smoke.py.
# Two isolated R/DSLite custodians; source custody and all secret keys stay local.
main <- function() {
args <- commandArgs(TRUE)
root <- normalizePath(if (length(args)) args[[1L]] else "..")
client <- file.path(root, "v3-dsVertClient")
server <- file.path(root, "v3-dsVert")
out <- file.path(root, ".tmp", "v3")
test_rows <- Sys.getenv("DSVERT_B3_TEST_ROWS", "")
stopifnot(test_rows %in% c("", "4"))
is_integration_test <- identical(test_rows, "4")
# Proof-of-concept mode: the admitted producer profile with real signal in the
# synthetic outcome, a 15-minute stage envelope and all three releases. It is
# not the B3 cost gate and records no gate verdict.
poc <- identical(Sys.getenv("DSVERT_B3_POC"), "1")
stopifnot(!(poc && is_integration_test))
n <- if (is_integration_test) 4L else 64L
stopifnot(normalizePath(Sys.getenv("DSVERT_SERVER_SOURCE")) == normalizePath(server),
          normalizePath(Sys.getenv("DS_VERT_MPC")) == normalizePath(file.path(out, "bin", "dsvert-mpc")))
campaign <- Sys.getenv("DSVERT_SMOKE_CAMPAIGN", "v3-b3")
stopifnot(grepl("^[A-Za-z0-9][A-Za-z0-9._-]{0,79}$", campaign))
source(file.path(client, "inst/validation/v1.3.0/worked_example_custodian.R"))
# Preserve the genuine allowlisted endpoint condition inside its custodian;
# this observer neither returns private values nor substitutes endpoint results.
assign(".we_worker_submit", function(worker, operation, args, async) {
  if (!identical(worker$get_state(), "idle")) stop("The isolated DSLite peer already has a pending command")
  job <- switch(operation,
    aggregate = function(expr) withCallingHandlers(
      get(".dsvert_e2e_server", .GlobalEnv)$aggregate(get(".dsvert_e2e_sid", .GlobalEnv), expr),
      error = function(e) assign(".irls_b3_endpoint_error", list(operation = deparse(expr[[1L]]),
        message = conditionMessage(e)), .GlobalEnv)),
    assign = function(symbol, expr) {
      get(".dsvert_e2e_server", .GlobalEnv)$assignExpr(get(".dsvert_e2e_sid", .GlobalEnv), symbol, expr)
      NULL
    })
  if (isTRUE(async)) {worker$call(job, args = args); return(NULL)}
  worker$run(job, args = args)
}, .GlobalEnv)
pkgload::load_all(client, quiet = TRUE)
cf <- function(name) get(name, asNamespace("dsVertClient"))
seconds <- function(start) as.numeric(difftime(Sys.time(), start, units = "secs"))
evidence_path <- file.path(out, paste0(campaign, "-evidence.json"))
evidence <- list(version = "dsvert-v3-irls-b3-gate-v1", profile = if (is_integration_test) "dsvert-dp-irls-producer-p3-q16-q12-v1" else "irls-smoke-k2-c64-p3-v1",
  integration_test_not_gate = is_integration_test, n = n, p = 3L, family = "binomial", stage_limit_seconds = 240,
  path_limit_seconds = 840, compile_in_stage_gate = TRUE,
  binary_sha256 = digest::digest(file = Sys.getenv("DS_VERT_MPC"), algo = "sha256"),
  setup = list(), stages = list(), gate = "not-run")
write_evidence <- function() {
  jsonlite::write_json(evidence, evidence_path, auto_unbox = TRUE, pretty = TRUE,
                        null = "null", digits = 15)
}
write_evidence()
peers <- list()
state <- tempfile(paste0(campaign, "-custodians-"), tmpdir = out)
dir.create(state, recursive = TRUE, mode = "0700")
# Retain private authority stores for crash/cold replay inspection; the launcher
# owns process cleanup, and this handler closes the ordinary process handles.
on.exit(we_close_peers(peers), add = TRUE)
setup_started <- Sys.time()
owners <- c("peer_a", "peer_b")
ids <- sprintf("person-%03d", seq_len(n))
x <- rep(c(-1, -.5, .5, 1), length.out = n)
z <- rep(c(-.75, .75, .25, -.25), each = max(1L, n %/% 4L), length.out = n)
y <- as.integer((seq_len(n) * 17L %% 23L) < 12L)
if (poc) {
  set.seed(20261015)
  x <- round(stats::runif(n, -1, 1), 3); z <- round(stats::runif(n, -1, 1), 3)
  y <- stats::rbinom(n, 1L, stats::plogis(-0.25 + 1.2 * x - 0.9 * z))
}
raw <- list(peer_a = data.frame(patient_id = ids, y = y),
            peer_b = data.frame(patient_id = rev(ids), x = rev(x), z = rev(z)))
specs <- stats::setNames(lapply(owners, function(owner) {
  path <- file.path(state, owner); dir.create(path, mode = "0700")
  list(name = owner, state_dir = path, dslite_home = file.path(path, "dslite"),
       mpc_binary = Sys.getenv("DS_VERT_MPC"), dataset_id = "irls-b3-synthetic",
       dataset_version = "v1")
}), owners)
for (owner in owners) peers[[owner]] <- we_boot_peer(specs[[owner]], server, raw = raw[[owner]])
# Administrative identity discovery is anchored by a second direct read inside
# each custodian process before the name-bound pinset is installed.
identities <- ds.getIdentityPks(we_connections(peers))
pins <- vapply(peers, function(peer) peer$worker$run(function()
  dsVert:::.dsvert_relay_normalize_identity_pk(dsVert:::.get_identity_keypair()$identity_pk)), character(1L))
stopifnot(identical(unname(unlist(identities)), unname(pins)))
we_close_peers(peers); peers <- list()
for (owner in owners) peers[[owner]] <- we_boot_peer(specs[[owner]], server,
  pins = pins[setdiff(owners, owner)], raw = raw[[owner]])
conns <- we_connections(peers)
if (poc) {
  # Custodian transport policy for this DSLite run only: 4 MiB exact-GC frames
  # (the package default is 480 KiB; the admitted range is 16 KiB..8 MiB).
  for (peer in peers) peer$worker$run(function() options(dsvert.exact_gc.chunk_bytes = 4L * 1024L^2))
}
options(dsvert.dp_irls.trusted_pins = pins)
evidence$pids <- vapply(peers, function(peer) peer$worker$get_pid(), integer(1L))
stopifnot(length(unique(evidence$pids)) == 2L, !Sys.getpid() %in% evidence$pids)
evidence$setup$boot_and_pin_seconds <- seconds(setup_started)
psi_started <- Sys.time()
attestation <- ds.psiAlign("D", "patient_id", "DA", datasources = conns, verbose = FALSE)
evidence$setup$psi_seconds <- seconds(psi_started)
evidence$psi_attestation <- attestation
evidence$pins <- as.list(pins)
write_evidence()

provision_started <- Sys.time()
contract <- peers[[1L]]$worker$run(function(server, pins, attestation, n, poc) {
  helper <- new.env(parent = asNamespace("dsVert"))
  sys.source(file.path(server, "tests/testthat/helper-dp-irls-contract.R"), helper)
  sys.source(file.path(server, "tests/testthat/helper-dp-irls-profile.R"), helper)
  spec <- if (n == 64L && !poc) helper$.irls_smoke_fixture()$spec else helper$.irls_producer_fixture("binomial", as.character(n))$spec
  if (poc) { spec$operations$max_stage_seconds <- "900"; spec$operations$max_path_seconds <- "2700" }
  spec$privacy_source$peer_pins <- as.list(pins)
  spec$privacy_source$pinset_sha256 <- dsVert:::.dsvert_dp_synopsis_pinset_hash_v1(pins)
  spec$privacy_source$alignment_purpose <- attestation$alignment_purpose
  source_names <- c(peer_a = "outcomes", peer_b = "predictors")
  deps <- lapply(names(pins), function(peer) {
    epoch <- Filter(function(e) e$peer_identity == pins[[peer]], attestation$authorization_epochs)[[1L]]
    list(source_identity = unname(pins[[peer]]), dataset = unname(source_names[[peer]]),
      dependency = list(version = dsVert:::.DSVERT_DP_SYNOPSIS_DEPENDENCY_VERSION,
        dataset_id = attestation$dataset_id, dataset_version = attestation$dataset_version,
        authorization_epoch = epoch$authorization_epoch, source_contract_id = epoch$source_contract_id,
        derivation_policy_id = attestation$derivation_policy_id, provenance_kind = "padded_psi",
        authorization_epochs = attestation$authorization_epochs))
  })
  deps <- deps[order(vapply(deps, function(d) paste(d$source_identity, d$dataset, sep = "|"),
                            character(1L)), method = "radix")]
  spec$privacy_source$closure$dependencies <- deps
  spec$privacy_source$source_ids <- lapply(deps, `[[`, "dataset")
  spec$privacy_source$schema_sha256 <- dsVert:::.dsvert_dp_irls_hash(
    "dsvert-v3-irls-b3-public-schema-v1", list(columns = spec$estimand$columns,
      source_abi = dsVert:::.DSVERT_DP_IRLS_SOURCE_ABI,
      conversion = dsVert:::.dsvert_dp_irls_conversion_certificate(),
      ordering_rule = spec$privacy_source$ordering_rule,
      attestation = attestation))
  dsVert:::.dsvert_dp_irls_contract_normalize(spec)
}, args = list(server, pins, attestation, n, poc))
unsigned <- list(version = "dsvert-dp-irls-contract-v1", contract = contract)
# Use server's actual version rather than assume fixture envelope spelling.
unsigned$version <- peers[[1L]]$worker$run(function() dsVert:::.DSVERT_DP_IRLS_CONTRACT_VERSION)
signed <- c(unsigned, list(signatures = lapply(peers, function(peer)
  peer$worker$run(function(unsigned) {
    dsVert:::.dsvert_relay_sign_message(dsVert:::.dsvert_dp_irls_contract_message(unsigned),
                                       dsVert:::.get_identity_keypair()$identity_sk)
  }, args = list(unsigned)))))
transport <- lapply(peers, function(peer) peer$worker$run(function() {
  keys <- dsVert:::.callMpcTool("transport-keygen", list())
  assign(".irls_b3_transport", keys, .GlobalEnv)
  list(public_key = keys$public_key,
    signature = dsVert:::.sign_transport_pk(keys$public_key, dsVert:::.get_identity_keypair()$identity_sk))
}))
for (owner in owners) {
  other <- setdiff(owners, owner)
  dir.create(file.path(specs[[owner]]$state_dir, "irls"), mode = "0700")
  peers[[owner]]$worker$run(function(owner, signed, pins, directory, remote) {
    server <- get(".dsvert_e2e_server", .GlobalEnv); sid <- get(".dsvert_e2e_sid", .GlobalEnv)
    aligned <- server$getSessionData(sid, "DA")
    source <- if (owner == "peer_a") "outcomes" else "predictors"
    dsVert:::.dsvert_dp_irls_provision("irls_b3_smoke", signed, pins, owner,
      dsVert:::.get_identity_keypair()$identity_sk,
      sources = stats::setNames(list(aligned), source), directory = directory,
      transport = get(".irls_b3_transport", .GlobalEnv),
      peer_transport_pk = remote$public_key, peer_transport_signature = remote$signature)
    TRUE
  }, args = list(owner, signed, pins, file.path(specs[[owner]]$state_dir, "irls"), transport[[other]]))
}
evidence$setup$provision_seconds <- seconds(provision_started)
evidence$setup$total_seconds <- seconds(setup_started)
evidence$signed_contract <- signed
write_evidence()
cat("B3 SETUP", jsonlite::toJSON(evidence$setup, auto_unbox = TRUE), "\n"); flush.console()

# These failures are tested on the real allowlisted R/DSLite path. Extra
# arguments are rejected before materialization or native worker startup.
injections <- c("beta", "predecessor")
for (argument in injections) {
  rejected <- try(cf(".dsvert_aggregate_strict")(conns = conns[1L],
    expr = as.call(c(list(as.name("dpIRLSStartDS"), analysis_id = "irls_b3_smoke"),
                stats::setNames(list("injected"), argument))), operation = "IRLS injection rejection"), silent = TRUE)
  stopifnot(inherits(rejected, "try-error"))
}
evidence$dsi_injection_rejections <- as.list(injections)
if (identical(Sys.getenv("DSVERT_B3_SETUP_ONLY"), "true")) {
  evidence$gate <- "setup-only"; write_evidence(); return(invisible(NULL))
}
# Measure encrypted transport roundtrips directly. They overlap native MPC
# spans because native workers wait for these same frames; never add them.
transport_clock <- new.env(parent = emptyenv()); transport_clock$seconds <- 0
transport_clock$calls <- 0L
trace(".dsvert_aggregate_strict", where = asNamespace("dsVertClient"),
  tracer = quote({.b3_transport_started <- proc.time()[["elapsed"]]}),
  exit = substitute({
    if (identical(operation, "DP IRLS encrypted transport")) {
      # CLOCK is substituted as an environment value, so it cannot be an
      # assignment target; update its bindings explicitly.
      assign("seconds", get("seconds", envir = CLOCK) +
        proc.time()[["elapsed"]] - .b3_transport_started, envir = CLOCK)
      assign("calls", get("calls", envir = CLOCK) + 1L, envir = CLOCK)
    }
  }, list(CLOCK = transport_clock)), print = FALSE)
on.exit(untrace(".dsvert_aggregate_strict", where = asNamespace("dsVertClient")), add = TRUE)

path_started <- Sys.time()
resolved <- tryCatch(cf(".dsvert_dp_irls_client_metadata")("irls_b3_smoke", conns, as.list(pins)),
                     error = identity)
if (inherits(resolved, "error")) {
  evidence$gate <- "metadata-failed-before-stage"
  evidence$metadata_error <- conditionMessage(resolved)
  evidence$endpoint_errors <- lapply(peers, function(peer) peer$worker$run(function()
    get0(".irls_b3_endpoint_error", .GlobalEnv, ifnotfound = NULL)))
  write_evidence()
  cat("B3 METADATA ERRORS", jsonlite::toJSON(evidence$endpoint_errors, auto_unbox = TRUE), "\n")
  stop(resolved)
}
metadata <- resolved$metadata
for (index in if (is_integration_test) 0L else 0:2) {
  stage_started <- Sys.time()
  prior_transport <- transport_clock$seconds; prior_calls <- transport_clock$calls
  publication <- tryCatch(cf(".dsvert_dp_irls_client_execute_stage")("irls_b3_smoke", conns, metadata, as.list(pins)),
                          error = identity)
  stage_wall <- seconds(stage_started)
  metrics <- lapply(peers, function(peer) tryCatch(peer$worker$run(function(index) {
    dsVert:::.dsvert_dp_irls_timings("irls_b3_smoke")[[as.character(index)]]
  }, args = list(index)), error = function(e) list(error = conditionMessage(e))))
  if (inherits(publication, "error")) {
    evidence$stages[[index + 1L]] <- list(index = index, wall_seconds = stage_wall,
      authority_timings = metrics,
      custodian_diagnostics = lapply(peers, function(peer) tryCatch(peer$worker$run(function()
        dsVert:::.dsvert_dp_irls_diagnostics("irls_b3_smoke")),
        error = function(e) list(error = conditionMessage(e)))),
      dsi_transport_roundtrip_seconds = transport_clock$seconds - prior_transport,
      dsi_transport_calls = transport_clock$calls - prior_calls, pass = FALSE,
      error = conditionMessage(publication),
      endpoint_errors = lapply(peers, function(peer) tryCatch(peer$worker$run(function()
        get0(".irls_b3_endpoint_error", .GlobalEnv, ifnotfound = NULL)),
        error = function(e) list(error = conditionMessage(e)))))
    evidence$gate <- if (is_integration_test) "integration-failed-not-cost-gate" else "stage-fail-stop"; evidence$path_seconds <- seconds(path_started)
    write_evidence(); stop(publication)
  }
  components <- c("compile", "source", "protected_compute", "noise", "publication")
  # Both native roles overlap. Component maxima and measured DSI roundtrips
  # overlap while native workers wait for frames; their sum is not stage wall.
  breakdown <- stats::setNames(lapply(components, function(name)
    max(vapply(metrics, function(m) as.numeric(m[[name]]), numeric(1L)))), components)
  breakdown$dsi_transport_roundtrips <- transport_clock$seconds - prior_transport
  breakdown$spans_overlap <- TRUE
  evidence$stages[[index + 1L]] <- list(index = index, wall_seconds = stage_wall,
    component_seconds = breakdown, authority_timings = metrics,
    dsi_transport_calls = transport_clock$calls - prior_calls,
    pass = stage_wall <= 240, published_stage_id = publication$stages[[index + 1L]]$body$stage_id)
  metadata$stages <- publication$stages
  metadata$resume_required <- NULL
  evidence$path_seconds <- seconds(path_started)
  evidence$gate <- if (poc) "poc-running" else if (stage_wall <= 240) "stage-pass" else "stage-fail-stop"
  write_evidence()
  cat("B3 STAGE", jsonlite::toJSON(evidence$stages[[index + 1L]], auto_unbox = TRUE), "\n"); flush.console()
  if (!poc && stage_wall > 240) stop("B3 stage cost gate failed; stop for maintainer decision.")
  if (!poc && evidence$path_seconds > 840) {
    evidence$gate <- "path-fail-stop"; write_evidence()
    stop("B3 three-release path exceeded 840 seconds; stop for maintainer decision.")
  }
}
if (is_integration_test) {
  evidence$gate <- "integration-pass-not-cost-gate"; write_evidence()
  cat("B3 INTEGRATION PASS rows", n, "stage_seconds", evidence$stages[[1L]]$wall_seconds, "\n")
  return(invisible(NULL))
}
evidence$gate <- if (poc) "poc-complete" else if (evidence$path_seconds <= 840) "pass" else "path-fail-stop"
if (poc) {
  # Everything the analyst received, plus the synthetic source table that only
  # this offline report uses for its centralised glm() comparison.
  saveRDS(list(signed = signed, pins = as.list(pins), root = metadata$root,
    stages = metadata$stages, data = data.frame(y = y, x = x, z = z),
    evidence = evidence), file.path(out, "poc-run.rds"))
}
# Repeated exported client calls read the retained chain and do no fresh work.
first <- ds.vertDPIRLS("irls_b3_smoke", datasources = conns, verbose = FALSE)
second <- ds.vertDPIRLS("irls_b3_smoke", datasources = conns, verbose = FALSE)
stopifnot(identical(serialize(first$stages, NULL), serialize(second$stages, NULL)))
evidence$replay_byte_identical <- TRUE
write_evidence()
cat("B3 GATE", evidence$gate, "path_seconds", evidence$path_seconds, "evidence", evidence_path, "\n")
if (!evidence$gate %in% c("pass", "poc-complete")) quit(status = 1L)
}
tryCatch(main(), error = function(e) {
  message("B3 stopped: ", conditionMessage(e)); quit(status = 1L)
})
