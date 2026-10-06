test_that("mixed provenance protocols stop before any private endpoint", {
  sources <- list(alpha = list(), beta = list())
  calls <- character()
  aggregate <- function(conns, expr, ...) {
    method <- as.character(expr[[1L]])
    calls <<- c(calls, method)
    if (!identical(method, "dsvertProvenanceCapabilitiesDS")) {
      stop("Private materialization was reached")
    }
    values <- lapply(names(conns), function(peer) list(
      protocol = if (peer == "alpha") "dsvert-provenance-epochs-v1" else
        "dsvert-provenance-epochs-v0",
      version = "1.4.1",
      families = c("psi", "count", "frequency", "synopsis")))
    stats::setNames(values, names(conns))
  }
  expect_error(.dsvert_provenance_preflight(sources, aggregate),
               "All peers require")
  expect_error(.dsvert_psi_padded_align(
    "D", "id", "DA", sources, .aggregate = aggregate,
    .assign = function(...) stop("private assignment"),
    .exact_run = function(...) stop("private execution")),
    "All peers require")
  expect_error(.dsvert_dp_count_compile_authorize_bind_v1(
    "DA", sources, .aggregate = aggregate), "All peers require")
  expect_error(.dsvert_dp_frequency_prepare_v1(
    "DA", "category", "alpha", sources, .aggregate = aggregate),
    "All peers require")
  expect_error(.dsvert_dp_synopsis_bootstrap_build_v1(
    sources, .aggregate = aggregate), "All peers require")
  expect_true(length(calls) >= 5L)
  expect_true(all(calls == "dsvertProvenanceCapabilitiesDS"))
})

test_that("epoch vectors require canonical complete independent authorities", {
  identities <- sort(vapply(1:3, function(index) {
    .dsvert_exact_gc_b64url_encode(as.raw(rep(index, 32L)))
  }, character(1L)), method = "radix")
  epochs <- lapply(seq_along(identities), function(index) list(
    peer_identity = identities[[index]],
    authorization_epoch = sprintf("%064x", index),
    source_contract_id = strrep("c", 64L)))
  expect_identical(.dsvert_provenance_epochs_validate(epochs, identities), epochs)
  expect_error(.dsvert_provenance_epochs_validate(rev(epochs), identities),
               "authorization epoch")
  expect_error(.dsvert_provenance_epochs_validate(epochs[-3L], identities),
               "authorization epoch")
  duplicate <- epochs
  duplicate[[3L]]$authorization_epoch <- duplicate[[1L]]$authorization_epoch
  expect_error(.dsvert_provenance_epochs_validate(duplicate, identities),
               "authorization epoch")
  malformed <- epochs
  malformed[[2L]]$content_digest <- strrep("d", 64L)
  expect_error(.dsvert_provenance_epochs_validate(malformed, identities),
               "authorization epoch")
})

.provenance_fixture_helpers <- function(file) {
  environment <- new.env(parent = asNamespace("dsVertClient"))
  for (expression in parse(testthat::test_path(file))) {
    if (is.call(expression) &&
        identical(as.character(expression[[1L]]), "test_that")) next
    eval(expression, envir = environment)
  }
  environment
}

test_that("Count stable identity binds every epoch and ignores new run tokens", {
  helpers <- .provenance_fixture_helpers("test-dp-count-analysis.R")
  fixture <- helpers$.count_client_fixture(3L)
  compile <- function(envelopes) .dsvert_dp_count_client_compile_v1(
    envelopes, names(envelopes))
  original <- compile(fixture$envelopes)
  update <- function(field = NULL, value = NULL, run = FALSE) {
    envelopes <- fixture$envelopes
    for (peer in names(envelopes)) {
      receipt <- envelopes[[peer]]$receipt
      if (!is.null(field)) receipt[[field]] <- value
      if (run) receipt$psi_run_sha256 <- strrep("f", 64L)
      receipt$snapshot_commitment <- .dsvert_provenance_snapshot_commitment(
        receipt, fixture$config, peer,
        list(family = "count", statistic = "aligned_privacy_unit_count"))
      envelopes[[peer]]$receipt <- helpers$.count_client_resign_receipt(
        receipt, fixture$keys$private[[peer]])
    }
    compile(envelopes)
  }
  replay <- update(run = TRUE)
  expect_identical(.dsvert_provenance_stable_receipt(replay$contract),
                   .dsvert_provenance_stable_receipt(original$contract))
  expect_false(identical(replay$receipt_set_sha256,
                         original$receipt_set_sha256))
  epochs <- original$receipts[[1L]]$authorization_epochs
  epochs[[3L]]$authorization_epoch <- strrep("e", 64L)
  rotated <- update("authorization_epochs", epochs)
  expect_false(identical(rotated$contract$artifact_key,
                         original$contract$artifact_key))
  expect_false(identical(rotated$worker_static_sha256,
                         original$worker_static_sha256))
  policy <- update("derivation_policy_id", strrep("a", 64L))
  expect_false(identical(policy$contract$artifact_key,
                         original$contract$artifact_key))
  stale <- fixture$envelopes
  for (peer in names(stale)) {
    stale[[peer]]$receipt$authorization_epochs <- epochs
    stale[[peer]]$receipt <- helpers$.count_client_resign_receipt(
      stale[[peer]]$receipt, fixture$keys$private[[peer]])
  }
  expect_error(compile(stale), "complete authorization provenance")
  mixed <- fixture$envelopes
  peer <- names(mixed)[[1L]]
  receipt <- mixed[[peer]]$receipt
  receipt$authorization_epochs <- epochs
  receipt$snapshot_commitment <- .dsvert_provenance_snapshot_commitment(
    receipt, fixture$config, peer,
    list(family = "count", statistic = "aligned_privacy_unit_count"))
  mixed[[peer]]$receipt <- helpers$.count_client_resign_receipt(
    receipt, fixture$keys$private[[peer]])
  expect_error(compile(mixed), "disagree on complete authorization provenance")
  stable <- .dsvert_provenance_stable_receipt(original$contract)
  stable$semantic$cohort_id <- "tampered"
  expect_error(.dsvert_provenance_stable_receipt(stable), "artifact binding")
})

test_that("Synopsis closure binds membership-only peers and every dataset", {
  identities <- sort(vapply(1:3, function(index) {
    .dsvert_exact_gc_b64url_encode(as.raw(rep(index, 32L)))
  }, character(1L)), method = "radix")
  peers <- c("alpha", "beta", "gamma")
  context <- list(servers = peers, pinset = stats::setNames(identities, peers))
  epochs <- lapply(seq_along(identities), function(index) list(
    peer_identity = identities[[index]],
    authorization_epoch = sprintf("%064x", index),
    source_contract_id = strrep("c", 64L)))
  dependency <- function(name) list(
    version = .DSVERT_CLIENT_SYNOPSIS_DEPENDENCY_VERSION,
    dataset_id = name, dataset_version = "v1",
    authorization_epoch = strrep("a", 64L),
    source_contract_id = strrep("b", 64L),
    derivation_policy_id = strrep("d", 64L),
    provenance_kind = "padded_psi", authorization_epochs = epochs)
  bootstraps <- stats::setNames(lapply(seq_along(peers), function(index) list(
    peer_identity_pk = identities[[index]],
    draft = list(datasets = list(primary = list(dataset_id = "primary",
      dataset_version = "v1"))),
    authorization_dependencies = list(primary = dependency("primary"),
      secondary = dependency("secondary")))), peers)
  closure <- .dsvert_dp_synopsis_client_closure_v2(bootstraps, context)
  expect_length(closure$dependencies, 6L)
  changed <- bootstraps
  for (peer in peers) changed[[peer]]$authorization_dependencies$secondary$
    authorization_epochs[[3L]]$authorization_epoch <- strrep("e", 64L)
  expect_false(identical(closure,
    .dsvert_dp_synopsis_client_closure_v2(changed, context)))
  incomplete <- bootstraps
  incomplete[[1L]]$authorization_dependencies$primary$authorization_epochs <-
    epochs[-3L]
  expect_error(.dsvert_dp_synopsis_client_closure_v2(incomplete, context),
               "authorization epoch")
  missing <- bootstraps
  missing[[1L]]$authorization_dependencies$primary <- NULL
  expect_error(.dsvert_dp_synopsis_client_closure_v2(missing, context),
               "dependency coverage")
  root <- bootstraps
  root[[1L]]$authorization_dependencies$primary$provenance_kind <-
    "frozen_materialized_root"
  expect_error(.dsvert_dp_synopsis_client_closure_v2(root, context),
               "authorization dependency")
})

test_that("Frequency binds non-source epoch rotations and retains stable run replay", {
  helpers <- .provenance_fixture_helpers("test-dp-frequency-analysis.R")
  fixture <- helpers$.frequency_client_fixture(3L, "convolution")
  sign <- function(value, peer, domain) {
    unsigned <- value[setdiff(names(value), "signature")]
    value$signature <- helpers$.frequency_client_sign(charToRaw(paste0(
      domain, .dsvert_joint_dp_client_json(unsigned))),
      fixture$keys$private[[peer]])
    value
  }
  change <- function(rotate = FALSE, pin = TRUE) {
    claim <- fixture$claim
    if (rotate) {
      other <- which(vapply(claim$authorization_epochs, function(entry) {
        !identical(entry$peer_identity, claim$source_identity_pk)
      }, logical(1L)))[[1L]]
      claim$authorization_epochs[[other]]$authorization_epoch <- strrep("e", 64L)
    } else claim$psi_run_sha256 <- strrep("f", 64L)
    claim <- sign(claim, fixture$source_peer,
      "dsVert/dp-frequency/factor-claim/v3|")
    envelopes <- fixture$envelopes
    for (peer in names(envelopes)) {
      receipt <- envelopes[[peer]]$receipt
      receipt$authorization_epochs <- claim$authorization_epochs
      receipt$psi_run_sha256 <- claim$psi_run_sha256
      receipt$source_claim_sha256 <- .dsvert_dp_frequency_client_claim_hash_v1(
        claim, fixture$config)
      if (pin) receipt$snapshot_commitment <-
        .dsvert_provenance_snapshot_commitment(claim, fixture$config, peer,
          list(family = "frequency", factor_domain = fixture$config$factor_domain,
            source_owner = fixture$config$source_owner,
            missingness_policy = fixture$config$missingness_policy,
            repeated_record_policy = fixture$config$repeated_record_policy))
      envelopes[[peer]]$receipt <- sign(receipt, peer,
        "dsVert/dp-frequency/receipt/v2|")
    }
    .dsvert_dp_frequency_client_compile_v1(envelopes, names(envelopes), claim)
  }
  original <- fixture$compiled
  replay <- change()
  expect_identical(.dsvert_provenance_stable_receipt(replay$contract),
                   .dsvert_provenance_stable_receipt(original$contract))
  expect_false(identical(replay$receipt_set_sha256,
                         original$receipt_set_sha256))
  rotated <- change(rotate = TRUE)
  expect_false(identical(rotated$contract$artifact_key,
                         original$contract$artifact_key))
  expect_false(identical(rotated$worker_static_sha256,
                         original$worker_static_sha256))
  expect_error(change(rotate = TRUE, pin = FALSE),
               "complete authorization provenance")
})

test_that("Synopsis source Claims authenticate provenance-only commitments", {
  key <- openssl::read_ed25519_key(as.raw(seq.int(0L, 31L)))
  identity <- .dsvert_exact_gc_b64url_encode(as.list(key)$pubkey$data)
  context <- list(pinset = c(alpha = identity))
  projection <- list(sha256 = strrep("a", 64L),
    catalog = list(authorization_dependencies = list()))
  claim <- list(
    version = "dsvert-stateless-catalog-synopsis-source-vector-claim-v2",
    source_peer_name = "alpha", source_identity_pk = identity,
    catalog_sha256 = projection$sha256,
    source_vector_commitment = .dsvert_dp_synopsis_client_hash(
      "dsVert/synopsis/public-source/v2|", list(
        version = "dsvert-synopsis-provenance-source-commitment-v1",
        catalog_sha256 = projection$sha256, source_identity = identity)))
  sign <- function(value) {
    unsigned <- value[setdiff(names(value), "signature")]
    value$signature <- .dsvert_exact_gc_b64url_encode(openssl::ed25519_sign(
      charToRaw(paste0("dsVert/stateless-catalog-synopsis/source-vector-claim/v2|",
        .dsvert_joint_dp_client_json(unsigned))), key))
    value
  }
  claim <- sign(claim)
  expect_identical(.dsvert_dp_synopsis_client_source_claim_v2(
    claim, projection, "alpha", context), claim)
  content_commitment <- claim
  content_commitment$source_vector_commitment <- strrep("b", 64L)
  expect_error(.dsvert_dp_synopsis_client_source_claim_v2(
    sign(content_commitment), projection, "alpha", context), "provenance-only")
  old <- claim
  old$version <- "dsvert-stateless-catalog-synopsis-source-vector-claim-v1"
  expect_error(.dsvert_dp_synopsis_client_source_claim_v2(
    sign(old), projection, "alpha", context), "signed Synopsis source Claim")
  claim$signature <- .dsvert_exact_gc_b64url_encode(as.raw(rep(0L, 64L)))
  expect_error(.dsvert_dp_synopsis_client_source_claim_v2(
    claim, projection, "alpha", context), "signature")
})
