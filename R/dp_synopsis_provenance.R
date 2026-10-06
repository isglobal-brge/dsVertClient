.DSVERT_CLIENT_SYNOPSIS_DEPENDENCY_VERSION <-
  "dsvert-synopsis-authorization-dependency-v1"
.DSVERT_CLIENT_SYNOPSIS_CLOSURE_VERSION <-
  "dsvert-synopsis-authorization-closure-v1"

.dsvert_dp_synopsis_client_dependency_v2 <- function(value) {
  fields <- c("version", "dataset_id", "dataset_version",
    "authorization_epoch", "source_contract_id", "derivation_policy_id",
    "provenance_kind", "authorization_epochs")
  fail <- function() stop("Invalid Synopsis authorization dependency.",
                          call. = FALSE)
  if (!.dsvert_dp_has_exact_names(value, fields) ||
      !identical(value$version, .DSVERT_CLIENT_SYNOPSIS_DEPENDENCY_VERSION) ||
      !is.character(value$provenance_kind) ||
      length(value$provenance_kind) != 1L ||
      !value$provenance_kind %in% c("padded_psi", "frozen_materialized_root") ||
      !all(vapply(value[c("dataset_id", "dataset_version")],
        .dsvert_dp_analysis_client_scalar_id, logical(1L))) ||
      !all(vapply(value[c("authorization_epoch", "source_contract_id",
        "derivation_policy_id")], .dsvert_vector_hex, logical(1L)))) fail()
  epochs <- .dsvert_provenance_epochs_validate(value$authorization_epochs)
  if (identical(value$provenance_kind, "frozen_materialized_root") &&
      (length(epochs) != 1L ||
       !identical(epochs[[1L]]$authorization_epoch, value$authorization_epoch) ||
       !identical(epochs[[1L]]$source_contract_id, value$source_contract_id))) fail()
  .dsvert_joint_dp_client_canonical(value)
}

.dsvert_dp_synopsis_client_closure_validate_v2 <- function(value) {
  fail <- function() stop("Invalid Synopsis authorization closure.",
                          call. = FALSE)
  if (!.dsvert_dp_has_exact_names(value, c("version", "dependencies")) ||
      !identical(value$version, .DSVERT_CLIENT_SYNOPSIS_CLOSURE_VERSION) ||
      !is.list(value$dependencies) || !is.null(names(value$dependencies))) fail()
  keys <- vapply(value$dependencies, function(entry) {
    if (!.dsvert_dp_has_exact_names(entry,
        c("source_identity", "dataset", "dependency")) ||
        !.dsvert_dp_analysis_client_scalar_id(entry$dataset)) fail()
    identity <- .dsvert_dp_analysis_client_identity_pk(entry$source_identity)
    if (!identical(identity, entry$source_identity)) fail()
    .dsvert_dp_synopsis_client_dependency_v2(entry$dependency)
    paste(identity, entry$dataset, sep = "|")
  }, character(1L))
  if (anyDuplicated(keys) || !identical(keys, sort(keys, method = "radix"))) fail()
  .dsvert_joint_dp_client_canonical(value)
}

.dsvert_dp_synopsis_client_closure_v2 <- function(bootstraps, context) {
  entries <- list()
  for (peer in context$servers) {
    bootstrap <- bootstraps[[peer]]
    dependencies <- bootstrap$authorization_dependencies
    if (!is.list(dependencies) ||
        (length(dependencies) && (is.null(names(dependencies)) ||
         anyDuplicated(names(dependencies)) || !identical(names(dependencies),
           sort(names(dependencies), method = "radix")))) ||
        !all(names(bootstrap$draft$datasets) %in% names(dependencies))) {
      stop("The signed Synopsis dependency coverage is incomplete.",
           call. = FALSE)
    }
    for (dataset in names(dependencies)) {
      dependency <- .dsvert_dp_synopsis_client_dependency_v2(
        dependencies[[dataset]])
      draft <- bootstrap$draft$datasets[[dataset]]
      if (!is.null(draft) &&
          (!identical(draft$dataset_id, dependency$dataset_id) ||
           !identical(draft$dataset_version, dependency$dataset_version))) {
        stop("The signed Synopsis dependency conflicts with its dataset.",
             call. = FALSE)
      }
      if (identical(dependency$provenance_kind, "padded_psi")) {
        .dsvert_provenance_epochs_validate(
          dependency$authorization_epochs, unname(context$pinset))
      } else if (!identical(
          dependency$authorization_epochs[[1L]]$peer_identity,
          bootstrap$peer_identity_pk)) {
        stop("The signed Synopsis import names another custodian.",
             call. = FALSE)
      }
      entries[[length(entries) + 1L]] <- list(
        source_identity = bootstrap$peer_identity_pk,
        dataset = dataset, dependency = dependency)
    }
  }
  if (!length(entries)) stop(
    "The signed Synopsis authorization closure is empty.", call. = FALSE)
  keys <- vapply(entries, function(entry) paste(
    entry$source_identity, entry$dataset, sep = "|"), character(1L))
  .dsvert_dp_synopsis_client_closure_validate_v2(list(
    version = .DSVERT_CLIENT_SYNOPSIS_CLOSURE_VERSION,
    dependencies = unname(entries[order(keys, method = "radix")])))
}

.dsvert_dp_synopsis_client_source_claim_v2 <- function(claim, projection, peer,
                                                       context) {
  fields <- c("version", "source_peer_name", "source_identity_pk",
    "catalog_sha256", "source_vector_commitment", "signature")
  if (!.dsvert_dp_has_exact_names(claim, fields) ||
      !identical(claim$version,
        "dsvert-stateless-catalog-synopsis-source-vector-claim-v2") ||
      !identical(claim$source_peer_name, peer) ||
      !identical(claim$source_identity_pk, unname(context$pinset[[peer]])) ||
      !identical(claim$catalog_sha256, projection$sha256) ||
      !.dsvert_vector_hex(claim$source_vector_commitment)) {
    stop("Invalid signed Synopsis source Claim.", call. = FALSE)
  }
  if (!is.null(projection$catalog$authorization_dependencies)) {
    commitment <- .dsvert_dp_synopsis_client_hash(
      "dsVert/synopsis/public-source/v2|", list(
        version = "dsvert-synopsis-provenance-source-commitment-v1",
        catalog_sha256 = claim$catalog_sha256,
        source_identity = claim$source_identity_pk))
    if (!identical(commitment, claim$source_vector_commitment)) {
      stop("The Synopsis source commitment is not provenance-only.",
           call. = FALSE)
    }
  }
  unsigned <- .dsvert_joint_dp_client_canonical(
    claim[setdiff(names(claim), "signature")])
  .dsvert_dp_synopsis_client_signature(unsigned, claim$signature,
    claim$source_identity_pk,
    "dsVert/stateless-catalog-synopsis/source-vector-claim/v2|",
    "Synopsis source Claim")
  claim
}

.dsvert_dp_synopsis_client_effective_dependencies_v2 <- function(
    context, manifest) {
  closure <- .dsvert_dp_synopsis_client_closure_validate_v2(
    context$authorization_dependencies)
  release <- .dsvert_dp_capsule_vector_layout(manifest)
  cross <- .dsvert_dp_gaussian_cross_layout_client(manifest, release)
  required <- unique(vapply(c(release$blocks, cross$blocks), function(block) {
    if (is.null(block$owner_peer) ||
        !block$owner_peer %in% names(context$pinset) || is.null(block$dataset)) {
      stop("The public Synopsis layout lacks source dependency provenance.",
           call. = FALSE)
    }
    paste(unname(context$pinset[[block$owner_peer]]), block$dataset, sep = "|")
  }, character(1L)))
  keys <- vapply(closure$dependencies, function(entry) paste(
    entry$source_identity, entry$dataset, sep = "|"), character(1L))
  if (!all(required %in% keys)) stop(
    "The Synopsis authorization closure omits an effective source.",
    call. = FALSE)
  closure$dependencies <- closure$dependencies[keys %in% required]
  .dsvert_dp_synopsis_client_closure_validate_v2(closure)
}
