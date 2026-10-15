# Offline B1 chain validation over public synthetic stage certificates.
# Wire reconstruction mirrors the server; update() below is independent R/gmp.
.dsvert_dp_irls_client_state_json <- function(value) {
  .dsvert_joint_dp_client_json(.dsvert_joint_dp_client_canonical(value))
}

.dsvert_dp_irls_client_state_decode <- function(bytes) {
  if (!is.character(bytes) || length(bytes) != 1L || is.na(bytes) ||
      nchar(bytes, type = "bytes") > 16L * 1024L^2)
    stop("Invalid IRLS canonical record.", call. = FALSE)
  value <- jsonlite::fromJSON(bytes, simplifyVector = FALSE)
  if (!identical(.dsvert_dp_irls_client_state_json(value), bytes))
    stop("Noncanonical IRLS record.", call. = FALSE)
  value
}

.dsvert_dp_irls_client_stage_message <- function(body) {
  charToRaw(paste0("dsvert-v3-irls-successor-v1|",
                  .dsvert_dp_irls_client_state_json(body)))
}

.dsvert_dp_irls_client_stage_signatures <- function(stage, contract, complete = TRUE) {
  if (!is.list(stage) || !setequal(names(stage), c("body", "signatures")) ||
      anyDuplicated(names(stage)) || !is.list(stage$signatures) ||
      is.null(names(stage$signatures)) || anyDuplicated(names(stage$signatures)))
    stop("Invalid IRLS stage envelope.", call. = FALSE)
  pair <- unlist(contract$privacy_source$computation_peers, use.names = FALSE)
  pins <- contract$privacy_source$peer_pins
  peers <- names(stage$signatures)
  if ((!length(peers)) || !all(peers %in% pair) ||
      (complete && !setequal(peers, pair)))
    stop("IRLS stage requires designated authority signatures.", call. = FALSE)
  message <- .dsvert_dp_irls_client_stage_message(stage$body)
  for (peer in peers) {
    if (!isTRUE(.dsvert_dp_glm_grid_cross_verify(message, pins[[peer]],
                                            stage$signatures[[peer]])))
      stop("Invalid IRLS stage signature or pin.", call. = FALSE)
  }
  invisible(TRUE)
}

.dsvert_dp_irls_client_stage_body <- function(contract, t, predecessor, payload,
                                      chain, candidate_beta = NULL,
                                      update = .dsvert_dp_irls_verify_update) {
  id <- .dsvert_dp_irls_client_stage_id(contract, t, predecessor)
  index <- as.integer(t)
  if (length(chain) != index)
    stop("IRLS predecessor chain length mismatch.", call. = FALSE)
  expected_predecessor <- if (index == 0L) {
    .dsvert_dp_irls_client_start_digest(contract)
  } else chain[[index]]$body$payload_digest
  if (!identical(predecessor, expected_predecessor))
    stop("Alternate IRLS predecessor rejected.", call. = FALSE)
  beta <- if (index == 0L) contract$path$start$beta else
    chain[[index]]$body$update_result$beta
  if (!is.null(candidate_beta) && !identical(
      .dsvert_dp_irls_client_state_json(beta), .dsvert_dp_irls_client_state_json(candidate_beta)))
    stop("Alternate IRLS beta rejected.", call. = FALSE)
  payload_digest <- .dsvert_dp_irls_client_payload_digest(contract, t, payload)
  H <- payload[["H"]]
  if (is.null(H)) {
    ref <- payload[["H_reference"]]
    matches <- which(vapply(chain, function(stage) {
      identical(stage$body$stage_id, ref$stage_id) &&
        identical(stage$body$payload_digest, ref$payload_digest) &&
        !is.null(stage$body$payload[["H"]])
    }, logical(1L)))
    # The refresh schedule identifies the most recent signed curvature release.
    refreshed <- which(vapply(chain, function(stage)
      !is.null(stage$body$payload[["H"]]), logical(1L)))
    if (length(matches) != 1L || !length(refreshed) ||
        matches != tail(refreshed, 1L))
      stop("Invalid IRLS signed prior-H reference.", call. = FALSE)
    H <- chain[[matches]]$body$payload[["H"]]
  }
  purpose <- if (index == as.integer(contract$path$T))
    "terminal-audit" else "optimizer"
  input <- result <- NULL
  if (purpose == "optimizer") {
    numeric <- contract$numeric
    path <- contract$path
    input <- list(beta = beta, box_lower = path$box$lower,
      box_upper = path$box$upper, damping = path$damping[[index + 1L]],
      floor = path$floor, floor_quantum = path$floor_quantum,
      g = payload$g, g_translation = numeric$translation$g,
      g_unit = numeric$lattice$g, h = H,
      h_translation = numeric$translation$H, h_unit = numeric$lattice$H,
      lattice = numeric$lattice$beta, limits = numeric$certificate,
      ridge = path$ridge, trust = path$trust[[index + 1L]],
      version = path$update_version)
    result <- update(input)
  }
  list(root = .dsvert_dp_irls_client_root(contract), t = t, purpose = purpose,
    stage_id = id, predecessor = predecessor, payload = payload,
    payload_digest = payload_digest, beta = beta, update_input = input,
    update_result = result, update_version_digest = .dsvert_dp_irls_client_hash(
      "dsvert-v3-irls-update-version-v1", contract$path$update_version))
}


.dsvert_dp_irls_verify_chain <- function(signed_contract, trusted_pins, stages,
                                        expected_root = NULL) {
  contract <- .dsvert_dp_irls_client_contract_validate(signed_contract, trusted_pins)
  root <- .dsvert_dp_irls_client_root(contract)
  if (!is.null(expected_root) && !identical(root, expected_root))
    stop("IRLS root does not match the externally expected path.", call. = FALSE)
  if (!is.list(stages) || !is.null(names(stages)) ||
      length(stages) > as.integer(contract$path$T) + as.integer(contract$path$terminal_audit))
    stop("Invalid IRLS chain length.", call. = FALSE)
  verified <- observed <- list()
  beta <- contract$path$start$beta
  for (i in seq_along(stages)) {
    stage <- stages[[i]]
    .dsvert_dp_irls_client_stage_signatures(stage, contract)
    expected <- .dsvert_dp_irls_client_stage_body(contract, as.character(i - 1L),
      stage$body$predecessor, stage$body$payload, verified,
      candidate_beta = stage$body$beta)
    if (!identical(.dsvert_dp_irls_client_state_json(expected),
                   .dsvert_dp_irls_client_state_json(stage$body)))
      stop("IRLS chain contains an incorrect update or identity.", call. = FALSE)
    verified[[i]] <- stage
    if (!is.null(expected$update_result)) beta <- expected$update_result$beta
    observed[[i]] <- list(t = expected$t, lanes = if (is.null(expected$payload[["H"]]))
      list("g") else list("g", "H"))
  }
  list(root = root, beta = beta,
    chain_complete = length(stages) == length(contract$accounting$stages),
    complete = length(stages) == length(contract$accounting$stages) &&
      length(contract$accounting$auxiliary) == 0L,
    accounting = .dsvert_dp_irls_client_accounting(contract, observed))
}
