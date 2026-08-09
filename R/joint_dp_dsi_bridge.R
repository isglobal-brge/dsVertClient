# Canonical message helpers and the complete pinned-peer context shared by the
# promoted manifest, source-transport and vector protocols.  The retired
# scalar joint-DP allocation/delivery bridge deliberately does not live here.

.dsvert_joint_dp_client_canonical <- function(value) {
  if (is.null(value)) return(NULL)
  if (is.object(value)) {
    stop("Unsupported joint-DP message value", call. = FALSE)
  }
  if (is.list(value)) {
    fields <- names(value)
    if (!is.null(fields)) {
      if (anyNA(fields) || any(!nzchar(fields)) || anyDuplicated(fields)) {
        stop("Invalid joint-DP message fields", call. = FALSE)
      }
      value <- value[order(fields, method = "radix")]
    }
    return(lapply(value, .dsvert_joint_dp_client_canonical))
  }
  if (!typeof(value) %in% c("logical", "integer", "double", "character") ||
      !is.null(names(value)) || anyNA(value) ||
      (is.numeric(value) && any(!is.finite(value)))) {
    stop("Unsupported joint-DP message value", call. = FALSE)
  }
  if (is.character(value)) return(enc2utf8(unname(value)))
  if (is.numeric(value)) {
    value <- unname(as.numeric(value))
    value[value == 0] <- 0
    return(value)
  }
  unname(value)
}

.dsvert_joint_dp_client_json <- function(value) {
  as.character(jsonlite::toJSON(
    .dsvert_joint_dp_client_canonical(value), auto_unbox = TRUE,
    null = "null", na = "null", digits = 17, pretty = FALSE))
}

.dsvert_joint_dp_client_decode <- function(value, what, maximum_bytes) {
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !nzchar(value) || nchar(value, type = "bytes") > maximum_bytes) {
    stop("Invalid joint-DP ", what, " returned by a peer", call. = FALSE)
  }
  decoded <- tryCatch(
    jsonlite::fromJSON(value, simplifyVector = FALSE),
    error = function(error) NULL)
  canonical <- tryCatch(
    .dsvert_joint_dp_client_json(decoded), error = function(error) NULL)
  if (is.null(canonical) || !identical(canonical, value) ||
      !is.list(decoded)) {
    stop("A peer returned a non-canonical joint-DP ", what, call. = FALSE)
  }
  decoded
}

.dsvert_joint_dp_client_b64url <- function(value, bytes, what) {
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !grepl("^[A-Za-z0-9_-]+$", value) || nchar(value) %% 4L == 1L) {
    stop("Invalid ", what, " in a joint-DP receipt", call. = FALSE)
  }
  padding <- (4L - nchar(value) %% 4L) %% 4L
  standard <- paste0(chartr("-_", "+/", value), strrep("=", padding))
  decoded <- tryCatch(jsonlite::base64_dec(standard),
                      error = function(error) raw(0L))
  canonical <- if (is.raw(decoded)) {
    chartr("+/", "-_", sub(
      "=+$", "", gsub("[\r\n]", "", jsonlite::base64_enc(decoded)),
      perl = TRUE))
  } else {
    ""
  }
  if (!is.raw(decoded) || length(decoded) != bytes ||
      !identical(canonical, value)) {
    stop("Invalid ", what, " in a joint-DP receipt", call. = FALSE)
  }
  decoded
}

.dsvert_joint_dp_client_context <- function(
    datasources, status = NULL, .aggregate = DSI::datashield.aggregate) {
  if (!is.list(datasources) || length(datasources) < 2L ||
      is.null(names(datasources)) || anyNA(names(datasources)) ||
      any(!nzchar(names(datasources))) || anyDuplicated(names(datasources))) {
    stop("Joint DP requires the complete named datasource set", call. = FALSE)
  }
  if (is.null(status)) {
    status <- .dsvert_joint_dp_capsule_status_impl(datasources, .aggregate)
  } else {
    status <- .dsvert_joint_dp_capsule_status_consensus(status, datasources)
  }
  servers <- sort(names(datasources), method = "radix")
  if (!is.list(status) || is.null(names(status)) || anyNA(names(status)) ||
      anyDuplicated(names(status)) || !setequal(names(status), servers)) {
    stop("Joint-DP status does not cover the complete datasource set",
         call. = FALSE)
  }
  status <- status[servers]
  policies <- lapply(status, `[[`, "policy")
  if (any(!vapply(policies, is.list, logical(1L)))) {
    stop("Joint-DP status lacks a policy", call. = FALSE)
  }
  reference <- policies[[1L]]
  pinset <- reference$peer_pinset
  designated <- reference$designated_noise_peers
  valid <- is.character(pinset) && !is.null(names(pinset)) &&
    identical(names(pinset), servers) && !anyNA(pinset) &&
    !anyDuplicated(pinset) &&
    all(vapply(pinset, function(value) {
      !is.null(.dsvert_dp_normalize_identity_pk(value))
    }, logical(1L))) &&
    .dsvert_dp_is_integer(reference$peer_count, 2, .Machine$integer.max) &&
    reference$peer_count == length(servers) &&
    is.character(designated) && length(designated) == 2L &&
    !anyNA(designated) && !anyDuplicated(designated) &&
    identical(designated, sort(designated, method = "radix")) &&
    all(designated %in% servers)
  if (!isTRUE(valid)) {
    stop("Joint-DP status has an invalid full pinset or designated pair",
         call. = FALSE)
  }
  for (server in servers) {
    policy <- policies[[server]]
    if (!identical(policy$peer_name, server) ||
        !identical(policy$peer_pinset, pinset) ||
        !identical(policy$peer_count, reference$peer_count) ||
        !identical(policy$designated_noise_peers, designated) ||
        !identical(policy$peer_pinset_sha256,
                   reference$peer_pinset_sha256) ||
        !identical(policy$own_identity_pk, unname(pinset[[server]]))) {
      stop("Connected peers disagree on the custodian-owned joint-DP pinset",
           call. = FALSE)
    }
  }
  list(
    status = status, servers = servers, pinset = pinset,
    designated = designated, conns = datasources[designated],
    all_conns = datasources[servers])
}
