# Fail-closed validation of the outer DataSHIELD transport.  Peer encryption
# protects purpose-bound MPC payloads, but it does not protect DataSHIELD
# credentials or ordinary DP releases from a downgraded connector.

.dsvert_dsi_member <- function(value, name) {
  if (isS4(value) && name %in% methods::slotNames(value)) {
    return(methods::slot(value, name))
  }
  if (is.environment(value)) {
    if (exists(name, envir = value, inherits = FALSE)) {
      return(get(name, envir = value, inherits = FALSE))
    }
    return(NULL)
  }
  if (is.list(value) && name %in% names(value)) return(value[[name]])
  NULL
}

.dsvert_dsi_scalar_string <- function(value) {
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !nzchar(trimws(value))) return(NULL)
  trimws(value)
}

# `format(externalptr)` is stable while one libcurl handle is alive, but an
# allocator may reuse that address after the handle is finalized.  Attach a
# process-local random lifecycle id to the exact external pointer and retain
# the address as a second discriminator. Token refresh and ordinary S4 copies
# share the pointer (and therefore the id); a fresh login receives a new id
# even in the pathological address-reuse case.
.dsvert_dsi_armadillo_handle_identity <- function(connection) {
  if (!inherits(connection, "ArmadilloConnection")) return(NULL)
  handle <- .dsvert_dsi_member(connection, "handle")
  curl_handle <- .dsvert_dsi_member(handle, "handle")
  if (!identical(typeof(curl_handle), "externalptr")) return(NULL)
  address <- format(curl_handle)
  if (!is.character(address) || length(address) != 1L || is.na(address) ||
      !nzchar(address)) return(NULL)
  attribute <- "dsvert.dsi.handle.lifecycle.v1"
  lifecycle <- attr(curl_handle, attribute, exact = TRUE)
  if (!is.character(lifecycle) || length(lifecycle) != 1L ||
      is.na(lifecycle) || !grepl("^h_[0-9a-f]{32}$", lifecycle)) {
    lifecycle <- paste0(
      "h_", paste(sprintf("%02x", as.integer(openssl::rand_bytes(16L))),
                   collapse = ""))
    attr(curl_handle, attribute) <- lifecycle
    if (!identical(attr(curl_handle, attribute, exact = TRUE), lifecycle)) {
      return(NULL)
    }
  }
  paste0("curl-handle:", lifecycle, "@", address)
}

.dsvert_dsi_endpoint <- function(url) {
  url <- .dsvert_dsi_scalar_string(url)
  if (is.null(url) ||
      !grepl("^[A-Za-z][A-Za-z0-9+.-]*://", url)) return(NULL)
  scheme <- tolower(sub(
    "^([A-Za-z][A-Za-z0-9+.-]*)://.*$", "\\1", url))
  remainder <- sub("^[A-Za-z][A-Za-z0-9+.-]*://", "", url)
  authority_full <- sub("[/?#].*$", "", remainder)
  tail <- substr(remainder, nchar(authority_full, type = "chars") + 1L,
                 nchar(remainder, type = "chars"))
  path <- sub("[?#].*$", "", tail)
  if (!nzchar(path)) path <- "/"
  authority <- authority_full
  authority <- sub("^.*@", "", authority)
  if (!nzchar(authority)) return(NULL)
  port <- NULL
  if (startsWith(authority, "[")) {
    close <- regexpr("]", authority, fixed = TRUE)[[1L]]
    if (close < 2L) return(NULL)
    host <- substr(authority, 2L, close - 1L)
    suffix <- substr(authority, close + 1L, nchar(authority))
    if (nzchar(suffix)) {
      if (!grepl("^:[0-9]{1,5}$", suffix)) return(NULL)
      port <- substring(suffix, 2L)
    }
  } else {
    if (lengths(regmatches(authority, gregexpr(":", authority,
                                                fixed = TRUE))) > 1L) {
      return(NULL)
    }
    if (grepl(":", authority, fixed = TRUE)) {
      host <- sub(":.*$", "", authority)
      port <- sub("^.*:", "", authority)
      if (!grepl("^[0-9]{1,5}$", port)) return(NULL)
    } else {
      host <- authority
    }
  }
  host <- sub("\\.$", "", tolower(host))
  if (!nzchar(host)) return(NULL)
  if (is.null(port)) {
    port <- if (identical(scheme, "https")) "443" else
      if (identical(scheme, "http")) "80" else ""
  }
  if (nzchar(port)) {
    port_number <- suppressWarnings(as.integer(port))
    if (length(port_number) != 1L || is.na(port_number) ||
        port_number < 1L || port_number > 65535L) return(NULL)
    port <- as.character(port_number)
  }
  display_host <- if (grepl(":", host, fixed = TRUE)) {
    paste0("[", host, "]")
  } else {
    host
  }
  origin <- paste0(
    scheme, "://", display_host,
    if (nzchar(port)) paste0(":", port) else "")
  list(scheme = scheme, host = host, port = port,
       scope = paste0(origin, path))
}

.dsvert_dsi_loopback_host <- function(host) {
  host <- .dsvert_dsi_scalar_string(host)
  if (is.null(host)) return(FALSE)
  host <- sub("%.*$", "", tolower(host))
  if (host %in% c("localhost", "::1", "0:0:0:0:0:0:0:1")) return(TRUE)
  parts <- strsplit(host, ".", fixed = TRUE)[[1L]]
  length(parts) == 4L && identical(parts[[1L]], "127") &&
    all(grepl("^[0-9]{1,3}$", parts)) &&
    all(as.integer(parts) >= 0L & as.integer(parts) <= 255L)
}

.dsvert_dsi_tls_verification_state <- function(options) {
  if (!is.list(options)) return("unknown")
  peer <- options[["ssl_verifypeer"]]
  if (!is.null(peer)) {
    peer <- suppressWarnings(as.numeric(peer))
    if (length(peer) != 1L || is.na(peer) || peer <= 0) return("disabled")
  }
  host <- options[["ssl_verifyhost"]]
  if (!is.null(host)) {
    if (is.logical(host) && length(host) == 1L && !is.na(host)) {
      return(if (isTRUE(host)) "verified" else "disabled")
    }
    host <- suppressWarnings(as.numeric(host))
    if (length(host) != 1L || is.na(host) || host < 2) return("disabled")
  }
  "verified"
}

.dsvert_dsi_httr_global_options <- function() {
  config <- getOption("httr_config", NULL)
  if (is.null(config)) return(list())
  options <- .dsvert_dsi_member(config, "options")
  # A real empty httr request may retain options=NULL, which means curl's
  # verification defaults. Any other present but uninspectable global config
  # is ambiguous and must not be treated like the option being absent.
  if (is.null(options) && inherits(config, "request")) return(list())
  if (!is.list(options)) {
    stop(
      "The process-wide httr_config is present but its TLS options cannot be inspected",
      call. = FALSE)
  }
  options
}

.dsvert_dsi_effective_opal_tls_state <- function(local_options) {
  global_options <- .dsvert_dsi_httr_global_options()
  if (!is.list(global_options) || !is.list(local_options)) {
    if (identical(
          .dsvert_dsi_tls_verification_state(global_options), "disabled")) {
      return("disabled")
    }
    return("unknown")
  }
  local_options <- local_options[
    !vapply(local_options, is.null, logical(1L))]
  .dsvert_dsi_tls_verification_state(utils::modifyList(
    global_options, local_options, keep.null = TRUE))
}

.dsvert_dsi_inspect_connection <- function(connection) {
  if (inherits(connection, "DSLiteConnection")) {
    return(list(kind = "dslite", endpoint = NULL,
                verification = "in_process"))
  }
  if (inherits(connection, "OpalConnection")) {
    opal <- .dsvert_dsi_member(connection, "opal")
    config <- .dsvert_dsi_member(opal, "config")
    return(list(
      kind = "opal",
      endpoint = .dsvert_dsi_endpoint(.dsvert_dsi_member(opal, "url")),
      verification = .dsvert_dsi_effective_opal_tls_state(
        .dsvert_dsi_member(config, "options"))))
  }
  if (inherits(connection, "ArmadilloConnection")) {
    handle <- .dsvert_dsi_member(connection, "handle")
    global_state <- .dsvert_dsi_tls_verification_state(
      .dsvert_dsi_httr_global_options())
    return(list(
      kind = "armadillo",
      endpoint = .dsvert_dsi_endpoint(.dsvert_dsi_member(handle, "url")),
      verification = if (identical(global_state, "disabled")) {
        "disabled"
      } else {
        "unknown"
      }))
  }
  list(
    kind = "generic",
    endpoint = .dsvert_dsi_endpoint(
      .dsvert_dsi_member(connection, "url")),
    verification = "unknown")
}

.dsvert_dsi_test_loopback_allowed <- function() {
  # The installed client has no option or environment variable that permits a
  # plaintext connector. Tests replace this private binding locally while
  # exercising an in-process HTTP stub; that replacement is never installed.
  FALSE
}

.dsvert_validate_dsi_transport_security <- function(connections) {
  sites <- names(connections)
  if (!is.list(connections) || !length(connections) || is.null(sites) ||
      anyNA(sites) || any(!nzchar(sites)) || anyDuplicated(sites)) {
    stop("DSI connections must have unique logical site names",
         call. = FALSE)
  }
  allow_test_loopback <- .dsvert_dsi_test_loopback_allowed()
  for (index in seq_along(connections)) {
    site <- sites[[index]]
    details <- .dsvert_dsi_inspect_connection(connections[[index]])
    if (identical(details$kind, "dslite")) next
    endpoint <- details$endpoint
    if (is.null(endpoint)) {
      stop("Site '", site,
           "' has no inspectable DSI endpoint; verified HTTPS is required",
           call. = FALSE)
    }
    if (.dsvert_dsi_loopback_host(endpoint$host) &&
        isTRUE(allow_test_loopback)) next
    if (!identical(endpoint$scheme, "https")) {
      stop("Site '", site,
           "' does not use verified HTTPS for its outer DSI transport",
           call. = FALSE)
    }
    if (identical(details$verification, "disabled")) {
      stop("Site '", site,
           "' has TLS certificate verification disabled",
           call. = FALSE)
    }
    if (identical(details$kind, "opal") &&
        identical(details$verification, "verified")) next
    # DSMolgenisArmadillo retains the authenticated endpoint URL, but its
    # connection object does not retain inspectable curl verification flags.
    # Match dsFlowerClient: a recognized Armadillo connector with a valid
    # HTTPS endpoint is accepted automatically.  The connector/reverse proxy
    # remains responsible for certificate and hostname validation; an
    # inherited process-wide TLS downgrade was rejected above.
    if (identical(details$kind, "armadillo")) next
    # An unknown connector has no inspectable TLS-verification contract.  A
    # process-local analyst option cannot turn it into one, so generic
    # connectors remain rejected.
    if (identical(details$kind, "generic")) {
      stop("Site '", site,
           "' uses an unsupported generic DSI connector; an inspectable ",
           "provider-specific HTTPS contract is required", call. = FALSE)
    }
    stop("Site '", site,
         "' has no supported inspectable DSI security contract",
         call. = FALSE)
  }
  invisible(TRUE)
}

.dsvert_require_dsi_fanout_version <- function(
    .package_version = utils::packageVersion) {
  version <- tryCatch(
    as.character(.package_version("DSI")),
    error = function(e) NA_character_)
  if (length(version) != 1L || is.na(version) ||
      utils::compareVersion(version, "1.8.0") < 0L) {
    stop("dsVertClient requires DSI >= 1.8.0 for named asynchronous fan-out",
         call. = FALSE)
  }
  invisible(TRUE)
}

.dsvert_validate_real_dsi_transport <- function(connections, aggregate) {
  if (identical(aggregate, DSI::datashield.aggregate)) {
    .dsvert_require_dsi_fanout_version()
    .dsvert_validate_dsi_transport_security(connections)
  }
  invisible(TRUE)
}
