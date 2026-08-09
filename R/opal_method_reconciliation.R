.dsvert_validate_opal_profile <- function(profile_name) {
  if (!is.character(profile_name) || length(profile_name) != 1L ||
      is.na(profile_name) ||
      !grepl("^[A-Za-z0-9._-]{1,64}$", profile_name)) {
    stop("Invalid Opal DataSHIELD profile name.", call. = FALSE)
  }
  if (identical(profile_name, "default")) {
    stop(
      "dsVert requires a dedicated non-default Opal profile.",
      call. = FALSE)
  }
  invisible(profile_name)
}

.dsvert_opal_method_contract <- function(description_path) {
  if (!is.character(description_path) || length(description_path) != 1L ||
      is.na(description_path) || !file.exists(description_path)) {
    stop("The dsVert DESCRIPTION file is unavailable.", call. = FALSE)
  }
  description <- tryCatch(
    read.dcf(description_path),
    error = function(e) stop(
      "The dsVert DESCRIPTION file is invalid: ", conditionMessage(e),
      call. = FALSE))
  if (nrow(description) != 1L || !"Package" %in% colnames(description) ||
      is.na(description[1L, "Package"]) ||
      !identical(unname(trimws(description[1L, "Package"])), "dsVert")) {
    stop("The method contract does not belong to dsVert.", call. = FALSE)
  }

  parse_field <- function(field, type) {
    if (!field %in% colnames(description)) {
      return(data.frame(name = character(), type = character(),
                        value = character(), stringsAsFactors = FALSE))
    }
    raw <- description[1L, field]
    if (is.na(raw)) {
      stop("The dsVert remote-method contract contains a missing field.",
           call. = FALSE)
    }
    names <- trimws(strsplit(raw, ",", fixed = TRUE)[[1L]])
    names <- names[nzchar(names)]
    aliases <- names[grepl("=", names, fixed = TRUE)]
    if (length(aliases)) {
      stop(
        "Remote method aliases are forbidden: ",
        paste(aliases, collapse = ", "), call. = FALSE)
    }
    if (any(!grepl("^[A-Za-z][A-Za-z0-9.]*DS$", names))) {
      stop("The dsVert remote-method contract contains an invalid name.",
           call. = FALSE)
    }
    data.frame(
      name = names, type = rep(type, length(names)),
      value = paste0("dsVert::", names), stringsAsFactors = FALSE)
  }

  contract <- rbind(
    parse_field("AggregateMethods", "aggregate"),
    parse_field("AssignMethods", "assign"))
  if (!nrow(contract)) {
    stop("The dsVert remote-method allowlist is empty.", call. = FALSE)
  }
  if (anyDuplicated(paste(contract$type, contract$name, sep = "\r")) ||
      anyDuplicated(contract$name)) {
    stop("The dsVert remote-method contract contains duplicate names.",
         call. = FALSE)
  }
  contract
}

.dsvert_remote_surface_attestation <- function(contract) {
  if (!is.data.frame(contract) ||
      !identical(names(contract), c("name", "type", "value")) ||
      !nrow(contract) || anyNA(contract) ||
      any(!contract$type %in% c("aggregate", "assign")) ||
      any(!grepl("^[A-Za-z][A-Za-z0-9.]*DS$", contract$name)) ||
      anyDuplicated(contract$name) ||
      any(contract$value != paste0("dsVert::", contract$name))) {
    stop("Invalid dsVert remote-method contract.", call. = FALSE)
  }
  ordered <- contract[order(contract$type, contract$name, method = "radix"),
                      , drop = FALSE]
  canonical <- paste(
    paste(ordered$type, ordered$name, ordered$value, sep = "\x1f"),
    collapse = "\x1e")
  hash <- digest::digest(canonical, algo = "sha256", serialize = FALSE)
  paste0(
    "dsvert-custodian-surface-attestation-v1:",
    "dsvert-disclosure-safe-v1:", hash)
}

.dsvert_reconcile_opal_methods <- function(
    opal, description_path, profile_name = "dsvert",
    get_methods = NULL, remove_method = NULL, set_method = NULL) {
  .dsvert_validate_opal_profile(profile_name)
  contract <- .dsvert_opal_method_contract(description_path)
  if (is.null(get_methods) || is.null(remove_method) || is.null(set_method)) {
    if (!requireNamespace("opalr", quietly = TRUE)) {
      stop("The opalr package is required to reconcile Opal methods.",
           call. = FALSE)
    }
    get_methods <- opalr::dsadmin.get_methods
    remove_method <- opalr::dsadmin.rm_method
    set_method <- opalr::dsadmin.set_method
  }
  if (!is.function(get_methods) || !is.function(remove_method) ||
      !is.function(set_method)) {
    stop("Invalid Opal method reconciliation backend.", call. = FALSE)
  }

  read_inventory <- function(stage) {
    by_type <- lapply(c("aggregate", "assign"), function(type) {
      value <- get_methods(opal, type = type, profile = profile_name)
      if (!is.data.frame(value)) {
        stop("Opal returned an invalid ", stage, " method inventory.",
             call. = FALSE)
      }
      # opalr returns a zero-column data.frame for an empty method list.
      if (!nrow(value)) {
        return(data.frame(
          name = character(), type = character(), value = character(),
          stringsAsFactors = FALSE))
      }
      if (!all(c("name", "class", "value") %in% names(value))) {
        stop("Opal returned an invalid ", stage, " method inventory.",
             call. = FALSE)
      }
      value$name <- as.character(value$name)
      value$class <- as.character(value$class)
      value$value <- as.character(value$value)
      if (anyNA(value$name) || anyNA(value$class) || anyNA(value$value) ||
          any(!nzchar(value$name)) || any(!nzchar(value$value))) {
        stop("Opal returned an invalid ", stage, " method inventory.",
             call. = FALSE)
      }
      if (any(value$class != "function")) {
        stop("Opal returned a script-backed method in the dedicated profile.",
             call. = FALSE)
      }
      value$type <- type
      value[, c("name", "type", "value"), drop = FALSE]
    })
    do.call(rbind, by_type)
  }

  sorted <- function(value) {
    value <- value[order(value$type, value$name, value$value,
                         method = "radix"), , drop = FALSE]
    rownames(value) <- NULL
    value
  }
  expected <- sorted(contract)
  existing <- read_inventory("pre-registration")

  # This profile is exclusive to dsVert. A shared profile cannot support an
  # exact-surface attestation because any unrelated callable expands the
  # analyst-visible disclosure surface.
  if (!identical(sorted(existing), expected)) {
    removal_keys <- unique(existing[, c("name", "type"), drop = FALSE])
    for (row in seq_len(nrow(removal_keys))) {
      remove_method(
        opal, name = removal_keys$name[[row]],
        type = removal_keys$type[[row]], profile = profile_name)
    }
    for (row in seq_len(nrow(contract))) {
      set_method(
        opal, name = contract$name[[row]], func = contract$value[[row]],
        type = contract$type[[row]], profile = profile_name)
    }
  }

  after <- read_inventory("post-registration")
  if (!identical(sorted(after), expected)) {
    stop(
      "Opal did not converge to dsVert's exact exclusive allowlist.",
      call. = FALSE)
  }
  invisible(contract)
}

.dsvert_attest_opal_surface <- function(
    opal, contract, profile_name = "dsvert",
    set_option = NULL, get_options = NULL) {
  .dsvert_validate_opal_profile(profile_name)
  token <- .dsvert_remote_surface_attestation(contract)
  if (is.null(set_option) || is.null(get_options)) {
    if (!requireNamespace("opalr", quietly = TRUE)) {
      stop("The opalr package is required to attest the Opal surface.",
           call. = FALSE)
    }
    set_option <- opalr::dsadmin.set_option
    get_options <- opalr::dsadmin.get_options
  }
  if (!is.function(set_option) || !is.function(get_options)) {
    stop("Invalid Opal surface-attestation backend.", call. = FALSE)
  }
  option_name <- "dsvert.remote_surface_attestation"
  set_option(
    opal, option_name, token, profile = profile_name)
  options <- get_options(opal, profile = profile_name)
  if (!is.data.frame(options) ||
      !all(c("name", "value") %in% names(options))) {
    stop("Opal returned an invalid option inventory.", call. = FALSE)
  }
  matches <- which(as.character(options$name) == option_name)
  if (length(matches) != 1L || is.na(options$value[[matches]]) ||
      !identical(as.character(options$value[[matches]]), token)) {
    stop("Opal did not persist the dsVert surface attestation.",
         call. = FALSE)
  }
  invisible(token)
}
