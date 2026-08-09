.opal_inventory_frame <- function(rows = list()) {
  if (!length(rows)) {
    return(data.frame(
      name = character(), type = character(), class = character(),
      value = character(), package = character(), version = character(),
      stringsAsFactors = FALSE))
  }
  do.call(rbind, lapply(rows, function(row) {
    data.frame(
      name = row$name, type = row$type, class = "function",
      value = row$value, package = row$package, version = "1",
      stringsAsFactors = FALSE)
  }))
}

.opal_contract_description <- function(lines = character()) {
  description <- tempfile("dsvert-description-")
  writeLines(c(
    "Package: dsVert",
    "Version: 1.0.0",
    "AggregateMethods: safeAggregateDS",
    "AssignMethods: safeAssignDS",
    lines), description)
  description
}

test_that("Opal reconciliation makes the dedicated profile exactly dsVert", {
  description <- .opal_contract_description()
  withr::defer(unlink(description))

  state <- new.env(parent = emptyenv())
  state$aggregate <- .opal_inventory_frame(list(
    list(name = "c", type = "aggregate", value = "base::c",
         package = "base"),
    list(name = "list", type = "aggregate", value = "base::list",
         package = "base"),
    list(name = "legacyExactDS", type = "aggregate",
         value = "dsVert::legacyExactDS", package = "dsVert"),
    list(name = "otherAggregateDS", type = "aggregate",
         value = "dsBase::otherAggregateDS", package = "dsBase")))
  state$assign <- .opal_inventory_frame(list(
    list(name = "numeric", type = "assign", value = "base::numeric",
         package = "base"),
    list(name = "character", type = "assign", value = "base::character",
         package = "base"),
    list(name = "safeAggregateDS", type = "assign",
         value = "other::safeAggregateDS", package = "other"),
    list(name = "legacyAssignDS", type = "assign",
         value = "dsVert::legacyAssignDS", package = "dsVert")))
  state$profiles <- character()
  state$mutations <- 0L

  get_methods <- function(opal, type, profile) {
    state$profiles <- c(state$profiles, profile)
    state[[type]]
  }
  remove_method <- function(opal, name, type, profile) {
    state$profiles <- c(state$profiles, profile)
    state$mutations <- state$mutations + 1L
    state[[type]] <- state[[type]][state[[type]]$name != name, , drop = FALSE]
    invisible(NULL)
  }
  set_method <- function(opal, name, func, type, profile) {
    state$profiles <- c(state$profiles, profile)
    state$mutations <- state$mutations + 1L
    state[[type]] <- rbind(state[[type]], .opal_inventory_frame(list(
      list(name = name, type = type, value = func, package = "dsVert"))))
    invisible(NULL)
  }

  expect_invisible(.dsvert_reconcile_opal_methods(
    opal = list(), description_path = description,
    profile_name = "research", get_methods = get_methods,
    remove_method = remove_method, set_method = set_method))
  expect_identical(state$aggregate$name, "safeAggregateDS")
  expect_identical(state$aggregate$value, "dsVert::safeAggregateDS")
  expect_identical(state$assign$name, "safeAssignDS")
  expect_identical(state$assign$value, "dsVert::safeAssignDS")
  expect_true(all(state$profiles == "research"))

  mutations <- state$mutations
  expect_invisible(.dsvert_reconcile_opal_methods(
    opal = list(), description_path = description,
    profile_name = "research", get_methods = get_methods,
    remove_method = remove_method, set_method = set_method))
  expect_identical(state$mutations, mutations)
})

test_that("Opal reconciliation rejects unsafe contracts before server I/O", {
  touched <- FALSE
  backend <- function(...) {
    touched <<- TRUE
    data.frame()
  }

  alias_description <- tempfile("dsvert-alias-description-")
  writeLines(c(
    "Package: dsVert",
    "Version: 1.0.0",
    "AggregateMethods: c=base::c"), alias_description)
  withr::defer(unlink(alias_description))
  expect_error(
    .dsvert_reconcile_opal_methods(
      opal = list(), description_path = alias_description,
      get_methods = backend, remove_method = backend, set_method = backend),
    "Remote method aliases are forbidden")
  expect_false(touched)

  wrong_package <- tempfile("wrong-package-description-")
  writeLines(c(
    "Package: other",
    "Version: 1.0.0",
    "AggregateMethods: safeAggregateDS"), wrong_package)
  withr::defer(unlink(wrong_package))
  expect_error(
    .dsvert_reconcile_opal_methods(
      opal = list(), description_path = wrong_package,
      get_methods = backend, remove_method = backend, set_method = backend),
    "does not belong to dsVert")
  expect_false(touched)

  safe_description <- .opal_contract_description()
  withr::defer(unlink(safe_description))
  expect_error(
    .dsvert_reconcile_opal_methods(
      opal = list(), description_path = safe_description,
      profile_name = "default", get_methods = backend,
      remove_method = backend, set_method = backend),
    "dedicated non-default")
  expect_false(touched)
})

test_that("invalid Opal inventory fails before any mutation", {
  description <- .opal_contract_description()
  withr::defer(unlink(description))
  calls <- 0L
  get_methods <- function(opal, type, profile) {
    if (identical(type, "aggregate")) {
      data.frame(name = NA_character_, value = "dsVert::badDS")
    } else {
      data.frame(name = character(), value = character())
    }
  }
  mutate <- function(...) {
    calls <<- calls + 1L
  }
  expect_error(
    .dsvert_reconcile_opal_methods(
      list(), description, get_methods = get_methods,
      remove_method = mutate, set_method = mutate),
    "invalid pre-registration")
  expect_identical(calls, 0L)
})

test_that("fresh Opal zero-column inventories are normalized", {
  description <- .opal_contract_description()
  withr::defer(unlink(description))
  state <- new.env(parent = emptyenv())
  state$aggregate <- data.frame()
  state$assign <- data.frame()
  get_methods <- function(opal, type, profile) state[[type]]
  remove_method <- function(...) {
    stop("an empty profile must not require removal")
  }
  set_method <- function(opal, name, func, type, profile) {
    state[[type]] <- rbind(
      state[[type]],
      data.frame(
        name = name, class = "function", value = func,
        stringsAsFactors = FALSE))
  }
  expect_invisible(.dsvert_reconcile_opal_methods(
    list(), description, get_methods = get_methods,
    remove_method = remove_method, set_method = set_method))
  expect_identical(state$aggregate$name, "safeAggregateDS")
  expect_identical(state$assign$name, "safeAssignDS")
})

test_that("foreign method injected before final inventory check fails closed", {
  description <- .opal_contract_description()
  withr::defer(unlink(description))
  reads <- c(aggregate = 0L, assign = 0L)
  inventories <- list(
    aggregate = .opal_inventory_frame(list(list(
      name = "safeAggregateDS", type = "aggregate",
      value = "dsVert::safeAggregateDS", package = "dsVert"))),
    assign = .opal_inventory_frame(list(list(
      name = "safeAssignDS", type = "assign",
      value = "dsVert::safeAssignDS", package = "dsVert"))))
  get_methods <- function(opal, type, profile) {
    reads[[type]] <<- reads[[type]] + 1L
    value <- inventories[[type]]
    if (identical(type, "aggregate") && reads[[type]] > 1L) {
      value <- rbind(value, .opal_inventory_frame(list(list(
        name = "foreignExactDS", type = "aggregate",
        value = "other::foreignExactDS", package = "other"))))
    }
    value
  }
  mutate <- function(...) stop("an initially exact inventory must not mutate")

  expect_error(
    .dsvert_reconcile_opal_methods(
      list(), description, get_methods = get_methods,
      remove_method = mutate, set_method = mutate),
    "exact exclusive allowlist")
})

test_that("surface attestation is connector-neutral and Opal persistence is verified", {
  description <- .opal_contract_description()
  withr::defer(unlink(description))
  contract <- .dsvert_opal_method_contract(description)
  reversed <- contract[rev(seq_len(nrow(contract))), , drop = FALSE]
  expect_identical(
    .dsvert_remote_surface_attestation(contract),
    .dsvert_remote_surface_attestation(reversed))
  expect_match(
    .dsvert_remote_surface_attestation(contract),
    paste0("^dsvert-custodian-surface-attestation-v1:",
           "dsvert-disclosure-safe-v1:[0-9a-f]{64}$"))
  expect_identical(
    .dsvert_remote_surface_attestation(contract),
    paste0(
      "dsvert-custodian-surface-attestation-v1:",
      "dsvert-disclosure-safe-v1:",
      "dcca8dcf15580d69dc64893918cf8e78063438497c83d7c8203cd8284a1117b7"))

  state <- new.env(parent = emptyenv())
  state$options <- data.frame(
    name = character(), value = character(), stringsAsFactors = FALSE)
  set_option <- function(opal, name, value, profile) {
    state$profile <- profile
    state$options <- data.frame(
      name = name, value = value, stringsAsFactors = FALSE)
  }
  get_options <- function(opal, profile) state$options
  token <- .dsvert_attest_opal_surface(
    list(), contract, set_option = set_option, get_options = get_options)
  expect_identical(state$profile, "dsvert")
  expect_identical(state$options$value, token)

  state$options <- data.frame(
    name = character(), value = character(), stringsAsFactors = FALSE)
  token_research <- .dsvert_attest_opal_surface(
    list(), contract, profile_name = "research",
    set_option = set_option, get_options = get_options)
  expect_identical(state$profile, "research")
  expect_identical(token_research, token)

  expect_error(
    .dsvert_attest_opal_surface(
      list(), contract, set_option = function(...) NULL,
      get_options = function(...) data.frame(
        name = "dsvert.remote_surface_attestation", value = "tampered")),
    "did not persist")
})

test_that("installed scripts contain no legacy remote-disclosure shortcuts", {
  scripts <- system.file("scripts", package = "dsVertClient")
  if (!nzchar(scripts)) {
    scripts <- testthat::test_path("..", "..", "inst", "scripts")
  }
  files <- list.files(scripts, pattern = "[.]R$", full.names = TRUE)
  expect_true(length(files) > 0L)
  contents <- paste(unlist(lapply(
    files, readLines, warn = FALSE, encoding = "UTF-8")), collapse = "\n")
  forbidden <- c(
    "opal-demo",
    "ssl_verify(host|peer)\\s*=\\s*0",
    "opal\\.verifyssl\\s*=\\s*FALSE",
    "datashield\\.privacyLevel\\s*[,=]\\s*[\"']?0",
    "\\bn_common\\b",
    "(?:c|list|numeric|character)\\s*=\\s*base::",
    "\\buser\\s*=\\s*[\"'][^\"']+[\"']",
    "\\bpassword\\s*=\\s*[\"'][^\"']+[\"']",
    "permission\\s*=\\s*[\"'](?:view-values|edit|edit-values|administrate)[\"']")
  for (pattern in forbidden) {
    expect_false(
      grepl(pattern, contents, ignore.case = TRUE, perl = TRUE),
      info = paste("forbidden installed-script pattern:", pattern))
  }

  validation <- file.path(scripts, "validation.R")
  expressions <- parse(validation, keep.source = FALSE)
  expect_length(expressions, 1L)
  expect_identical(as.character(expressions[[1L]][[1L]]), "stop")
})

test_that("provisioner binds attestation to an exclusive staged surface", {
  script <- system.file(
    "scripts", "provision_opal.R", package = "dsVertClient")
  if (!nzchar(script)) {
    script <- testthat::test_path(
      "..", "..", "inst", "scripts", "provision_opal.R")
  }
  lines <- readLines(script, warn = FALSE, encoding = "UTF-8")
  text <- paste(lines, collapse = "\n")
  expect_match(text, 'PROFILE <- "dsvert"', fixed = TRUE)
  expect_match(
    text,
    'DSVERT_DESCRIPTION <- .tarball_description(TARBALLS$dsVert, "dsVert")',
    fixed = TRUE)
  expect_match(text, ".stage_tarballs(TARBALL_INPUTS)", fixed = TRUE)
  expect_match(text, ".reconcile_methods(", fixed = TRUE)
  expect_match(text, ".attest_surface(", fixed = TRUE)
  expect_match(text, "changed before enablement", fixed = TRUE)
  expect_match(text, "dsvert.remote_surface_attestation", fixed = TRUE)
  expect_match(text, "datashield.privacyLevel=5", fixed = TRUE)
  expect_match(text, "enabled = FALSE", fixed = TRUE)
  expect_match(text, "restricted = TRUE", fixed = TRUE)
  expect_match(text, "dsadmin.profile_perm_delete", fixed = TRUE)
  expect_match(text, "unexpected principals", fixed = TRUE)
  expect_match(text, "opal.table_perm_delete", fixed = TRUE)
  expect_match(text, 'permission = "view"', fixed = TRUE)
  expect_match(text, "unexpected or raw-value permissions", fixed = TRUE)
  expect_match(text, '.required_env("DSVERT_OPAL_RUNNER")', fixed = TRUE)
  expect_match(
    text, "DSVERT_OPAL_RUNNER must differ from the administrator account.",
    fixed = TRUE)
  expect_false(grepl('profile_name = "default"', text, fixed = TRUE))
  expect_false(grepl("try(dsadmin.profile_enable", text, fixed = TRUE))
  expect_false(grepl('permission = "view-values"', text, fixed = TRUE))

  profile <- grep("ensure_profile\\(o\\)", lines)
  project <- grep("ensure_project\\(o, PROJ\\)", lines)
  installation <- grep("install_packages\\(o\\)", lines)
  options <- grep("set_options\\(o\\)", lines)
  registration <- grep("register_methods\\(o\\)", lines)
  enabling <- grep("enable_profile\\(o\\)", lines)
  expect_length(profile, 1L)
  expect_length(project, 1L)
  expect_length(installation, 1L)
  expect_length(options, 1L)
  expect_length(registration, 1L)
  expect_length(enabling, 1L)
  expect_true(profile < project)
  expect_true(project < installation)
  expect_true(installation < options)
  expect_true(options < registration)
  expect_true(registration < enabling)
})

test_that("repository-level legacy entry points are fail-closed stubs", {
  repository <- normalizePath(
    testthat::test_path("..", "..", ".."), mustWork = FALSE)
  legacy <- file.path(repository, c(
    "register_methods.R",
    "scripts/opal_demo_provision_register.R",
    "demo_dsvert.R"))
  if (!all(file.exists(legacy))) {
    skip("repository-level scripts are not part of the installed package")
  }
  for (script in legacy) {
    expressions <- parse(script, keep.source = FALSE)
    expect_length(expressions, 1L)
    expect_identical(
      as.character(expressions[[1L]][[1L]]), "stop")
  }
})
