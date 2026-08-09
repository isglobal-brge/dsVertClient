#!/usr/bin/env Rscript
## provision_opal.R — out-of-band uploader for vert_validation_evidence.
##
## Idempotent and fail-closed. The target backend, administrator credentials
## and trusted HTTPS endpoints are supplied explicitly through environment
## variables; this script never disables certificate verification.
##
## Steps per server:
##   1. ensure project vert_demo and the restricted, dedicated dsvert profile
##   2. disable the profile and clear its previous surface attestation
##   3. install dsVert + dsVertClient from private copies of local tarballs
##   4. set/verify privacy level 5; peer pins remain administrator-owned
##   5. replace the full profile inventory with the exact dsVert allowlist
##   6. upload probe data and install exact non-value runner permissions
##   7. re-verify the Opal inventory, persist the connector-neutral custodian
##      surface token, then enable the profile
##
## Access:
##   Rscript "$(Rscript -e 'cat(system.file(\"scripts/provision_opal.R\", package=\"dsVertClient\"))')"
##
## Inputs:
##   DSVERT_PROVISION_BACKEND   "remote" / "local" (required)
##   DSVERT_PROVISION_TARBALLS  path:path of dsVert / dsVertClient
##                              tarballs (defaults to repo build dir)
##   DSVERT_OPAL_USER           administrator account (required)
##   DSVERT_OPAL_PASSWORD       administrator secret (required)
##   DSVERT_OPAL_RUNNER         distinct DataSHIELD runner account (required)
##   DSVERT_OPAL_URL            remote HTTPS URL when backend=remote
##
## This validation provisioner grants the runner only `use` on the restricted
## dsvert profile and Opal's `view` table permission (dictionary and summaries,
## DataSHIELD-capable, but no individual-value queries). It never grants
## `view-values`, edit or administrative table access.

suppressPackageStartupMessages({ library(opalr) })

PROJ <- "vert_demo"
PROFILE <- "dsvert"

## Tarball locations: prefer DSVERT_PROVISION_TARBALLS else latest repo builds.
.repo_root <- function() {
  r <- Sys.getenv("DSVERT_REPO_ROOT", unset = "")
  if (nzchar(r)) return(r)
  here <- normalizePath(getwd(), mustWork = FALSE)
  ## Walk up looking for the dsvert-paper repo signature.
  while (here != "/" && here != "") {
    if (dir.exists(file.path(here, "dsVert")) &&
        dir.exists(file.path(here, "dsVertClient"))) return(here)
    here <- dirname(here)
  }
  "."
}

.latest_tarball <- function(root, pkg) {
  hits <- list.files(root, pattern = paste0("^", pkg, "_.*[.]tar[.]gz$"),
                     full.names = TRUE)
  if (!length(hits)) {
    stop("Cannot locate a built ", pkg, " tarball under ", root,
         ". Run R CMD build first or set DSVERT_PROVISION_TARBALLS.",
         call. = FALSE)
  }
  hits[order(file.info(hits)$mtime, decreasing = TRUE)][1L]
}

TARBALL_INPUTS <- local({
  ovr <- Sys.getenv("DSVERT_PROVISION_TARBALLS", unset = "")
  if (nzchar(ovr)) {
    parts <- strsplit(ovr, ":", fixed = TRUE)[[1]]
    if (length(parts) != 2L || any(!nzchar(parts))) {
      stop(
        "DSVERT_PROVISION_TARBALLS must contain exactly two paths: ",
        "dsVert:dsVertClient.", call. = FALSE)
    }
    list(dsVert = parts[1], dsVertClient = parts[2])
  } else {
    r <- .repo_root()
    list(dsVert       = .latest_tarball(r, "dsVert"),
         dsVertClient = .latest_tarball(r, "dsVertClient"))
  }
})

.stage_tarballs <- function(tarballs) {
  directory <- tempfile("dsvert-provision-artifacts-")
  if (!dir.create(directory, mode = "0700")) {
    stop("Cannot create the private tarball staging directory.",
         call. = FALSE)
  }
  staged <- lapply(names(tarballs), function(package) {
    source <- tarballs[[package]]
    if (!is.character(source) || length(source) != 1L || is.na(source) ||
        !file.exists(source) || dir.exists(source)) {
      stop("Missing ", package, " tarball.", call. = FALSE)
    }
    target <- file.path(directory, paste0(package, ".tar.gz"))
    if (!file.copy(source, target, overwrite = FALSE, copy.mode = FALSE)) {
      stop("Cannot stage the ", package, " tarball.", call. = FALSE)
    }
    Sys.chmod(target, mode = "0600")
    target
  })
  names(staged) <- names(tarballs)
  staged
}

# Private copies close the inspection/upload time-of-check/time-of-use gap.
TARBALLS <- .stage_tarballs(TARBALL_INPUTS)

.tarball_description <- function(tarball, package) {
  if (!is.character(tarball) || length(tarball) != 1L || is.na(tarball) ||
      !file.exists(tarball)) {
    stop("Missing ", package, " tarball.", call. = FALSE)
  }
  members <- utils::untar(tarball, list = TRUE)
  expected <- paste0(package, "/DESCRIPTION")
  if (!identical(members[members == expected], expected)) {
    stop("The ", package, " tarball has no canonical DESCRIPTION entry.",
         call. = FALSE)
  }
  directory <- tempfile(paste0(package, "-description-"))
  if (!dir.create(directory, mode = "0700")) {
    stop("Cannot create the temporary tarball inspection directory.",
         call. = FALSE)
  }
  utils::untar(tarball, files = expected, exdir = directory)
  path <- file.path(directory, expected)
  if (!file.exists(path)) {
    stop("Cannot read DESCRIPTION from the uploaded tarball.",
         call. = FALSE)
  }
  description <- tryCatch(
    read.dcf(path),
    error = function(e) stop(
      "Cannot parse DESCRIPTION from the ", package, " tarball: ",
      conditionMessage(e), call. = FALSE))
  if (nrow(description) != 1L || !"Package" %in% colnames(description) ||
      is.na(description[1L, "Package"]) ||
      !identical(unname(trimws(description[1L, "Package"])), package)) {
    stop("The staged tarball does not contain package ", package, ".",
         call. = FALSE)
  }
  path
}

# Registration is derived from the same immutable server tarball that is
# uploaded below, never from an unrelated locally installed dsVert version.
DSVERT_DESCRIPTION <- .tarball_description(TARBALLS$dsVert, "dsVert")
invisible(.tarball_description(TARBALLS$dsVertClient, "dsVertClient"))

.description_import_names <- function(path) {
  description <- tryCatch(
    read.dcf(path, fields = "Imports"),
    error = function(e) stop(
      "Cannot parse server Imports: ", conditionMessage(e),
      call. = FALSE))
  if (nrow(description) != 1L || ncol(description) != 1L) {
    stop("The server tarball has an invalid Imports field.", call. = FALSE)
  }
  value <- description[1L, 1L]
  if (is.na(value) || !nzchar(trimws(value))) return(character())
  specifications <- trimws(strsplit(value, ",", fixed = TRUE)[[1L]])
  valid <- grepl(
    paste0(
      "^[[:alpha:]][[:alnum:].]*",
      "([[:space:]]*\\((>=|<=|==|>|<)[[:space:]]*",
      "[[:digit:]][[:alnum:].-]*\\))?$"),
    specifications)
  if (!length(specifications) || any(!nzchar(specifications)) || any(!valid)) {
    stop("The server tarball has an invalid Imports field.", call. = FALSE)
  }
  packages <- sub(
    "^([[:alpha:]][[:alnum:].]*).*$", "\\1", specifications)
  if (anyDuplicated(packages)) {
    stop("The server tarball has duplicate Imports entries.", call. = FALSE)
  }
  packages
}

DSVERT_IMPORTS <- .description_import_names(DSVERT_DESCRIPTION)

# Validate registration before opening a remote session. Registration is
# derived from the DESCRIPTION inside the staged dsVert tarball uploaded below.
.method_contract <- getFromNamespace(
  ".dsvert_opal_method_contract", "dsVertClient")
.reconcile_methods <- getFromNamespace(
  ".dsvert_reconcile_opal_methods", "dsVertClient")
.attest_surface <- getFromNamespace(
  ".dsvert_attest_opal_surface", "dsVertClient")
DSVERT_CONTRACT <- .method_contract(DSVERT_DESCRIPTION)

.required_env <- function(name) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) {
    stop("Missing required environment variable ", name, ".",
         call. = FALSE)
  }
  value
}

.validated_https_url <- function(value, what) {
  parsed <- if (is.character(value) && length(value) == 1L && !is.na(value)) {
    tryCatch(httr::parse_url(value), error = function(e) NULL)
  } else {
    NULL
  }
  if (is.null(parsed) || !identical(tolower(parsed$scheme), "https") ||
      is.null(parsed$hostname) || !nzchar(parsed$hostname) ||
      !is.null(parsed$username) || !is.null(parsed$password) ||
      grepl("[[:space:]]", value)) {
    stop(what, " must be an HTTPS URL.", call. = FALSE)
  }
  value
}

.opal_user <- .required_env("DSVERT_OPAL_USER")
.opal_password <- .required_env("DSVERT_OPAL_PASSWORD")
.opal_runner <- .required_env("DSVERT_OPAL_RUNNER")
if (identical(.opal_runner, .opal_user)) {
  stop(
    "DSVERT_OPAL_RUNNER must differ from the administrator account.",
    call. = FALSE)
}
.local_cluster <- list(
  s1 = list(url = Sys.getenv("DSVERT_OPAL_S1_URL",
                            "https://localhost:8443")),
  s2 = list(url = Sys.getenv("DSVERT_OPAL_S2_URL",
                            "https://localhost:8444")),
  s3 = list(url = Sys.getenv("DSVERT_OPAL_S3_URL",
                            "https://localhost:8445")))
.remote_cluster <- list(default = list(
  url = Sys.getenv("DSVERT_OPAL_URL", unset = "")))

login_for <- function(server_id, backend) {
  spec <- if (backend == "local") .local_cluster[[server_id]] else
    .remote_cluster$default
  url <- .validated_https_url(spec$url, paste0("Opal URL for ", server_id))
  opal.login(.opal_user, .opal_password, url = url)
}

with_opal_connection <- function(server_id, backend, code) {
  o <- login_for(server_id, backend)
  on.exit(try(opal.logout(o), silent = TRUE), add = TRUE)
  code(o)
}

## ============================================================
## Per-cluster operations (idempotent)
## ============================================================
ensure_project <- function(o, proj) {
  if (!proj %in% opal.projects(o)$name) {
    cat(sprintf("[create] project %s\n", proj))
    opal.project_create(o, project = proj, database = "mongodb",
      title = "dsVert validation cluster",
      description = "vert_validation_evidence K=2/K=3 probe datasets")
  } else cat(sprintf("[skip]  project %s exists\n", proj))
}

.validated_profile_permissions <- function(permissions) {
  required <- c("subject", "type", "permission")
  if (is.data.frame(permissions) && nrow(permissions) == 0L &&
      length(names(permissions)) == 0L) {
    permissions <- data.frame(
      subject = character(), type = character(), permission = character(),
      stringsAsFactors = FALSE)
  }
  if (!is.data.frame(permissions) ||
      !identical(names(permissions), required) || anyNA(permissions)) {
    stop("Opal returned an invalid profile-permission inventory.",
         call. = FALSE)
  }
  permissions
}

ensure_profile <- function(o, profile_name = PROFILE, runner = .opal_runner) {
  if (!dsadmin.profile_exists(o, profile_name)) {
    dsadmin.profile_create(o, profile_name)
    cat(sprintf("  [create] exclusive profile %s\n", profile_name))
  }
  # The profile remains unavailable while packages, options and its exact
  # method inventory are being reconciled.
  dsadmin.profile_enable(o, profile_name, enabled = FALSE)
  dsadmin.profile_access(o, profile_name, restricted = TRUE)
  permissions <- .validated_profile_permissions(
    dsadmin.profile_perm(o, profile_name))
  principals <- unique(permissions[, c("subject", "type"), drop = FALSE])
  for (row in seq_len(nrow(principals))) {
    dsadmin.profile_perm_delete(
      o, profile_name, subject = principals$subject[[row]],
      type = principals$type[[row]])
  }
  dsadmin.profile_perm_add(
    o, profile_name, subject = runner, type = "user", permission = "use")
  permissions <- dsadmin.profile_perm(o, profile_name)
  expected_permission <- data.frame(
    subject = runner, type = "user", permission = "use",
    stringsAsFactors = FALSE)
  rownames(permissions) <- NULL
  if (!identical(permissions, expected_permission)) {
    stop("The dedicated dsVert profile has unexpected principals.",
         call. = FALSE)
  }
  profile <- dsadmin.profile(o, profile_name)
  if (!identical(profile$enabled, FALSE) ||
      !identical(profile$restrictedAccess, TRUE)) {
    stop("The dedicated dsVert profile is not disabled and restricted.",
         call. = FALSE)
  }

  option_name <- "dsvert.remote_surface_attestation"
  current <- dsadmin.get_options(o, profile = profile_name)
  if (!is.data.frame(current) ||
      (nrow(current) && !all(c("name", "value") %in% names(current)))) {
    stop("Opal returned an invalid option inventory.", call. = FALSE)
  }
  if (nrow(current) && option_name %in% as.character(current$name)) {
    dsadmin.rm_option(o, option_name, profile = profile_name)
  }
  after <- dsadmin.get_options(o, profile = profile_name)
  if (!is.data.frame(after) ||
      (nrow(after) && (!"name" %in% names(after) ||
                       option_name %in% as.character(after$name)))) {
    stop("Cannot clear the previous dsVert surface attestation.",
         call. = FALSE)
  }
}

upload_table <- function(o, proj, table_name, df, runner = .opal_runner) {
  df$id <- sprintf("R%05d", seq_len(nrow(df)))
  opal.table_save(o, df, project = proj, table = table_name,
                  id.name = "id", overwrite = TRUE, force = TRUE)

  # Opal documents `view` as dictionary/summaries access with DataSHIELD but
  # without individual-value queries. Purge inherited explicit ACLs and grant
  # exactly that permission to the non-administrator runner.
  permissions <- opal.table_perm(o, project = proj, table = table_name)
  if (!is.data.frame(permissions) ||
      !all(c("subject", "type", "permission") %in% names(permissions)) ||
      anyNA(permissions)) {
    stop("Opal returned an invalid table-permission inventory.",
         call. = FALSE)
  }
  principals <- unique(permissions[, c("subject", "type"), drop = FALSE])
  for (row in seq_len(nrow(principals))) {
    opal.table_perm_delete(
      o, project = proj, table = table_name,
      subject = principals$subject[[row]], type = principals$type[[row]])
  }
  opal.table_perm_add(
    o, project = proj, table = table_name,
    subject = runner, type = "user", permission = "view")
  permissions <- opal.table_perm(o, project = proj, table = table_name)
  expected_permission <- data.frame(
    subject = runner, type = "user", permission = "view",
    stringsAsFactors = FALSE)
  rownames(permissions) <- NULL
  if (!identical(permissions, expected_permission)) {
    stop("The validation table has unexpected or raw-value permissions.",
         call. = FALSE)
  }
  cat(sprintf("  [up] %-18s rows=%d cols=%s\n",
              table_name, nrow(df),
              paste(setdiff(names(df), "id"), collapse = ",")))
}

install_server_imports <- function(o, packages, cluster_name) {
  for (pkg in packages) {
    installed <- tryCatch(
      oadmin.installed_package(o, pkg, profile = cluster_name),
      error = function(e) stop(
        "Cannot inspect required server dependency ", pkg, ": ",
        conditionMessage(e), call. = FALSE))
    if (!is.logical(installed) || length(installed) != 1L ||
        is.na(installed)) {
      stop("Invalid Rock inventory state for required server dependency ",
           pkg, ".", call. = FALSE)
    }
    if (installed) next
    # The repository used here is the CRAN repository configured by the
    # custodian for this Rock profile; the provisioner never overrides it.
    tryCatch(
      oadmin.install_package(o, pkg, profile = cluster_name),
      error = function(e) stop(
        "Cannot install required server dependency ", pkg, ": ",
        conditionMessage(e), call. = FALSE))
    verified <- tryCatch(
      oadmin.installed_package(o, pkg, profile = cluster_name),
      error = function(e) stop(
        "Cannot verify required server dependency ", pkg, ": ",
        conditionMessage(e), call. = FALSE))
    if (!identical(verified, TRUE)) {
      stop("Rock did not install required server dependency ", pkg, ".",
           call. = FALSE)
    }
    cat(sprintf("  [verify] dependency %s installed\n", pkg))
  }
  invisible(TRUE)
}

.profile_rock_cluster <- function(o, profile_name) {
  profile <- tryCatch(
    dsadmin.profile(o, profile_name),
    error = function(e) stop(
      "Cannot resolve the Rock cluster for profile ", profile_name, ": ",
      conditionMessage(e), call. = FALSE))
  if (!is.list(profile) ||
      !is.character(profile$cluster) || length(profile$cluster) != 1L ||
      is.na(profile$cluster) || !nzchar(profile$cluster) ||
      !grepl("^[A-Za-z0-9._-]{1,64}$", profile$cluster)) {
    stop("Opal returned an invalid Rock cluster for profile ", profile_name,
         ".", call. = FALSE)
  }
  profile$cluster
}

install_packages <- function(o, profile_name = PROFILE) {
  cluster_name <- .profile_rock_cluster(o, profile_name)
  # R packages are cluster-scoped. The exact callable isolation below is the
  # DataSHIELD profile method allowlist, not a private package library claim.
  for (pkg in c("dsVert", "dsVertClient")) {
    # The DataSHIELD package endpoint may answer HTTP 200 with an empty body
    # after removal. Query the profile's real Rock inventory instead.
    installed <- oadmin.installed_package(o, pkg, profile = cluster_name)
    if (!is.logical(installed) || length(installed) != 1L ||
        is.na(installed)) {
      stop("Cannot determine whether ", pkg, " is installed.",
           call. = FALSE)
    }
    if (!installed) next
    tryCatch(
      dsadmin.rm_package_methods(o, pkg, profile = profile_name),
      error = function(e) stop(
        "Cannot remove stale ", pkg, " methods: ", conditionMessage(e),
        call. = FALSE))
    tryCatch(
      dsadmin.remove_package(o, pkg, profile = profile_name),
      error = function(e) stop(
        "Cannot unpublish stale ", pkg, " from DataSHIELD profile ",
        profile_name, ": ", conditionMessage(e), call. = FALSE))
    tryCatch(
      oadmin.remove_package(o, pkg, profile = cluster_name),
      error = function(e) stop(
        "Cannot remove stale ", pkg, " from Rock cluster ", cluster_name,
        ": ", conditionMessage(e), call. = FALSE))
    if (!identical(
        oadmin.installed_package(o, pkg, profile = cluster_name), FALSE)) {
      stop("Opal did not remove the previous ", pkg, " package.",
           call. = FALSE)
    }
  }
  install_server_imports(o, DSVERT_IMPORTS, cluster_name)
  dsadmin.install_local_package(o, TARBALLS$dsVert, profile = profile_name)
  dsadmin.install_local_package(o, TARBALLS$dsVertClient, profile = profile_name)
  for (pkg in c("dsVert", "dsVertClient")) {
    if (!identical(
        oadmin.installed_package(o, pkg, profile = cluster_name), TRUE)) {
      stop("Opal did not install ", pkg, " in profile ", profile_name, ".",
           call. = FALSE)
    }
    cat(sprintf("  [verify] %s installed\n", pkg))
  }
}

register_methods <- function(o, profile_name = PROFILE) {
  contract <- .reconcile_methods(
    o, description_path = DSVERT_DESCRIPTION,
    profile_name = profile_name)
  if (!identical(contract, DSVERT_CONTRACT)) {
    stop("The staged dsVert method contract changed during provisioning.",
         call. = FALSE)
  }
  cat(sprintf("  [register] %d methods (%d agg + %d assign)\n",
              nrow(contract), sum(contract$type == "aggregate"),
              sum(contract$type == "assign")))
}

attest_methods <- function(o, profile_name = PROFILE) {
  contract <- .reconcile_methods(
    o, description_path = DSVERT_DESCRIPTION,
    profile_name = profile_name)
  if (!identical(contract, DSVERT_CONTRACT)) {
    stop("The staged dsVert method contract changed before attestation.",
         call. = FALSE)
  }
  token <- .attest_surface(
    o, contract = contract, profile_name = profile_name)
  verified <- .reconcile_methods(
    o, description_path = DSVERT_DESCRIPTION,
    profile_name = profile_name)
  if (!identical(verified, contract)) {
    stop("The attested dsVert method contract changed before enablement.",
         call. = FALSE)
  }
  cat(sprintf("  [attest] %s\n", token))
}

set_options <- function(o, profile_name = PROFILE) {
  dsadmin.set_option(
    o, "datashield.privacyLevel", "5", profile = profile_name)
  options <- dsadmin.get_options(o, profile = profile_name)
  matches <- if (is.data.frame(options) &&
                 all(c("name", "value") %in% names(options))) {
    which(as.character(options$name) == "datashield.privacyLevel")
  } else {
    integer()
  }
  if (length(matches) != 1L ||
      !identical(as.character(options$value[[matches]]), "5")) {
    stop("Opal did not persist datashield.privacyLevel=5.",
         call. = FALSE)
  }
  cat("  [opt] datashield.privacyLevel=5; mandatory peer pins unchanged\n")
}

enable_profile <- function(o, profile_name = PROFILE) {
  dsadmin.profile_enable(o, profile_name, enabled = TRUE)
  profile <- dsadmin.profile(o, profile_name)
  if (!identical(profile$enabled, TRUE) ||
      !identical(profile$restrictedAccess, TRUE)) {
    stop("The dedicated dsVert profile was not safely enabled.",
         call. = FALSE)
  }
  cat(sprintf("  [enable] exclusive profile %s\n", profile_name))
}

## ============================================================
## Builders — sourced from generate_demo_datasets.R sibling.
## ============================================================
.find_gen_script <- function() {
  ## (1) Sibling next to this script when run via Rscript.
  cargs <- commandArgs(trailingOnly = FALSE)
  hit <- grep("--file=", cargs, fixed = TRUE, value = TRUE)
  if (length(hit) > 0L) {
    self <- sub("^--file=", "", hit[1L])
    cand <- file.path(dirname(normalizePath(self)),
                       "generate_demo_datasets.R")
    if (file.exists(cand)) return(cand)
  }
  ## (2) Package-installed path.
  cand <- system.file("scripts/generate_demo_datasets.R",
                       package = "dsVertClient")
  if (nzchar(cand) && file.exists(cand)) return(cand)
  ## (3) Repo development path.
  cand <- file.path(.repo_root(), "dsVertClient", "inst", "scripts",
                     "generate_demo_datasets.R")
  if (file.exists(cand)) return(cand)
  stop("Cannot locate generate_demo_datasets.R sibling")
}
source(.find_gen_script(), local = TRUE, chdir = TRUE)

## ============================================================
## Backend dispatch
## ============================================================
PROVISION <- function(backend) {
  cat(sprintf("\n=== PROVISION backend=%s ===\n", backend))
  servers <- if (backend == "local") names(.local_cluster) else "default"
  bw     <- build_birthwt_k2()
  pim2   <- build_pima_k2()
  lng    <- build_lung_k2()
  glm_l  <- build_long_glmm_k2()
  lmm_l  <- build_long_lmm_k2()
  ipw_l  <- build_cgd_ipw_k2()
  pim3   <- build_pima_k3()
  warm3  <- build_pima_warm_k3()

  for (sid in servers) {
    cat(sprintf("\n[server %s]\n", sid))
    with_opal_connection(sid, backend, function(o) {
    ensure_profile(o)
    ensure_project(o, PROJ)
    install_packages(o)
    set_options(o)
    register_methods(o)
    if (backend == "local") {
      uploads <- switch(sid,
        s1 = list(bw_hospital_a = bw$bw_hospital_a,
                   na_data_s1   = pim2$na_data_s1,
                   lung_s1      = lng$lung_s1,
                   long_glmm_s1 = glm_l$long_glmm_s1,
                   long_lmm_s1  = lmm_l$long_lmm_s1,
                   cgd_ipw_s1   = ipw_l$cgd_ipw_s1,
                   pima_server1 = pim3$pima_server1,
                   pima_warm_s1 = warm3$pima_warm_s1),
        s2 = list(bw_hospital_b = bw$bw_hospital_b,
                   na_data_s2   = pim2$na_data_s2,
                   lung_s2      = lng$lung_s2,
                   long_glmm_s2 = glm_l$long_glmm_s2,
                   long_lmm_s2  = lmm_l$long_lmm_s2,
                   cgd_ipw_s2   = ipw_l$cgd_ipw_s2,
                   pima_server2 = pim3$pima_server2,
                   pima_warm_s2 = warm3$pima_warm_s2),
        s3 = list(pima_server3 = pim3$pima_server3,
                   pima_warm_s3 = warm3$pima_warm_s3))
      for (nm in names(uploads))
        upload_table(o, PROJ, nm, uploads[[nm]])
    } else {
      all_tabs <- list(bw_hospital_a = bw$bw_hospital_a,
                        bw_hospital_b = bw$bw_hospital_b,
                        na_data_s1    = pim2$na_data_s1,
                        na_data_s2    = pim2$na_data_s2,
                        lung_s1       = lng$lung_s1,
                        lung_s2       = lng$lung_s2,
                        long_glmm_s1  = glm_l$long_glmm_s1,
                        long_glmm_s2  = glm_l$long_glmm_s2,
                        long_lmm_s1   = lmm_l$long_lmm_s1,
                        long_lmm_s2   = lmm_l$long_lmm_s2,
                        cgd_ipw_s1    = ipw_l$cgd_ipw_s1,
                        cgd_ipw_s2    = ipw_l$cgd_ipw_s2,
                        pima_server1  = pim3$pima_server1,
                        pima_server2  = pim3$pima_server2,
                        pima_server3  = pim3$pima_server3,
                        pima_warm_s1  = warm3$pima_warm_s1,
                        pima_warm_s2  = warm3$pima_warm_s2,
                        pima_warm_s3  = warm3$pima_warm_s3)
      for (nm in names(all_tabs))
        upload_table(o, PROJ, nm, all_tabs[[nm]])
    }
    attest_methods(o)
    enable_profile(o)
    })
  }
  cat(paste0(
    "\n=== VALIDATION SURFACE PROVISIONED ===\n",
    "The runner has only DataSHIELD profile use plus dictionary/summary ",
    "table view. It cannot query, read or export individual row values.\n"))
}

## ============================================================
## main
## ============================================================
backend_env <- Sys.getenv("DSVERT_PROVISION_BACKEND", unset = "")
if (!backend_env %in% c("local", "remote")) {
  stop("DSVERT_PROVISION_BACKEND must be explicitly set to local or remote.",
       call. = FALSE)
}
backend <- backend_env
PROVISION(backend)
