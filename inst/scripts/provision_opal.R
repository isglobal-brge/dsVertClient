#!/usr/bin/env Rscript
## provision_opal.R — out-of-band uploader for vert_validation_evidence.
##
## Idempotent. Probes opal-demo first; auto-falls back to local docker
## cluster (opal1+opal2+opal3) if /ws/datashield/packages?manager=local
## returns 403 (public-sandbox lock-down).
##
## Steps per server:
##   1. ensure project vert_demo
##   2. set datashield.privacyLevel=0 + dsvert.require_trusted_peers=FALSE
##   3. install dsVert + dsVertClient from local tarballs
##   4. register 140 DataSHIELD methods
##   5. upload all probe datasets (from generate_demo_datasets.R)
##   6. grant view permission to the runner account
##
## Access:
##   Rscript "$(Rscript -e 'cat(system.file(\"scripts/provision_opal.R\", package=\"dsVertClient\"))')"
##
## Env overrides (optional):
##   DSVERT_PROVISION_BACKEND   "remote" / "local" / "" (auto)
##   DSVERT_PROVISION_TARBALLS  path:path of dsVert / dsVertClient
##                              tarballs (defaults to repo build dir)

suppressPackageStartupMessages({ library(opalr) })
options(opal.verifyssl = FALSE)

`%||%` <- function(a, b) if (is.null(a) || !nzchar(a)) b else a

PROJ <- "vert_demo"

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

TARBALLS <- local({
  ovr <- Sys.getenv("DSVERT_PROVISION_TARBALLS", unset = "")
  if (nzchar(ovr)) {
    parts <- strsplit(ovr, ":", fixed = TRUE)[[1]]
    list(dsVert = parts[1], dsVertClient = parts[2])
  } else {
    r <- .repo_root()
    list(dsVert       = .latest_tarball(r, "dsVert"),
         dsVertClient = .latest_tarball(r, "dsVertClient"))
  }
})

.local_cluster <- list(
  s1 = list(url = "https://localhost:8443", user = "administrator",
             pass = "admin123"),
  s2 = list(url = "https://localhost:8444", user = "administrator",
             pass = "admin123"),
  s3 = list(url = "https://localhost:8445", user = "administrator",
             pass = "admin123"))
.demo_cluster <- list(default = list(url = "https://opal-demo.obiba.org",
                                       user = "administrator",
                                       pass = "password"))

login_for <- function(server_id, backend) {
  spec <- if (backend == "local") .local_cluster[[server_id]] else
    .demo_cluster$default
  opal.login(spec$user, spec$pass, url = spec$url,
             opts = list(ssl_verifyhost = 0, ssl_verifypeer = 0))
}

## ============================================================
## Backend probe — does opal-demo accept local-manager installs?
## ============================================================
.try_remote_install <- function() {
  o <- tryCatch(opal.login("administrator", "password",
                            url = "https://opal-demo.obiba.org",
                            opts = list(ssl_verifyhost = 0,
                                         ssl_verifypeer = 0)),
                error = function(e) NULL)
  if (is.null(o)) return(FALSE)
  on.exit(try(opal.logout(o), silent = TRUE))
  res <- tryCatch(
    dsadmin.install_local_package(o, TARBALLS$dsVert, profile = "default"),
    error = function(e) e)
  if (inherits(res, "error") &&
      grepl("403", conditionMessage(res), fixed = TRUE)) {
    cat("[backend] opal-demo /ws/datashield/packages locked (403); ",
        "falling back to local docker cluster.\n", sep = "")
    return(FALSE)
  }
  cat("[backend] opal-demo install endpoint accepted; staying remote.\n")
  TRUE
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

upload_table <- function(o, proj, table_name, df, runner = "administrator") {
  df$id <- sprintf("R%05d", seq_len(nrow(df)))
  opal.table_save(o, df, project = proj, table = table_name,
                  id.name = "id", overwrite = TRUE, force = TRUE)
  try(opal.table_perm_add(o, project = proj, table = table_name,
                           subject = runner, type = "USER",
                           permission = "view"), silent = TRUE)
  cat(sprintf("  [up] %-18s rows=%d cols=%s\n",
              table_name, nrow(df),
              paste(setdiff(names(df), "id"), collapse = ",")))
}

install_packages <- function(o, profile_name = "default") {
  for (pkg in c("dsVert", "dsVertClient")) {
    tryCatch(dsadmin.remove_package(o, pkg, profile = profile_name),
             error = function(e) NULL)
  }
  Sys.sleep(2)
  dsadmin.install_local_package(o, TARBALLS$dsVert, profile = profile_name)
  Sys.sleep(3)
  dsadmin.install_local_package(o, TARBALLS$dsVertClient, profile = profile_name)
  for (pkg in c("dsVert", "dsVertClient"))
    cat(sprintf("  [verify] %s installed: %s\n", pkg,
                dsadmin.installed_package(o, pkg, profile = profile_name)))
}

register_methods <- function(o, profile_name = "default") {
  desc_path <- system.file("DESCRIPTION", package = "dsVert")
  if (!nzchar(desc_path) || !file.exists(desc_path))
    desc_path <- file.path(.repo_root(), "dsVert", "DESCRIPTION")
  desc <- read.dcf(desc_path, all = TRUE)
  agg <- trimws(strsplit(gsub("\n", ",", desc$AggregateMethods), ",")[[1]])
  agg <- agg[nzchar(agg)]
  asn <- if (!is.null(desc$AssignMethods))
    trimws(strsplit(gsub("\n", ",", desc$AssignMethods), ",")[[1]]) else
    character(0)
  asn <- asn[nzchar(asn)]
  n <- 0L
  for (fn in agg) {
    if (grepl("=", fn)) {
      parts <- strsplit(fn, "=")[[1]]
      name <- trimws(parts[1]); func <- trimws(parts[2])
    } else { name <- fn; func <- paste0("dsVert::", fn) }
    tryCatch({
      dsadmin.set_method(o, name, func = func, type = "aggregate")
      n <- n + 1L
    }, error = function(e) cat(sprintf("  WARN agg %s: %s\n", name,
                                        conditionMessage(e))))
  }
  for (fn in asn) {
    tryCatch({
      dsadmin.set_method(o, fn, func = paste0("dsVert::", fn),
                          type = "assign")
      n <- n + 1L
    }, error = function(e) cat(sprintf("  WARN assign %s: %s\n", fn,
                                        conditionMessage(e))))
  }
  try(dsadmin.profile_enable(o, profile_name), silent = TRUE)
  cat(sprintf("  [register] %d methods (%d agg + %d assign)\n",
              n, length(agg), length(asn)))
}

set_options <- function(o) {
  try(dsadmin.set_option(o, "dsvert.require_trusted_peers", "FALSE"),
      silent = TRUE)
  try(dsadmin.set_option(o, "datashield.privacyLevel", "0"), silent = TRUE)
  cat("  [opt] dsvert.require_trusted_peers=FALSE, datashield.privacyLevel=0\n")
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
    o <- login_for(sid, backend)
    on.exit(try(opal.logout(o), silent = TRUE), add = TRUE)
    ensure_project(o, PROJ)
    set_options(o)
    install_packages(o)
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
  }
  cat("\n=== PROVISION COMPLETE ===\n")
}

## ============================================================
## main
## ============================================================
backend_env <- Sys.getenv("DSVERT_PROVISION_BACKEND", unset = "")
backend <- if (nzchar(backend_env)) backend_env else
  (if (.try_remote_install()) "remote" else "local")
PROVISION(backend)
