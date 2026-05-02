#!/usr/bin/env Rscript
## generate_demo_datasets.R — seed-documented builders for vert_demo.
##
## Reproducible recipe for the K=2 / K=3 splits consumed by
## vert_validation_evidence. Every builder is a pure function of its
## seed; running it locally rebuilds identical data to what is loaded
## on the validation cluster. Sourced both from
## dsVertClient/vignettes/data-preparation.Rmd (documentation /
## pkgdown article) and from provision_opal.R (out-of-band uploader).
##
## Access from a user session (after dsVertClient install):
##   source(system.file("scripts/generate_demo_datasets.R",
##                       package = "dsVertClient"))
##   bw <- build_birthwt_k2()         # default seed 123
##
## All builders return a named list of data frames (one per server
## table). Synthetic longitudinal builders also return a `diagnostics`
## slot reporting the realised p̄ / σ_b² so users can verify the well-
## identified regime.

suppressPackageStartupMessages({
  library(MASS); library(survival)
})

## ---------------------------------------------------------------
## K=2 splits
## ---------------------------------------------------------------

#' birthwt K=2 split (`bw_hospital_a` / `bw_hospital_b`).
#' Used by the GLM K=2 probe.
build_birthwt_k2 <- function(seed = 123L) {
  data(birthwt, package = "MASS")
  bdf <- birthwt
  bdf$patient_id <- sprintf("BW%03d", seq_len(nrow(bdf)))
  set.seed(seed)
  core <- sample(bdf$patient_id, 160)
  rem  <- setdiff(bdf$patient_id, core)
  ids1 <- sample(c(core, sample(rem, 15)))
  ids2 <- sample(c(core, sample(rem, 14)))
  list(
    bw_hospital_a = bdf[bdf$patient_id %in% ids1,
                        c("patient_id", "age", "lwt", "race", "smoke")],
    bw_hospital_b = bdf[bdf$patient_id %in% ids2,
                        c("patient_id", "ptl", "ht", "ui", "ftv",
                          "bwt", "low")])
}

#' NHANES Pima K=2 split + bp-tertile factor + indicator + cumulative
#' columns on s2. Used by NB, ord_joint, mnl_joint, ordinal_warm K=2,
#' multinom_warm K=2, lasso K=2.
build_pima_k2 <- function(seed = 123L, bp_seed = 124L) {
  pdf <- rbind(MASS::Pima.tr, MASS::Pima.te)
  pdf$diabetes <- as.integer(pdf$type == "Yes")  # MASS uses `type` factor
  pdf$patient_id <- sprintf("P%03d", seq_len(nrow(pdf)))
  set.seed(seed)
  core <- sample(pdf$patient_id, 150)
  rem  <- setdiff(pdf$patient_id, core)
  ids1 <- sample(c(core, sample(rem, 25)))
  ids2 <- sample(c(core, sample(rem, 20)))
  s1 <- pdf[pdf$patient_id %in% ids1,
             c("patient_id", "age", "bmi", "ped", "npreg", "skin")]
  s2 <- pdf[pdf$patient_id %in% ids2,
             c("patient_id", "glu", "bp", "diabetes")]
  set.seed(bp_seed)
  brks <- quantile(s2$bp, c(0, 1/3, 2/3, 1), na.rm = TRUE)
  s2$bp_cls <- factor(
    cut(s2$bp, breaks = brks, include.lowest = TRUE,
        labels = c("low", "med", "high")),
    levels = c("low", "med", "high"), ordered = TRUE)
  s2$low_ind  <- as.integer(s2$bp_cls == "low")
  s2$med_ind  <- as.integer(s2$bp_cls == "med")
  s2$high_ind <- as.integer(s2$bp_cls == "high")
  s2$low_leq  <- as.integer(as.integer(s2$bp_cls) <= 1L)
  s2$med_leq  <- as.integer(as.integer(s2$bp_cls) <= 2L)
  s2$high_leq <- 1L
  list(na_data_s1 = s1, na_data_s2 = s2)
}

#' survival::lung K=2 split. Used by Cox K=2.
build_lung_k2 <- function(seed = 125L) {
  data(lung, package = "survival")
  lng <- get("lung")
  lng$patient_id <- sprintf("LNG%03d", seq_len(nrow(lng)))
  keep <- complete.cases(lng[, c("age", "sex", "ph.ecog", "ph.karno",
                                  "time", "status")])
  lng <- lng[keep, ]
  set.seed(seed)
  core <- sample(lng$patient_id, 150)
  rem  <- setdiff(lng$patient_id, core)
  ids1 <- sample(c(core, sample(rem, 20)))
  ids2 <- sample(c(core, sample(rem, 15)))
  s1 <- lng[lng$patient_id %in% ids1,
             c("patient_id", "age", "sex", "ph.ecog", "ph.karno")]
  s2 <- lng[lng$patient_id %in% ids2,
             c("patient_id", "time", "status", "meal.cal", "wt.loss")]
  s2$status <- as.integer(s2$status == 2)
  list(lung_s1 = s1, lung_s2 = s2)
}

#' Synthetic longitudinal binomial GLMM. K = n_obs observations per
#' cluster across n_clusters clusters. Intercept is tuned so p̄ ≈ 0.3
#' (well-identified regime — far from both 0.5 and the rare-event tail
#' where the simplified Laplace BLUP becomes biased per Breslow &
#' Clayton 1993 §3).
build_long_glmm_k2 <- function(seed = 200L,
                                n_clusters = 200L, n_obs = 10L,
                                sigma_b = sqrt(0.5),
                                beta_intercept = -0.85,
                                beta_x1 = 0.5, beta_x2 = 0.3) {
  set.seed(seed)
  cluster_id  <- rep(sprintf("C%04d", seq_len(n_clusters)), each = n_obs)
  obs_id      <- sprintf("O%05d", seq_along(cluster_id))
  b_i         <- rnorm(n_clusters, 0, sigma_b)
  b           <- rep(b_i, each = n_obs)
  x1          <- rnorm(length(cluster_id))
  x2          <- rbinom(length(cluster_id), 1, 0.5)
  eta         <- beta_intercept + beta_x1 * x1 + beta_x2 * x2 + b
  p           <- plogis(eta)
  y           <- rbinom(length(cluster_id), 1, p)
  list(long_glmm_s1 = data.frame(patient_id = obs_id,
                                   x1 = x1, anchor = 1L),
       long_glmm_s2 = data.frame(patient_id = obs_id,
                                   cluster_id = cluster_id,
                                   x2 = x2, y = y),
       diagnostics = list(p_bar = mean(p), y_rate = mean(y),
                          sigma_b2 = sigma_b^2,
                          n = length(cluster_id),
                          n_clusters = n_clusters))
}

#' Synthetic longitudinal Gaussian LMM/GEE. Same cluster structure as
#' long_glmm but linear DGP — used by both ds.vertLMM (REML) and
#' ds.vertGEE (sandwich) probes.
build_long_lmm_k2 <- function(seed = 201L,
                               n_clusters = 200L, n_obs = 10L,
                               sigma_b = sqrt(0.5), sigma_e = 1,
                               beta_intercept = 5,
                               beta_x1 = 0.5, beta_x2 = 0.3) {
  set.seed(seed)
  cluster_id <- rep(sprintf("C%04d", seq_len(n_clusters)), each = n_obs)
  obs_id     <- sprintf("O%05d", seq_along(cluster_id))
  b_i        <- rnorm(n_clusters, 0, sigma_b)
  b          <- rep(b_i, each = n_obs)
  x1         <- rnorm(length(cluster_id))
  x2         <- rbinom(length(cluster_id), 1, 0.5)
  y          <- beta_intercept + beta_x1 * x1 + beta_x2 * x2 + b +
                rnorm(length(cluster_id), 0, sigma_e)
  list(long_lmm_s1 = data.frame(patient_id = obs_id,
                                  x1 = x1, anchor = 1L),
       long_lmm_s2 = data.frame(patient_id = obs_id,
                                  cluster_id = cluster_id,
                                  x2 = x2, y = y),
       diagnostics = list(sigma_b2 = sigma_b^2, sigma_e = sigma_e,
                          n = length(cluster_id),
                          n_clusters = n_clusters))
}

#' CGD IPW K=2 — enum=1 single row per subject + propensity-derived
#' weights pre-computed centrally for stage-2 reproducibility.
build_cgd_ipw_k2 <- function(seed_s1 = 128L, seed_s2 = 129L) {
  data(cgd, package = "survival")
  cg <- cgd[cgd$enum == 1L, ]
  cg$patient_id   <- sprintf("CGD%03d", cg$id)
  cg$sex_num      <- as.integer(cg$sex == "male")
  cg$treat_num    <- as.integer(cg$treat == "rIFN-g")
  n_inf <- tapply(cgd$status, cgd$id, sum)
  cg$n_infections <- as.integer(n_inf[as.character(cg$id)])
  pooled <- cg[, c("patient_id", "age", "sex_num", "height", "weight",
                   "treat_num", "status", "n_infections")]
  pooled <- pooled[complete.cases(pooled), ]
  prop <- glm(treat_num ~ age + sex_num + height + weight,
              data = pooled, family = binomial())
  phat <- fitted(prop)
  pooled$ipw <- ifelse(pooled$treat_num == 1L, 1/phat, 1/(1 - phat))
  set.seed(seed_s1); ids1 <- sample(pooled$patient_id)
  set.seed(seed_s2); ids2 <- sample(pooled$patient_id)
  s1 <- pooled[match(ids1, pooled$patient_id),
                c("patient_id", "age", "sex_num", "height", "weight")]
  s2 <- pooled[match(ids2, pooled$patient_id),
                c("patient_id", "treat_num", "status",
                  "n_infections", "ipw")]
  list(cgd_ipw_s1 = s1, cgd_ipw_s2 = s2,
       propensity_coef_central = coef(prop))
}

## ---------------------------------------------------------------
## K=3 splits
## ---------------------------------------------------------------

#' Pima K=3 split (`pima_server1..3`). Used by GLM K=3 + lasso K=3.
build_pima_k3 <- function(seed = 123L) {
  pdf <- rbind(MASS::Pima.tr, MASS::Pima.te)
  pdf$diabetes <- as.integer(pdf$type == "Yes")
  pdf$patient_id <- sprintf("P%03d", seq_len(nrow(pdf)))
  set.seed(seed)
  core <- sample(pdf$patient_id, 150)
  rem  <- setdiff(pdf$patient_id, core)
  ids1 <- sample(c(core, sample(rem, 30)))
  ids2 <- sample(c(core, sample(rem, 25)))
  ids3 <- sample(c(core, sample(rem, 20)))
  list(pima_server1 = pdf[pdf$patient_id %in% ids1,
                           c("patient_id", "age", "bmi")],
       pima_server2 = pdf[pdf$patient_id %in% ids2,
                           c("patient_id", "ped", "glu")],
       pima_server3 = pdf[pdf$patient_id %in% ids3,
                           c("patient_id", "bp", "skin", "npreg",
                             "diabetes")])
}

#' Pima K=3 + bp-tertile augmentation on s3. Used by multinom_warm K=3
#' and ordinal_warm K=3.
build_pima_warm_k3 <- function(seed = 123L, bp_seed = 124L) {
  k3 <- build_pima_k3(seed)
  set.seed(bp_seed)
  brks <- quantile(k3$pima_server3$bp, c(0, 1/3, 2/3, 1), na.rm = TRUE)
  s3 <- k3$pima_server3
  s3$bp_cls <- factor(
    cut(s3$bp, breaks = brks, include.lowest = TRUE,
        labels = c("low", "med", "high")),
    levels = c("low", "med", "high"), ordered = TRUE)
  s3$low_ind  <- as.integer(s3$bp_cls == "low")
  s3$med_ind  <- as.integer(s3$bp_cls == "med")
  s3$high_ind <- as.integer(s3$bp_cls == "high")
  s3$low_leq  <- as.integer(as.integer(s3$bp_cls) <= 1L)
  s3$med_leq  <- as.integer(as.integer(s3$bp_cls) <= 2L)
  s3$high_leq <- 1L
  list(pima_warm_s1 = k3$pima_server1,
       pima_warm_s2 = k3$pima_server2,
       pima_warm_s3 = s3)
}

## When sourced as a script, print a manifest.
if (sys.nframe() == 0L) {
  manifest <- list(
    birthwt_k2 = build_birthwt_k2(),
    pima_k2    = build_pima_k2(),
    lung_k2    = build_lung_k2(),
    glmm_k2    = build_long_glmm_k2(),
    lmm_k2     = build_long_lmm_k2(),
    ipw_k2     = build_cgd_ipw_k2(),
    pima_k3    = build_pima_k3(),
    pima_warm_k3 = build_pima_warm_k3())
  for (nm in names(manifest)) {
    obj <- manifest[[nm]]
    cat(sprintf("== %s ==\n", nm))
    for (tnm in setdiff(names(obj), "diagnostics"))
      if (is.data.frame(obj[[tnm]]))
        cat(sprintf("  %-18s rows=%d cols=%d\n", tnm,
                    nrow(obj[[tnm]]), ncol(obj[[tnm]])))
    if (!is.null(obj$diagnostics)) {
      d <- obj$diagnostics
      diag_str <- paste(sprintf("%s=%g", names(d),
                                  unlist(d, use.names = FALSE)),
                          collapse = " ")
      cat(sprintf("  diagnostics: %s\n", diag_str))
    }
  }
}
