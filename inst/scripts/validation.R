## ----setup, message=FALSE-----------------------------------------------------
library(DSI)
library(DSOpal)
library(dsVertClient)
library(MASS)
library(survival)


## ----bw-ground-truth----------------------------------------------------------
data(birthwt)
df_bw <- birthwt
df_bw$patient_id <- sprintf("BW%03d", seq_len(nrow(df_bw)))

set.seed(123)
core <- sample(df_bw$patient_id, 160)
remaining <- setdiff(df_bw$patient_id, core)
ids1 <- sample(c(core, sample(remaining, 15)))
ids2 <- sample(c(core, sample(remaining, 10)))

bw_a <- df_bw[df_bw$patient_id %in% ids1, c("patient_id", "age", "lwt", "race", "smoke")]
bw_b <- df_bw[df_bw$patient_id %in% ids2, c("patient_id", "ptl", "ht", "ui", "ftv", "bwt", "low")]
bw_merged <- merge(bw_a, bw_b, by = "patient_id")

nrow(bw_merged)


## ----bw-centralised-----------------------------------------------------------
bw_formula <- ~ age + lwt + race + smoke + ptl + ht + ui

gt_bw_gaussian <- glm(update(bw_formula, bwt ~ .), data = bw_merged, family = gaussian())
gt_bw_binomial <- glm(update(bw_formula, low ~ .), data = bw_merged, family = binomial())
gt_bw_poisson  <- glm(update(bw_formula, ftv ~ .), data = bw_merged, family = poisson())


## ----cgd-ground-truth---------------------------------------------------------
data(cgd)
cgd_one <- cgd[cgd$enum == 1, ]
cgd_one$n_infections <- as.integer(tapply(cgd$status, cgd$id, sum)[as.character(cgd_one$id)])
cgd_one$patient_id <- sprintf("CGD%03d", seq_len(nrow(cgd_one)))
cgd_one$treat_num <- as.integer(cgd_one$treat == "rIFN-g")
cgd_one$sex_num <- as.integer(cgd_one$sex == "male")

nrow(cgd_one)


## ----cgd-centralised----------------------------------------------------------
gt_cgd_gaussian <- glm(weight ~ age + sex_num + height + treat_num + steroids, data = cgd_one, family = gaussian())
gt_cgd_binomial <- glm(status ~ age + sex_num + height + weight + treat_num + steroids, data = cgd_one, family = binomial())
gt_cgd_poisson  <- glm(n_infections ~ age + sex_num + height + weight + treat_num + steroids, data = cgd_one, family = poisson())


## ----bw-connect---------------------------------------------------------------
opal_url <- "https://opal-demo.obiba.org"

builder <- DSI::newDSLoginBuilder()
builder$append(server = "hospitalA", url = opal_url,
               user = "administrator", password = "password",
               table = "vert_demo.bw_hospital_a")
builder$append(server = "hospitalB", url = opal_url,
               user = "administrator", password = "password",
               table = "vert_demo.bw_hospital_b")
bw_conns <- DSI::datashield.login(logins = builder$build())
DSI::datashield.assign.table(bw_conns, "D",
  list(hospitalA = "vert_demo.bw_hospital_a",
       hospitalB = "vert_demo.bw_hospital_b"))


## ----bw-psi-------------------------------------------------------------------
bw_psi <- ds.psiAlign("D", "patient_id", "DA", verbose = FALSE, datasources = bw_conns)
bw_psi$n_common


## ----bw-glm-------------------------------------------------------------------
dv_bw_gaussian <- ds.vertGLM(bwt ~ age + lwt + race + smoke + ptl + ht + ui,
                              data = "DA", family = "gaussian",
                              verbose = FALSE, datasources = bw_conns)

dv_bw_binomial <- ds.vertGLM(low ~ age + lwt + race + smoke + ptl + ht + ui,
                              data = "DA", family = "binomial",
                              verbose = FALSE, datasources = bw_conns)

dv_bw_poisson  <- ds.vertGLM(ftv ~ age + lwt + race + smoke + ptl + ht + ui,
                              data = "DA", family = "poisson",
                              verbose = FALSE, datasources = bw_conns)


## ----bw-cor-pca---------------------------------------------------------------
bw_vars <- c("age", "lwt", "race", "smoke", "ptl", "ht", "ui")
dv_bw_cor <- ds.vertCor("DA", bw_vars, verbose = FALSE, datasources = bw_conns)
dv_bw_pca <- ds.vertPCA(cor_result = dv_bw_cor)


## ----bw-disconnect------------------------------------------------------------
DSI::datashield.logout(bw_conns)


## ----cgd-connect--------------------------------------------------------------
builder <- DSI::newDSLoginBuilder()
builder$append(server = "demographics", url = opal_url,
               user = "administrator", password = "password",
               table = "vert_demo.cgd_demographics")
builder$append(server = "treatment", url = opal_url,
               user = "administrator", password = "password",
               table = "vert_demo.cgd_treatment")
builder$append(server = "outcomes", url = opal_url,
               user = "administrator", password = "password",
               table = "vert_demo.cgd_outcomes")
cgd_conns <- DSI::datashield.login(logins = builder$build())
DSI::datashield.assign.table(cgd_conns, "D",
  list(demographics = "vert_demo.cgd_demographics",
       treatment    = "vert_demo.cgd_treatment",
       outcomes     = "vert_demo.cgd_outcomes"))


## ----cgd-psi------------------------------------------------------------------
cgd_psi <- ds.psiAlign("D", "patient_id", "DA", verbose = FALSE, datasources = cgd_conns)
cgd_psi$n_common


## ----cgd-glm------------------------------------------------------------------
dv_cgd_gaussian <- ds.vertGLM(weight ~ age + sex_num + height + treat_num + steroids,
                               data = "DA", family = "gaussian",
                               verbose = FALSE, datasources = cgd_conns)

dv_cgd_binomial <- ds.vertGLM(status ~ age + sex_num + height + weight + treat_num + steroids,
                               data = "DA", family = "binomial",
                               verbose = FALSE, datasources = cgd_conns)

dv_cgd_poisson  <- ds.vertGLM(n_infections ~ age + sex_num + height + weight + treat_num + steroids,
                               data = "DA", family = "poisson",
                               verbose = FALSE, datasources = cgd_conns)


## ----cgd-cor-pca--------------------------------------------------------------
cgd_vars <- c("age", "sex_num", "height", "weight", "treat_num", "steroids")
dv_cgd_cor <- ds.vertCor("DA", cgd_vars, verbose = FALSE, datasources = cgd_conns)
dv_cgd_pca <- ds.vertPCA(cor_result = dv_cgd_cor)


## ----cgd-disconnect-----------------------------------------------------------
DSI::datashield.logout(cgd_conns)


## ----coef-comparison----------------------------------------------------------
compare_glm <- function(gt, dv, label) {
  nms <- names(dv$coefficients)
  data.frame(
    variable    = nms,
    centralised = round(coef(gt)[nms], 4),
    dsVert      = round(dv$coefficients[nms], 4),
    delta       = round(abs(coef(gt)[nms] - dv$coefficients[nms]), 4),
    SE_central  = round(summary(gt)$coefficients[nms, "Std. Error"], 4),
    SE_dsVert   = round(dv$std_errors[nms], 4),
    row.names   = NULL
  )
}


## ----bw-g-compare-------------------------------------------------------------
compare_glm(gt_bw_gaussian, dv_bw_gaussian, "birthwt Gaussian")


## ----bw-b-compare-------------------------------------------------------------
compare_glm(gt_bw_binomial, dv_bw_binomial, "birthwt Binomial")


## ----bw-p-compare-------------------------------------------------------------
compare_glm(gt_bw_poisson, dv_bw_poisson, "birthwt Poisson")


## ----cgd-g-compare------------------------------------------------------------
compare_glm(gt_cgd_gaussian, dv_cgd_gaussian, "CGD Gaussian")


## ----cgd-b-compare------------------------------------------------------------
compare_glm(gt_cgd_binomial, dv_cgd_binomial, "CGD Binomial")


## ----cgd-p-compare------------------------------------------------------------
compare_glm(gt_cgd_poisson, dv_cgd_poisson, "CGD Poisson")


## ----deviance-comparison------------------------------------------------------
deviance_table <- data.frame(
  dataset = rep(c("birthwt", "CGD"), each = 3),
  K       = rep(c(2, 3), each = 3),
  family  = rep(c("Gaussian", "Binomial", "Poisson"), 2),
  centralised = c(
    deviance(gt_bw_gaussian), deviance(gt_bw_binomial), deviance(gt_bw_poisson),
    deviance(gt_cgd_gaussian), deviance(gt_cgd_binomial), deviance(gt_cgd_poisson)),
  dsVert = c(
    dv_bw_gaussian$deviance, dv_bw_binomial$deviance, dv_bw_poisson$deviance,
    dv_cgd_gaussian$deviance, dv_cgd_binomial$deviance, dv_cgd_poisson$deviance)
)
deviance_table$error_pct <- round(
  100 * abs(deviance_table$dsVert - deviance_table$centralised) / deviance_table$centralised, 2)
deviance_table$centralised <- round(deviance_table$centralised, 2)
deviance_table$dsVert <- round(deviance_table$dsVert, 2)
deviance_table


## ----bw-cor-compare-----------------------------------------------------------
gt_bw_cor <- cor(bw_merged[, bw_vars])
max(abs(dv_bw_cor$correlation - gt_bw_cor))


## ----cgd-cor-compare----------------------------------------------------------
gt_cgd_cor <- cor(cgd_one[, cgd_vars])
max(abs(dv_cgd_cor$correlation - gt_cgd_cor))


## ----bw-pca-compare-----------------------------------------------------------
gt_bw_pca <- eigen(gt_bw_cor)
data.frame(
  component   = paste0("PC", 1:length(gt_bw_pca$values)),
  centralised = round(gt_bw_pca$values, 4),
  dsVert      = round(dv_bw_pca$eigenvalues, 4),
  delta       = round(abs(gt_bw_pca$values - dv_bw_pca$eigenvalues), 6)
)


## ----cgd-pca-compare----------------------------------------------------------
gt_cgd_pca <- eigen(gt_cgd_cor)
data.frame(
  component   = paste0("PC", 1:length(gt_cgd_pca$values)),
  centralised = round(gt_cgd_pca$values, 4),
  dsVert      = round(dv_cgd_pca$eigenvalues, 4),
  delta       = round(abs(gt_cgd_pca$values - dv_cgd_pca$eigenvalues), 6)
)

