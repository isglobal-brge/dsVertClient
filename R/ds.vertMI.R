#' @title Federated multiple imputation with server-local imputations (skeleton)
#' @description Multiple imputation with imputations that stay on the
#'   server holding the variable with missingness. Outer loop:
#'
#'     for m in 1..M:
#'       1. server imputes its own missing values using its own local
#'          sampler (Bayesian ridge / PMM style, parameters derived from
#'          the previous round's aggregate moments)
#'       2. client runs ds.vertGLM on the imputed dataset, collecting
#'          beta_m and Cov(beta_m)
#'     Rubin-pool (beta_m, Cov(beta_m)) client-side -> pooled_beta,
#'     total variance with within/between decomposition.
#'
#'   Privacy: the M imputations stay server-local; only the M pooled
#'   coefficient vectors cross back to the client. This is a stronger
#'   privacy story than standard federated MI (which often relays
#'   imputed values).
#'
#'   Requires:
#'     - Server-side imputation function (new: dsvertImputeColumnDS)
#'       that takes the aggregate moments as a prior and draws local
#'       imputations from a Bayesian ridge posterior.
#'     - Rubin-pool helpers client-side (already trivial).
#'
#'   Scheduled for Month 3; signature reserved.
#' @export
ds.vertMI <- function(formula, data = NULL,
                      impute_columns = NULL,
                      m = 20L, max_iter_per_impute = 30L,
                      verbose = TRUE, datasources = NULL, ...) {
  stop("ds.vertMI: server-local multiple imputation is scheduled for ",
       "Month 3. The design keeps imputations on the server holding the ",
       "missing variable; client only sees pooled coefficient vectors ",
       "(Rubin's rules). See V2_PROGRESS.md for the server-side ",
       "dsvertImputeColumnDS helper plan.",
       call. = FALSE)
}
