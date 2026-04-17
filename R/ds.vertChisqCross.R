#' @title Cross-server chi-square via Beaver cross-products (skeleton)
#' @description Placeholder for a chi-square test where the two
#'   categorical variables are on DIFFERENT servers. The plan:
#'
#'     1. Each server one-hot encodes its variable into K (resp. L)
#'        binary columns (new helper: dsvertOneHotDS).
#'     2. The binary indicator columns are secret-shared between the
#'        two DCF parties via the existing k2ShareInputDS machinery.
#'     3. For each (k, l) pair: compute sum_i X_{ki} * Y_{li} via a
#'        Beaver dot product on the shared indicator vectors. This is
#'        exactly the n_{kl} cell count.
#'     4. Client assembles the K x L contingency table and computes
#'        Pearson chi-square / Fisher exact as usual.
#'
#'   Requires:
#'     - dsvertOneHotDS server helper (new).
#'     - K*L Beaver dot-product invocations (reuse localCorDS-style
#'       infrastructure).
#'
#'   Same-server chi-square is already supported via ds.vertChisq /
#'   ds.vertFisher. Scheduled for Month 3; signature reserved.
#' @export
ds.vertChisqCross <- function(data, var1, var2, var1_server = NULL,
                              var2_server = NULL,
                              correct = TRUE, datasources = NULL) {
  stop("ds.vertChisqCross: cross-server chi-square is scheduled for ",
       "Month 3. For same-server chi-square use ds.vertChisq. See ",
       "V2_PROGRESS.md for the one-hot sharing + Beaver cross-product ",
       "implementation plan.",
       call. = FALSE)
}
