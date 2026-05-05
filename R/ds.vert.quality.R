.dsvert_quality <- function(status = c("ok", "approximate", "degraded",
                                       "failed"),
                            warnings = character(0),
                            metrics = list()) {
  status <- match.arg(status)
  list(status = status,
       warnings = as.character(warnings),
       metrics = metrics)
}

.dsvert_quality_from_convergence <- function(converged, metric = NA_real_,
                                             tolerance = NA_real_,
                                             label = "optimizer") {
  warnings <- character(0)
  status <- "ok"
  if (!isTRUE(converged)) {
    status <- "degraded"
    warnings <- c(warnings, paste0(label, " did not meet its convergence criterion."))
  }
  if (is.finite(metric) && is.finite(tolerance) && metric > 10 * tolerance) {
    status <- "degraded"
    warnings <- c(warnings, sprintf(
      "%s final diagnostic %.4g is more than 10x tolerance %.4g.",
      label, metric, tolerance))
  }
  .dsvert_quality(status = status, warnings = unique(warnings),
                  metrics = list(metric = metric, tolerance = tolerance))
}
