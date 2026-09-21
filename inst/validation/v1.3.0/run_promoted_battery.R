# Evidence-only adaptation of tools/validate_promoted_historical_e2e.R.
# Preserves rc1 source; applies the v1.2.1 manifest's K split and continues
# after failed expectations. Each completed gate is checkpointed to CSV.
root <- normalizePath(Sys.getenv("DSVERT_VALIDATION_ROOT", "/workspace/rc1"))
client <- file.path(root, "dsVertClient")
out <- file.path(client, "inst/validation/v1.3.0")
Sys.setenv(DSVERT_SERVER_SOURCE = file.path(root, "dsVert"))
plan <- data.frame(
  name = c("count_k235", "synopsis_full_k2", "lasso_focal_k235", "describe_k35", "gaussian_k35", "glm_grid_k35", "survival_k35", "cross_owner_tamper_k3"),
  k = c("2,3,5", "2", "2,3,5", "3,5", "3,5", "3,5", "3,5", "3"),
  family = c("", "", "gaussian_lasso_focal", "describe", "gaussian", "glm_grid", "survival", "cross_owner_tamper"))
source_text <- readLines(file.path(client, "tools/validate_promoted_historical_e2e.R"))
# Reuse the committed driver's actual devtools call, changing only the
# required failure policy and recording its returned structured results.
expr <- parse(text = source_text)
last <- expr[[length(expr)]]
call <- last[[4]][[5]]
stopifnot(identical(call[[1]], quote(devtools::test)))
call$stop_on_failure <- FALSE
results <- list()
for (i in seq_len(nrow(plan))) {
  gate <- list(name = plan$name[i], filter = if (i == 1L) "dp-count-synopsis-rock-e2e" else "dp-synopsis-describe-rock-e2e")
  client_root <- client
  Sys.setenv(DSVERT_TEST_SYNOPSIS_E2E_K = plan$k[i], DSVERT_TEST_SYNOPSIS_E2E_FAMILY = plan$family[i])
  cat("\n==", gate$name, "K=", plan$k[i], "==\n")
  result <- tryCatch(eval(call), error = identity)
  if (inherits(result, "error")) {
    writeLines(conditionMessage(result), file.path(out, paste0(gate$name, "_error.txt")))
    next
  }
  d <- as.data.frame(result)
  d <- d[, intersect(c("file", "test", "nb", "passed", "failed", "skipped", "error", "real"), names(d)), drop = FALSE]
  d <- cbind(run = gate$name, k = plan$k[i], d)
  results[[gate$name]] <- d
  write.csv(do.call(rbind, results), file.path(out, "harnessB_v130_results.csv"), row.names = FALSE)
}
