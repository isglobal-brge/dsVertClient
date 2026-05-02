test_that("ds.vertOrdinalJointNewton is diagnostic-only by default", {
  old_opt <- getOption("dsvert.allow_patient_level_ordinal_joint", NULL)
  old_env <- Sys.getenv("DSVERT_ALLOW_PATIENT_LEVEL_ORDINAL_JOINT", unset = NA)
  on.exit({
    if (is.null(old_opt)) {
      options(dsvert.allow_patient_level_ordinal_joint = NULL)
    } else {
      options(dsvert.allow_patient_level_ordinal_joint = old_opt)
    }
    if (is.na(old_env)) {
      Sys.unsetenv("DSVERT_ALLOW_PATIENT_LEVEL_ORDINAL_JOINT")
    } else {
      Sys.setenv(DSVERT_ALLOW_PATIENT_LEVEL_ORDINAL_JOINT = old_env)
    }
  }, add = TRUE)
  options(dsvert.allow_patient_level_ordinal_joint = FALSE)
  Sys.unsetenv("DSVERT_ALLOW_PATIENT_LEVEL_ORDINAL_JOINT")

  expect_error(
    ds.vertOrdinalJointNewton(
      y ~ x,
      levels_ordered = c("low", "mid", "high"),
      datasources = list()),
    "disabled under strict non-disclosure")
})
