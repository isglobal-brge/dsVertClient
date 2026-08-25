test_that("local categorical manifest comparison canonicalizes JSON scalars", {
  references <- c("site_b$disease", "site_b$exposure")
  selector <- list(
    dataset = "cohort", owner_peer = "site_b",
    columns = lapply(references, function(reference) list(
      reference = reference, column = sub("^site_b\\$", "", reference))),
    parent = list(logical_snapshot = list(
      logical_snapshot_id = "cohort", version = "schema-v1-parent",
      alignment_protocol_version = 1L)))
  empty_artifacts <- list(artifacts = list())
  manifest <- list(
    logical_snapshot = list(
      logical_snapshot_id = "cohort", version = "schema-v1-parent",
      alignment_protocol_version = 1),
    workload = list(
      primitive_scope = list(
        mode = "catalog_v1",
        selection = list(explicit_catalog = list(
          categorical_pairs = list(references)))),
      families = list(
        numeric_moments = empty_artifacts,
        numeric_pair_moments = empty_artifacts,
        gaussian_models = empty_artifacts,
        fixed_numeric_histograms = empty_artifacts,
        correlation_artifacts = list(),
        describe_artifacts = list(),
        survival_artifacts = list(),
        categorical_pairs = list(
          sets = list(list(
            dataset = "cohort", owner_peer = "site_b",
            included_pairs = list(references))),
          cross_artifacts = list()))))

  expect_no_error(.dsvert_dp_synopsis_local_pair_manifest_v1(
    manifest, selector))
  manifest$logical_snapshot$version <- "other-snapshot"
  expect_error(.dsvert_dp_synopsis_local_pair_manifest_v1(
    manifest, selector), "local categorical Synopsis projection")
})

test_that("local categorical projection commits strict-missing scope", {
  schema <- list(unsigned = list(
    datasets = list(cohort = list(
      dataset_id = "cohort-v1", dataset_version = "v1",
      schema_version = "test-v1", alignment_group = "cohort",
      patient_keys = list(site_a = "patient_id"),
      columns = list(
        disease = list(kind = "categorical", owner_peer = "site_a",
                       levels = c("no", "yes")),
        exposure = list(kind = "categorical", owner_peer = "site_a",
                        levels = c("unexposed", "exposed")))))),
    logical_snapshot = list(
      logical_snapshot_id = "cohort", version = "schema-v1-parent",
      alignment_protocol_version = 1L))
  context <- list(policy = list(
    domain = "test-domain", cohort_id = "cohort",
    peer_pinset_sha256 = strrep("a", 64L), primitive_scope = list()))

  projection <- .dsvert_dp_synopsis_local_pair_project_components_v1(
    schema, "cohort", "site_a", c("disease", "exposure"), context)
  expect_identical(projection$primitive_scope$strict_missing_categorical,
                   character())

  context$policy$primitive_scope$strict_missing_categorical <-
    list("disease", "exposure")
  strict_projection <- .dsvert_dp_synopsis_local_pair_project_components_v1(
    schema, "cohort", "site_a", c("disease", "exposure"), context)
  expect_identical(
    strict_projection$primitive_scope$strict_missing_categorical,
    c("disease", "exposure"))
})
