test_that("cross-grid bounds agree with the existing client finite grids", {
  for (family in c("binomial", "poisson")) {
    for (g in c(8L, 16L, 18L)) {
      for (adjacency in c("add_remove_patient", "replace_one_fixed_cohort")) {
        beta <- list(c(0, 0), c(0, -1), c(8, 8), c(-8, -8))
        maximum <- if (family == "binomial") 1L else 1024L
        bounds <- .dsvert_dp_glm_grid_loss_bounds(family, beta, maximum)
        cross <- .dsvert_dp_glm_grid_cross_sensitivity(
          beta, family, maximum, g, 100, adjacency)
        caps <- ceiling(2^g * bounds)
        multiplier <- if (adjacency == "add_remove_patient") 1 else 2
        expect_equal(vapply(cross$candidate_bounds, `[[`, numeric(1L), "loss_bound"),
                     bounds)
        expect_equal(vapply(cross$candidate_bounds, `[[`, numeric(1L), "per_patient_cap"),
                     caps)
        expect_equal(unlist(cross$maximum_coordinates), 100 * caps)
        expect_equal(cross$raw_l1_sensitivity, multiplier * sum(caps))
        expect_equal(cross$raw_l2_sensitivity, multiplier * sqrt(sum(caps^2)))
        expect_equal(cross$natural_l1_sensitivity, multiplier * sum(caps) / 2^g)
        expect_equal(cross$natural_l2_sensitivity,
                     multiplier * sqrt(sum(caps^2)) / 2^g)
      }
    }
  }
})
