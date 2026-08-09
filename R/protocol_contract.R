.dsvert_select_eta_privacy <- function(eta_privacy, n_partitions) {
  if (!is.character(eta_privacy) || length(eta_privacy) != 1L ||
      is.na(eta_privacy) ||
      !eta_privacy %in% c("auto", "k2_beaver", "secure_agg")) {
    stop("eta_privacy must be 'auto', 'k2_beaver', or 'secure_agg'",
         call. = FALSE)
  }
  if (!is.numeric(n_partitions) || length(n_partitions) != 1L ||
      !is.finite(n_partitions) || n_partitions != floor(n_partitions) ||
      n_partitions < 2L) {
    stop("At least two vertical partitions are required", call. = FALSE)
  }
  n_partitions <- as.integer(n_partitions)
  selected <- if (identical(eta_privacy, "auto")) {
    if (n_partitions == 2L) "k2_beaver" else "secure_agg"
  } else {
    eta_privacy
  }
  if (identical(selected, "k2_beaver") && n_partitions != 2L) {
    stop("k2_beaver requires exactly 2 servers", call. = FALSE)
  }
  if (identical(selected, "secure_agg") && n_partitions < 3L) {
    stop("secure_agg requires at least 3 servers", call. = FALSE)
  }
  selected
}

.dsvert_label_joint_warm_start <- function(out, estimator) {
  if (!is.list(out) || !is.character(estimator) || length(estimator) != 1L ||
      is.na(estimator)) {
    stop("invalid joint-estimator output contract", call. = FALSE)
  }

  if (identical(estimator, "multinomial")) {
    warm_fields <- c("fits", "coefficients", "coefficients_ovr",
                     "std_errors", "family")
  } else if (identical(estimator, "ordinal")) {
    warm_fields <- c("fits", "thresholds", "thresholds_ovr", "beta",
                     "beta_po", "covariance_po", "po_test",
                     "std_errors_po", "joint_mle", "family")
  } else {
    stop("unknown joint estimator", call. = FALSE)
  }

  present <- intersect(warm_fields, names(out))
  out$warm_start <- out[present]
  out[present] <- NULL
  out$inference_status <- "unavailable_for_joint_estimator"
  out$inference_reason <- paste(
    "Warm-start uncertainty is not valid for the final joint",
    estimator, "estimator and is retained only under warm_start."
  )
  out
}
