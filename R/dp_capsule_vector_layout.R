# Public-coordinate map for the signed biomedical capsule manifest.  This must
# remain byte-for-byte compatible with dsVert's server-side layout: the client
# uses it only after the exact final DP vector has been jointly released.

.DSVERT_CLIENT_DP_CAPSULE_LAYOUT_VERSION <-
  "dsvert-biomedical-capsule-coordinate-layout-v4"
.DSVERT_CLIENT_DP_CAPSULE_RELEASE_CACHE_VERSION <-
  "dsvert-biomedical-capsule-public-release-cache-v1"
.DSVERT_CLIENT_DP_CAPSULE_RELEASE_CACHE_MAX_ENTRIES <- 4L
.DSVERT_CLIENT_DP_CAPSULE_RELEASE_CACHE_MAX_BYTES <- 64L * 1024L^2
.dsvert_dp_capsule_release_cache <- new.env(parent = emptyenv())
.dsvert_dp_capsule_release_cache$entries <- list()

.dsvert_dp_capsule_release_cache_clear <- function() {
  .dsvert_dp_capsule_release_cache$entries <- list()
  invisible(TRUE)
}

.dsvert_dp_capsule_release_cache_key <- function(
    datasources, status, manifest_bundle) {
  if (!is.list(datasources) || !length(datasources) ||
      is.null(names(datasources)) || anyNA(names(datasources)) ||
      any(!nzchar(names(datasources))) || anyDuplicated(names(datasources)) ||
      !is.list(status) || is.null(names(status)) || anyNA(names(status)) ||
      any(!nzchar(names(status))) || anyDuplicated(names(status)) ||
      !setequal(names(status), names(datasources)) ||
      !is.list(manifest_bundle) ||
      !.dsvert_vector_hex(manifest_bundle$manifest_sha256) ||
      !.dsvert_vector_hex(manifest_bundle$capsule_id)) {
    return(NULL)
  }
  servers <- sort(names(datasources), method = "radix")
  stable_fields <- c(
    "version", "enabled", "privacy_contract", "policy", "noise_root",
    "release_domain", "role")
  status_fields <- c(
    stable_fields, "composition_telemetry", "release_instance_telemetry")
  status <- unclass(status)[servers]
  if (!all(vapply(status, function(entry) {
        .dsvert_dp_has_exact_names(unclass(entry), status_fields)
      }, logical(1L)))) return(NULL)
  stable_status <- lapply(status, function(entry) {
    unclass(entry)[stable_fields]
  })
  tryCatch(.dsvert_vector_hash(list(
    version = .DSVERT_CLIENT_DP_CAPSULE_RELEASE_CACHE_VERSION,
    servers = as.list(servers),
    stable_control_plane_status = stable_status,
    manifest_sha256 = manifest_bundle$manifest_sha256,
    capsule_id = manifest_bundle$capsule_id)),
    error = function(error) NULL)
}

.dsvert_dp_capsule_release_cache_get <- function(key) {
  if (!is.character(key) || length(key) != 1L || is.na(key) ||
      !grepl("^[0-9a-f]{64}$", key)) return(NULL)
  entries <- .dsvert_dp_capsule_release_cache$entries
  index <- which(vapply(entries, function(entry) {
    identical(entry$key, key)
  }, logical(1L)))
  if (length(index) != 1L) return(NULL)
  entry <- entries[[index]]
  .dsvert_dp_capsule_release_cache$entries <- c(
    list(entry), entries[-index])
  entry$value
}

.dsvert_dp_capsule_release_cache_put <- function(
    key, value,
    .max_entries = .DSVERT_CLIENT_DP_CAPSULE_RELEASE_CACHE_MAX_ENTRIES,
    .max_bytes = .DSVERT_CLIENT_DP_CAPSULE_RELEASE_CACHE_MAX_BYTES) {
  valid_limits <- is.numeric(.max_entries) && length(.max_entries) == 1L &&
    !is.na(.max_entries) && is.finite(.max_entries) && .max_entries >= 1 &&
    .max_entries == floor(.max_entries) &&
    is.numeric(.max_bytes) && length(.max_bytes) == 1L &&
    !is.na(.max_bytes) && is.finite(.max_bytes) && .max_bytes >= 1
  if (!isTRUE(valid_limits) || !is.character(key) || length(key) != 1L ||
      is.na(key) || !grepl("^[0-9a-f]{64}$", key) || !is.list(value) ||
      !is.list(value$manifest_bundle)) return(invisible(FALSE))
  cached <- value
  cached$manifest_bundle$context <- NULL
  bytes <- as.numeric(utils::object.size(cached))
  if (!is.finite(bytes) || bytes > .max_bytes) return(invisible(FALSE))
  entries <- .dsvert_dp_capsule_release_cache$entries
  entries <- entries[!vapply(entries, function(entry) {
    identical(entry$key, key)
  }, logical(1L))]
  entries <- c(list(list(key = key, bytes = bytes, value = cached)), entries)
  while (length(entries) > .max_entries ||
         sum(vapply(entries, `[[`, numeric(1L), "bytes")) > .max_bytes) {
    entries <- entries[-length(entries)]
  }
  .dsvert_dp_capsule_release_cache$entries <- entries
  invisible(TRUE)
}

.dsvert_dp_capsule_sorted_artifact_names <- function(value) {
  if (!length(value)) character() else sort(names(value), method = "radix")
}

.dsvert_dp_capsule_manifest_strings <- function(
    value, what, sorted = FALSE, unique = TRUE) {
  if (is.list(value) && is.null(names(value))) {
    valid <- all(vapply(value, function(item) {
      is.character(item) && length(item) == 1L && !is.na(item) && nzchar(item)
    }, logical(1L)))
    if (isTRUE(valid)) value <- unname(unlist(value, use.names = FALSE))
  }
  valid <- is.character(value) && length(value) > 0L && is.null(names(value)) &&
    !anyNA(value) && all(nzchar(value)) &&
    (!isTRUE(unique) || !anyDuplicated(value)) &&
    (!isTRUE(sorted) || identical(value, sort(value, method = "radix")))
  if (!isTRUE(valid)) stop("Invalid signed capsule ", what, call. = FALSE)
  enc2utf8(unname(value))
}

.dsvert_dp_capsule_manifest_numbers <- function(value, what) {
  if (is.list(value) && is.null(names(value))) {
    valid <- all(vapply(value, function(item) {
      is.numeric(item) && length(item) == 1L && !is.na(item) &&
        is.finite(item)
    }, logical(1L)))
    if (isTRUE(valid)) value <- unname(unlist(value, use.names = FALSE))
  }
  if (!is.numeric(value) || !length(value) || !is.null(names(value)) ||
      anyNA(value) || any(!is.finite(value))) {
    stop("Invalid signed capsule ", what, call. = FALSE)
  }
  unname(as.numeric(value))
}

.dsvert_dp_capsule_vector_layout <- function(manifest) {
  families <- tryCatch(manifest$workload$families,
                       error = function(error) NULL)
  required <- c(
    "admitted_count", "numeric_moments", "numeric_pair_moments",
    "gaussian_models",
    "fixed_numeric_histograms", "categorical_marginals",
    "categorical_pairs", "correlation_artifacts", "describe_artifacts",
    "survival_artifacts")
  if (!is.list(families) || !all(required %in% names(families))) {
    stop("The biomedical capsule coordinate contract is invalid",
         call. = FALSE)
  }

  cursor <- 1L
  blocks <- list()
  add <- function(family, key, coordinate_length, owner_peer, dataset,
                  descriptor) {
    valid_id <- function(value) {
      is.character(value) && length(value) == 1L && !is.na(value) &&
        nzchar(value)
    }
    if (!valid_id(family) || !valid_id(key) || !valid_id(owner_peer) ||
        !valid_id(dataset) || !is.list(descriptor) ||
        !is.numeric(coordinate_length) || length(coordinate_length) != 1L ||
        is.na(coordinate_length) || !is.finite(coordinate_length) ||
        coordinate_length < 1 || coordinate_length != floor(coordinate_length) ||
        coordinate_length > .DSVERT_DP_MAX_COORDINATES ||
        cursor > .DSVERT_DP_MAX_COORDINATES - coordinate_length + 1L) {
      stop("The biomedical capsule coordinate layout is invalid",
           call. = FALSE)
    }
    start <- cursor
    end <- cursor + as.integer(coordinate_length) - 1L
    block_id <- paste(family, key, sep = "::")
    if (block_id %in% names(blocks)) {
      stop("The biomedical capsule coordinate layout is ambiguous",
           call. = FALSE)
    }
    blocks[[block_id]] <<- list(
      family = family, key = key, start = start, end = end,
      length = as.integer(coordinate_length), owner_peer = owner_peer,
      dataset = dataset, descriptor = descriptor)
    cursor <<- end + 1L
  }

  count <- families$admitted_count
  add("admitted_count", "canonical", 1L, count$owner_peer,
      count$dataset, count)

  numeric <- families$numeric_moments$artifacts
  for (key in .dsvert_dp_capsule_sorted_artifact_names(numeric)) {
    artifact <- numeric[[key]]
    add("numeric_moments", key, 3L, artifact$owner_peer,
        artifact$dataset, artifact)
  }

  numeric_pairs <- families$numeric_pair_moments$artifacts
  for (key in .dsvert_dp_capsule_sorted_artifact_names(numeric_pairs)) {
    artifact <- numeric_pairs[[key]]
    add("numeric_pair_moments", key, 6L, artifact$owner_peer,
        artifact$dataset, artifact)
  }

  gaussian <- families$gaussian_models$artifacts
  for (key in .dsvert_dp_capsule_sorted_artifact_names(gaussian)) {
    artifact <- gaussian[[key]]
    add("gaussian_models", key, artifact$coordinate_count,
        artifact$owner_peer, artifact$dataset, artifact)
  }

  histograms <- families$fixed_numeric_histograms$artifacts
  for (key in .dsvert_dp_capsule_sorted_artifact_names(histograms)) {
    artifact <- histograms[[key]]
    add("fixed_numeric_histograms", key, artifact$coordinate_count,
        artifact$owner_peer, artifact$dataset, artifact)
  }

  marginals <- families$categorical_marginals$artifacts
  for (key in .dsvert_dp_capsule_sorted_artifact_names(marginals)) {
    artifact <- marginals[[key]]
    add("categorical_marginals", key, length(artifact$levels),
        artifact$owner_peer, artifact$dataset, artifact)
  }

  pair_sets <- families$categorical_pairs$sets
  for (set_key in .dsvert_dp_capsule_sorted_artifact_names(pair_sets)) {
    set <- pair_sets[[set_key]]
    columns <- set$columns
    if (!is.list(columns) || !length(columns)) {
      stop("The biomedical categorical-pair layout is invalid",
           call. = FALSE)
    }
    column_names <- tryCatch(
      vapply(columns, `[[`, character(1L), "column"),
      error = function(error) character())
    if (length(column_names) != length(columns) || anyNA(column_names) ||
        any(!nzchar(column_names)) || anyDuplicated(column_names)) {
      stop("The biomedical categorical-pair layout is invalid",
           call. = FALSE)
    }
    names(columns) <- column_names
    included <- set$included_pairs
    if (!is.list(included) || !length(included)) {
      stop("The biomedical categorical-pair layout is invalid",
           call. = FALSE)
    }
    included <- lapply(included, function(pair) {
      if (is.list(pair) && is.null(names(pair)) && length(pair) == 2L &&
          all(vapply(pair, function(value) {
            is.character(value) && length(value) == 1L && !is.na(value) &&
              nzchar(value)
          }, logical(1L)))) {
        pair <- unname(unlist(pair, use.names = FALSE))
      }
      if (!is.character(pair) || length(pair) != 2L || !is.null(names(pair)) ||
          anyNA(pair) || any(!nzchar(pair)) || anyDuplicated(pair) ||
          !identical(pair, sort(pair, method = "radix")) ||
          any(!pair %in% column_names)) {
        stop("The biomedical categorical-pair layout is invalid",
             call. = FALSE)
      }
      unname(pair)
    })
    pair_keys <- vapply(
      included, paste, character(1L), collapse = "::")
    if (anyDuplicated(pair_keys) ||
        !identical(pair_keys, sort(pair_keys, method = "radix")) ||
        !.dsvert_dp_is_integer(set$pair_count, length(included),
                               length(included))) {
      stop("The biomedical categorical-pair layout is invalid",
           call. = FALSE)
    }
    pair_coordinates <- 0
    for (pair in included) {
      left_column <- columns[[pair[[1L]]]]
      right_column <- columns[[pair[[2L]]]]
      coordinate_length <- as.double(length(left_column$levels)) *
        as.double(length(right_column$levels))
      if (!is.finite(coordinate_length) || coordinate_length < 1 ||
          coordinate_length > .DSVERT_DP_MAX_COORDINATES - pair_coordinates) {
        stop("The biomedical categorical-pair layout is invalid",
             call. = FALSE)
      }
      pair_coordinates <- pair_coordinates + coordinate_length
      key <- paste(left_column$column, right_column$column, sep = "::")
      descriptor <- list(
        left = left_column, right = right_column,
        repeated_record_policy = set$repeated_record_policy,
        missingness_policy = set$missingness_policy)
      add(
        "categorical_pairs", paste(set_key, key, sep = "::"),
        coordinate_length, set$owner_peer, set$dataset, descriptor)
    }
    if (!.dsvert_dp_is_integer(set$coordinate_count, pair_coordinates,
                               pair_coordinates)) {
      stop("The biomedical categorical-pair layout is invalid",
           call. = FALSE)
    }
  }

  cross_pairs <- families$categorical_pairs$cross_artifacts %||% list()
  for (analysis_id in
       .dsvert_dp_capsule_sorted_artifact_names(cross_pairs)) {
    artifact <- cross_pairs[[analysis_id]]
    add(
      "categorical_pairs", paste("cross", analysis_id, sep = "::"),
      artifact$coordinate_count, artifact$left$owner_peer,
      artifact$left$dataset, artifact)
  }

  survival <- families$survival_artifacts
  for (key in .dsvert_dp_capsule_sorted_artifact_names(survival)) {
    artifact <- survival[[key]]
    add("survival_artifacts", key, artifact$coordinate_count,
        artifact$owner_peer, artifact$dataset, artifact)
  }

  coordinate_count <- cursor - 1L
  expected <- manifest$workload$coordinate_count
  if (!is.numeric(expected) || length(expected) != 1L || is.na(expected) ||
      !is.finite(expected) || expected != floor(expected) ||
      !identical(as.numeric(coordinate_count), as.numeric(expected))) {
    stop("The biomedical capsule coordinate layout does not match its manifest",
         call. = FALSE)
  }
  list(
    version = .DSVERT_CLIENT_DP_CAPSULE_LAYOUT_VERSION,
    coordinate_count = as.integer(coordinate_count), blocks = blocks,
    sha256 = .dsvert_vector_hash(list(
      version = .DSVERT_CLIENT_DP_CAPSULE_LAYOUT_VERSION,
      coordinate_count = as.integer(coordinate_count),
      admitted_count = count,
      numeric_moments = families$numeric_moments,
      numeric_pair_moments = families$numeric_pair_moments,
      gaussian_models = families$gaussian_models,
      fixed_numeric_histograms = families$fixed_numeric_histograms,
      categorical_marginals = families$categorical_marginals,
      categorical_pairs = families$categorical_pairs,
      correlation_artifacts = families$correlation_artifacts,
      survival_artifacts = families$survival_artifacts)))
}

.dsvert_dp_capsule_vector_blocks <- function(
    layout, family, dataset = NULL, owner_peer = NULL) {
  blocks <- layout$blocks[vapply(
    layout$blocks, function(block) identical(block$family, family),
    logical(1L))]
  if (!is.null(dataset)) {
    blocks <- blocks[vapply(
      blocks, function(block) identical(block$dataset, dataset), logical(1L))]
  }
  if (!is.null(owner_peer)) {
    blocks <- blocks[vapply(
      blocks, function(block) identical(block$owner_peer, owner_peer),
      logical(1L))]
  }
  blocks
}

.dsvert_dp_capsule_vector_values <- function(release, block) {
  if (!inherits(release, "dsvert_joint_dp_vector") ||
      !is.numeric(release$values) || anyNA(release$values) ||
      any(!is.finite(release$values)) ||
      length(release$values) != release$coordinate_count ||
      !is.list(block) || !all(c("start", "end", "length") %in% names(block)) ||
      block$start < 1L || block$end > length(release$values) ||
      block$end - block$start + 1L != block$length) {
    stop("The final DP vector cannot be mapped to its signed coordinate layout",
         call. = FALSE)
  }
  unname(release$values[block$start:block$end])
}

.dsvert_dp_capsule_vector_run <- function(
    datasources, status = NULL, .aggregate = DSI::datashield.aggregate) {
  datasources <- .dsvert_dp_datasources(datasources)
  if (is.null(status)) {
    status <- .dsvert_joint_dp_capsule_status_impl(datasources, .aggregate)
  }
  manifest_bundle <- .dsvert_dp_capsule_manifest_build(
    datasources, status = status, .aggregate = .aggregate)
  cache_key <- .dsvert_dp_capsule_release_cache_key(
    datasources, status, manifest_bundle)
  cached <- .dsvert_dp_capsule_release_cache_get(cache_key)
  if (!is.null(cached)) {
    cached$status <- status
    cached$manifest_bundle <- manifest_bundle
    return(cached)
  }
  release <- .dsvert_joint_dp_vector_capsule(
    datasources, status = status, manifest_bundle = manifest_bundle,
    .aggregate = .aggregate)
  layout <- .dsvert_dp_capsule_vector_layout(release$manifest)
  if (!identical(as.numeric(layout$coordinate_count),
                 as.numeric(release$coordinate_count)) ||
      length(release$values) != layout$coordinate_count ||
      (!is.null(release$coordinate_order_sha256) &&
       !identical(release$coordinate_order_sha256, layout$sha256))) {
    stop("The released DP vector does not match its signed coordinate order",
         call. = FALSE)
  }
  result <- list(
    release = release, layout = layout, status = status,
    manifest_bundle = manifest_bundle)
  .dsvert_dp_capsule_release_cache_put(cache_key, result)
  result
}

.dsvert_dp_capsule_single_block <- function(
    layout, family, dataset = NULL, owner_peer = NULL, predicate = NULL,
    description = family) {
  blocks <- .dsvert_dp_capsule_vector_blocks(
    layout, family, dataset = dataset, owner_peer = owner_peer)
  if (!is.null(predicate)) {
    if (!is.function(predicate)) stop("Invalid capsule block predicate",
                                      call. = FALSE)
    blocks <- blocks[vapply(blocks, predicate, logical(1L))]
  }
  if (length(blocks) != 1L) {
    stop("The signed biomedical capsule does not contain exactly one ",
         description, call. = FALSE)
  }
  blocks[[1L]]
}

.dsvert_dp_vector_fraction <- function(value) {
  if (!.dsvert_vector_string(
      value, "^(0|[1-9][0-9]*)/[1-9][0-9]*\\z",
      .DSVERT_CLIENT_VECTOR_MAX_RECEIPT_BYTES)) {
    stop("Invalid vector implementation-delta certificate", call. = FALSE)
  }
  fields <- strsplit(value, "/", fixed = TRUE)[[1L]]
  numerator <- fields[[1L]]
  denominator <- fields[[2L]]
  numerator_bytes <- nchar(numerator, type = "bytes")
  denominator_bytes <- nchar(denominator, type = "bytes")
  same_width_less <- numerator_bytes == denominator_bytes &&
    !identical(numerator, denominator) &&
    identical(order(c(numerator, denominator), method = "radix"), 1:2)
  if (numerator_bytes > denominator_bytes ||
      (numerator_bytes == denominator_bytes && !same_width_less)) {
    stop("Invalid vector implementation-delta certificate", call. = FALSE)
  }
  min(1, .dsvert_dp_vector_fraction_upper(
    numerator, denominator,
    maximum_bytes = .DSVERT_CLIENT_VECTOR_MAX_RECEIPT_BYTES))
}

# Exact ideal tail for Z = X + Y where X and Y are independent symmetric
# two-sided-geometric variables with P(X=k) proportional to exp(-eps*|k|/S).
# The productive finite sampler is within the separately signed TV bound.
.dsvert_dp_vector_convolution_log_tail <- function(
    radius_steps, epsilon, sensitivity_steps) {
  if (!is.numeric(radius_steps) || length(radius_steps) != 1L ||
      is.na(radius_steps) || !is.finite(radius_steps) || radius_steps < 0 ||
      !is.numeric(epsilon) || length(epsilon) != 1L || is.na(epsilon) ||
      !is.finite(epsilon) || epsilon <= 0 ||
      !is.numeric(sensitivity_steps) || length(sensitivity_steps) != 1L ||
      is.na(sensitivity_steps) || !is.finite(sensitivity_steps) ||
      sensitivity_steps <= 0) {
    stop("Invalid vector accuracy inputs", call. = FALSE)
  }
  n <- floor(radius_steps) + 1
  log_probability <- -epsilon / sensitivity_steps
  one_minus_probability <- -expm1(log_probability)
  one_plus_probability <- 2 - one_minus_probability
  bracket <- 1 + (n - 2) * one_minus_probability +
    2 / one_plus_probability
  log(2) + n * log_probability - 2 * log(one_plus_probability) +
    log(bracket)
}

# Exact ideal tail for one symmetric two-sided-geometric variable:
# P(|Z| > r) = 2 p^(floor(r)+1) / (1+p), p=exp(-epsilon/S).
.dsvert_dp_vector_laplace_log_tail <- function(
    radius_steps, epsilon, sensitivity_steps) {
  if (!is.numeric(radius_steps) || length(radius_steps) != 1L ||
      is.na(radius_steps) || !is.finite(radius_steps) || radius_steps < 0 ||
      !is.numeric(epsilon) || length(epsilon) != 1L || is.na(epsilon) ||
      !is.finite(epsilon) || epsilon <= 0 ||
      !is.numeric(sensitivity_steps) || length(sensitivity_steps) != 1L ||
      is.na(sensitivity_steps) || !is.finite(sensitivity_steps) ||
      sensitivity_steps <= 0) {
    stop("Invalid vector accuracy inputs", call. = FALSE)
  }
  n <- floor(radius_steps) + 1
  log_p <- -epsilon / sensitivity_steps
  log(2) + n * log_p - log1p(exp(log_p))
}

.dsvert_dp_vector_adjacent_double <- function(value, upward) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      !is.logical(upward) || length(upward) != 1L || is.na(upward)) {
    stop("Invalid outward floating-point interval", call. = FALSE)
  }
  if (is.infinite(value)) {
    if (value > 0) return(if (upward) value else .Machine$double.xmax)
    return(if (upward) -.Machine$double.xmax else value)
  }
  smallest <- .Machine$double.xmin * .Machine$double.eps
  if (value == 0) return(if (upward) smallest else -smallest)
  bytes <- writeBin(value, raw(), size = 8L, endian = "little")
  increment <- xor(value < 0, upward)
  indices <- seq_len(8L)
  if (isTRUE(increment)) {
    for (index in indices) {
      byte <- as.integer(bytes[[index]])
      if (byte < 255L) {
        bytes[[index]] <- as.raw(byte + 1L)
        break
      }
      bytes[[index]] <- as.raw(0L)
    }
  } else {
    for (index in indices) {
      byte <- as.integer(bytes[[index]])
      if (byte > 0L) {
        bytes[[index]] <- as.raw(byte - 1L)
        break
      }
      bytes[[index]] <- as.raw(255L)
    }
  }
  readBin(bytes, what = "double", n = 1L, size = 8L,
          endian = "little")
}

.dsvert_dp_vector_next_up <- function(value) {
  .dsvert_dp_vector_adjacent_double(value, TRUE)
}

.dsvert_dp_vector_next_down <- function(value) {
  .dsvert_dp_vector_adjacent_double(value, FALSE)
}

.dsvert_dp_vector_fraction_upper <- function(
    numerator, denominator, maximum_bytes = 128L) {
  if (!.dsvert_vector_integer_text(numerator,
                                   maximum_bytes = maximum_bytes) ||
      !.dsvert_vector_integer_text(denominator, positive = TRUE,
                                   maximum_bytes = maximum_bytes)) {
    stop("Invalid signed rational accuracy bound", call. = FALSE)
  }
  if (identical(numerator, "0")) return(0)
  significand <- function(value, upward) {
    width <- nchar(value, type = "bytes")
    prefix <- substr(value, 1L, min(width, 32L))
    normalized <- if (nchar(prefix, type = "bytes") == 1L) prefix else
      paste0(substr(prefix, 1L, 1L), ".",
             substr(prefix, 2L, nchar(prefix, type = "bytes")))
    parsed <- suppressWarnings(as.numeric(normalized))
    if (isTRUE(upward)) .dsvert_dp_vector_next_up(parsed) else
      .dsvert_dp_vector_next_down(parsed)
  }
  top <- significand(numerator, TRUE)
  bottom <- significand(denominator, FALSE)
  exponent <- nchar(numerator, type = "bytes") -
    nchar(denominator, type = "bytes")
  scale <- .dsvert_dp_vector_next_up(
    suppressWarnings(as.numeric(paste0("1e", exponent))))
  value <- .dsvert_dp_vector_next_up(
    .dsvert_dp_vector_next_up(top / bottom) * scale)
  if (!is.finite(value) || value < 0) {
    stop("The signed rational accuracy bound is not representable",
         call. = FALSE)
  }
  for (index in seq_len(8L)) value <- .dsvert_dp_vector_next_up(value)
  value
}

.dsvert_dp_vector_dyadic_fraction_interval <- function(numerator, bits) {
  if (!.dsvert_vector_integer_text(numerator, positive = TRUE) ||
      !.dsvert_vector_whole(bits, 1, 1023)) {
    stop("The signed dyadic tail plan is invalid", call. = FALSE)
  }
  denominator_integer <- openssl::bignum("2") ^ as.integer(bits)
  numerator_integer <- tryCatch(
    openssl::bignum(numerator), error = function(error) NULL)
  if (is.null(numerator_integer) || numerator_integer >= denominator_integer) {
    stop("The signed dyadic tail plan is not representable",
         call. = FALSE)
  }
  outward <- function(value, upward) {
    for (index in seq_len(8L)) {
      value <- if (isTRUE(upward)) {
        .dsvert_dp_vector_next_up(value)
      } else {
        .dsvert_dp_vector_next_down(value)
      }
    }
    value
  }
  width <- nchar(numerator, type = "bytes")
  first_width <- width %% 9L
  if (first_width == 0L) first_width <- 9L
  starts <- c(1L, if (width > first_width) {
    seq.int(first_width + 1L, width, by = 9L)
  } else integer())
  widths <- c(first_width, rep.int(9L, length(starts) - 1L))
  first <- suppressWarnings(as.numeric(substr(
    numerator, starts[[1L]], starts[[1L]] + widths[[1L]] - 1L)))
  value <- lower <- upper <- first
  if (length(starts) > 1L) {
    for (index in seq.int(2L, length(starts))) {
      chunk <- suppressWarnings(as.numeric(substr(
        numerator, starts[[index]],
        starts[[index]] + widths[[index]] - 1L)))
      base <- 10^widths[[index]]
      value <- value * base + chunk
      lower <- outward(outward(lower * base, FALSE) + chunk, FALSE)
      upper <- outward(outward(upper * base, TRUE) + chunk, TRUE)
    }
  }
  denominator <- 2^as.numeric(bits)
  value <- value / denominator
  lower <- outward(lower / denominator, FALSE)
  upper <- outward(upper / denominator, TRUE)
  if (!is.finite(value) || !is.finite(lower) || !is.finite(upper) ||
      value <= 0 || value >= 1 || lower <= 0 || upper >= 1 ||
      lower > value || upper < value) {
    stop("The signed dyadic tail interval is not representable",
         call. = FALSE)
  }
  list(q = value, q_lower = lower, q_upper = upper)
}

.dsvert_dp_vector_dyadic_tail_context <- function(plan) {
  if (!is.list(plan) ||
      !.dsvert_vector_whole(plan$stop_bits, 1, 1023) ||
      !.dsvert_vector_integer_text(plan$stop_numerator)) {
    stop("The signed dyadic tail plan is invalid", call. = FALSE)
  }
  interval <- .dsvert_dp_vector_dyadic_fraction_interval(
    plan$stop_numerator, plan$stop_bits)
  q_lower <- interval$q_lower
  q_upper <- interval$q_upper
  log_p_upper <- .dsvert_dp_vector_next_up(
    .dsvert_dp_vector_next_up(log1p(-q_lower)))
  log_one_plus_p_lower <- .dsvert_dp_vector_next_down(
    .dsvert_dp_vector_next_down(log(2 - q_upper)))
  list(q = interval$q, q_lower = q_lower, q_upper = q_upper,
       log_p_upper = log_p_upper,
       log_one_plus_p_lower = log_one_plus_p_lower)
}

.dsvert_dp_vector_plan_log_tail_upper <- function(
    radius_steps, tail, convolution = FALSE) {
  if (!.dsvert_dp_is_integer(radius_steps, 0, 2^53 - 1) ||
      !is.list(tail)) {
    stop("Invalid certified dyadic tail input", call. = FALSE)
  }
  n <- as.numeric(radius_steps) + 1
  log_two_upper <- .dsvert_dp_vector_next_up(log(2))
  power_upper <- .dsvert_dp_vector_next_up(
    n * tail$log_p_upper)
  if (!isTRUE(convolution)) {
    return(.dsvert_dp_vector_next_up(
      .dsvert_dp_vector_next_up(
        log_two_upper + power_upper) - tail$log_one_plus_p_lower))
  }
  coefficient <- n - 2
  q_for_upper <- if (coefficient >= 0) tail$q_upper else tail$q_lower
  linear_upper <- .dsvert_dp_vector_next_up(coefficient * q_for_upper)
  denominator_lower <- .dsvert_dp_vector_next_down(2 - tail$q_upper)
  reciprocal_upper <- .dsvert_dp_vector_next_up(2 / denominator_lower)
  bracket_upper <- .dsvert_dp_vector_next_up(
    .dsvert_dp_vector_next_up(1 + linear_upper) + reciprocal_upper)
  log_bracket_upper <- .dsvert_dp_vector_next_up(log(bracket_upper))
  value <- .dsvert_dp_vector_next_up(log_two_upper + power_upper)
  twice_denominator_lower <- .dsvert_dp_vector_next_down(
    2 * tail$log_one_plus_p_lower)
  value <- .dsvert_dp_vector_next_up(
    value - twice_denominator_lower)
  .dsvert_dp_vector_next_up(value + log_bracket_upper)
}

.dsvert_dp_vector_sampler_tv_upper <- function(plan, exact_gc) {
  if (!is.list(plan) ||
      !.dsvert_vector_whole(plan$total_coordinate_count, 1,
                            .DSVERT_DP_MAX_COORDINATES)) {
    stop("The signed vector sampler-TV plan is invalid", call. = FALSE)
  }
  one <- .dsvert_dp_vector_fraction_upper(
    plan$one_geometric_tv_numerator,
    plan$one_geometric_tv_denominator)
  factor <- if (isTRUE(exact_gc)) 2 else 4
  total <- .dsvert_dp_vector_next_up(
    factor * as.numeric(plan$total_coordinate_count) * one)
  if (!is.finite(total) || total < 0 || total >= 1) {
    stop("The signed vector sampler-TV bound cannot certify accuracy",
         call. = FALSE)
  }
  total
}

.dsvert_dp_vector_release_profile <- function(release, manifest) {
  selection <- if (is.list(release$backend_selection)) {
    release$backend_selection
  } else NULL
  backend <- if (is.list(selection)) selection$backend else release$backend
  if (is.null(backend)) {
    backend <- if (identical(
        release$mechanism,
        .DSVERT_CLIENT_VECTOR_EXACT_RELEASE_MECHANISM)) {
      .DSVERT_CLIENT_VECTOR_EXACT_BACKEND
    } else if (identical(
        release$mechanism,
        .DSVERT_CLIENT_VECTOR_RELEASE_MECHANISM)) {
      .DSVERT_CLIENT_VECTOR_BACKEND
    } else if (identical(
        release$mechanism,
        .DSVERT_CLIENT_VECTOR_GAUSSIAN_RELEASE_MECHANISM)) {
      .DSVERT_CLIENT_VECTOR_GAUSSIAN_BACKEND
    } else NULL
  }
  profile <- .dsvert_vector_profile(
    manifest$workload$capsule_mechanism,
    manifest$workload$mechanism_selection, backend = backend)
  if (isTRUE(profile$selection_bound) && !is.null(selection)) {
    manifest_sha256 <- release$manifest_sha256
    if (!.dsvert_vector_hex(manifest_sha256)) {
      stop("The vector backend selection has no signed manifest binding",
           call. = FALSE)
    }
    .dsvert_joint_dp_vector_exact_gc_client_selection(
      selection, manifest_sha256,
      require_exact = isTRUE(profile$exact_gc))
  }
  profile
}

.dsvert_dp_vector_gaussian_accuracy_steps <- function(
    plan, coordinate_count, confidence) {
  required <- c(
    "sigma_squared_numerator", "sigma_squared_denominator",
    "maximum_noise_magnitude_two_peers")
  if (!is.list(plan) || !all(required %in% names(plan)) ||
      !.dsvert_vector_integer_text(
        plan$sigma_squared_numerator, positive = TRUE) ||
      !.dsvert_vector_integer_text(
        plan$sigma_squared_denominator, positive = TRUE) ||
      !.dsvert_vector_integer_text(
        plan$maximum_noise_magnitude_two_peers) ||
      !.dsvert_dp_is_integer(coordinate_count, 1,
                             .DSVERT_DP_MAX_COORDINATES) ||
      !.dsvert_dp_is_number(confidence, 0, 1, lower_open = TRUE) ||
      confidence >= 1) {
    stop("The signed Gaussian accuracy certificate is invalid",
         call. = FALSE)
  }

  support <- .dsvert_dp_vector_fraction_upper(
    plan$maximum_noise_magnitude_two_peers, "1")
  support_result <- function() {
    list(steps = support, finite_support = TRUE)
  }
  sigma_squared <- .dsvert_dp_vector_fraction_upper(
    plan$sigma_squared_numerator, plan$sigma_squared_denominator)
  vector_tv <- .dsvert_dp_vector_fraction_upper(
    plan$vector_total_tv_upper_numerator,
    plan$vector_total_tv_upper_denominator)

  # The signed plan proves, for the sum of its two complete independent
  # draws, P(max_j |Z_j| >= r) <=
  #   2*d*exp(-r^2/(4*sigma^2)) + 2*vector_TV.
  # Mirror the server's exact-rational radius proof: if
  # k=ceil(log2(target)), ln(target) < 7*k/10. Every floating operation below
  # is rounded outwards; an exact integer boundary therefore moves to the
  # next, conservative radius. The signed finite support is the unconditional
  # fallback when the requested alpha is smaller than the transfer term.
  alpha_lower <- 1 - as.numeric(confidence)
  for (index in seq_len(8L)) {
    alpha_lower <- .dsvert_dp_vector_next_down(alpha_lower)
  }
  transfer_upper <- 2 * vector_tv
  for (index in seq_len(8L)) {
    transfer_upper <- .dsvert_dp_vector_next_up(transfer_upper)
  }
  remaining_lower <- alpha_lower - transfer_upper
  for (index in seq_len(8L)) {
    remaining_lower <- .dsvert_dp_vector_next_down(remaining_lower)
  }
  if (!is.finite(remaining_lower) || remaining_lower <= 0) {
    return(support_result())
  }

  target_upper <- (2 * as.numeric(coordinate_count)) / remaining_lower
  for (index in seq_len(8L)) {
    target_upper <- .dsvert_dp_vector_next_up(target_upper)
  }
  log2_upper <- log2(target_upper)
  for (index in seq_len(8L)) {
    log2_upper <- .dsvert_dp_vector_next_up(log2_upper)
  }
  k <- ceiling(log2_upper)
  if (!is.finite(k) || k < 1) return(support_result())

  exponent_upper <- 7 * k / 10
  for (index in seq_len(8L)) {
    exponent_upper <- .dsvert_dp_vector_next_up(exponent_upper)
  }
  threshold_upper <- 4 * sigma_squared
  for (index in seq_len(8L)) {
    threshold_upper <- .dsvert_dp_vector_next_up(threshold_upper)
  }
  threshold_upper <- threshold_upper * exponent_upper
  for (index in seq_len(8L)) {
    threshold_upper <- .dsvert_dp_vector_next_up(threshold_upper)
  }
  radius_upper <- sqrt(threshold_upper)
  for (index in seq_len(8L)) {
    radius_upper <- .dsvert_dp_vector_next_up(radius_upper)
  }
  if (!is.finite(radius_upper) || radius_upper < 0) {
    return(support_result())
  }
  steps <- ceiling(radius_upper)
  if (!is.finite(steps) || steps >= support) return(support_result())
  list(steps = unname(steps), finite_support = FALSE)
}

.dsvert_dp_vector_accuracy_radius <- function(
    release, manifest, coordinate_count = 1L, confidence = 0.95,
    maximum_error = Inf) {
  lattice <- manifest$workload$release_lattice
  profile <- .dsvert_dp_vector_release_profile(release, manifest)
  sensitivity <- as.numeric(if (isTRUE(profile$gaussian)) {
    lattice$integer_l2_sensitivity_steps
  } else {
    lattice$integer_l1_sensitivity_steps
  })
  natural_sensitivity <- as.numeric(if (isTRUE(profile$gaussian)) {
    lattice$natural_l2_sensitivity
  } else {
    lattice$natural_l1_sensitivity
  })
  scale <- as.numeric(lattice$output_lattice_scale)
  epsilon <- as.numeric(release$epsilon)
  if (!.dsvert_dp_is_integer(coordinate_count, 1,
                             .DSVERT_DP_MAX_COORDINATES) ||
      !.dsvert_dp_is_number(confidence, 0, 1, lower_open = TRUE) ||
      confidence >= 1 || !is.numeric(maximum_error) ||
      length(maximum_error) != 1L || is.na(maximum_error) ||
      maximum_error < 0 ||
      !.dsvert_dp_is_number(epsilon, 0, .DSVERT_DP_MAXIMUM_EPSILON,
                            lower_open = TRUE) ||
      !is.finite(sensitivity) || sensitivity <= 0 ||
      !is.finite(scale) || scale <= 0 ||
      !identical(as.numeric(lattice$output_lattice_bits), log2(scale)) ||
      !identical(as.numeric(sensitivity),
                 natural_sensitivity * scale) ||
      !identical(release$mechanism, profile$release_mechanism)) {
    stop("Invalid signed vector accuracy contract", call. = FALSE)
  }
  implementation_delta <- .dsvert_dp_vector_fraction(
    release$implementation_delta)
  if (isTRUE(profile$gaussian)) {
    selection <- profile$manifest_selection
    plan <- if (is.list(selection)) selection$gaussian_plan else NULL
    request <- if (is.list(selection)) {
      selection$gaussian_calibration_request
    } else NULL
    plan_hash <- if (is.list(selection)) {
      selection$gaussian_plan_sha256
    } else NULL
    sensitivity_text <- if (is.list(request)) {
      request$l2_sensitivity_steps
    } else NULL
    .dsvert_vector_plan_validate(
      plan, plan_hash, profile,
      manifest$workload$coordinate_count, sensitivity_text)
    implementation_tv <- .dsvert_dp_vector_fraction_upper(
      plan$vector_total_tv_upper_numerator,
      plan$vector_total_tv_upper_denominator)
    if (identical(as.numeric(confidence), 0.95)) {
      steps <- suppressWarnings(as.numeric(plan$simultaneous_95_abs))
      if (!is.finite(steps) || steps < 0 || steps > 2^53 - 1) {
        stop("The fixed-work Gaussian v2 accuracy radius is not representable",
             call. = FALSE)
      }
      finite_support <- FALSE
      method <- paste(
        "signed fixed-work dyadic discrete-Gaussian plan v2 simultaneous",
        "95% bound; tail and CDF TV transfers already charged;",
        "fixed-clamp range applied")
    } else {
      derived <- .dsvert_dp_vector_gaussian_accuracy_steps(
        plan, coordinate_count, confidence)
      steps <- derived$steps
      finite_support <- isTRUE(derived$finite_support)
      if (confidence > 0.95) {
        signed_95_steps <- suppressWarnings(
          as.numeric(plan$simultaneous_95_abs))
        if (!is.finite(signed_95_steps) || signed_95_steps < 0 ||
            signed_95_steps > 2^53 - 1) {
          stop("The fixed-work Gaussian v2 accuracy radius is not representable",
               call. = FALSE)
        }
        steps <- max(steps, signed_95_steps)
      }
      method <- if (isTRUE(finite_support)) {
        paste(
          "signed fixed-work dyadic discrete-Gaussian plan v2 absolute",
          "finite-support bound; valid at every confidence;",
          "fixed-clamp range applied")
      } else {
        paste(
          "confidence-specific signed fixed-work dyadic discrete-Gaussian",
          "plan v2 subgaussian bound; signed rational variance and total",
          "variation transfer with conservative outward rounding;",
          "fixed-clamp range applied")
      }
    }
    radius <- steps / scale
    if (!identical(as.numeric(confidence), 0.95)) {
      radius <- .dsvert_dp_vector_next_up(radius)
    }
    if (is.finite(maximum_error)) radius <- min(radius, maximum_error)
    return(list(
      radius = unname(radius), confidence = confidence,
      coordinate_count = as.integer(coordinate_count),
      method = method,
      implementation_delta_bound = implementation_delta,
      implementation_tv_upper_bound = implementation_tv,
      sampler_tv_upper_bound = implementation_tv,
      finite_support_fallback = finite_support,
      accuracy_plan_certified = TRUE,
      additional_privacy_cost = c(epsilon = 0, delta = 0)))
  }
  plan <- release$mechanism_plan
  if (!is.list(plan)) {
    prepare <- tryCatch(
      release$signed_provenance$prepare_receipts[[1L]],
      error = function(error) NULL)
    plan <- if (is.list(prepare)) prepare$mechanism_plan else NULL
  }
  certified_plan <- is.list(plan)
  if (isTRUE(certified_plan)) {
    plan_sha256 <- release$plan_sha256
    if (!.dsvert_vector_hex(plan_sha256)) {
      plan_sha256 <- .dsvert_vector_hash(plan)
    }
    .dsvert_vector_plan_validate(
      plan, plan_sha256, profile,
      manifest$workload$coordinate_count,
      as.character(lattice$integer_l1_sensitivity_steps))
    delta_fields <- if (isTRUE(profile$exact_gc)) {
      c("implementation_delta_numerator",
        "implementation_delta_denominator")
    } else {
      c("per_peer_implementation_delta_numerator",
        "per_peer_implementation_delta_denominator")
    }
    signed_delta <- paste0(plan[[delta_fields[[1L]]]], "/",
                           plan[[delta_fields[[2L]]]])
    if (!identical(signed_delta, release$implementation_delta) ||
        (is.list(release$backend_selection) &&
         !identical(release$backend_selection$exact_gc_plan_sha256,
                    if (isTRUE(profile$exact_gc)) plan_sha256 else
                      release$backend_assessment$plan_sha256))) {
      stop("The vector accuracy law is detached from its signed plan",
           call. = FALSE)
    }
    if (as.numeric(coordinate_count) >
        as.numeric(plan$total_coordinate_count)) {
      stop("The requested accuracy dimension exceeds the signed vector",
           call. = FALSE)
    }
  }
  implementation_tv <- if (isTRUE(certified_plan)) {
    .dsvert_dp_vector_sampler_tv_upper(plan, profile$exact_gc)
  } else {
    # Compatibility for pre-v5 public records, which did not carry the full
    # sampler plan. New releases never enter this branch.
    min(.dsvert_dp_vector_next_up(2 * implementation_delta),
        1 - .Machine$double.eps)
  }
  # Utility transfers the signed finite sampler to its ideal law. The exact-GC
  # route emits one joint draw and therefore charges its vector TV once; the
  # certified convolution fallback charges the two independent peer TVs.
  alpha <- .dsvert_dp_vector_next_down(
    .dsvert_dp_vector_next_down(1 - confidence - implementation_tv))
  if (!is.finite(alpha) || alpha <= 0) {
    stop("The vector sampler cannot certify the requested confidence",
         call. = FALSE)
  }
  target_probability <- .dsvert_dp_vector_next_down(
    alpha / as.numeric(coordinate_count))
  target <- .dsvert_dp_vector_next_down(log(target_probability))
  natural_scale <- sensitivity / epsilon
  high <- max(0, ceiling(
    2 * natural_scale * log(4 * as.numeric(coordinate_count) / alpha)))
  tail <- if (isTRUE(certified_plan)) {
    .dsvert_dp_vector_dyadic_tail_context(plan)
  } else NULL
  log_tail <- if (isTRUE(certified_plan)) {
    function(radius, epsilon, sensitivity) {
      .dsvert_dp_vector_plan_log_tail_upper(
        radius, tail, convolution = !isTRUE(profile$exact_gc))
    }
  } else if (isTRUE(profile$exact_gc)) {
    .dsvert_dp_vector_laplace_log_tail
  } else {
    .dsvert_dp_vector_convolution_log_tail
  }
  while (log_tail(
      high, epsilon, sensitivity) > target) {
    high <- ceiling(2 * high + 1)
    if (!is.finite(high)) {
      stop("The vector accuracy radius is not representable", call. = FALSE)
    }
  }
  low <- -1
  while (high - low > 1) {
    middle <- low + floor((high - low) / 2)
    if (log_tail(
        middle, epsilon, sensitivity) <= target) {
      high <- middle
    } else {
      low <- middle
    }
  }
  radius <- high / scale
  if (is.finite(maximum_error)) radius <- min(radius, maximum_error)
  list(
    radius = unname(radius), confidence = confidence,
    coordinate_count = as.integer(coordinate_count),
    method = if (isTRUE(profile$exact_gc)) paste(
      "exact ideal one-draw two-sided-geometric tail with union bound;",
      "signed vector sampler TV deducted once; clamp inside exact GC applied")
    else paste(
      "exact ideal two-sided-geometric convolution tail with union bound;",
      "two-peer finite-sampler TV deducted; fixed-clamp range applied"),
    implementation_delta_bound = implementation_delta,
    implementation_tv_upper_bound = implementation_tv,
    sampler_tv_upper_bound = implementation_tv,
    accuracy_plan_certified = isTRUE(certified_plan),
    additional_privacy_cost = c(epsilon = 0, delta = 0))
}
