# Canonical public projection of one signed Ring128 lattice coordinate.
#
# The integer text is the authoritative public value.  The adjacent double is
# display-only and must agree with its deterministic IEEE-754 projection; it
# is never trusted as an independent numeric source.
.dsvert_formal_lattice_value <- function(
    steps, fraction_bits, reported_value, label = "formal coefficient") {
  if (!is.character(steps) || length(steps) != 1L || is.na(steps) ||
      !grepl("^(0|-[1-9][0-9]*|[1-9][0-9]*)$", steps) ||
      !is.numeric(fraction_bits) || length(fraction_bits) != 1L ||
      is.na(fraction_bits) || !is.finite(fraction_bits) ||
      fraction_bits < 0 || fraction_bits > 127 ||
      fraction_bits != floor(fraction_bits) ||
      !is.numeric(reported_value) || length(reported_value) != 1L ||
      is.na(reported_value) || !is.finite(reported_value)) {
    stop("Invalid ", label, " lattice projection", call. = FALSE)
  }
  magnitude <- sub("^-", "", steps)
  ring128_limit <- if (startsWith(steps, "-")) {
    "170141183460469231731687303715884105728"
  } else {
    "170141183460469231731687303715884105727"
  }
  if (nchar(magnitude, type = "bytes") > nchar(ring128_limit) ||
      (nchar(magnitude, type = "bytes") == nchar(ring128_limit) &&
       magnitude > ring128_limit)) {
    stop("Invalid ", label, " lattice coordinate", call. = FALSE)
  }
  projected <- suppressWarnings(as.numeric(steps)) * 2^(-as.integer(fraction_bits))
  if (!is.finite(projected)) {
    stop("Invalid ", label, " lattice projection", call. = FALSE)
  }
  tolerance <- 8 * .Machine$double.eps * max(
    abs(projected), abs(reported_value), .Machine$double.xmin)
  if (abs(reported_value - projected) > tolerance) {
    stop("The ", label,
         " value does not match its signed lattice steps", call. = FALSE)
  }
  projected
}
