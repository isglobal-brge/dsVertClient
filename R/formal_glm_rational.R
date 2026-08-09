# Exact signed rational arithmetic for the formal GLM reference compiler.
#
# OpenSSL's BIGNUM type gives us arbitrary-precision unsigned integers without
# adding another package dependency.  Signs are kept separately because the R
# binding deliberately rejects negative BIGNUM values.  These helpers are
# internal and intentionally small; they are a specification oracle, not the
# production MPC arithmetic backend.

.dsvert_glm_bn <- function(value) {
  openssl::bignum(as.character(value))
}

.dsvert_glm_bn_zero <- function() .dsvert_glm_bn("0")
.dsvert_glm_bn_one <- function() .dsvert_glm_bn("1")

.dsvert_glm_bn_gcd <- function(left, right) {
  zero <- .dsvert_glm_bn_zero()
  while (right != zero) {
    remainder <- left %% right
    left <- right
    right <- remainder
  }
  left
}

.dsvert_glm_rat_new <- function(sign, numerator, denominator) {
  zero <- .dsvert_glm_bn_zero()
  if (denominator == zero) {
    stop("A formal GLM rational denominator cannot be zero", call. = FALSE)
  }
  if (numerator == zero) {
    return(structure(list(
      sign = 0L, numerator = zero, denominator = .dsvert_glm_bn_one()),
      class = "dsvert_glm_rational"))
  }
  divisor <- .dsvert_glm_bn_gcd(numerator, denominator)
  structure(list(
    sign = if (sign < 0L) -1L else 1L,
    numerator = numerator %/% divisor,
    denominator = denominator %/% divisor),
    class = "dsvert_glm_rational")
}

.dsvert_glm_rat_from_decimal <- function(value) {
  if (is.numeric(value)) {
    if (length(value) != 1L || !is.finite(value)) {
      stop("Formal GLM numeric values must be finite scalars", call. = FALSE)
    }
    value <- sprintf("%.17g", value)
  }
  if (!is.character(value) || length(value) != 1L || is.na(value)) {
    stop("Formal GLM rational values must be decimal scalars", call. = FALSE)
  }
  value <- trimws(value)
  if (!nzchar(value) || nchar(value, type = "bytes") > 4096L) {
    stop("Formal GLM decimal rational is too large", call. = FALSE)
  }
  match <- regexec(
    "^([+-]?)([0-9]+)(?:\\.([0-9]*))?(?:[eE]([+-]?[0-9]+))?$",
    value, perl = TRUE)
  fields <- regmatches(value, match)[[1L]]
  if (!length(fields)) {
    stop("Invalid formal GLM decimal rational", call. = FALSE)
  }
  sign <- if (identical(fields[[2L]], "-")) -1L else 1L
  fraction <- fields[[4L]]
  if (is.na(fraction)) fraction <- ""
  exponent <- fields[[5L]]
  if (is.na(exponent) || !nzchar(exponent)) exponent <- "0"
  exponent <- suppressWarnings(as.integer(exponent))
  if (is.na(exponent) || abs(exponent) > 4096L) {
    stop("Formal GLM decimal exponent is out of range", call. = FALSE)
  }
  digits <- sub("^0+(?=[0-9])", "", paste0(fields[[3L]], fraction),
                perl = TRUE)
  numerator <- .dsvert_glm_bn(digits)
  power <- exponent - nchar(fraction, type = "bytes")
  ten <- .dsvert_glm_bn("10")
  if (power >= 0L) {
    numerator <- numerator * (ten ^ power)
    denominator <- .dsvert_glm_bn_one()
  } else {
    denominator <- ten ^ (-power)
  }
  .dsvert_glm_rat_new(sign, numerator, denominator)
}

.dsvert_glm_rat <- function(value) {
  if (inherits(value, "dsvert_glm_rational")) return(value)
  if (is.list(value) && identical(sort(names(value)),
                                  c("denominator", "numerator"))) {
    numerator <- value$numerator
    denominator <- value$denominator
    if (!is.character(numerator) || length(numerator) != 1L ||
        !grepl("^[+-]?[0-9]+$", numerator) ||
        !is.character(denominator) || length(denominator) != 1L ||
        !grepl("^[0-9]+$", denominator)) {
      stop("Invalid canonical formal GLM rational", call. = FALSE)
    }
    sign <- if (startsWith(numerator, "-")) -1L else 1L
    numerator <- sub("^[+-]", "", numerator)
    return(.dsvert_glm_rat_new(
      sign, .dsvert_glm_bn(numerator), .dsvert_glm_bn(denominator)))
  }
  .dsvert_glm_rat_from_decimal(value)
}

.dsvert_glm_rat_json <- function(value) {
  value <- .dsvert_glm_rat(value)
  numerator <- as.character(value$numerator)
  if (value$sign < 0L) numerator <- paste0("-", numerator)
  list(numerator = numerator, denominator = as.character(value$denominator))
}

.dsvert_glm_rat_neg <- function(value) {
  value <- .dsvert_glm_rat(value)
  .dsvert_glm_rat_new(-value$sign, value$numerator, value$denominator)
}

.dsvert_glm_rat_abs <- function(value) {
  value <- .dsvert_glm_rat(value)
  .dsvert_glm_rat_new(abs(value$sign), value$numerator, value$denominator)
}

.dsvert_glm_rat_add <- function(left, right) {
  left <- .dsvert_glm_rat(left)
  right <- .dsvert_glm_rat(right)
  if (left$sign == 0L) return(right)
  if (right$sign == 0L) return(left)
  left_scaled <- left$numerator * right$denominator
  right_scaled <- right$numerator * left$denominator
  denominator <- left$denominator * right$denominator
  if (left$sign == right$sign) {
    return(.dsvert_glm_rat_new(
      left$sign, left_scaled + right_scaled, denominator))
  }
  if (left_scaled == right_scaled) return(.dsvert_glm_rat("0"))
  if (left_scaled > right_scaled) {
    .dsvert_glm_rat_new(left$sign, left_scaled - right_scaled, denominator)
  } else {
    .dsvert_glm_rat_new(right$sign, right_scaled - left_scaled, denominator)
  }
}

.dsvert_glm_rat_sub <- function(left, right) {
  .dsvert_glm_rat_add(left, .dsvert_glm_rat_neg(right))
}

.dsvert_glm_rat_mul <- function(left, right) {
  left <- .dsvert_glm_rat(left)
  right <- .dsvert_glm_rat(right)
  .dsvert_glm_rat_new(
    left$sign * right$sign,
    left$numerator * right$numerator,
    left$denominator * right$denominator)
}

.dsvert_glm_rat_div <- function(left, right) {
  left <- .dsvert_glm_rat(left)
  right <- .dsvert_glm_rat(right)
  if (right$sign == 0L) {
    stop("Division by zero in the formal GLM rational oracle", call. = FALSE)
  }
  .dsvert_glm_rat_new(
    left$sign * right$sign,
    left$numerator * right$denominator,
    left$denominator * right$numerator)
}

.dsvert_glm_rat_cmp <- function(left, right) {
  difference <- .dsvert_glm_rat_sub(left, right)
  difference$sign
}

.dsvert_glm_rat_min <- function(left, right) {
  if (.dsvert_glm_rat_cmp(left, right) <= 0L) left else right
}

.dsvert_glm_rat_max <- function(left, right) {
  if (.dsvert_glm_rat_cmp(left, right) >= 0L) left else right
}

.dsvert_glm_rat_pow <- function(value, exponent) {
  value <- .dsvert_glm_rat(value)
  if (length(exponent) != 1L || is.na(exponent) || exponent < 0 ||
      exponent != floor(exponent) || exponent > .Machine$integer.max) {
    stop("Invalid formal GLM rational exponent", call. = FALSE)
  }
  exponent <- as.integer(exponent)
  if (exponent == 0L) return(.dsvert_glm_rat("1"))
  .dsvert_glm_rat_new(
    if (value$sign < 0L && exponent %% 2L) -1L else 1L,
    value$numerator ^ exponent,
    value$denominator ^ exponent)
}

.dsvert_glm_rat_double <- function(value) {
  value <- .dsvert_glm_rat(value)
  if (value$sign == 0L) return(0)
  numerator <- as.character(value$numerator)
  denominator <- as.character(value$denominator)
  digits <- 16L
  numerator_digits <- min(digits, nchar(numerator, type = "bytes"))
  denominator_digits <- min(digits, nchar(denominator, type = "bytes"))
  leading_numerator <- as.numeric(substr(
    numerator, 1L, numerator_digits)) / (10 ^ (numerator_digits - 1L))
  leading_denominator <- as.numeric(substr(
    denominator, 1L, denominator_digits)) /
    (10 ^ (denominator_digits - 1L))
  exponent <- nchar(numerator, type = "bytes") -
    nchar(denominator, type = "bytes")
  value$sign * (leading_numerator / leading_denominator) * (10 ^ exponent)
}

.dsvert_glm_rat_round_dyadic <- function(value, bits) {
  value <- .dsvert_glm_rat(value)
  if (length(bits) != 1L || is.na(bits) || bits < 0L || bits > 256L ||
      bits != floor(bits)) {
    stop("Invalid formal GLM dyadic precision", call. = FALSE)
  }
  scale <- .dsvert_glm_bn("2") ^ as.integer(bits)
  scaled <- value$numerator * scale
  quotient <- scaled %/% value$denominator
  remainder <- scaled %% value$denominator
  twice <- remainder * .dsvert_glm_bn("2")
  odd <- (quotient %% .dsvert_glm_bn("2")) == .dsvert_glm_bn_one()
  if (twice > value$denominator ||
      (twice == value$denominator && isTRUE(odd))) {
    quotient <- quotient + .dsvert_glm_bn_one()
  }
  .dsvert_glm_rat_new(value$sign, quotient, scale)
}

.dsvert_glm_rat_ceil_abs_scaled <- function(value, bits) {
  value <- .dsvert_glm_rat_abs(value)
  scale <- .dsvert_glm_bn("2") ^ as.integer(bits)
  numerator <- value$numerator * scale
  quotient <- numerator %/% value$denominator
  if (numerator %% value$denominator != .dsvert_glm_bn_zero()) {
    quotient <- quotient + .dsvert_glm_bn_one()
  }
  quotient
}

.dsvert_glm_rat_interval_midpoint <- function(interval) {
  .dsvert_glm_rat_div(
    .dsvert_glm_rat_add(interval$lower, interval$upper), "2")
}

.dsvert_glm_rat_exp_interval <- function(value, precision_bits = 112L) {
  value <- .dsvert_glm_rat(value)
  if (value$sign < 0L) {
    positive <- .dsvert_glm_rat_exp_interval(
      .dsvert_glm_rat_neg(value), precision_bits)
    return(list(
      lower = .dsvert_glm_rat_div("1", positive$upper),
      upper = .dsvert_glm_rat_div("1", positive$lower),
      terms = positive$terms))
  }
  target <- .dsvert_glm_rat_div(
    "1", .dsvert_glm_rat_pow("2", precision_bits))
  total <- .dsvert_glm_rat("1")
  term <- .dsvert_glm_rat("1")
  maximum <- max(64L, as.integer(ceiling(.dsvert_glm_rat_double(value))) +
                   8L * precision_bits)
  for (index in seq_len(maximum)) {
    term <- .dsvert_glm_rat_div(
      .dsvert_glm_rat_mul(term, value), as.character(index))
    total <- .dsvert_glm_rat_add(total, term)
    ratio <- .dsvert_glm_rat_div(value, as.character(index + 1L))
    if (.dsvert_glm_rat_cmp(ratio, "1") < 0L) {
      next_term <- .dsvert_glm_rat_mul(term, ratio)
      remainder <- .dsvert_glm_rat_div(
        next_term, .dsvert_glm_rat_sub("1", ratio))
      if (.dsvert_glm_rat_cmp(remainder, target) <= 0L) {
        return(list(
          lower = total,
          upper = .dsvert_glm_rat_add(total, remainder),
          terms = as.integer(index)))
      }
    }
  }
  stop("High-precision exponential oracle did not converge", call. = FALSE)
}

.dsvert_glm_rat_log_interval <- function(value, precision_bits = 112L) {
  value <- .dsvert_glm_rat(value)
  if (value$sign <= 0L) {
    stop("The formal GLM logarithm requires a positive rational", call. = FALSE)
  }
  one <- .dsvert_glm_rat("1")
  two <- .dsvert_glm_rat("2")
  mantissa <- value
  exponent <- 0L
  while (.dsvert_glm_rat_cmp(mantissa, two) >= 0L) {
    mantissa <- .dsvert_glm_rat_div(mantissa, two)
    exponent <- exponent + 1L
    if (exponent > 4096L) stop("Logarithm range is unsupported", call. = FALSE)
  }
  while (.dsvert_glm_rat_cmp(mantissa, one) < 0L) {
    mantissa <- .dsvert_glm_rat_mul(mantissa, two)
    exponent <- exponent - 1L
    if (exponent < -4096L) stop("Logarithm range is unsupported", call. = FALSE)
  }
  series <- function(argument) {
    z <- .dsvert_glm_rat_div(
      .dsvert_glm_rat_sub(argument, one),
      .dsvert_glm_rat_add(argument, one))
    z2 <- .dsvert_glm_rat_mul(z, z)
    power <- z
    total <- z
    target <- .dsvert_glm_rat_div(
      "1", .dsvert_glm_rat_pow("2", precision_bits + 4L))
    for (index in seq_len(8L * precision_bits)) {
      power <- .dsvert_glm_rat_mul(power, z2)
      denominator <- 2L * index + 1L
      addend <- .dsvert_glm_rat_div(power, as.character(denominator))
      total <- .dsvert_glm_rat_add(total, addend)
      remainder <- .dsvert_glm_rat_div(
        .dsvert_glm_rat_mul(
          .dsvert_glm_rat_mul("2", power), z2),
        .dsvert_glm_rat_mul(
          as.character(denominator + 2L),
          .dsvert_glm_rat_sub(one, z2)))
      if (.dsvert_glm_rat_cmp(remainder, target) <= 0L) {
        centre <- .dsvert_glm_rat_mul("2", total)
        return(list(
          lower = centre,
          upper = .dsvert_glm_rat_add(centre, remainder)))
      }
    }
    stop("High-precision logarithm oracle did not converge", call. = FALSE)
  }
  log_mantissa <- series(mantissa)
  log_two <- series(two)
  if (exponent >= 0L) {
    lower <- .dsvert_glm_rat_add(
      log_mantissa$lower,
      .dsvert_glm_rat_mul(as.character(exponent), log_two$lower))
    upper <- .dsvert_glm_rat_add(
      log_mantissa$upper,
      .dsvert_glm_rat_mul(as.character(exponent), log_two$upper))
  } else {
    lower <- .dsvert_glm_rat_sub(
      log_mantissa$lower,
      .dsvert_glm_rat_mul(as.character(-exponent), log_two$upper))
    upper <- .dsvert_glm_rat_sub(
      log_mantissa$upper,
      .dsvert_glm_rat_mul(as.character(-exponent), log_two$lower))
  }
  list(lower = lower, upper = upper)
}
