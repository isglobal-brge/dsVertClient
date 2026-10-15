# Independent public arithmetic for dp_irls_v1. This file has no process,
# network, Go, BLAS, or floating-point arithmetic dependency. Machine integers
# are used only for bounded dimensions and bit counts, never signed values.

.dsvert_dp_irls_arithmetic_fail <- function(message) {
  stop(paste0("Invalid dp_irls_update_v1: ", message), call. = FALSE)
}

.dsvert_dp_irls_arithmetic_json <- function(value) {
  order_fields <- function(x) {
    if (!is.list(x)) return(x)
    if (!is.null(names(x))) x <- x[order(names(x), method = "radix")]
    lapply(x, order_fields)
  }
  as.character(jsonlite::toJSON(order_fields(value), auto_unbox = TRUE,
                               null = "null", digits = NA))
}

.dsvert_dp_irls_arithmetic_fields <- function(x, expected, label) {
  if (!is.list(x) || is.null(names(x)) || anyDuplicated(names(x)) ||
      !setequal(names(x), expected)) {
    .dsvert_dp_irls_arithmetic_fail(paste(label, "fields"))
  }
  invisible(TRUE)
}

.dsvert_dp_irls_arithmetic_integer <- function(x, nonnegative = FALSE) {
  if (!is.character(x) || length(x) != 1L || is.na(x) ||
      nchar(x, type = "bytes") > 1500L ||
      !grepl("^(0|[1-9][0-9]*|-[1-9][0-9]*)$", x)) {
    .dsvert_dp_irls_arithmetic_fail("noncanonical integer")
  }
  out <- gmp::as.bigz(x)
  if (nonnegative && out < 0) .dsvert_dp_irls_arithmetic_fail("negative integer")
  out
}

.dsvert_dp_irls_arithmetic_rational <- function(x) {
  .dsvert_dp_irls_arithmetic_fields(x, c("d", "n"), "rational")
  n <- .dsvert_dp_irls_arithmetic_integer(x$n)
  d <- .dsvert_dp_irls_arithmetic_integer(x$d, TRUE)
  if (d <= 0 || gmp::gcd.bigz(abs(n), d) != 1) {
    .dsvert_dp_irls_arithmetic_fail("nonreduced rational")
  }
  gmp::as.bigq(n, d)
}

.dsvert_dp_irls_arithmetic_pair <- function(x) {
  list(d = as.character(gmp::denominator(x)),
       n = as.character(gmp::numerator(x)))
}

.dsvert_dp_irls_arithmetic_bits <- function(x) {
  if (x == 0) return(0L)
  nchar(as.character(abs(x), b = 2L), type = "bytes")
}

.dsvert_dp_irls_arithmetic_math <- function(bits) {
  check <- function(x) {
    if (.dsvert_dp_irls_arithmetic_bits(x) > bits) {
      .dsvert_dp_irls_arithmetic_fail("temporary bit ceiling exceeded")
    }
    x
  }
  add <- function(a, b) {
    an <- gmp::numerator(a); ad <- gmp::denominator(a)
    bn <- gmp::numerator(b); bd <- gmp::denominator(b)
    gmp::as.bigq(check(check(an * bd) + check(bn * ad)), check(ad * bd))
  }
  mul <- function(a, b) gmp::as.bigq(
    check(gmp::numerator(a) * gmp::numerator(b)),
    check(gmp::denominator(a) * gmp::denominator(b)))
  div <- function(a, b) {
    if (b == 0) .dsvert_dp_irls_arithmetic_fail("division by zero")
    gmp::as.bigq(check(gmp::numerator(a) * gmp::denominator(b)),
                check(gmp::denominator(a) * gmp::numerator(b)))
  }
  compare <- function(a, b) {
    left <- check(gmp::numerator(a) * gmp::denominator(b))
    right <- check(gmp::numerator(b) * gmp::denominator(a))
    if (left < right) -1L else if (left > right) 1L else 0L
  }
  list(check = check, add = add, mul = mul, div = div, compare = compare)
}

.dsvert_dp_irls_arithmetic_parse <- function(request, with_limits = TRUE) {
  if (!requireNamespace("gmp", quietly = TRUE)) {
    .dsvert_dp_irls_arithmetic_fail("gmp is required for independent verification")
  }
  fields <- c("beta", "box_lower", "box_upper", "damping", "floor",
              "floor_quantum", "g", "g_translation", "g_unit", "h",
              "h_translation", "h_unit", "lattice", "ridge", "trust", "version")
  if (with_limits) fields <- c(fields, "limits")
  .dsvert_dp_irls_arithmetic_fields(request, fields, "request")
  ceiling <- 5000000L
  if (with_limits) {
    .dsvert_dp_irls_arithmetic_fields(request$limits,
      c("common_denominator", "input_bits", "k_max", "temporary_bits"), "limits")
    supplied <- .dsvert_dp_irls_arithmetic_integer(request$limits$temporary_bits, TRUE)
    if (supplied <= 0 || supplied > ceiling) {
      .dsvert_dp_irls_arithmetic_fail("temporary bit ceiling")
    }
    ceiling <- as.integer(supplied)
  }
  math <- .dsvert_dp_irls_arithmetic_math(ceiling)
  parse_q <- function(x) {
    value <- .dsvert_dp_irls_arithmetic_rational(x)
    math$check(gmp::numerator(value)); math$check(gmp::denominator(value))
    value
  }
  if (!identical(request$version, "dp_irls_update_v1")) {
    .dsvert_dp_irls_arithmetic_fail("version")
  }
  p <- length(request$beta)
  if (p < 1L || p > 128L) .dsvert_dp_irls_arithmetic_fail("dimension")
  vector <- function(x, size, parse) {
    if (!is.list(x) || !is.null(names(x)) || length(x) != size) {
      .dsvert_dp_irls_arithmetic_fail("vector shape")
    }
    lapply(x, parse)
  }
  r <- list(p = p, math = math)
  for (name in c("beta", "box_lower", "box_upper")) {
    r[[name]] <- vector(request[[name]], p, parse_q)
  }
  for (name in c("damping", "floor", "floor_quantum", "g_unit", "h_unit",
                 "lattice", "ridge", "trust")) {
    r[[name]] <- parse_q(request[[name]])
  }
  for (name in c("g", "g_translation", "h", "h_translation")) {
    size <- if (substr(name, 1L, 1L) == "g") p else (p * (p + 1L)) %/% 2L
    r[[name]] <- vector(request[[name]], size, function(x) {
      math$check(.dsvert_dp_irls_arithmetic_integer(x, TRUE))
    })
  }
  for (name in c("floor", "floor_quantum", "g_unit", "h_unit", "lattice")) {
    if (r[[name]] <= 0) .dsvert_dp_irls_arithmetic_fail("nonpositive constant")
  }
  if (r$damping <= 0 || r$damping > 1 || r$ridge < 0 || r$trust < 0) {
    .dsvert_dp_irls_arithmetic_fail("damping, ridge or trust range")
  }
  for (i in seq_len(p)) {
    if (math$compare(r$box_lower[[i]], r$box_upper[[i]]) > 0L ||
        math$compare(r$beta[[i]], r$box_lower[[i]]) < 0L ||
        math$compare(r$beta[[i]], r$box_upper[[i]]) > 0L) {
      .dsvert_dp_irls_arithmetic_fail("coefficient box")
    }
    for (name in c("beta", "box_lower", "box_upper")) {
      if (gmp::denominator(math$div(r[[name]][[i]], r$lattice)) != 1) {
        .dsvert_dp_irls_arithmetic_fail("coefficient off lattice")
      }
    }
  }
  r$saturated <- FALSE
  for (lane in c("g", "h")) {
    for (i in seq_along(r[[lane]])) {
      a <- r[[paste0(lane, "_translation")]][[i]]
      v <- r[[lane]][[i]]
      if (a <= 0) .dsvert_dp_irls_arithmetic_fail("nonpositive translation")
      upper <- math$check(2L * a)
      if (v > upper) .dsvert_dp_irls_arithmetic_fail("finalizer bound")
      r$saturated <- r$saturated || v == 0 || v == upper
    }
  }
  r
}

# Derive the public numeric envelope independently of the signed limits.
.dsvert_dp_irls_arithmetic_envelope <- function(r, signed_w = NULL) {
  z <- gmp::as.bigz
  q <- gmp::as.bigq
  math <- r$math
  maximum <- function(a, b) if (math$compare(a, b) > 0L) a else b
  d <- z(1L)
  for (v in c(list(r$h_unit, r$floor_quantum, r$ridge, r$floor, r$g_unit),
              list(math$mul(r$ridge, r$lattice)))) {
    vd <- gmp::denominator(v)
    d <- math$check((d %/% gmp::gcd.bigz(d, vd)) * vd)
  }
  row_bounds <- rep(list(q(0L)), r$p)
  hmax <- q(0L)
  index <- 0L
  for (i in seq_len(r$p)) for (j in i:r$p) {
    index <- index + 1L
    bound <- math$mul(q(r$h_translation[[index]]), r$h_unit)
    row_bounds[[i]] <- math$add(row_bounds[[i]], bound)
    if (i != j) row_bounds[[j]] <- math$add(row_bounds[[j]], bound)
    hmax <- maximum(hmax, bound)
  }
  rh <- Reduce(maximum, row_bounds)
  ceiling_q <- function(v) {
    n <- gmp::numerator(v); dd <- gmp::denominator(v)
    math$check(math$check(n + dd) - 1L) %/% dd
  }
  kmax <- ceiling_q(math$div(math$add(r$floor, rh), r$floor_quantum))
  boxmax <- Reduce(maximum, lapply(c(r$box_lower, r$box_upper), abs))
  gmax <- Reduce(maximum, lapply(r$g_translation, function(a) math$mul(q(a), r$g_unit)))
  hbound <- math$add(math$add(math$add(hmax,
    math$mul(q(kmax), r$floor_quantum)), r$ridge), r$floor)
  gbound <- math$add(gmax, math$mul(r$ridge, boxmax))
  bound <- maximum(hbound, gbound)
  raw_q <- c(r$box_lower, r$box_upper,
             list(r$damping, r$floor, r$floor_quantum, r$g_unit, r$h_unit,
                  r$lattice, r$ridge, r$trust))
  raw_z <- c(lapply(raw_q, gmp::numerator), lapply(raw_q, gmp::denominator),
             lapply(c(r$g_translation, r$h_translation), function(a) math$check(2L * a)),
             list(d, ceiling_q(math$mul(q(d), bound)),
                  ceiling_q(math$mul(boxmax, q(gmp::denominator(r$lattice))))))
  w <- max(1L, vapply(raw_z, .dsvert_dp_irls_arithmetic_bits, integer(1L)))
  if (!is.null(signed_w)) {
    supplied <- .dsvert_dp_irls_arithmetic_integer(signed_w, TRUE)
    if (supplied < w || supplied > 4096L) {
      .dsvert_dp_irls_arithmetic_fail("input bit certificate")
    }
    w <- as.integer(supplied)
  }
  if (w > 4096L) .dsvert_dp_irls_arithmetic_fail("input bit certificate")
  m <- 0L
  while (z(2L) ^ (2L * m) < z(r$p) ^ r$p) m <- m + 1L
  ell <- r$p * w + m + 2L
  temporary <- 8L * ell + 16L * w + 64L
  list(common_denominator = as.character(d), input_bits = as.character(w),
       k_max = as.character(kmax), temporary_bits = as.character(temporary))
}

.dsvert_dp_irls_numeric_limits <- function(request) {
  request$limits <- NULL
  .dsvert_dp_irls_arithmetic_envelope(
    .dsvert_dp_irls_arithmetic_parse(request, with_limits = FALSE))
}

.dsvert_dp_irls_verify_update <- function(request) {
  r <- .dsvert_dp_irls_arithmetic_parse(request)
  .dsvert_dp_irls_arithmetic_fields(request$limits,
    c("common_denominator", "input_bits", "k_max", "temporary_bits"), "limits")
  for (x in request$limits) .dsvert_dp_irls_arithmetic_integer(x, TRUE)
  expected <- .dsvert_dp_irls_arithmetic_envelope(r, request$limits$input_bits)
  if (!identical(.dsvert_dp_irls_arithmetic_json(expected),
                 .dsvert_dp_irls_arithmetic_json(request$limits))) {
    .dsvert_dp_irls_arithmetic_fail("numeric certificate limits")
  }
  ceiling <- as.integer(expected$temporary_bits)
  if (is.na(ceiling) || ceiling > 5000000L) {
    .dsvert_dp_irls_arithmetic_fail("temporary bit ceiling")
  }
  check_z <- r$math$check
  add <- r$math$add
  mul <- r$math$mul
  div <- r$math$div
  compare <- r$math$compare
  d <- gmp::as.bigz(expected$common_denominator)
  kmax <- gmp::as.bigz(expected$k_max)
  w <- as.integer(expected$input_bits)
  minor_bits <- (ceiling - 16L * w - 64L) %/% 8L
  minor_check <- function(x, bits = minor_bits) {
    if (.dsvert_dp_irls_arithmetic_bits(x) > bits) {
      .dsvert_dp_irls_arithmetic_fail("Bareiss bit ceiling")
    }
    check_z(x)
  }
  h <- matrix(rep(list(gmp::as.bigq(0L)), r$p * r$p), r$p, r$p)
  index <- 0L
  for (i in seq_len(r$p)) for (j in i:r$p) {
    index <- index + 1L
    h[[i, j]] <- mul(gmp::as.bigq(r$h[[index]] - r$h_translation[[index]]),
                      r$h_unit)
    h[[j, i]] <- h[[i, j]]
  }
  g <- lapply(seq_len(r$p), function(i) {
    mul(gmp::as.bigq(r$g[[i]] - r$g_translation[[i]]), r$g_unit)
  })
  clear <- function(v) {
    value <- mul(v, gmp::as.bigq(d))
    if (gmp::denominator(value) != 1) {
      .dsvert_dp_irls_arithmetic_fail("common denominator")
    }
    n <- gmp::numerator(value)
    if (.dsvert_dp_irls_arithmetic_bits(n) > w) {
      .dsvert_dp_irls_arithmetic_fail("input bit ceiling")
    }
    check_z(n)
  }
  shifted <- function(k, solve = FALSE) {
    a <- h
    delta <- add(mul(gmp::as.bigq(k), r$floor_quantum),
                 if (solve) r$ridge else -r$floor)
    for (i in seq_len(r$p)) a[[i, i]] <- add(a[[i, i]], delta)
    for (i in seq_len(r$p)) for (j in seq_len(r$p)) a[[i, j]] <- clear(a[[i, j]])
    a
  }
  bareiss_entry <- function(pivot, entry, column, row, previous) {
    n <- minor_check(minor_check(pivot * entry, 2L * minor_bits + 2L) -
                       minor_check(column * row, 2L * minor_bits + 2L),
                     2L * minor_bits + 2L)
    if (n %% previous != 0) .dsvert_dp_irls_arithmetic_fail("inexact Bareiss division")
    minor_check(n %/% previous)
  }
  psd <- function(k) {
    a <- shifted(k)
    previous <- gmp::as.bigz(1L)
    for (i in seq_len(r$p)) {
      pivot <- a[[i, i]]
      if (pivot < 0) return(FALSE)
      if (pivot == 0) {
        for (j in i:r$p) if (a[[i, j]] != 0) return(FALSE)
        # Erratum 7: skipping a zero residual row does not change the
        # previous successful pivot; in particular it never becomes zero.
        next
      }
      if (i < r$p) for (j in (i + 1L):r$p) for (l in j:r$p) {
        a[[j, l]] <- bareiss_entry(pivot, a[[j, l]], a[[j, i]], a[[i, l]], previous)
        a[[l, j]] <- a[[j, l]]
      }
      previous <- pivot
    }
    TRUE
  }
  k <- gmp::as.bigz(0L)
  if (!psd(k)) {
    lower <- k
    upper <- gmp::as.bigz(1L)
    repeat {
      if (upper > kmax) upper <- kmax
      if (psd(upper)) break
      if (upper == kmax) .dsvert_dp_irls_arithmetic_fail("k_max exhausted")
      lower <- upper
      upper <- check_z(2L * upper)
    }
    while (upper - lower > 1L) {
      middle <- (lower + upper) %/% 2L
      if (psd(middle)) upper <- middle else lower <- middle
    }
    k <- upper
  }
  a <- shifted(k, solve = TRUE)
  rhs <- lapply(seq_len(r$p), function(i) clear(add(g[[i]], -mul(r$ridge, r$beta[[i]]))))
  a <- cbind(a, matrix(rhs, r$p, 1L))
  previous <- gmp::as.bigz(1L)
  for (i in seq_len(r$p)) {
    pivot <- a[[i, i]]
    if (pivot <= 0) .dsvert_dp_irls_arithmetic_fail("nonpositive solve pivot")
    if (i < r$p) for (j in (i + 1L):r$p) {
      for (l in (i + 1L):(r$p + 1L)) {
        a[[j, l]] <- bareiss_entry(pivot, a[[j, l]], a[[j, i]], a[[i, l]], previous)
      }
      a[[j, i]] <- gmp::as.bigz(0L)
    }
    previous <- pivot
  }
  numerators <- rep(list(gmp::as.bigz(0L)), r$p)
  determinant <- a[[r$p, r$p]]
  for (i in rev(seq_len(r$p))) {
    v <- check_z(a[[i, r$p + 1L]] * determinant)
    if (i < r$p) for (j in (i + 1L):r$p) {
      v <- check_z(v - check_z(a[[i, j]] * numerators[[j]]))
    }
    if (v %% a[[i, i]] != 0) .dsvert_dp_irls_arithmetic_fail("inexact back substitution")
    numerators[[i]] <- minor_check(v %/% a[[i, i]])
  }
  step <- lapply(numerators, function(n) gmp::as.bigq(n, determinant))
  step <- lapply(step, function(v) mul(r$damping, v))
  norm <- Reduce(function(a, b) if (compare(a, b) > 0L) a else b, lapply(step, abs))
  trusted <- compare(norm, r$trust) > 0L
  if (trusted) {
    scale <- div(r$trust, norm)
    step <- lapply(step, function(v) mul(scale, v))
  }
  projected <- rep(list(FALSE), r$p)
  beta <- lapply(seq_len(r$p), function(i) {
    candidate <- add(r$beta[[i]], step[[i]])
    if (compare(candidate, r$box_lower[[i]]) < 0L) {
      candidate <- r$box_lower[[i]]; projected[[i]] <<- TRUE
    }
    if (compare(candidate, r$box_upper[[i]]) > 0L) {
      candidate <- r$box_upper[[i]]; projected[[i]] <<- TRUE
    }
    value <- div(candidate, r$lattice)
    n <- gmp::numerator(value); dd <- gmp::denominator(value)
    integer <- n %/% dd
    rem <- n %% dd
    twice <- check_z(2L * rem)
    if (twice > dd || (twice == dd && integer %% 2L != 0)) integer <- integer + 1L
    .dsvert_dp_irls_arithmetic_pair(mul(gmp::as.bigq(integer), r$lattice))
  })
  list(beta = beta, k = as.character(k), projected = projected,
       saturated = r$saturated, trust_scaled = trusted, version = "dp_irls_update_v1")
}

.dsvert_dp_irls_verify_update_json <- function(canonical_json) {
  if (!is.character(canonical_json) || length(canonical_json) != 1L ||
      is.na(canonical_json) || nchar(canonical_json, type = "bytes") > 4194304L) {
    .dsvert_dp_irls_arithmetic_fail("JSON input")
  }
  request <- tryCatch(jsonlite::fromJSON(canonical_json, simplifyVector = FALSE),
    error = function(e) .dsvert_dp_irls_arithmetic_fail("JSON decoding"))
  if (!identical(.dsvert_dp_irls_arithmetic_json(request), canonical_json)) {
    .dsvert_dp_irls_arithmetic_fail("noncanonical JSON")
  }
  .dsvert_dp_irls_arithmetic_json(.dsvert_dp_irls_verify_update(request))
}
