# Reproducible synthetic public arithmetic fixtures; no private data or noise.
# Run from the client repository root. Both implementations consume the exact
# same JSONL file, copied to the server's inst/dsvert-mpc/testdata directory.
source("R/dp_irls_verify.R")
q <- function(n, d = 1L) gmp::as.bigq(as.character(n), as.character(d))
pair <- .dsvert_dp_irls_arithmetic_pair
request <- function(h, g, beta = rep(0L, length(g)), lattice = q(1L),
                    floor = q(1L), quantum = q(1L), ridge = q(0L),
                    damping = q(1L), trust = q(1000L),
                    lower = rep(-1000L, length(g)), upper = -lower,
                    g_unit = q(1L), h_unit = q(1L),
                    g_translation = NULL, h_translation = NULL) {
  p <- length(g)
  hi <- unlist(lapply(seq_len(p), function(i) h[i, i:p]), use.names = FALSE)
  if (is.null(g_translation)) g_translation <- abs(g) + 1L
  if (is.null(h_translation)) h_translation <- abs(hi) + 1L
  ints <- function(x) lapply(x, as.character)
  out <- list(beta = lapply(beta, function(b) pair(q(b) * lattice)),
    box_lower = lapply(lower, function(b) pair(q(b) * lattice)),
    box_upper = lapply(upper, function(b) pair(q(b) * lattice)),
    damping = pair(damping), floor = pair(floor), floor_quantum = pair(quantum),
    g = ints(g + g_translation), g_translation = ints(g_translation),
    g_unit = pair(g_unit), h = ints(hi + h_translation),
    h_translation = ints(h_translation), h_unit = pair(h_unit),
    lattice = pair(lattice), ridge = pair(ridge), trust = pair(trust),
    version = "dp_irls_update_v1")
  out$limits <- .dsvert_dp_irls_numeric_limits(out)
  out
}
fixtures <- list()
add <- function(name, x, beta = NULL, k = NULL, projected = NULL, trust = NULL,
                saturated = NULL) {
  value <- .dsvert_dp_irls_verify_update(x)
  # Known-answer cases assert independently specified answers before writing.
  if (!is.null(beta)) stopifnot(identical(value$beta, lapply(beta, function(b) pair(q(b)))))
  if (!is.null(k)) stopifnot(identical(value$k, as.character(k)))
  if (!is.null(projected)) stopifnot(identical(value$projected, as.list(projected)))
  if (!is.null(trust)) stopifnot(identical(value$trust_scaled, trust))
  if (!is.null(saturated)) stopifnot(identical(value$saturated, saturated))
  fixtures[[length(fixtures) + 1L]] <<- list(expected = value, name = name, request = x)
}
add("positive-definite", request(matrix(c(2, 1, 1, 2), 2), c(3, 3)),
    beta = c(1, 1), k = 0, trust = FALSE)
add("semidefinite", request(diag(c(0, 2)), c(1, 3)), beta = c(1, 1), k = 1)
add("zero-matrix", request(matrix(0L, 3, 3), c(1, 2, -3)), beta = c(1, 2, -3), k = 1)
skip <- matrix(c(3, 0, 1, 0, 0, 1, 0, 0, 1, 0, 3, 1, 0, 0, 1, 3), 4)
add("zero-row-keeps-previous-pivot", request(skip, rowSums(skip)), beta = rep(1, 4), k = 0)
add("zero-pivot-nonzero-row", request(matrix(c(1, 1, 1, 2), 2), c(3, 4)),
    beta = c(1, 1), k = 1)
add("negative-pivot", request(diag(c(-2, 1)), c(2, 3)), beta = c(2, 1), k = 3)
for (shift in c(1L, 2L, 3L, 4L, 5L, 7L, 8L, 9L)) {
  add(paste0("bracket-", shift), request(matrix(1L - shift, 1), 1L), beta = 1L, k = shift)
}
add("kmax-endpoint", request(matrix(-4L, 1), 1L, h_translation = 4L),
    beta = 1L, k = 5L, saturated = TRUE)
add("damping-before-trust", request(diag(2), c(8, 4), damping = q(1, 4), trust = q(1)),
    beta = c(1, 0), k = 0, trust = TRUE)
add("damping-avoids-trust", request(diag(2), c(2, 0), damping = q(1, 4), trust = q(1)),
    beta = c(0, 0), k = 0, trust = FALSE)
add("zero-trust", request(diag(2), c(2, -2), trust = q(0)),
    beta = c(0, 0), k = 0, trust = TRUE)
add("zero-step-zero-trust", request(diag(2), c(0, 0), trust = q(0)),
    beta = c(0, 0), k = 0, trust = FALSE)
add("box-both-directions", request(diag(2), c(4, -4), lower = c(-1, -1), upper = c(1, 1)),
    beta = c(1, -1), projected = c(TRUE, TRUE), k = 0)
for (score in c(1L, 3L, -1L, -3L)) {
  expected <- switch(as.character(score), "1" = 0L, "3" = 2L, "-1" = 0L, "-3" = -2L)
  add(paste0("nearest-even-", score), request(matrix(2L, 1), score), beta = expected, k = 0)
}
add("ridge-right-hand-side", request(matrix(2L, 1), 5L, beta = 1L, ridge = q(1)),
    beta = 2L, k = 0)
add("rational-scaled-translation", request(diag(c(4, 8)), c(4, -8),
    g_unit = q(1, 2), h_unit = q(1, 4), lattice = q(1, 8), floor = q(1, 2),
    quantum = q(1, 4)), beta = c(2, -2), k = 0)
# Include exact integers well beyond binary64's precision.
large <- request(matrix(1L, 1), 1L)
large$g <- list("9007199254740994")
large$g_translation <- list("9007199254740993")
large$limits <- .dsvert_dp_irls_numeric_limits(large)
add("beyond-double-integer", large, beta = 1L, k = 0)
dense <- matrix(1L, 10L, 10L)
diag(dense) <- 11L
add("dense-p10-positive-definite", request(dense, rowSums(dense)),
    beta = rep(1L, 10L), k = 0L)
wide <- request(matrix(c(2L, 1L, 1L, 4L), 2L), c(3L, 5L),
                h_unit = q(1L, 3L), g_unit = q(1L, 3L), floor = q(1L, 3L))
huge <- gmp::as.bigz(2L)^255L + 19L
hvalues <- list(huge, gmp::as.bigz(1L), huge + 2L)
gvalues <- list(huge + 1L, huge + 3L)
wide$h_translation <- lapply(hvalues, function(z) as.character(z + 1L))
wide$h <- lapply(hvalues, function(z) as.character(2L * z + 1L))
wide$g_translation <- lapply(gvalues, function(z) as.character(z + 1L))
wide$g <- lapply(gvalues, function(z) as.character(2L * z + 1L))
wide$limits <- .dsvert_dp_irls_numeric_limits(wide)
add("rational-256-bit-bareiss", wide, beta = c(1L, 1L), k = 0L)
set.seed(20261015L)
for (i in seq_len(240L)) {
  p <- sample.int(5L, 1L)
  h <- matrix(sample(-7L:7L, p * p, replace = TRUE), p)
  h[lower.tri(h)] <- t(h)[lower.tri(h)]
  if (i %% 4L == 0L) h <- crossprod(h)
  if (i %% 17L == 0L) h[1L, ] <- h[, 1L] <- 0L
  g <- sample(-15L:15L, p, replace = TRUE)
  lattice <- q(1L, sample(c(1L, 2L, 4L, 8L), 1L))
  x <- request(h, g, beta = sample(-5L:5L, p, replace = TRUE), lattice = lattice,
    lower = rep(-8L, p), upper = rep(8L, p),
    floor = q(sample.int(4L, 1L), sample(c(1L, 2L, 4L), 1L)),
    quantum = q(1L, sample(c(1L, 2L, 3L, 4L), 1L)),
    ridge = q(sample(0L:2L, 1L), sample(c(1L, 2L), 1L)),
    damping = q(1L, sample(c(1L, 2L, 3L, 4L), 1L)),
    trust = q(sample(0L:7L, 1L), sample(c(1L, 2L, 3L), 1L)),
    g_unit = q(sample.int(3L, 1L), sample(c(1L, 2L, 3L), 1L)),
    h_unit = q(sample.int(3L, 1L), sample(c(1L, 2L, 3L), 1L)))
  add(sprintf("random-%03d", i), x)
}
writeLines(vapply(fixtures, .dsvert_dp_irls_arithmetic_json, character(1L)),
           "inst/testdata/dp_irls_update_v1.jsonl", useBytes = TRUE)
cat(sprintf("Wrote %d known-answer + 240 random = %d fixtures\n", length(fixtures) - 240L, length(fixtures)))
