# Cryptographically secure identifiers for MPC protocol state.  There is no
# pseudo-random fallback: predictable identifiers can cause session collisions
# and state confusion across concurrent analyses.
.dsvert_random_bytes <- function(n, .urandom_path = "/dev/urandom",
                                 .openssl_rand_bytes = if (
                                   requireNamespace("openssl", quietly = TRUE)
                                 ) openssl::rand_bytes else NULL) {
  if (length(n) != 1L || !is.numeric(n) || is.na(n) || !is.finite(n) ||
      n < 1 || n != floor(n) || n > .Machine$integer.max) {
    stop("n must be a positive integer", call. = FALSE)
  }
  n <- as.integer(n)

  if (is.character(.urandom_path) && length(.urandom_path) == 1L &&
      !is.na(.urandom_path) && file.exists(.urandom_path)) {
    bytes <- tryCatch({
      con <- file(.urandom_path, open = "rb", raw = TRUE)
      on.exit(close(con), add = TRUE)
      readBin(con, what = "raw", n = n)
    }, error = function(e) NULL)
    if (is.raw(bytes) && length(bytes) == n) return(bytes)
  }

  if (is.function(.openssl_rand_bytes)) {
    bytes <- tryCatch(.openssl_rand_bytes(n), error = function(e) NULL)
    if (is.raw(bytes) && length(bytes) == n) return(bytes)
  }

  stop("cryptographically secure random source is unavailable", call. = FALSE)
}

.dsvert_uuid4 <- function(.random_bytes = .dsvert_random_bytes) {
  if (!is.function(.random_bytes)) {
    stop(".random_bytes must be a function", call. = FALSE)
  }
  bytes <- .random_bytes(16L)
  if (!is.raw(bytes) || length(bytes) != 16L) {
    stop("UUIDv4 generation requires exactly 16 random bytes", call. = FALSE)
  }

  # RFC 9562 UUIDv4: version 4 in byte 7 and variant 10xx in byte 9.
  bytes[[7L]] <- as.raw(bitwOr(bitwAnd(as.integer(bytes[[7L]]), 0x0fL), 0x40L))
  bytes[[9L]] <- as.raw(bitwOr(bitwAnd(as.integer(bytes[[9L]]), 0x3fL), 0x80L))
  hex <- sprintf("%02x", as.integer(bytes))

  paste0(
    paste0(hex[1:4], collapse = ""), "-",
    paste0(hex[5:6], collapse = ""), "-",
    paste0(hex[7:8], collapse = ""), "-",
    paste0(hex[9:10], collapse = ""), "-",
    paste0(hex[11:16], collapse = "")
  )
}

# Historical internal name retained by LMM and compatibility routes.  Keep the
# implementation centralized so every package-generated session obeys the
# same CSPRNG-backed canonical UUID contract.
.mpc_session_id <- function() .dsvert_uuid4()
