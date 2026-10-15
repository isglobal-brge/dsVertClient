# Public, balanced synthetic fixtures for the process-isolated federation smoke.
# K=3 includes a within-owner pair and pairs spanning every data owner.
# K=5 uses the same requested variables, with two additional signed members.
vertcor142_cor_fixture <- function(k, n = 64L) {
  stopifnot(length(k) == 1L, k %in% c(2L, 3L, 5L),
            length(n) == 1L, is.finite(n), n >= 8L, n <= 1024L,
            n == as.integer(n), n %% 8L == 0L)
  owners <- paste0("site_", letters[seq_len(k)])
  refs <- if (k == 2L) c("site_a$a", "site_b$x") else
    c("site_a$a", "site_a$b", "site_b$x", "site_c$x")
  masks <- if (k == 2L) c(1L, 2L) else c(1L, 2L, 4L, 3L)
  raw <- stats::setNames(lapply(owners, function(owner)
    data.frame(patient_id = sprintf("u%04d", seq_len(n)))), owners)
  for (i in seq_along(refs)) {
    parts <- strsplit(refs[[i]], "$", fixed = TRUE)[[1L]]
    raw[[parts[[1L]]]][[parts[[2L]]]] <- vapply(seq_len(n) - 1L,
      function(row) sum(as.integer(intToBits(bitwAnd(row, masks[[i]])))) %% 2,
      numeric(1L))
  }
  # Bootstrap requires every member to register a bounded local column.
  # The two extra members contribute no variable to the requested Gaussian block.
  for (owner in owners) if (ncol(raw[[owner]]) == 1L)
    raw[[owner]]$membership_marker <- rep(c(0, 1), length.out = n)
  physical <- sub("^[^$]+[$]", "", refs)
  duplicate <- duplicated(physical) | duplicated(physical, fromLast = TRUE)
  list(owners = owners, refs = refs, raw = raw,
       registered = ifelse(duplicate, refs, physical),
       n = as.integer(n), designated = owners[1:2])
}
