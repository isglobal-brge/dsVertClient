# Leak-free fixed-iteration control.
#
# The iterative fit methods (IRLS / Newton / PQL) emit one batch of Beaver /
# Ring127 round-trips per iteration. Stopping the loop on a data-derived
# convergence test (`if (max_diff < tol) break`) shortens that round sequence by
# an amount that depends on the data's conditioning/separability, so a peer
# server counting rounds learns the iteration-to-convergence — a weak,
# NON-reconstructive, but data-dependent side channel.
#
# To close it, dsVert runs a FIXED, PUBLIC, data-INDEPENDENT number of
# iterations by DEFAULT: the per-family count below is a function of the model
# family only (public DataSHIELD metadata), never of any shared value, so every
# fit of a given family emits the same number of rounds regardless of the data.
# Convergence is still detected (the `converged` flag is recorded), but the loop
# does not exit early — post-convergence iterations are numerically harmless
# no-ops (the update scales with the near-zero gradient; L-BFGS history freezes
# via its curvature gate). The final beta is bit-identical to the early-stopped
# fit within the ring noise floor.
#
# Opt into the faster, data-dependent early stop with
# `options(dsvert.early_stop = TRUE)` — this restores the data-dependent
# early stop and re-exposes the (weak) iteration-count side channel. Custodians may tune a
# per-family fixed count with `options(dsvert.fixed_iters_<family> = N)`.

#' @keywords internal
.dsvert_early_stop <- function() {
  isTRUE(getOption("dsvert.early_stop", FALSE))
}

#' @keywords internal
# Public tol-adaptive growth: a stricter convergence tolerance legitimately
# needs more iterations, and `tol` is a PUBLIC parameter, so growing the fixed
# count with tol leaks nothing while keeping default-tol fits fast. ~5 extra
# iterations per order of magnitude below the 1e-4 reference (empirically the
# L-BFGS/Newton cost of tightening), capped.
.dsvert_tol_growth <- function(tol) {
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0) return(0L)
  g <- ceiling(5 * log10(1e-4 / tol))
  as.integer(min(max(0L, g), 40L))
}

#' @keywords internal
.dsvert_fixed_iters <- function(family, base_n, tol = 1e-4) {
  fam <- tolower(as.character(family)[1L])
  opt <- getOption(paste0("dsvert.fixed_iters_", fam))
  if (!is.null(opt)) {                    # custodian override = exact total, no growth
    n <- suppressWarnings(as.integer(opt)[1L])
    if (length(n) == 1L && !is.na(n) && n >= 1L) return(n)
  }
  as.integer(base_n) + .dsvert_tol_growth(tol)
}

#' @keywords internal
# Resolve the loop bound: leak-free fixed-N (public, family + public-tol driven)
# by default; the data-dependent early stop (up to max_iter) only when opted in.
# `max_iter` always acts as an upper safety clamp and can only LOWER the count.
.dsvert_loop_n <- function(family, base_n, max_iter, tol = 1e-4) {
  if (.dsvert_early_stop()) {
    return(max(1L, as.integer(max_iter)))
  }
  min(as.integer(max_iter), .dsvert_fixed_iters(family, base_n, tol))
}
