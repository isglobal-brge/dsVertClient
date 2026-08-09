#!/usr/bin/env Rscript
# Kept as a fail-closed migration marker for callers of the former installed
# validation runner. It must never open a session or invoke quarantined routes.
stop(
  paste(
    "This historical validation runner is retired.",
    "Use validate_dp_statistical_methods.R for offline statistical checks",
    "and the current DP capsule examples for supported remote workflows."
  ),
  call. = FALSE
)
