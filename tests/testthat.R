# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.

library(testthat)
library(dsVertClient)

# Each installed-check process owns its state; never inherit a service identity.
.dsvert_test_state_dir <- tempfile("dsvert-check-state-")
stopifnot(dir.create(.dsvert_test_state_dir, mode = "0700"))
Sys.setenv(DSVERT_STATE_DIR = .dsvert_test_state_dir)
testthat::set_max_fails(Inf)

test_check("dsVertClient")
