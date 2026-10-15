# Differential harness only. The R verifier itself never calls this program,
# the Go binary, or any process. Run from the client repository root.
source("R/dp_irls_verify.R")
binary <- Sys.getenv("DS_VERT_MPC")
stopifnot(nzchar(binary), file.exists(binary))
fixtures <- lapply(readLines("inst/testdata/dp_irls_update_v1.jsonl", warn = FALSE),
                   jsonlite::fromJSON, simplifyVector = FALSE)
input <- tempfile(); output <- tempfile(); errors <- tempfile()
cleanup <- function() unlink(c(input, output, errors))
count <- 0L
tryCatch({
  for (fixture in fixtures) {
    wire <- .dsvert_dp_irls_arithmetic_json(fixture$request)
    writeChar(wire, input, eos = NULL, useBytes = TRUE)
    status <- system2(binary, "dp-irls-update", stdin = input, stdout = output, stderr = errors)
    if (status != 0L) stop(paste(fixture$name, paste(readLines(errors), collapse = "\n")))
    go <- rawToChar(readBin(output, "raw", n = file.info(output)$size))
    independent <- .dsvert_dp_irls_verify_update_json(wire)
    expected <- .dsvert_dp_irls_arithmetic_json(fixture$expected)
    if (!identical(go, independent) || !identical(go, expected)) {
      stop(paste("Byte disagreement:", fixture$name))
    }
    count <- count + 1L
  }
  cat(sprintf("Go/R byte agreement: %d/%d (29 known answers + 240 random); zero mismatches\n",
              count, length(fixtures)))
}, finally = cleanup())
