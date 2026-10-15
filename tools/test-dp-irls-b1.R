# Run only B1 offline arithmetic and signed-chain tests from this worktree.
pkgload::load_all('.',quiet=TRUE)
files <- c('test-dp-irls-verify.R','test-dp-irls-chain.R')
counts <- lapply(files,function(file) {
 r <- testthat::test_file(file.path('tests/testthat',file),reporter='summary')
 d <- as.data.frame(r)
 data.frame(file=file,tests=nrow(d),assertions=sum(d$passed),
   failed=sum(d$failed),errors=sum(d$error),warnings=sum(d$warning),skipped=sum(d$skipped))
})
counts <- do.call(rbind,counts); print(counts,row.names=FALSE)
if(any(counts$failed>0|counts$errors>0|counts$warnings>0|counts$skipped>0)) quit(status=1L)
