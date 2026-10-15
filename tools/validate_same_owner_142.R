# Run from the dsVertClient source root with the supported R library active:
# timeout 900 Rscript tools/validate_same_owner_142.R [output_dir] [baseline_commit]
# Optional third argument is the pre-edit Gaussian consumer RDS captured in phase A/B.
# Inputs are synthetic. Bundle/compilation/accuracy fixture seams are documented
# in helper-gaussian-current-certificate.R; signed RELEASE/REPLAY, cross attestations,
# certificate builders/readers and consumer arithmetic are real production code.
library(testthat)
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(getwd(), mustWork = TRUE)
stopifnot(file.exists(file.path(root, 'DESCRIPTION')),
          file.exists(file.path(root, 'tests/testthat/helper-gaussian-current-certificate.R')))
output <- if (length(args)) args[[1L]] else file.path(tempdir(), 'dsvert-same-owner-142')
dir.create(output, recursive = TRUE, showWarnings = FALSE)
base_ref <- if (length(args) >= 2L) args[[2L]] else '04a6903'
stopifnot(grepl('^[0-9a-f]{7,40}$', base_ref))
pkgload::load_all(root, quiet = TRUE)
cat('Reference client commit:', base_ref, '\n')
current <- asNamespace('dsVertClient')
helper <- new.env(parent=current)
sys.source('tests/testthat/helper-gaussian-current-certificate.R',helper)
legacy <- new.env(parent=current)
for (node in parse('tests/testthat/test-dp-gaussian.R')) {
 if (is.call(node) && identical(node[[1L]],as.name('<-'))) eval(node,legacy)
}
reference <- new.env(parent=current)
reference_files <- c('ds.vertDPGaussian.R','dp_gaussian_certificate.R','dp_synopsis_public_vector.R',
 'ds.vertLASSO.R','ds.vertLASSOCapsule.R','ds.vertLASSOProximal.R','ds.vertLASSOCV.R','ds.vertInference.R')
for(file in reference_files) {
 lines <- system2('git',c('show',paste0(base_ref,':R/',file)),stdout=TRUE)
 stopifnot(is.null(attr(lines,'status')))
 eval(parse(text=lines),reference)
}
# Authenticated same-owner inputs, keys and signed wire frames are fixed once.
f <- helper$.gaussian_current_fixture(2L,artifact=legacy$.dp_gaussian_artifact(),
 coordinates=c(80,80,40,30,40,24,30))
reference$.dsvert_dp_synopsis_client_bundle <- function(...) f$trusted
reference$.dsvert_dp_synopsis_client_compile <- function(...) f$compiled
reference$.dsvert_dp_gaussian_synopsis_trusted <- function(...) list(trusted=f$trusted,bundle=f$bundle)
reference$.dsvert_dp_vector_accuracy_radius <- function(...) list(radius=.01,confidence=.95,
 method='fixture_public_calibration',implementation_tv_upper_bound=0,sampler_tv_upper_bound=0,
 additional_privacy_cost=c(epsilon=0,delta=0))
reference$.gaussian_current_source <- helper$.gaussian_current_source
environment(reference$.gaussian_current_source) <- reference
helper$.with_gaussian_current_fixture(f, {
 old_source <- reference$.gaussian_current_source(f)
 new_source <- helper$.gaussian_current_source(f)
 reference$.dsvert_dp_gaussian_synopsis_release <- function(...) old_source
 old_fit <- reference$.dsvert_dp_gaussian_impl('cohort','gaussian_primary',0,NULL,f$conns,
  function(...) stop('Unexpected server call'))
 new_fit <- with_mocked_bindings(ds.vertDPGaussian('cohort','gaussian_primary',datasources=f$conns),
  .dsvert_dp_gaussian_synopsis_release=function(...) new_source,.package='dsVertClient')
 consume <- function(ns,fit) list(
  path=ns$ds.vertLASSO(fit,.02,alpha_grid=c(1,.5),trusted_pinset=f$pins),
  one_step=ns$ds.vertLASSO1Step(fit,c(.02,.01),trusted_pinset=f$pins),
  proximal=ns$ds.vertLASSOProximal(fit,.02,trusted_pinset=f$pins),
  cv=ns$ds.vertLASSOCV(fit,c(.02,.01),trusted_pinset=f$pins),
  confint=ns$ds.vertConfint(fit,type='mechanism'),
  wald=ns$ds.vertWald(fit,'x',type='mechanism'),
  contrast=ns$ds.vertContrast(fit,'x',type='mechanism'),
  lr=ns$ds.vertLR('(Intercept)',fit,type='mechanism'))
 old <- list(artifact=old_source$artifact,release=old_source$context$release,
  certificate=old_source$certificate,fit=old_fit,consumers=consume(reference,old_fit))
 new <- list(artifact=new_source$artifact,release=new_source$context$release,
  certificate=new_source$certificate,fit=new_fit,consumers=consume(current,new_fit))
 for(name in names(old)) {
  x <- serialize(old[[name]],NULL,version=3); y <- serialize(new[[name]],NULL,version=3)
  stopifnot(identical(x,y))
  cat(name,'byte-identical',length(x),'sha256',digest::digest(x,'sha256',serialize=FALSE),'\n')
 }
 saveRDS(list(fixture=f,reference=old,current=new),
  file.path(output, 'same-owner-current-byte-comparison.rds'),version=3)
})
if (length(args) >= 3L) {
baseline <- readRDS(args[[3L]])
fit <- baseline$fit; pins <- baseline$pins
again <- list(path=ds.vertLASSO(fit,.1,trusted_pinset=pins),
 onestep=ds.vertLASSO1Step(fit,c(.1,.05),trusted_pinset=pins),
 proximal=ds.vertLASSOProximal(fit,.05,trusted_pinset=pins),
 cv=ds.vertLASSOCV(fit,c(.1,.05),trusted_pinset=pins))
stopifnot(identical(serialize(again,NULL,version=3),serialize(baseline$result,NULL,version=3)))
cat('Pre-edit legacy same-owner LASSO golden',digest::digest(serialize(again,NULL,version=3),'sha256',serialize=FALSE),'identical\n')
}
cat('Saved fixed fixture and both result sets:', file.path(output, 'same-owner-current-byte-comparison.rds'), '\n')
