# PUBLIC SYNTHETIC DATA ONLY. The tagged Go executable and these DSLite methods
# are test fixtures, not a production producer, DP sampler or deployment gate.
test_that("two DSLite custodians exercise tagged Cox oracle and DP grid selection", {
  skip_if_not_installed("DSLite")
  skip_if_not_installed("survival")
  binary <- Sys.getenv("DSVERT_COX_TEST_BINARY")
  if (!nzchar(binary)) skip("requires explicitly compiled dsvert_cox_plaintext_test binary")
  expect_true(file.exists(binary))
  set.seed(91027)
  n <- 192L
  # Exact dyadic normalized inputs/coefs make the independent eta oracle exact.
  pooled <- data.frame(x = sample(0:16, n, TRUE) / 16,
                       z = sample(0:16, n, TRUE) / 16)
  eta <- pooled$x - pooled$z / 2
  failure <- rexp(n, exp(eta)); censor <- rexp(n, 0.25)
  pooled$time <- pmin(20, ceiling(pmin(failure, censor) * 4) / 4)
  pooled$event <- as.integer(failure <= censor & failure <= 20)
  grid <- lapply(seq(-1, 2, by=0.5), function(x) c(x, -0.5))
  f <- .cox_cross_client_fixture(n, grid)
  spec <- f$contract$spec
  beta <- do.call(rbind, lapply(spec$beta_grid, unlist))
  caps <- unlist(spec$sensitivity$maximum_coordinates)
  # Each peer contributes independent Gamma-Poisson r=1/2 differences.
  # Their sum is a discrete Laplace reference draw at the signed L1 scale.
  # This checks utility plumbing; production uses the existing exact GC sampler.
  refs <- character()
  for (index in 1:2) {
    local_data <- if (index == 1L) pooled[c("x", "time", "event")] else pooled["z"]
    server <- DSLite::newDSLiteServer(tables=list(t=local_data))
    # DSLite intentionally rebinds aggregate environments; freeze only PUBLIC
    # fixture constants in the method body and read the assigned local table.
    endpoint <- eval(substitute(function(epsilon) {
      d <- get("coxdata", envir=parent.frame())
      public_beta <- matrix(BETA, ncol=2)
      partial <- if (OWNER == 1L) outer(d$x, public_beta[,1]) else outer(d$z, public_beta[,2])
      probability <- -expm1(-epsilon / DELTA)
      noise <- stats::rnbinom(nrow(public_beta), size=0.5, prob=probability) -
        stats::rnbinom(nrow(public_beta), size=0.5, prob=probability)
      list(partial_q16=partial*65536, noise=noise,
        time=if (OWNER==1L) d$time else NULL,
        event=if (OWNER==1L) d$event else NULL)
    }, list(BETA=beta, OWNER=index, DELTA=spec$sensitivity$raw_l1_sensitivity)))
    server$aggregateMethod("coxPublicSyntheticFixture", endpoint)
    name <- paste0("cox_fixture_", Sys.getpid(), "_", index)
    assign(name, server, envir=.GlobalEnv); refs[[c("site_a", "site_b")[[index]]]] <- name
  }
  on.exit(rm(list=unname(refs), envir=.GlobalEnv), add=TRUE)
  login <- DSI::newDSLoginBuilder()
  for (site in names(refs)) login$append(server=site, url=refs[[site]], table="t", driver="DSLiteDriver")
  conns <- DSI::datashield.login(login$build(), assign=FALSE)
  DSI::datashield.assign.table(conns, symbol="coxdata", table="t")
  on.exit(DSI::datashield.logout(conns), add=TRUE)
  exact_loss <- function(b) {
    e <- drop(as.matrix(pooled[c("x","z")]) %*% b)
    sum(vapply(which(pooled$event==1), function(i) {
      log(sum(exp(e[pooled$time >= pooled$time[[i]]]))) - e[[i]]
    }, numeric(1L)))
  }
  exact <- apply(beta, 1, exact_loss)
  fit <- survival::coxph(survival::Surv(time,event) ~ x+z, pooled, ties="breslow")
  expect_lte(exact_loss(stats::coef(fit)), min(exact)+1e-7)
  input <- tempfile(fileext=".json"); output <- tempfile(); errors <- tempfile()
  on.exit(unlink(c(input,output,errors)), add=TRUE)
  report <- list()
  for (epsilon in c(1,4,8)) {
    peers <- DSI::datashield.aggregate(conns, as.call(list(
      as.name("coxPublicSyntheticFixture"), epsilon=epsilon)), errors.print=TRUE)
    eta_q16 <- peers$site_a$partial_q16 + peers$site_b$partial_q16
    request <- list(eta_q16=unname(eta_q16), time=peers$site_a$time,
      event=peers$site_a$event, valid=rep(TRUE,n), capacity=n,
      numeric_grid_bits=spec$numeric_grid_bits, coefficient_grid=unname(beta))
    writeLines(jsonlite::toJSON(request, auto_unbox=TRUE, digits=NA), input)
    status <- system2(binary, "-test.run=^TestCoxGridCrossPlaintextCommand$",
      stdin=input, stdout=output, stderr=errors, env="DSVERT_COX_PLAINTEXT_TEST=1")
    expect_equal(status,0)
    answer <- jsonlite::fromJSON(output)
    expect_true(answer$test_only)
    reference <- vapply(seq_len(nrow(beta)), function(j) .cox_integer_loss(
      eta_q16[,j], pooled$time, pooled$event, rep(TRUE,n), caps[[j]], 8), numeric(1L))
    expect_equal(answer$loss_coordinates, reference, tolerance=0)
    expect_true(all(abs(reference/256-exact) <= n/64+1/512))
    noisy <- pmin(caps, pmax(0, reference+peers$site_a$noise+peers$site_b$noise))
    release <- function(...) c(f[c("contract","policy","schema_manifest")],list(coordinates=noisy))
    result <- .dsvert_dp_cox_grid_cross_impl(
      Surv(site_a$time,site_a$event)~site_a$x+site_b$z,
      "aligned","cox_grid",conns,.release=release)
    selected <- result$selected_candidate
    expect_equal(selected, which.min(noisy))
    expect_null(result$std_errors)
    expect_false(result$production_ready)
    report[[as.character(epsilon)]] <- list(epsilon=epsilon,
      selected=selected, exact_best=which.min(exact), selected_beta=unname(beta[selected,]),
      exact_best_beta=unname(beta[which.min(exact),]), coxph_beta=unname(stats::coef(fit)),
      exact_excess_loss=exact[[selected]]-min(exact), coxph_loss=exact_loss(stats::coef(fit)),
      grid_best_loss=min(exact), loss_error_max=max(abs(reference/256-exact)))
  }
  path <- Sys.getenv("DSVERT_COX_TEST_REPORT")
  if (nzchar(path)) writeLines(jsonlite::toJSON(list(public_synthetic=TRUE,
    production_protocol=FALSE,n=n,ties="breslow",results=report),auto_unbox=TRUE,pretty=TRUE,digits=16),path)
})
