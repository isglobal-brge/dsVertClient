# Offline proof-of-concept report for one completed DP-IRLS run.
# Input: .tmp/v3/poc-run.rds from `DSVERT_B3_POC=1 tools/validate_dp_irls_b3.R`.
# Uses only released stages, the signed contract and public noise laws, plus
# the synthetic source table for a centralised glm() comparison. The inference
# shown here is the experimental replay candidate, not qualified inference.
args <- commandArgs(TRUE)
root <- normalizePath(if (length(args)) args[[1L]] else "..")
out <- file.path(root, ".tmp", "v3")
pkgload::load_all(file.path(root, "v3-dsVertClient"), quiet = TRUE)
cf <- function(name) get(name, asNamespace("dsVertClient"))
run <- readRDS(file.path(out, "poc-run.rds"))
contract <- run$signed$contract
num <- function(q) {
  if (is.list(q) && !is.null(q$n)) return(as.numeric(gmp::as.bigq(q$n, q$d)))
  as.numeric(gmp::as.bigq(as.character(q)))
}

# 1. Independent verification of every signature, identity and public update.
verified <- cf(".dsvert_dp_irls_verify_chain")(run$signed, run$pins, run$stages,
                                               expected_root = run$root)
stopifnot(isTRUE(verified$complete))
beta <- vapply(verified$beta, num, numeric(1L))
terms <- unlist(contract$estimand$coefficient_order)
p <- length(beta)

# 2. Released terminal-audit statistics at beta_T.
audit <- cf(".dsvert_dp_irls_client_decode")(run$stages[[length(run$stages)]], contract)
g_T <- vapply(audit$released$g, num, numeric(1L))
upper <- vapply(audit$released$H, num, numeric(1L))
sym <- function(u) {  # signed layout is the row-major upper triangle H:i:j, i <= j
  H <- matrix(0, p, p); k <- 0L
  for (i in seq_len(p)) for (j in i:p) { k <- k + 1L; H[i, j] <- H[j, i] <- u[[k]] }
  H
}
stopifnot(identical(unlist(contract$numeric$layout)[-seq_len(p)],
  unlist(lapply(seq_len(p), function(i) sprintf("H:%d:%d", i, i:p)))))
H_T <- sym(upper)

# 3. Public path and noise-law constants from the signed contract.
path <- contract$path; acc <- contract$accounting
floor_ <- num(path$floor); quantum <- num(path$floor_quantum); ridge <- num(path$ridge)
damping <- vapply(path$damping, num, numeric(1L)); trust <- vapply(path$trust, num, numeric(1L))
lower <- vapply(path$box$lower, num, numeric(1L)); upperb <- vapply(path$box$upper, num, numeric(1L))
lat <- c(g = num(contract$numeric$lattice$g), H = num(contract$numeric$lattice$H))
S1 <- c(g = as.numeric(acc$sensitivity$g), H = as.numeric(acc$sensitivity$H))
lane_eps <- function(i, lane) num(acc$stages[[i]][[lane]])
dlap2 <- function(m, eps, lane) {  # two independent full discrete-Laplace draws
  q <- exp(-eps / S1[[lane]])
  k <- (stats::rgeom(m, 1 - q) - stats::rgeom(m, 1 - q)) +
       (stats::rgeom(m, 1 - q) - stats::rgeom(m, 1 - q))
  k * lat[[lane]]
}
update <- function(b, g, H, t) {
  shift <- max(0, ceiling((floor_ - min(eigen(H, TRUE, only.values = TRUE)$values)) / quantum))
  v <- solve(H + (shift * quantum + ridge) * diag(p), g - ridge * b) * damping[[t]]
  if (max(abs(v)) > trust[[t]]) v <- v * trust[[t]] / max(abs(v))
  pmin(upperb, pmax(lower, b + v))
}

# 4. Experimental replay candidate on the released local quadratic (Claude review, #26).
set.seed(20261015); R <- 2000L
T_ <- as.integer(path$T)
start <- vapply(path$start$beta, num, numeric(1L))
ends <- t(vapply(seq_len(R), function(r) {
  b <- start
  for (t in seq_len(T_)) {
    E <- sym(dlap2(p * (p + 1) / 2, lane_eps(t, "H"), "H"))
    g <- drop(H_T %*% (beta - b)) + dlap2(p, lane_eps(t, "g"), "g")
    b <- update(b, g, H_T + E, t)
  }
  b
}, numeric(p)))
V_dp <- stats::cov(ends); replay_bias <- colMeans(ends) - beta
pd <- min(eigen(H_T, TRUE, only.values = TRUE)$values) > 0
V_s <- if (pd) solve(H_T) else matrix(NA_real_, p, p)
se <- sqrt(diag(V_s) + diag(V_dp))

# 5. Public refusal reasons (PoC subset of design §8).
flags <- lapply(run$stages, function(s) s$body$update_result)
reasons <- character()
for (i in seq_along(flags)) if (!is.null(f <- flags[[i]])) {
  if (!identical(as.character(f$k), "0")) reasons <- c(reasons, sprintf("stage %d floor shift", i - 1L))
  if (isTRUE(f$trust_scaled)) reasons <- c(reasons, sprintf("stage %d trust truncation", i - 1L))
  if (any(unlist(f$projected))) reasons <- c(reasons, sprintf("stage %d box projection", i - 1L))
}
if (isTRUE(audit$endpoint_hit)) reasons <- c(reasons, "finaliser endpoint")
if (!pd) reasons <- c(reasons, "terminal curvature not positive definite")
if (pd && max(abs(replay_bias) / sqrt(diag(V_s))) > 0.1) reasons <- c(reasons, "replay bias > 0.1 sampling SD")

# 6. Centralised comparison (synthetic data; never available to a real analyst).
fit <- stats::glm(y ~ x + z, family = stats::binomial(), data = run$data)
table <- data.frame(term = terms, dp_estimate = round(beta, 3),
  dp_total_se = round(se, 3), sampling_sd = round(sqrt(diag(V_s)), 3),
  dp_sd = round(sqrt(diag(V_dp)), 3),
  ci_low = round(beta - 1.96 * se, 3), ci_high = round(beta + 1.96 * se, 3),
  glm_estimate = round(unname(stats::coef(fit)), 3),
  glm_se = round(unname(sqrt(diag(stats::vcov(fit)))), 3))
stage_seconds <- vapply(run$evidence$stages, function(s) s$wall_seconds, numeric(1L))
report <- list(version = "dsvert-v3-irls-poc-report-v1", n = nrow(run$data),
  epsilon_total = num(acc$epsilon_path), stages = length(run$stages),
  chain_verified_offline = TRUE, accounting = verified$accounting,
  inference_status = if (length(reasons)) "refused (experimental values shown for diagnosis only)" else "experimental candidate",
  refusal_reasons = as.list(reasons), table = table,
  stage_seconds = stage_seconds, setup_seconds = run$evidence$setup$total_seconds)
jsonlite::write_json(report, file.path(out, "poc-report.json"), auto_unbox = TRUE, pretty = TRUE, digits = 8)
cat("Signed continuous-DP binomial analysis (proof of concept)\n")
cat(sprintf("n = %d, epsilon = %g over %d signed releases, chain verified offline\n",
            report$n, report$epsilon_total, report$stages))
print(table, row.names = FALSE)
cat("Inference:", report$inference_status, "\n")
if (length(reasons)) cat("Reasons:", paste(reasons, collapse = "; "), "\n")
cat("Stage seconds:", paste(round(stage_seconds, 1), collapse = ", "), "\n")
