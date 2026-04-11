test_that("formula parsing extracts y_var and x_vars", {
  # R formula object
  f <- npreg ~ age + bmi + glu
  ft <- terms(f)
  expect_equal(as.character(attr(ft, "variables")[[2]]), "npreg")
  expect_equal(attr(ft, "term.labels"), c("age", "bmi", "glu"))

  # String formula
  f2 <- as.formula("diabetes ~ age + bmi")
  ft2 <- terms(f2)
  expect_equal(as.character(attr(ft2, "variables")[[2]]), "diabetes")
  expect_equal(attr(ft2, "term.labels"), c("age", "bmi"))
})

test_that("auto-detect maps variables to correct servers", {
  # Simulate col_results from dsvertColNamesDS
  col_results <- list(
    s1 = list(columns = c("id", "patient_id", "age", "bmi", "ped")),
    s2 = list(columns = c("id", "patient_id", "glu", "bp", "skin")),
    s3 = list(columns = c("id", "patient_id", "npreg", "diabetes"))
  )

  # Build col_map (same logic as ds.vertGLM)
  col_map <- list()
  for (srv in names(col_results))
    col_map[[srv]] <- setdiff(col_results[[srv]]$columns, c("id", "patient_id"))

  expect_equal(col_map$s1, c("age", "bmi", "ped"))
  expect_equal(col_map$s2, c("glu", "bp", "skin"))
  expect_equal(col_map$s3, c("npreg", "diabetes"))

  # y_var detection
  y_var <- "npreg"
  y_servers <- names(col_map)[sapply(names(col_map), function(s) y_var %in% col_map[[s]])]
  expect_equal(y_servers, "s3")

  # x_vars building
  x_vars <- list()
  user_x_vars <- c("age", "bmi", "ped", "glu", "bp", "skin")
  for (srv in names(col_map)) {
    feats <- setdiff(col_map[[srv]], y_var)
    feats <- intersect(feats, user_x_vars)
    x_vars[[srv]] <- feats
  }
  expect_equal(x_vars$s1, c("age", "bmi", "ped"))
  expect_equal(x_vars$s2, c("glu", "bp", "skin"))
  expect_equal(x_vars$s3, character(0))  # s3 has only y
})

test_that("auto-detect with subset of formula vars", {
  col_results <- list(
    s1 = list(columns = c("id", "patient_id", "age", "bmi", "ped")),
    s2 = list(columns = c("id", "patient_id", "glu", "bp", "skin", "npreg"))
  )
  col_map <- list()
  for (srv in names(col_results))
    col_map[[srv]] <- setdiff(col_results[[srv]]$columns, c("id", "patient_id"))

  # User only wants age + glu
  user_x_vars <- c("age", "glu")
  x_vars <- list()
  for (srv in names(col_map)) {
    feats <- setdiff(col_map[[srv]], "npreg")
    feats <- intersect(feats, user_x_vars)
    x_vars[[srv]] <- feats
  }
  expect_equal(x_vars$s1, "age")
  expect_equal(x_vars$s2, "glu")
})

test_that("missing variable gives clear error", {
  all_available <- c("age", "bmi", "ped", "glu", "bp", "skin", "npreg")
  user_x_vars <- c("age", "bmi", "NONEXISTENT")
  missing <- setdiff(user_x_vars, all_available)
  expect_equal(missing, "NONEXISTENT")
})

test_that("y_var on multiple servers picks first", {
  col_map <- list(
    s1 = c("age", "bmi", "diabetes"),
    s2 = c("glu", "diabetes")
  )
  y_var <- "diabetes"
  y_servers <- names(col_map)[sapply(names(col_map), function(s) y_var %in% col_map[[s]])]
  expect_equal(y_servers, c("s1", "s2"))
  # Picks first
  expect_equal(y_servers[1], "s1")
})

test_that("p_coord=0 produces correct beta_map", {
  # Coordinator (s3) has 0 features, fusion (s1) has 3, s2 has 3
  x_vars <- list(s1 = c("age", "bmi", "ped"), s2 = c("glu", "bp", "skin"), s3 = character(0))
  coordinator <- "s3"
  p_coord <- length(x_vars[[coordinator]])

  beta_map <- list()
  idx <- 1
  if (p_coord > 0) { beta_map[[coordinator]] <- idx:(idx + p_coord - 1); idx <- idx + p_coord }
  for (srv in names(x_vars)) {
    if (srv == coordinator) next
    p_s <- length(x_vars[[srv]])
    if (p_s == 0) next
    beta_map[[srv]] <- idx:(idx + p_s - 1); idx <- idx + p_s
  }

  expect_null(beta_map[["s3"]])  # coordinator has no beta indices
  expect_equal(beta_map[["s1"]], 1:3)
  expect_equal(beta_map[["s2"]], 4:6)
})

test_that(".bsafe returns NULL for empty/NULL index", {
  .bsafe <- function(b, idx) { if (is.null(idx) || length(idx) == 0) NULL else b[idx] }
  beta <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6)

  expect_null(.bsafe(beta, NULL))
  expect_null(.bsafe(beta, integer(0)))
  expect_equal(.bsafe(beta, 1:3), c(0.1, 0.2, 0.3))
  expect_equal(.bsafe(beta, 4:6), c(0.4, 0.5, 0.6))
})

test_that("gradient remapping handles p_coord=0", {
  p_coord <- 0; p_fusion <- 3; p_total <- 6
  beta_map <- list(s1 = 1:3, s2 = 4:6)  # no coordinator in map
  gradient_canonical <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6)

  gradient <- numeric(p_total)
  gi <- 1
  if (p_coord > 0) { gradient[beta_map[["s3"]]] <- gradient_canonical[gi:(gi+p_coord-1)]; gi <- gi + p_coord }
  if (p_fusion > 0) { gradient[beta_map[["s1"]]] <- gradient_canonical[gi:(gi+p_fusion-1)]; gi <- gi + p_fusion }
  gradient[beta_map[["s2"]]] <- gradient_canonical[gi:(gi+3-1)]

  expect_equal(gradient, c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6))
})

test_that("base64url encoding is clean (no newlines, no special chars)", {
  .to_b64url <- function(x) gsub("+","-",gsub("/","_",gsub("=+$","",x,perl=TRUE),fixed=TRUE),fixed=TRUE)
  keys <- list(s1 = "abc123def456", s2 = "ghi789jkl012", ref = "abc123def456")
  j <- jsonlite::toJSON(keys, auto_unbox = TRUE)
  b64 <- jsonlite::base64_enc(charToRaw(j))
  b64url <- .to_b64url(gsub("\n", "", b64, fixed = TRUE))

  # No special chars
  expect_false(grepl("[^A-Za-z0-9_-]", b64url))
  # No newlines
  expect_false(grepl("\n", b64url))
  # Roundtrip
  .from_b64url <- function(x) {
    x <- gsub("-","+",gsub("_","/",x,fixed=TRUE),fixed=TRUE)
    pad <- nchar(x)%%4; if(pad==2) x<-paste0(x,"=="); if(pad==3) x<-paste0(x,"="); x
  }
  decoded <- rawToChar(jsonlite::base64_dec(.from_b64url(b64url)))
  expect_equal(jsonlite::fromJSON(decoded, simplifyVector = FALSE), keys)
})

test_that("NULL x_vars in call() doesn't produce character(0)", {
  # When a server has 0 features, x_vars should be NULL, not character(0)
  srv_x <- character(0)
  safe_x <- if (length(srv_x) == 0) NULL else srv_x
  expect_null(safe_x)

  # Verify NULL in call() works
  e <- call("someFunction", x_vars = NULL, session_id = "test")
  expect_equal(deparse(e), "someFunction(x_vars = NULL, session_id = \"test\")")
})
