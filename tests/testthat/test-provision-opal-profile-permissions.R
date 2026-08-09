.provision_profile_permissions_validator <- function() {
  source_script <- testthat::test_path(
    "..", "..", "inst", "scripts", "provision_opal.R")
  script <- if (file.exists(source_script)) source_script else system.file(
    "scripts", "provision_opal.R", package = "dsVertClient")
  expressions <- parse(script, keep.source = FALSE)
  match <- vapply(expressions, function(expression) {
    is.call(expression) && identical(expression[[1L]], as.name("<-")) &&
      identical(expression[[2L]], as.name(".validated_profile_permissions"))
  }, logical(1L))
  stopifnot(sum(match) == 1L)
  evaluation_env <- new.env(parent = baseenv())
  eval(expressions[[which(match)]], envir = evaluation_env)
  evaluation_env$.validated_profile_permissions
}

.provision_install_packages <- function(bindings) {
  source_script <- testthat::test_path(
    "..", "..", "inst", "scripts", "provision_opal.R")
  script <- if (file.exists(source_script)) source_script else system.file(
    "scripts", "provision_opal.R", package = "dsVertClient")
  expressions <- parse(script, keep.source = FALSE)
  match <- vapply(expressions, function(expression) {
    is.call(expression) && identical(expression[[1L]], as.name("<-")) &&
      identical(expression[[2L]], as.name("install_packages"))
  }, logical(1L))
  stopifnot(sum(match) == 1L)
  evaluation_env <- list2env(bindings, parent = baseenv())
  evaluation_env$PROFILE <- "dsvert"
  evaluation_env$DSVERT_IMPORTS <- character()
  evaluation_env$.profile_rock_cluster <- function(...) "default"
  evaluation_env$install_server_imports <- function(...) invisible(TRUE)
  evaluation_env$TARBALLS <- list(
    dsVert = "dsVert.tar.gz", dsVertClient = "dsVertClient.tar.gz")
  eval(expressions[[which(match)]], envir = evaluation_env)
  evaluation_env$install_packages
}

.provision_install_server_imports <- function(bindings) {
  source_script <- testthat::test_path(
    "..", "..", "inst", "scripts", "provision_opal.R")
  script <- if (file.exists(source_script)) source_script else system.file(
    "scripts", "provision_opal.R", package = "dsVertClient")
  expressions <- parse(script, keep.source = FALSE)
  match <- vapply(expressions, function(expression) {
    is.call(expression) && identical(expression[[1L]], as.name("<-")) &&
      identical(expression[[2L]], as.name("install_server_imports"))
  }, logical(1L))
  stopifnot(sum(match) == 1L)
  evaluation_env <- list2env(bindings, parent = baseenv())
  evaluation_env$PROFILE <- "dsvert"
  evaluation_env$DSVERT_IMPORTS <- "filelock"
  eval(expressions[[which(match)]], envir = evaluation_env)
  evaluation_env$install_server_imports
}

.provision_profile_rock_cluster <- function(bindings = list()) {
  source_script <- testthat::test_path(
    "..", "..", "inst", "scripts", "provision_opal.R")
  script <- if (file.exists(source_script)) source_script else system.file(
    "scripts", "provision_opal.R", package = "dsVertClient")
  expressions <- parse(script, keep.source = FALSE)
  match <- vapply(expressions, function(expression) {
    is.call(expression) && identical(expression[[1L]], as.name("<-")) &&
      identical(expression[[2L]], as.name(".profile_rock_cluster"))
  }, logical(1L))
  stopifnot(sum(match) == 1L)
  evaluation_env <- list2env(bindings, parent = baseenv())
  eval(expressions[[which(match)]], envir = evaluation_env)
  evaluation_env$.profile_rock_cluster
}

test_that("Opal provisioner accepts only canonical permission inventories", {
  validate <- .provision_profile_permissions_validator()

  empty <- validate(data.frame())
  expect_identical(
    empty,
    data.frame(
      subject = character(), type = character(), permission = character(),
      stringsAsFactors = FALSE))

  valid <- data.frame(
    subject = "runner", type = "user", permission = "use",
    stringsAsFactors = FALSE)
  expect_identical(validate(valid), valid)

  expect_error(
    validate(data.frame(unexpected = "value", stringsAsFactors = FALSE)),
    "invalid profile-permission inventory", fixed = TRUE)
  expect_error(
    validate(data.frame(unexpected = character())),
    "invalid profile-permission inventory", fixed = TRUE)
})

test_that("Opal provisioner uses the real Rock package inventory", {
  state <- new.env(parent = emptyenv())
  state$installed <- c(dsVert = FALSE, dsVertClient = FALSE)
  state$installed_paths <- character()
  state$removed <- character()
  state$removed_methods <- character()
  state$unpublished <- character()
  state$events <- character()
  state$dsadmin_probe_calls <- 0L

  install <- .provision_install_packages(list(
    oadmin.installed_package = function(o, pkg, profile) {
      stopifnot(identical(profile, "default"))
      unname(state$installed[[pkg]])
    },
    dsadmin.installed_package = function(...) {
      state$dsadmin_probe_calls <- state$dsadmin_probe_calls + 1L
      TRUE # Equivalent to Opal's post-delete HTTP 200 with an empty body.
    },
    dsadmin.rm_package_methods = function(o, pkg, profile) {
      stopifnot(identical(profile, "dsvert"))
      state$removed_methods <- c(state$removed_methods, pkg)
      state$events <- c(state$events, paste("methods", pkg, sep = ":"))
    },
    dsadmin.remove_package = function(o, pkg, profile) {
      stopifnot(identical(profile, "dsvert"))
      state$unpublished <- c(state$unpublished, pkg)
      state$events <- c(state$events, paste("unpublish", pkg, sep = ":"))
    },
    oadmin.remove_package = function(o, pkg, profile) {
      stopifnot(identical(profile, "default"))
      state$removed <- c(state$removed, pkg)
      state$events <- c(state$events, paste("remove", pkg, sep = ":"))
      state$installed[[pkg]] <- FALSE
    },
    dsadmin.install_local_package = function(o, path, profile) {
      stopifnot(identical(profile, "dsvert"))
      pkg <- sub("[.]tar[.]gz$", "", path)
      state$installed_paths <- c(state$installed_paths, path)
      state$events <- c(state$events, paste("install", pkg, sep = ":"))
      state$installed[[pkg]] <- TRUE
    }))

  install("opal")
  expect_identical(state$installed, c(dsVert = TRUE, dsVertClient = TRUE))
  expect_identical(
    state$installed_paths, c("dsVert.tar.gz", "dsVertClient.tar.gz"))
  expect_identical(state$removed, character())
  expect_identical(state$removed_methods, character())
  expect_identical(state$unpublished, character())
  expect_identical(state$dsadmin_probe_calls, 0L)

  state$installed_paths <- character()
  state$events <- character()
  install("opal")
  expect_identical(state$removed, c("dsVert", "dsVertClient"))
  expect_identical(state$removed_methods, c("dsVert", "dsVertClient"))
  expect_identical(state$unpublished, c("dsVert", "dsVertClient"))
  expect_identical(state$installed, c(dsVert = TRUE, dsVertClient = TRUE))
  expect_identical(state$dsadmin_probe_calls, 0L)
  expect_identical(state$events, c(
    "methods:dsVert", "unpublish:dsVert", "remove:dsVert",
    "methods:dsVertClient", "unpublish:dsVertClient",
    "remove:dsVertClient", "install:dsVert", "install:dsVertClient"))
})

test_that("Opal provisioner fails closed on package inventory anomalies", {
  bindings <- list(
    dsadmin.rm_package_methods = function(...) NULL,
    dsadmin.remove_package = function(...) NULL,
    oadmin.remove_package = function(...) NULL,
    dsadmin.install_local_package = function(...) NULL)

  malformed <- bindings
  malformed$oadmin.installed_package <- function(...) NA
  expect_error(
    .provision_install_packages(malformed)("opal"),
    "Cannot determine whether dsVert is installed", fixed = TRUE)

  unavailable <- bindings
  unavailable$oadmin.installed_package <- function(...) {
    stop("Rock inventory unavailable", call. = FALSE)
  }
  expect_error(
    .provision_install_packages(unavailable)("opal"),
    "Rock inventory unavailable", fixed = TRUE)

  removal_error <- bindings
  removal_error$oadmin.installed_package <- function(...) TRUE
  removal_error$oadmin.remove_package <- function(...) {
    stop("cluster unavailable", call. = FALSE)
  }
  expect_error(
    .provision_install_packages(removal_error)("opal"),
    paste0("Cannot remove stale dsVert from Rock cluster default: ",
           "cluster unavailable"),
    fixed = TRUE)

  not_removed <- bindings
  not_removed$oadmin.installed_package <- local({
    calls <- 0L
    function(...) {
      calls <<- calls + 1L
      TRUE
    }
  })
  expect_error(
    .provision_install_packages(not_removed)("opal"),
    "Opal did not remove the previous dsVert package", fixed = TRUE)

  not_installed <- bindings
  not_installed$oadmin.installed_package <- function(...) FALSE
  expect_error(
    .provision_install_packages(not_installed)("opal"),
    "Opal did not install dsVert in profile dsvert", fixed = TRUE)
})

test_that("Opal provisioner installs only missing server Imports", {
  state <- new.env(parent = emptyenv())
  state$installed <- TRUE
  state$install_calls <- character()
  install <- .provision_install_server_imports(list(
    oadmin.installed_package = function(o, pkg, profile) {
      stopifnot(identical(profile, "default"))
      state$installed
    },
    oadmin.install_package = function(o, pkg, profile) {
      stopifnot(identical(profile, "default"))
      state$install_calls <- c(state$install_calls, pkg)
      state$installed <- TRUE
    }))

  install("opal", "filelock", "default")
  expect_identical(state$install_calls, character())

  state$installed <- FALSE
  install("opal", "filelock", "default")
  expect_identical(state$install_calls, "filelock")
  expect_true(state$installed)
})

test_that("Opal provisioner resolves and validates the profile Rock cluster", {
  resolve <- .provision_profile_rock_cluster(list(
    dsadmin.profile = function(...) list(cluster = "rock-prod")))
  expect_identical(resolve("opal", "dsvert"), "rock-prod")

  for (bad in list(list(), list(cluster = NA_character_),
                   list(cluster = character()), list(cluster = "bad name"),
                   list(cluster = "../x"),
                   list(cluster = paste(rep("a", 65L), collapse = "")),
                   list(cluster = 1L))) {
    invalid <- .provision_profile_rock_cluster(list(
      dsadmin.profile = local({ value <- bad; function(...) value })))
    expect_error(
      invalid("opal", "dsvert"),
      "invalid Rock cluster for profile dsvert", fixed = TRUE)
  }

  unavailable <- .provision_profile_rock_cluster(list(
    dsadmin.profile = function(...) {
      stop("profile unavailable", call. = FALSE)
    }))
  expect_error(
    unavailable("opal", "dsvert"),
    "Cannot resolve the Rock cluster for profile dsvert: profile unavailable",
    fixed = TRUE)
})

test_that("Opal provisioner fails closed when an Import is not verified", {
  invalid <- .provision_install_server_imports(list(
    oadmin.installed_package = function(...) NA,
    oadmin.install_package = function(...) NULL))
  expect_error(
    invalid("opal", "filelock", "default"),
    "Invalid Rock inventory state", fixed = TRUE)

  unavailable <- .provision_install_server_imports(list(
    oadmin.installed_package = function(...) {
      stop("inventory unavailable", call. = FALSE)
    },
    oadmin.install_package = function(...) NULL))
  expect_error(
    unavailable("opal", "filelock", "default"),
    "Cannot inspect required server dependency filelock: inventory unavailable",
    fixed = TRUE)

  install_error <- .provision_install_server_imports(list(
    oadmin.installed_package = function(...) FALSE,
    oadmin.install_package = function(...) {
      stop("repository unavailable", call. = FALSE)
    }))
  expect_error(
    install_error("opal", "filelock", "default"),
    "Cannot install required server dependency filelock: repository unavailable",
    fixed = TRUE)

  post_false <- local({
    calls <- 0L
    .provision_install_server_imports(list(
      oadmin.installed_package = function(...) {
        calls <<- calls + 1L
        FALSE
      },
      oadmin.install_package = function(...) NULL))
  })
  expect_error(
    post_false("opal", "filelock", "default"),
    "Rock did not install required server dependency filelock", fixed = TRUE)

  post_error <- local({
    calls <- 0L
    .provision_install_server_imports(list(
      oadmin.installed_package = function(...) {
        calls <<- calls + 1L
        if (calls == 1L) FALSE else
          stop("postcondition unavailable", call. = FALSE)
      },
      oadmin.install_package = function(...) NULL))
  })
  expect_error(
    post_error("opal", "filelock", "default"),
    "Cannot verify required server dependency filelock: postcondition unavailable",
    fixed = TRUE)
})
