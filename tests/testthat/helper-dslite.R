# DSLite integration tests for dsVertClient
# Tests the complete workflow with multiple simulated DataSHIELD servers

# Skip if DSLite or dsVert not available
skip_if_not_installed("DSLite")
skip_if_not_installed("dsVert")

library(DSLite)
library(dsVert)


# =============================================================================
# Helper function to set up DSLite environment
# =============================================================================

setup_dslite_env <- function(seed = 123) {

  set.seed(seed)
  n <- 100

  # Common IDs across all servers (shuffled differently on each)
  ids <- paste0("PAT", sprintf("%04d", 1:n))

  # Generate correlated predictors first (in original order, linked to IDs)
  age <- round(rnorm(n, 45, 15))
  weight <- round(rnorm(n, 70, 12), 1)
  height <- round(rnorm(n, 170, 10), 1)
  bmi <- round(weight / (height/100)^2, 1)  # Realistic BMI
  glucose <- round(rnorm(n, 100, 20), 1)
  cholesterol <- round(rnorm(n, 200, 40), 1)

  # Generate outcomes as functions of predictors (true relationships)
  # Continuous: linear combination
  outcome_cont <- round(2 + 0.05 * age + 0.02 * weight + 0.1 * bmi + rnorm(n, 0, 1), 2)

  # Binary: logistic function
  logit_p <- -5 + 0.05 * age + 0.03 * bmi + 0.01 * glucose
  outcome_binary <- rbinom(n, 1, plogis(logit_p))

  # Count: Poisson with log link
  log_rate <- 0.5 + 0.01 * age + 0.02 * bmi
  outcome_count <- rpois(n, exp(pmin(log_rate, 3)))  # Cap to avoid huge values

  # Positive continuous: for Gamma/IG
  outcome_positive <- round(pmax(0.1, exp(1 + 0.01 * age + 0.02 * bmi + rnorm(n, 0, 0.3))), 2)

  # Server 1: Demographics + outcomes (shuffled)
  order1 <- sample(n)
  data_server1 <- data.frame(
    patient_id = ids[order1],
    age = age[order1],
    weight = weight[order1],
    outcome_cont = outcome_cont[order1],
    outcome_binary = outcome_binary[order1],
    outcome_count = outcome_count[order1],
    outcome_positive = outcome_positive[order1],
    stringsAsFactors = FALSE
  )

  # Server 2: Clinical measurements + outcomes (shuffled differently)
  order2 <- sample(n)
  data_server2 <- data.frame(
    patient_id = ids[order2],
    height = height[order2],
    bmi = bmi[order2],
    outcome_cont = outcome_cont[order2],
    outcome_binary = outcome_binary[order2],
    outcome_count = outcome_count[order2],
    outcome_positive = outcome_positive[order2],
    stringsAsFactors = FALSE
  )

  # Server 3: Lab results + outcomes (shuffled differently)
  order3 <- sample(n)
  data_server3 <- data.frame(
    patient_id = ids[order3],
    glucose = glucose[order3],
    cholesterol = cholesterol[order3],
    outcome_cont = outcome_cont[order3],
    outcome_binary = outcome_binary[order3],
    outcome_count = outcome_count[order3],
    outcome_positive = outcome_positive[order3],
    stringsAsFactors = FALSE
  )

  # Create DSLite server with all tables
  dslite_server <- DSLite::newDSLiteServer(
    tables = list(
      server1 = data_server1,
      server2 = data_server2,
      server3 = data_server3
    )
  )

  # Configure with dsVert package
dslite_server$config(DSLite::defaultDSConfiguration(include = c("dsVert")))

  # Build login credentials
  builder <- DSI::newDSLoginBuilder()
  builder$append(server = "server1", url = "dslite_server",
                 table = "server1", driver = "DSLiteDriver")
  builder$append(server = "server2", url = "dslite_server",
                 table = "server2", driver = "DSLiteDriver")
  builder$append(server = "server3", url = "dslite_server",
                 table = "server3", driver = "DSLiteDriver")

  login_data <- builder$build()

  # Store server in global env for DSLite to find
  assign("dslite_server", dslite_server, envir = globalenv())

  list(
    login_data = login_data,
    n = n,
    ids = ids
  )
}

teardown_dslite_env <- function(conns) {
  if (!is.null(conns)) {
    try(DSI::datashield.logout(conns), silent = TRUE)
  }
  if (exists("dslite_server", envir = globalenv())) {
    rm("dslite_server", envir = globalenv())
  }
}
