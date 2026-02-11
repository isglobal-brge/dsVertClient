#' =============================================================================
#' Complete Workflow: Vertically Partitioned Data Analysis with dsVert
#' =============================================================================
#'
#' This script demonstrates the complete workflow for analyzing vertically
#' partitioned data using dsVert and dsVertClient packages with DSLite.
#'
#' Author: David Sarrat González
#' Date: 2024
#' =============================================================================

# Load packages
library(DSI)
library(DSLite)
library(dsVert)
library(dsVertClient)

# =============================================================================
# STEP 1: Create Simulated Vertically Partitioned Data
# =============================================================================

message("STEP 1: Creating simulated vertically partitioned data\n")

set.seed(2024)
n <- 500  # Number of patients

# Generate patient identifiers
patient_ids <- paste0("PATIENT_", sprintf("%05d", 1:n))

# Generate correlated clinical variables
age <- rnorm(n, mean = 55, sd = 12)
sex <- rbinom(n, 1, 0.48)  # 0 = female, 1 = male
weight <- rnorm(n, mean = 75 + 10 * sex, sd = 15)
height <- rnorm(n, mean = 165 + 12 * sex, sd = 8)
bmi <- weight / ((height / 100)^2)
sbp <- rnorm(n, mean = 120 + 0.5 * age, sd = 15)  # Systolic BP
glucose <- rnorm(n, mean = 90 + 0.3 * age + 5 * bmi/25, sd = 20)
cholesterol <- rnorm(n, mean = 180 + 0.8 * age, sd = 35)

# Generate outcome (cardiovascular risk score)
# True coefficients
true_beta <- c(
  age = 0.03,
  sex = 0.5,
  weight = 0.01,
  height = -0.02,
  bmi = 0.1,
  sbp = 0.02,
  glucose = 0.01,
  cholesterol = 0.005
)

linear_pred <- true_beta["age"] * scale(age) +
               true_beta["sex"] * sex +
               true_beta["weight"] * scale(weight) +
               true_beta["height"] * scale(height) +
               true_beta["bmi"] * scale(bmi) +
               true_beta["sbp"] * scale(sbp) +
               true_beta["glucose"] * scale(glucose) +
               true_beta["cholesterol"] * scale(cholesterol)

# Continuous outcome (risk score)
risk_score <- as.vector(linear_pred) + rnorm(n, sd = 0.5)

# Binary outcome (high risk: yes/no)
prob_high_risk <- 1 / (1 + exp(-linear_pred))
high_risk <- rbinom(n, 1, as.vector(prob_high_risk))

# Create three hospital datasets (vertically partitioned)
# Hospital A: Demographics
hospital_A <- data.frame(
  patient_id = patient_ids,
  age = age,
  sex = sex,
  risk_score = risk_score,
  high_risk = high_risk
)

# Hospital B: Physical measurements (SHUFFLED to test alignment)
hospital_B <- data.frame(
  patient_id = patient_ids,
  weight = weight,
  height = height,
  bmi = bmi,
  risk_score = risk_score,
  high_risk = high_risk
)[sample(n), ]

# Hospital C: Lab results (SHUFFLED)
hospital_C <- data.frame(
  patient_id = patient_ids,
  sbp = sbp,
  glucose = glucose,
  cholesterol = cholesterol,
  risk_score = risk_score,
  high_risk = high_risk
)[sample(n), ]

message("Created 3 hospital datasets:")
message(sprintf("  Hospital A (Demographics): %d patients, %d variables",
                nrow(hospital_A), ncol(hospital_A) - 3))
message(sprintf("  Hospital B (Physical): %d patients, %d variables (SHUFFLED)",
                nrow(hospital_B), ncol(hospital_B) - 3))
message(sprintf("  Hospital C (Lab): %d patients, %d variables (SHUFFLED)",
                nrow(hospital_C), ncol(hospital_C) - 3))

# =============================================================================
# STEP 2: Set Up DSLite Environment
# =============================================================================

message("\nSTEP 2: Setting up DSLite environment\n")

# Create DSLite server
server <- newDSLiteServer(
  tables = list(
    hospitalA = hospital_A,
    hospitalB = hospital_B,
    hospitalC = hospital_C
  )
)

# Configure with dsVert methods
server$config(defaultDSConfiguration(include = c("dsVert")))

# Build login credentials
builder <- newDSLoginBuilder()
builder$append(server = "hospitalA", url = "server", table = "hospitalA",
               driver = "DSLiteDriver")
builder$append(server = "hospitalB", url = "server", table = "hospitalB",
               driver = "DSLiteDriver")
builder$append(server = "hospitalC", url = "server", table = "hospitalC",
               driver = "DSLiteDriver")

# Connect
conns <- datashield.login(logins = builder$build(), assign = TRUE, symbol = "D")
message("Connected to: ", paste(names(conns), collapse = ", "))

# =============================================================================
# STEP 3: Record Alignment
# =============================================================================

message("\nSTEP 3: Aligning records across hospitals\n")

# Get reference hashes from Hospital A
ref_hashes <- ds.hashId("D", "patient_id", algo = "sha256",
                        datasource = conns["hospitalA"])
message(sprintf("Generated %d reference hashes from Hospital A", ref_hashes$n))

# Align all hospitals
alignment <- ds.alignRecords(
  data_name = "D",
  id_col = "patient_id",
  reference_hashes = ref_hashes$hashes,
  newobj = "D_aligned",
  datasources = conns
)

# =============================================================================
# STEP 4: Correlation Analysis
# =============================================================================

message("\nSTEP 4: Computing correlation matrix\n")

# Define variable mapping
variables <- list(
  hospitalA = c("age", "sex"),
  hospitalB = c("weight", "height", "bmi"),
  hospitalC = c("sbp", "glucose", "cholesterol")
)

# Compute correlation matrix
cor_matrix <- ds.vertCor("D_aligned", variables, datasources = conns)

message("Correlation matrix (all 8 variables across 3 hospitals):")
print(round(cor_matrix, 2))

# =============================================================================
# STEP 5: Principal Component Analysis
# =============================================================================

message("\nSTEP 5: Performing PCA\n")

# Perform PCA
pca_result <- ds.vertPCA("D_aligned", variables, n_components = 4,
                         datasources = conns)
print(pca_result)

# =============================================================================
# STEP 6: Fit GLM (Linear Regression)
# =============================================================================

message("\nSTEP 6: Fitting Linear Regression (Gaussian GLM)\n")

# Define predictor variables per server
x_vars <- list(
  hospitalA = c("age", "sex"),
  hospitalB = c("weight", "height", "bmi"),
  hospitalC = c("sbp", "glucose", "cholesterol")
)

# Fit model
model_gaussian <- ds.vertGLM(
  data_name = "D_aligned",
  y_var = "risk_score",
  x_vars = x_vars,
  family = "gaussian",
  max_iter = 200,
  tol = 1e-6,
  lambda = 1e-4,
  verbose = TRUE,
  datasources = conns
)

print(model_gaussian)

# =============================================================================
# STEP 7: Fit GLM (Logistic Regression)
# =============================================================================

message("\nSTEP 7: Fitting Logistic Regression (Binomial GLM)\n")

model_binomial <- ds.vertGLM(
  data_name = "D_aligned",
  y_var = "high_risk",
  x_vars = x_vars,
  family = "binomial",
  max_iter = 200,
  tol = 1e-6,
  lambda = 1e-4,
  verbose = TRUE,
  datasources = conns
)

print(model_binomial)

# =============================================================================
# STEP 8: Compare with Centralized Analysis
# =============================================================================

message("\nSTEP 8: Comparison with centralized analysis\n")

# Merge data (pretending we could centralize - for validation only)
full_data <- data.frame(
  age = age,
  sex = sex,
  weight = weight,
  height = height,
  bmi = bmi,
  sbp = sbp,
  glucose = glucose,
  cholesterol = cholesterol,
  risk_score = risk_score,
  high_risk = high_risk
)

# Centralized correlation
cor_centralized <- cor(full_data[, 1:8])
message("Max correlation difference (distributed vs centralized): ",
        round(max(abs(cor_matrix - cor_centralized)), 6))

# Centralized GLM
glm_centralized <- glm(risk_score ~ age + sex + weight + height + bmi +
                       sbp + glucose + cholesterol - 1,
                       data = full_data, family = gaussian())

message("\nCoefficient comparison (Gaussian GLM):")
comparison <- data.frame(
  Variable = names(coef(model_gaussian)),
  Distributed = round(coef(model_gaussian), 4),
  Centralized = round(coef(glm_centralized), 4)
)
comparison$Ratio <- round(comparison$Distributed / comparison$Centralized, 3)
print(comparison)

# =============================================================================
# CLEANUP
# =============================================================================

message("\n=============================================================================")
message("WORKFLOW COMPLETE")
message("=============================================================================")
message("\nSummary:")
message(sprintf("- Processed %d patients across %d hospitals", n, 3))
message(sprintf("- Aligned records with 100%% match rate"))
message(sprintf("- Computed %dx%d correlation matrix", 8, 8))
message(sprintf("- Performed PCA with %d components (%.1f%% variance explained)",
                4, pca_result$cumulative_pct[4]))
message(sprintf("- Fitted Gaussian GLM (%s after %d iterations)",
                ifelse(model_gaussian$converged, "converged", "did not converge"),
                model_gaussian$iterations))
message(sprintf("- Fitted Binomial GLM (%s after %d iterations)",
                ifelse(model_binomial$converged, "converged", "did not converge"),
                model_binomial$iterations))

datashield.logout(conns)
message("\nDisconnected from all servers.")
