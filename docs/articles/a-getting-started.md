# Getting Started with dsVertClient

## What is Vertical Data Partitioning?

In traditional **horizontal partitioning**, different institutions hold
data for *different patients* with the *same variables*. For example,
Hospital A has 1000 patients and Hospital B has 2000 different patients,
but both measure the same things (age, blood pressure, etc.).

**Vertical partitioning** is different: multiple institutions hold data
for the *same patients* but with *different variables*. This commonly
occurs when:

- A hospital has clinical data (diagnoses, treatments)
- A laboratory has biomarker measurements
- A research center has genomic data
- A government agency has demographic data

All these institutions have records for the same individuals, identified
by a common ID (e.g., national health number), but each holds different
pieces of information.

### Example: Three Institutions with Shared Patients

**Hospital A**

| ID  | Age | Weight |
|:---:|:---:|:------:|
| P1  | 45  |   70   |
| P2  | 52  |   85   |
| P3  | 38  |   62   |

**Laboratory B**

| ID  | Glucose | HDL |
|:---:|:-------:|:---:|
| P1  |   95    | 55  |
| P2  |   110   | 42  |
| P3  |   88    | 61  |

**Research Center C**

| ID  | Gene1 | Gene2 |
|:---:|:-----:|:-----:|
| P1  |  0.2  |  0.8  |
| P2  |  0.5  |  0.3  |
| P3  |  0.1  |  0.9  |

*Same patients (P1, P2, P3) but different variables at each
institution.*

## The Privacy Challenge

To analyze relationships between variables held by different
institutions (e.g., “Does glucose level correlate with genetic
markers?”), traditionally you would need to:

1.  Send all data to a central location
2.  Merge the datasets
3.  Perform the analysis

This approach has serious privacy and legal concerns:

- **Patient privacy**: Sensitive data leaves protected environments
- **Legal barriers**: GDPR and other regulations restrict data transfer
- **Trust issues**: Institutions may not trust each other with their
  data

## The DataSHIELD Solution

**DataSHIELD** is a framework that enables privacy-preserving federated
analysis. The key principle is:

> *“Bring the analysis to the data, not the data to the analysis”*

Instead of moving data, DataSHIELD:

1.  Sends analysis commands to each server
2.  Each server executes computations locally on its data
3.  Only **aggregate results** (never individual-level data) are
    returned
4.  The client combines these aggregates to produce final results

**dsVertClient** extends DataSHIELD specifically for **vertically
partitioned data**, implementing:

- Privacy-preserving record alignment via cryptographic hashing
- Distributed correlation and PCA using Block SVD
- Distributed GLM fitting using Block Coordinate Descent

------------------------------------------------------------------------

## Our Test Environment

For this tutorial, we have a simulated environment with **3
institutions** and **200 patients**. Each institution holds different
variables for the same patients:

| Institution | Role       | Variables            |
|-------------|------------|----------------------|
| **inst_A**  | Hospital   | age, weight          |
| **inst_B**  | Clinic     | height, bmi          |
| **inst_C**  | Laboratory | glucose, cholesterol |

All institutions also store outcome variables (blood pressure, diabetes
status, hospital visits, costs) which will be used for GLM examples.

### Connecting to Servers

We connect to the DataSHIELD servers and assign the data to a symbol
`D`:

``` r

# Build login credentials for each server
builder <- newDSLoginBuilder()
builder$append(server = "inst_A", url = "dslite_server",
               table = "inst_A", driver = "DSLiteDriver")
builder$append(server = "inst_B", url = "dslite_server",
               table = "inst_B", driver = "DSLiteDriver")
builder$append(server = "inst_C", url = "dslite_server",
               table = "inst_C", driver = "DSLiteDriver")

# Connect and assign data to symbol "D" on each server
connections <- datashield.login(builder$build(), assign = TRUE, symbol = "D")
```

### What Each Institution Sees

**Institution A (Hospital)**

| patient_id    | age | weight |
|:--------------|----:|-------:|
| PATIENT_00093 |  49 |   94.3 |
| PATIENT_00091 |  31 |   92.8 |
| PATIENT_00013 |  52 |   65.8 |
| PATIENT_00178 |  59 |   54.6 |
| PATIENT_00066 |  52 |   60.5 |

**Institution B (Clinic)**

| patient_id    | height |  bmi |
|:--------------|-------:|-----:|
| PATIENT_00003 |  159.4 | 28.2 |
| PATIENT_00167 |  167.0 | 22.4 |
| PATIENT_00105 |  169.5 | 21.7 |
| PATIENT_00045 |  174.3 | 24.7 |
| PATIENT_00132 |  151.4 | 32.8 |

**Institution C (Laboratory)**

| patient_id    | glucose | cholesterol |
|:--------------|--------:|------------:|
| PATIENT_00134 |   136.6 |       180.9 |
| PATIENT_00099 |   131.5 |       216.2 |
| PATIENT_00184 |   101.7 |       237.5 |
| PATIENT_00145 |   125.6 |       223.6 |
| PATIENT_00041 |   119.1 |       223.1 |

Notice how:

- Each institution has **different variables**
- The **patient order is different** at each institution
- The same patient (e.g., PATIENT_00093) appears at different row
  positions

This is realistic: institutions store data independently and don’t
coordinate their internal ordering.

------------------------------------------------------------------------

## Step 1: Validating Identifier Formats

Before aligning records, we verify that patient identifiers have
consistent formats across all institutions. Common problems include:

- Different ID formats (e.g., “001” vs “1” vs “P001”)
- Leading/trailing whitespace
- Missing values or duplicates

``` r

validation_result <- ds.validateIdFormat(
  data_name = "D",
  id_col = "patient_id",
  datasources = connections
)

validation_result
```

``` bg-light
#> 
#> Identifier Format Validation
#> ============================
#> 
#> Overall Status: VALID 
#> Format Consistency: Yes 
#> 
#> Server Details:
#>  server n_obs n_unique n_missing  id_class
#>  inst_A   200      200         0 character
#>  inst_B   200      200         0 character
#>  inst_C   200      200         0 character
#>                                                  format_signature
#>  8d404dbb1d0a33a5bf77b9f34bfabf2e0930831a8dfd5dfa57302c6bbb345ac1
#>  8d404dbb1d0a33a5bf77b9f34bfabf2e0930831a8dfd5dfa57302c6bbb345ac1
#>  8d404dbb1d0a33a5bf77b9f34bfabf2e0930831a8dfd5dfa57302c6bbb345ac1
```

**What we’re looking for:**

- **n_obs**: Same number of observations across all institutions
- **n_unique**: Should equal n_obs (no duplicates)
- **n_missing**: Should be 0
- **format_signature**: Should be identical across institutions

------------------------------------------------------------------------

## Step 2: Privacy-Preserving Record Alignment

### The Problem

Even though all institutions have data for the same patients, the
records are in different order. To analyze relationships across
variables (e.g., correlation between age and glucose), we need
corresponding records to be in the same position.

**But we cannot simply share patient IDs!** That would violate privacy.

### The Solution: Cryptographic Hashing

Instead of sharing actual IDs, we:

1.  **Hash all IDs** using SHA-256 (a one-way cryptographic function)
2.  **Share only the hashes** - these cannot be reversed to reveal
    original IDs
3.  **Reorder records** at each institution to match a reference order

| Original ID     | →   | SHA-256 Hash                         |
|-----------------|-----|--------------------------------------|
| `PATIENT_00042` | →   | `a7f3b9c2d8e1f4a5...` (64 hex chars) |

The hash is shared safely - the original ID cannot be recovered from it.

### Step 2a: Get Reference Hashes

First, we obtain hashed identifiers from one institution (which becomes
our reference):

``` r

reference_hashes <- ds.hashId(
  data_name = "D",
  id_col = "patient_id",
  algo = "sha256",
  datasource = connections["inst_A"]
)
```

We retrieved 200 hashed identifiers.

### Step 2b: Align All Institutions

Now we send the reference hashes to all institutions. Each institution:

1.  Computes hashes of its own IDs
2.  Reorders its records to match the reference hash order
3.  Removes any records not in the reference set

``` r

ds.alignRecords(
  data_name = "D",
  id_col = "patient_id",
  reference_hashes = reference_hashes$hashes,
  newobj = "D_aligned",
  datasources = connections
)
```

### Verifying Alignment

``` r

alignment_counts <- sapply(names(connections), function(inst) {
  count_result <- DSI::datashield.aggregate(
    connections[inst],
    call("getObsCountDS", "D_aligned")
  )
  count_result[[1]]$n_obs
})

data.frame(Institution = names(alignment_counts), Records = alignment_counts, row.names = NULL)
```

``` bg-light
#>   Institution Records
#> 1      inst_A     200
#> 2      inst_B     200
#> 3      inst_C     200
```

All institutions now have:

- The **same number of records**
- Records in the **same order** (corresponding to the same patients)
- Ready for **cross-institutional analysis**

------------------------------------------------------------------------

## Next Steps

Now that records are aligned, we can perform statistical analyses:

- [**Statistical
  Analysis**](https://isglobal-brge.github.io/dsVertClient/articles/b-statistical-analysis.md):
  Correlation, PCA, and GLMs
- [**Methodology**](https://isglobal-brge.github.io/dsVertClient/articles/c-methodology.md):
  Mathematical details behind the algorithms
