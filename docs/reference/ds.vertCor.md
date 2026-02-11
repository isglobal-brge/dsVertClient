# Privacy-Preserving Correlation for Vertically Partitioned Data

Computes the Pearson correlation matrix for vertically partitioned data
using Multiparty Homomorphic Encryption (MHE) with threshold decryption.

## Usage

``` r
ds.vertCor(
  data_name,
  variables,
  log_n = 12,
  log_scale = 40,
  datasources = NULL
)
```

## Arguments

- data_name:

  Character string. Name of the (aligned) data frame on each server.

- variables:

  A named list where each name corresponds to a server name and each
  element is a character vector of variable names from that server.

- log_n:

  Integer. CKKS ring dimension parameter (12, 13, or 14). Controls the
  number of slots: N/2 = 2^(logN-1). Default is 12 (2048 slots, fast).
  Use 13 or 14 for larger datasets.

- log_scale:

  Integer. CKKS scale parameter controlling precision. Default is 40
  (approximately 12 decimal digits of precision).

- datasources:

  DataSHIELD connection object or list of connections. If NULL, uses all
  available connections.

## Value

A list with class `"ds.cor"` containing:

- `correlation`: The full correlation matrix (p x p)

- `var_names`: Variable names in order

- `n_obs`: Number of observations

- `method`: `"MHE-CKKS-Threshold"` indicating the method used

- `servers`: Names of servers involved

- `local_correlations`: List of within-server correlation matrices

- `cross_correlations`: List of cross-server correlation matrices

- `mhe_params`: List with `log_n` and `log_scale` used

## Details

### Threshold MHE Protocol (6 phases)

This function implements a full multiparty homomorphic encryption
protocol with threshold decryption:

1.  **Key Generation**: Each server generates its own secret key share
    and a public key share. Party 0 also generates the Common Reference
    Polynomial (CRP) shared by all parties.

2.  **Key Combination**: Public key shares are combined into a
    Collective Public Key (CPK). Data encrypted under the CPK can only
    be decrypted with ALL servers cooperating.

3.  **Encryption**: Each server standardizes its columns (Z-scores) and
    encrypts them column-by-column under the CPK.

4.  **Local Correlation**: Within-server correlations are computed in
    plaintext (no encryption needed for data that stays on-server).

5.  **Cross-Server Correlation**: For each pair of servers (A, B),
    server A receives Enc(Z_B) and computes the element-wise product Z_A
    \* Enc(Z_B) homomorphically. The result is still encrypted and
    requires threshold decryption. Each server produces a partial
    decryption share, and the client fuses all shares to recover the
    inner product.

6.  **Assembly**: The client assembles the full p x p correlation matrix
    from local correlations (diagonal blocks) and cross-server
    correlations (off-diagonal blocks).

### Security Guarantees

- Client privacy:

  The client (researcher) CANNOT decrypt any individual data. It only
  receives partial decryption shares that are useless alone. Only the
  final aggregate statistics (correlation coefficients) are revealed
  after fusing ALL shares.

- Server privacy:

  Each server's raw data never leaves the server. Other servers only see
  encrypted columns (opaque ciphertexts).

- Collusion resistance:

  Even K-1 colluding servers cannot decrypt without the K-th server's
  key share. Full decryption requires cooperation from ALL K servers.

## Performance Notes

- `log_n = 12`: Up to 2048 observations, fastest

- `log_n = 13`: Up to 4096 observations

- `log_n = 14`: Up to 8192 observations, slowest

- Precision: approximately 10^-3 to 10^-4 error due to CKKS
  approximation

- The dominant cost is the threshold decryption loop (Phase 5), which
  requires one round-trip per server per cross-correlation element.

## References

Mouchet, C. et al. (2021). "Multiparty Homomorphic Encryption from
Ring-Learning-With-Errors". *Proceedings on Privacy Enhancing
Technologies (PETS)*.

Cheon, J.H. et al. (2017). "Homomorphic Encryption for Arithmetic of
Approximate Numbers". *ASIACRYPT 2017*.

## See also

[`ds.vertPCA`](https://isglobal-brge.github.io/dsVertClient/reference/ds.vertPCA.md)
for PCA analysis built on this function

## Examples

``` r
if (FALSE) { # \dontrun{
# Connect to Opal/DataSHIELD servers
connections <- DSI::datashield.login(builder$build())

# Align records across servers
ref_hashes <- ds.hashId("D", "patient_id", datasource = connections["server1"])
ds.alignRecords("D", "patient_id", ref_hashes$hashes,
                newobj = "D_aligned", datasources = connections)

# Define which variables are on which server
vars <- list(
  hospital_A = c("age", "bmi"),
  hospital_B = c("glucose", "systolic_bp")
)

# Compute privacy-preserving correlation
result <- ds.vertCor("D_aligned", vars)
print(result)
} # }
```
