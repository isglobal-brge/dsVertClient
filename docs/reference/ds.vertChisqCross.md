# Cross-server chi-square (one-hot + Beaver cross-products)

Pearson chi-square and Fisher-exact tests on a two-way contingency table
where the row variable is held on one server and the column variable is
held on another. Uses the existing Ring63 Beaver cross-product
infrastructure: one-hot indicator matrices are constructed server-side
via
[`dsvertOneHotDS`](https://rdrr.io/pkg/dsVert/man/dsvertOneHotDS.html),
then the \\K \times L\\ cell counts \\n\_{kl} = \sum_i X\_{ik} Y\_{il}\\
are obtained by Beaver dot products on the shared indicator vectors. The
client never sees any \\n\\-length indicator vector; only the \\K \times
L\\ aggregate table.

## Usage

``` r
ds.vertChisqCross(
  data,
  var1,
  var2,
  correct = TRUE,
  fisher = FALSE,
  datasources = NULL,
  verbose = TRUE
)
```

## Arguments

- data:

  Aligned data-frame name.

- var1:

  Row variable (categorical or numeric treated as factor).

- var2:

  Column variable.

- correct:

  Apply Yates continuity correction (2x2 only).

- fisher:

  Also return a Fisher-exact p-value (via R's `fisher.test` applied to
  the reconstructed aggregate table).

- datasources:

  DataSHIELD connection object.

- verbose:

  Print progress.

## Value

An object of class `ds.vertChisq` (same print method as the same-server
helper) with components `observed`, `expected`, `chisq`, `df`,
`p_value`, `fisher_p` (if requested), `n`, `row_levels`, `col_levels`.
