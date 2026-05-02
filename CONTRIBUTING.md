# Contributing to dsVertClient

Thank you for considering a contribution. dsVertClient is the
client-side half of a tightly-paired R + Go codebase (with the
server-side companion [`dsVert`](https://github.com/isglobal-brge/dsVert))
that implements federated MPC primitives for vertically partitioned
DataSHIELD analyses. Contributions need to clear three quality bars:
**correctness** (matches centralised R within documented bounds),
**non-disclosure** (no observation-level data ever leaves a server),
and **reviewer-grade hygiene** (`R CMD check` clean, `testthat` green).

## Repository layout

- `R/` — exported `ds.vert*` user-facing functions. Auto-generated
  `NAMESPACE` via `roxygen2` (regenerable with
  `roxygen2::roxygenize(".", roclets = c("rd", "namespace"))`).
- `tests/testthat/` — client-side unit tests including DSLite mocks
  and input-validation probes.
- `vignettes/` — knitr R-Markdown vignettes (`Rbuildignored` from
  the source tarball; rendered for the pkgdown site).
- `_pkgdown.yml` — pkgdown site config; reference index must list
  every exported `ds.vert*` topic.

## Workflow

1. **Fork** and create a feature branch off `main`.
2. **Match existing style**: don't reformat adjacent code. New `R/`
   files use `roxygen2` markdown mode (set in `DESCRIPTION`); wrap
   math in `\eqn{}` / `\preformatted{}`, never bare LaTeX macros.
   Use `call(name = "fn", ...)` (the explicit `name = ...` form) in
   any `DSI::datashield.aggregate` / `assign` call expressions.
3. **Document every exported argument** — `R CMD check --as-cran`
   gates on this. Use `\code{}` (not bracket-link prose) inside
   markdown roxygen.
4. **Run the local validation gate** before pushing:
   ```sh
   bash scripts/quick_impl_check.sh        # ~3 min  (8 L1 probes + go test)
   bash scripts/continuous_validation.sh medium   # ~33 min  (L1 + L2 local-distributed)
   ```
5. **R CMD check**:
   ```sh
   R CMD build --no-build-vignettes dsVertClient
   R CMD check dsVertClient_*.tar.gz --no-manual --as-cran
   ```
   Acceptance: **`Status: OK`** (0 ERRORs / 0 WARNINGs / 0 NOTEs).
6. **`testthat`**: `Rscript -e 'devtools::test()'` — full suite must
   PASS. Add tests for any new method.
7. **pkgdown**: any new exported `ds.vert*` function MUST be added to
   the appropriate `reference:` section in `_pkgdown.yml` —
   `pkgdown::build_site()` errors hard if the index is incomplete.
8. **Submit a pull request** referencing any related issue. The CI
   workflow (`.github/workflows/R-CMD-check.yaml`) will reproduce the
   above check on clean Linux + macOS runners.

## Disclosure invariants (do not violate)

1. The client must only see p-dimensional aggregate sums, scalar
   deviances / log-likelihoods, coefficient vectors, scalar test
   statistics, or already-agreed-aggregate correlation / histogram
   tables. **No `n`-dimensional vectors ever reconstructed at the
   client.**
2. Beyond the existing Ring63 / Ring127 shares + transport-encrypted
   blobs, permissible new inter-server leakage is restricted to the
   patterns already listed in the disclosure ledger (paper §VI Table
   3). Any new leakage tier requires a ledger row + reviewer sign-off.
3. K=2 must work with exactly two servers, both acting as DCF parties.
   No method may require a third party for correctness.

## Acceptance bands

A method is reviewer-shippable in one of three bands:

| Band | Bound | Where it belongs |
|---|---|---|
| STRICT | max\|Δβ\| < 1e-3 vs centralised reference | Paper §V.A row, no caveat |
| SUB-NOISE | σ-probe ratio ≥ 100× the per-fit Wald SE | Paper §V.B row + (H10) sub-noise margin |
| PRACTICAL | Empirical max\|Δβ\| inside a peer-reviewed theoretical floor | Paper §V.A row with formal-bound citation |

Anything below PRACTICAL needs further work before merge.

## Reporting issues

Use the GitHub issue tracker:
<https://github.com/isglobal-brge/dsVertClient/issues>

For the server-side companion package, see
<https://github.com/isglobal-brge/dsVert>.
