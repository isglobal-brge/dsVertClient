## Summary

<!-- One paragraph: what does this PR do? -->

## Scope check

- [ ] R CMD check passes — `Status: OK` (0 ERRORs / 0 WARNINGs / 0 NOTEs).
- [ ] `testthat` passes — full suite green, new tests added for any
      new method.
- [ ] `bash scripts/quick_impl_check.sh` passes (8/8 L1 + go).
- [ ] If the change touches federated estimators,
      `bash scripts/continuous_validation.sh medium` passes.
- [ ] Roxygen docs cover every new exported argument; markdown
      bracket-link traps avoided (`\code{}` for inline code, no bare
      `[token]`); `call(name = "fn", ...)` style preserved.
- [ ] Any new exported `ds.vert*` function is added to the
      appropriate `reference:` section in `_pkgdown.yml`
      (pkgdown errors hard if the index is incomplete).

## Disclosure check

- [ ] Client-side output remains aggregate-only (p-vectors, scalars,
      coefficient matrices, per-cluster BLUPs, χ² statistics, etc.).
- [ ] No new observation-level reveal at the client.
- [ ] No new inter-server leakage tier beyond the disclosure ledger.
- [ ] K=2 path works.

## Validation evidence

| Method | max\|Δβ\| | Reference | Acceptance band |
|---|---|---|---|
| <method> | <e.g. 3.18e-04> | <`stats::glm`> | STRICT / SUB-NOISE / PRACTICAL |

## Related

Closes #<issue>.
