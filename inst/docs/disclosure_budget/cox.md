# P3 disclosure budget - Cox PH

Current product route: `ds.vertCox()` dispatches to
`ds.vertCoxProfileNonDisclosive()`. The historical rank/permutation Cox path
and K>=3 Poisson-trick prototype are removed from the client API.

## Summary

| channel | size | to | content | tier |
|---|---|---|---|---|
| client | p-vector | analyst | aggregate score | scalar/model aggregate |
| client | p x p matrix | analyst | aggregate observed information | scalar/model aggregate |
| client | scalar | analyst | partial log-likelihood diagnostics | scalar/model aggregate |
| inter-server | encrypted/share-domain state | selected DCF parties | risk-set/profile working values | no plaintext row vector |

No per-observation time, event, rank, risk-set membership, eta, probability, or
residual vector is returned to the analyst. Event/risk-set work stays local or
in Ring127 share-domain arithmetic.

## Accepted Disclosure

The product route reveals the same model-scale objects expected from a Cox PH
fit: coefficients, optional covariance/standard errors, convergence metadata,
and scalar likelihood diagnostics. These are O(p^2) in the number of model
parameters and independent of n as returned values.

This is in the same accepted aggregate tier as `ds.vertCor()` and the GLM
Fisher/Hessian summaries already used by the package. It does not expose the
rank/permutation metadata that made the older Cox path unsuitable as a product
route.

## Removed Routes

The following routes are not product routes and are no longer reachable through
the client API:

- Cox rank/permutation orchestration based on `k2SetCoxTimesDS`,
  `k2ReceiveCoxMetaDS`, `k2ApplyCoxPermutationDS`, and the old Newton/Path-B
  helpers.
- The K>=3 Poisson/person-time Cox prototype.

Server-side primitives from those experiments are not listed in dsVert
`AggregateMethods`; the client no longer calls them.
