---
name: Bug report
about: Report a problem in dsVertClient (the user-facing client package)
labels: bug
---

## What happened

<!-- A clear description of the unexpected behaviour. -->

## What did you expect

<!-- What you thought would happen. -->

## Reproducer

```r
library(DSI); library(DSOpal); library(dsVertClient)
# Minimal client-side call that reproduces the issue:
# ds.vertGLM(...)
```

If you ran against DSLite / a local harness rather than Opal, please
note that: some PSI-related issues only surface in DSLite due to
shared-process state.

## Environment

- dsVertClient version: <output of `packageVersion("dsVertClient")`>
- dsVert version (server-side): <output of `packageVersion("dsVert")` on the server>
- R version: <output of `R.version.string`>
- Platform: <macOS / Linux / Windows + arch>
- DataSHIELD backend: <Opal x.y.z, DSLite, or local-harness>
- K (number of servers): <2 or 3>

## Disclosure check

- [ ] Output you receive does not include any per-observation values
      (only p-dimensional aggregates / scalars / coefficient vectors
      / per-cluster BLUPs).

If a per-observation value DID leak, this is a **security issue** and
should also be reported privately to **david.sarrat@isglobal.org**.

## Logs / output

```
<paste relevant error / verbose output here, with `verbose = TRUE`>
```
