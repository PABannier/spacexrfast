# Fast RCTD Doublet Profile Notes

This package implements the first production-oriented fast path from
`IMPLEMENTATION_PLAN.md`.

Implemented:

- `run.RCTD.fast.doublet()` drop-in API.
- Adaptive chunked dense extraction over the RCTD regression gene list.
- Upstream `decompose_full()` retained for exact all-type weights and candidate
  selection.
- C++ likelihood interpolation, one- and two-type IRWLS, candidate pair loop,
  singlet scoring, final selected-pair refit, and classification.
- Parallel processing over cells via standard C++ threads.
- Optional compact diagnostics through `return_diagnostics`.

Validation:

- `tests/testthat/test-fast-doublet.R` compares categorical outputs and selected
  types against upstream `process_bead_doublet()` and `fitPixels()`.
- `benchmarks/xenium_smoke.R` samples a local Xenium H5 matrix from
  `~/Documents/10x-examples` and compares vanilla and fast outputs on a small
  smoke slice.

Known exactness boundary:

- The full all-type fit is still upstream R code.
- Two-type QP updates use a closed-form active-set solver instead of calling
  `quadprog` in workers. Categorical labels and selected types match the
  synthetic validation cases; floating-point weights show small drift around
  `1e-5` to `1e-4` in smoke tests.
